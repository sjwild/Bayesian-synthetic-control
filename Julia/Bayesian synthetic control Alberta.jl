using Pkg

Pkg.activate("COVID_19")

ENV["CMDSTAN"] = expanduser("~/cmdstan/")

using Plots, StatsPlots
using CSV, DataFrames, HTTP
using Dates, Measures
using Turing, ReverseDiff, Memoization
using StanSample


# For adjustments
cox_adjustment(x::Float64, σ::Float64) = exp(x + (σ^2)/2)
function cox_adjustment!(X::DataFrame, σ::Vector{Float64})
  I, J = size(X)

  for i in 1:I
    for j in 1:J
      X[i, j] = cox_adjustment(X[i, j], sigma[i])
    end
  end

  return X

end


# Load data for Canada
df_canada = CSV.read(HTTP.get(url_canada).body, DataFrame)
df_nyt = CSV.read(HTTP.get(url_nyt).body, DataFrame; normalizenames = true)

rename!(df_canada, :prname => :Province_State,
                   :numtoday => :cases,
                   :numconf => :total_cases)
df_canada = df_canada[ in([1, 99]).(df_canada.pruid) .== false, :]



# Load data for US
url_cdc = "https://data.cdc.gov/api/views/9mfq-cb36/rows.csv?accessType=DOWNLOAD"
df_cdc = CSV.read(HTTP.get(url_cdc).body, DataFrame; normalizenames = true)
rename!(df_cdc, :submission_date => :date,
                :new_case => :cases,
                :tot_cases => :total_cases)
df_cdc.date = replace.(df_cdc.date, "/" => "-")
df_cdc.date = Date.(df_cdc.date, dateformat"mm-dd-yyyy")


# get full state names and combine with CDC data
state_url = "https://github.com/jasonong/List-of-US-States/raw/master/states.csv"
state_names = CSV.read(HTTP.get(state_url).body, DataFrame; normalizenames = true)
rename!(state_names, :Abbreviation => :state,
                     :State => :Province_State)
df_cdc = leftjoin(df_cdc, state_names, on = :state)
dropmissing!(df_cdc, :Province_State)



# Load stringency index
url_stringency = "https://github.com/OxCGRT/covid-policy-tracker/raw/master/data/OxCGRT_latest_withnotes.csv"
df_stringency = CSV.read(HTTP.get(url_stringency).body, DataFrame; normalizenames = true)
dropmissing!(df_stringency, :RegionName)
filter!(row -> in(["Canada", "United States"]).(row.CountryName), df_stringency)

df_stringency.Date = Date.(string.(df_stringency.Date), dateformat"yyyymmdd")
filter!(row -> row.Date .== Date(2021, 07, 01), df_stringency)


# get names of all provinces and states that have a higher stringency on July 1,
# which is the date AB "opened for summer"
ab = df_stringency.StringencyIndex[df_stringency.RegionName .== "Alberta"]

stringency_idx = df_stringency.RegionName[df_stringency.StringencyIndex .> ab]



# Build dataframe combining Canada and CDC data, then subset
# so that we only get dates 90 days before or after July 1
df = append!(df_canada[:, [:date, :Province_State, :total_cases]], df_cdc[:, [:date, :Province_State, :total_cases]])
start_date = Date(2021, 07, 01) - Day(97)
end_date = Date(2021, 07, 01) + Day(90)
filter!(row -> start_date .≤ row.date .≤ end_date, df)

# Set up newdate so that it can become column names
df.newdate = [join(["X", df.date[i]]) for i in 1:size(df, 1)]
df.newdate = replace.(df.newdate, "-" => "_")
df = unstack(df[:, [:newdate, :Province_State, :total_cases]], :newdate, :total_cases)

# add in population totals to compute cases per 100k
df = innerjoin(df, pop, on = :Province_State)


# get 7 day moving average
y_raw = diff(Matrix{Float64}(df[df.Province_State .== "Alberta", 2:(end-1)]), dims = 2) ./ df[df.Province_State .== "Alberta", end] * scale
X_raw =  diff(Matrix{Float64}(df[df.Province_State .!= "Alberta", 2:(end-1)]), dims = 2) ./ df[df.Province_State .!= "Alberta", end] * scale


# Transform into 7-day moving averages. The problem otherwise is that there is a lot of zeros
# owing to non-reporting on holidays/weekends.
# This makes it difficult for the models to find appropriate weights
y = Vector{Float64}()
X = Matrix{Float64}(undef, (size(X_raw, 1), size(X_raw, 2) - 6))
for j in 7:length(y_raw)
  push!(y, sum(y_raw[(j-6):j])) / 7
  
  for i in 1:size(X_raw, 1)
    X[i, (j-6)] = sum(X_raw[i, (j-6):j]) / 7
  end

end

p = plot(legend = :false,
             bottommargin = 5mm,
             rightmargin = 5mm,
             ylabel = "Cases per 100k population")
for i in 1:63
  plot!(p, X[i, :], lc = :grey75)
end
plot!(p, y, lc = :red, linewidth = 2,
      title = "7-day moving avergage of cases in Alberta,\nbefore and after July 1, 2021",
      titleposition = :left)



indx = [sum(X[i, :] .≤ 0) for i in 1:63]
X = X[indx .== 0, :]
Province_State = df.Province_State[df.Province_State .!= "Alberta"]
Province_State = Province_State[indx .== 0]
Province_State = in(stringency_idx).(Province_State)

X = X[Province_State, :]

y_pre = y[1:90]
y_post = y[91:end]
X_pre = X[:, 1:90]
X_post = X[:, 91:end]


stan_data = Dict(
  
  "N_pre" => size(X_pre, 2), # Numbers of pre observations
  "N_post" => size(X_post, 2), # Number of post observations
  "p" => size(X_pre, 1), # Size of donor pool
  "y_pre" => log.(y_pre),
  "X_pre" => log.(X_pre),
  "X_post" => log.(X_post),
  
  "scale_alpha" => 1,
  "scale_g" => 1,
  
  "slab_scale" => 1,
  "slab_df" => 5,
  
  "nu_g" => 3,
  "nu_l" => 3
  
)


# load model
tmpdir = joinpath(@__DIR__, "tmp")

# run model
sm = SampleModel("bscm_alberta", read("Stan/bscm_horseshoe_modified.stan", String), tmpdir)
rc = stan_sample(sm; 
                 data = stan_data,
                 num_samples = 1000,
                 num_warmups = 1000,
                 num_chains = 4,
                 delta = 0.999,
                 max_depth = 15,
                 seed = 33294)





# check diagnostics
if success(rc)
    diagnose(sm)
end
  
# Extract parameters and transform them
params_out = DataFrame(read_samples(sm))
sigma = params_out.sigma
gq = stan_generate_quantities(sm)
  
sigma = params_out.sigma
y_fit = gq[:, 1:length(y_pre)]
y_post_est = gq[:, (length(y_pre) + 1):(length(y_pre) + length(y_post))]
  
  
# transform to correct scale
cox_adjustment!(y_fit, sigma)
y_post_est = exp.(y_post_est)
  

# Build series for plotting.
# Start with pre-period fitted values
y_fit_m = [mean(y_fit[:, i]) for i in 1:size(y_fit, 2)]
y_fit_ll = [quantile(y_fit[:, i], 0.025) for i in 1:size(y_fit, 2)]
y_fit_ul = [quantile(y_fit[:, i], 0.975) for i in 1:size(y_fit, 2)]

# Next build July 1 and after
y_post_m = [mean(y_post_est[:, i]) for i in 1:size(y_post_est, 2)]
y_post_ll = [quantile(y_post_est[:, i], 0.025) for i in 1:size(y_post_est, 2)]
y_post_ul = [quantile(y_post_est[:, i], 0.975) for i in 1:size(y_post_est, 2)]

# combine values into one vector each
synth_m = vcat(y_fit_m, y_post_m)
synth_ll = vcat(y_fit_ll, y_post_ll)
synth_ul = vcat(y_fit_ul, y_post_ul)
  

# dates
date_list = [(start_date + Day(7)):Day(1):end_date;]
  
AB_trend = plot(date_list, synth_m, 
                ribbon = (synth_m - synth_ll, synth_ul - synth_m), 
                lc = :red, fill = :red,
                label = "Estimated trend",
                title = "Synthetic control estimates: Alberta 7-day moving average",
                ylabel = "Num cases (per 100k)",
                lw = 3,
                title_position = :left,
                size = (800, 500),
                xrotation = 45,
                legend = :bottomleft,
                legend_foreground_color = nothing,
                bottom_margin = 8mm,
                left_margin = 5mm,
                tickfontsize = 12,
                titlefontsize = 16,
                legendfontsize = 10,
                guidefontsize = 14)
plot!(AB_trend, date_list, vcat(y_pre, y_post), 
      label = "Actual", 
      lc = :black,
      lw = 3)
vline!(AB_trend, [Date(2021, 07, 01)],
       lc = :grey, linestyle = :dash,
       linealpha = 0.5,
       label = "")
annotate!(AB_trend, Date(2021, 07, 03), 325, 
          StatsPlots.text("Restrictions dropped",
          10, :left))
png(AB_trend, "Images/Alberta/stronger_stringency_donor_pool")
  
AB_trend


