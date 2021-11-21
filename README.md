This repo is an attempt to use the code of [Lee, Kim, and Gupta (2020)](https://journals.sagepub.com/doi/full/10.1177/0022243720936230) to learn Bayesian synthetic control. I attempt to use their model to duplicate the results from Scott Cunningham's [Mixtape](https://mixtape.scunning.com/synthetic-control.html). You can find the code for Lee, Kim, and Gupta's (2020) in their supplementary materials.

A second Stan script uses the code of [Piironen and Vehtari (2017)](https://projecteuclid.org/journals/electronic-journal-of-statistics/volume-11/issue-2/Sparsity-information-and-regularization-in-the-horseshoe-and-other-shrinkage/10.1214/17-EJS1337SI.full), which is also used by brms. This version runs without transitions, unlike the code of Lee, Kim, and Gupta (2020). 

At this point I have not changed the code of either Cunningham or Lee, Kim, and Gupta, though I have modified the code of Piironen and Vehtari. You can find a modified version of their code in bscm_horseshoe_modified.stan. 

A third approach is included in `sc_spike_slab.Rmd`, which is a spike and slab formulation for synthetic controls.

This is a work in progress, and contributions are welcome. In particular, I am looking for convenient ways of using auxillary variables to help select the weights.
