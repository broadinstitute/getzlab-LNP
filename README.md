# Log-normal-Poisson regression model
Code to run the log-normal-Poisson regression model from Hess _et al._ 2019, "Passenger Hotspot Mutations in Cancer" (https://doi.org/10.1016/j.ccell.2019.08.002)

### What's in this repo?

This repository contains MATLAB functions for running a Bayesian log-normal-Poisson (LNP) regression. We include code to both sample from the posterior distribution of the model parameters via MCMC, and also compute statistical significance of counts, given a set of posterior samples (i.e., compute posterior predictive _p_-values).

The aforementioned functionality is suitable for modeling generic count data, and is not specific to modeling somatic mutation counts. For that purpose, we also include a function to process the specific somatic mutation calls analyzed in Hess _et al._

### Included tools

* To run the regression on raw count/covariate data, use `src/regr/pois_LN_reg.m` (see `demo.m`, or the [inline example](#Regression-demo) in this readme)
  * This tool generates samples from the posterior distribution of the LNP parameters (μ, τ, β); it does not identify significant outlying counts.
  * To compute significance (i.e., LNP posterior predictive _p_-values), use `src/regr/regr_post_pred.m`
* To run the model on the somatic mutation calls analyzed in the manuscript, use the wrapper script `src/pois_LN_reg_wrapper.m`. Currently, it is only intended for processing calls in a format specific to this manuscript, output by the [analysis notebooks](https://github.com/broadinstitute/getzlab-PHS) used to generate the publication. We are planning to release a user-friendly wrapper for running on generic `.maf` files soon.

### Regression demo

Here is a simple demo of running the regression model on simulated counts, sans covariates.

First, we generate 50,000 random samples from a log-normal-Poisson distribution with μ = -3, σ = 0.9:

``` MATLAB
x = poissrnd(exp(-3 + 0.9*randn(50000, 1)));
```
We thus expect samples drawn from the posterior distribution p(μ,σ|x) will be centered around (-3, 0.9).

Next, we set some basic parameters for running the MCMC:

``` MATLAB
P = [];
P.niter = 2000; % total iterations
P.burnin = 50; % burn-in
P.m0 = log(mean(x)); % initial guess for mu

% normal-gamma hyperparameters
P.mumu = -3;
P.taua = 10;
P.taub = 0.2;
```

Actually run the MCMC:

``` MATLAB
[~, ~, mu, tau] = pois_LN_reg(x, zeros(size(x)), P);
```

2,000 samples from μ and σ will be placed into `mu` and `tau`, respectively. Note that τ = 1/σ^2.

Finally, we plot samples from the posterior:

``` MATLAB
figure(1); clf
hold on
scatter(mu(P.burnin:end), 1./sqrt(tau(1, P.burnin:end)), 'marker', '.', ...
  'markeredgealpha', 0.6, 'markeredgecolor', 'k')
scatter(-3, 0.9, 70, 'marker', 'x', 'markeredgecolor', 'm', 'linewidth', 2)

title('MCMC draws from LNP posterior p(\mu, \sigma|x)')

xlabel('\mu')
ylabel('\sigma')

ax = gca;
ax.Box = 'on';

ax.XLim = [-3.15 -2.85];
ax.XTick = -3.15:0.05:-2.70;

ax.YLim = [0.75 1.05];
ax.YTick = 0.75:0.05:1.05;

grid on
```

<img src="https://github.com/broadinstitute/getzlab-LNP/blob/master/demo.png" alt="MCMC samples from posterior" width="600">

confirming that samples from the posterior distribution look as we would expect.
