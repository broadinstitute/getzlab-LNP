%
% sample from log-normal-Poisson with mu = -3, sigma = 0.9

x = poissrnd(exp(-3 + 0.9*randn(50000, 1)));

%
% set MCMC parameters

P = [];
P.niter = 5000; % total iterations
P.burnin = 50; % burn-in 
P.m0 = log(mean(x)); % initial guess for mu

% normal-gamma hyperparameters
P.mumu = -3;
P.taua = 10;
P.taub = 0.2;

%
% run MCMC

[~, ~, mu, tau] = pois_LN_reg(x, zeros(size(x)), P);

%
% plot posterior samples and ground truth mu, sigma

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
