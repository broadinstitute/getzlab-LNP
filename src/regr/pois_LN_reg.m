function [Beta, epsi, mu, tau, logmarg_lik, full_lik] = pois_LN_reg(varargin)
usage = sprintf(['pois_LN_reg(<counts>, <covariate matrix>, <P struct>) or\n' ...
         'pois_LN_reg(<counts>, <exposures>, <covariate matrix>, <P struct>)']);

if nargin == 3,
  x = varargin{1};
  c = varargin{2};
  P = varargin{3};
elseif nargin == 4,
  x = varargin{1};
  y = varargin{2};
  c = varargin{3};
  P = varargin{4};
else
  error(usage);
end

X = [];

fprintf('\n');

X.x = as_column(x);
X.len = length(x);

if any(size(c) == X.len),
  if find(size(c) == X.len) ~= 1,
    X.c = c';
  end

  X.c = c;

  %TODO: more robust handling of missing covariate information
  X.c(isnan(X.c)) = 0;
  X.c(isinf(X.c)) = 0;

  %add exposures, or replace with default of zero (log(1)) if missing
  if exist('y', 'var'),
    X.y = as_column(y);
  else
    X.y = 0;
  end

  %handle random effect categories, if present
  if isfield(P, 'random_effect_cat_covar_index'),
    X.reidx = X.c(:, P.random_effect_cat_covar_index);
    if size(X.c, 2) > 1,
      X.c(:, P.random_effect_cat_covar_index) = [];
    else
      X.c(:, P.random_effect_cat_covar_index) = 0;
    end
  else
    X.reidx = ones(size(X.x));
  end
  X.ncat = max(X.reidx);

  %memoize X sums
  X.sums = accumarray(X.reidx, X.x);

  %binary expansion of X.reidx for quick slicing
  X.reidx_bin = sparse(1:length(X.reidx), X.reidx, 1);

  %necessary for computing gradient of p(beta|-)
  X.cx = X.x'*X.c;

  %necessary for computing Hessian of p(beta|-)
  X.cc = NaN(size(X.c, 2), size(X.c, 2), X.len);

  for i = 1:X.len,
    X.cc(:, :, i) = X.c(i, :)'*X.c(i, :);
  end
  X.cc = reshape(X.cc, size(X.c, 2)^2, X.len);
else
  error('Covariate matrix must be of data dimension!');
end
clear x y c

%
%lookup Hermite polynomial roots/quadrature weights
[X.h_r X.q_w] = hermite();

%
% Parameters
%

if ~exist('P', 'var'), P = []; end
P = impose_default_value(P, 'fixed_beta', []); 

%whether to use hierarchical model.  otherwise, categories are treated independently.
P = impose_default_value(P, 'ignore_categories', false);

% total number of MCMC iterations/burnin
P = impose_default_value(P, 'niter', 300);
P = impose_default_value(P, 'burnin', 100);

%initial parameter values
P = impose_default_value(P, 't0', ones(X.ncat + 1, 1));
P = impose_default_value(P, 'm0', zeros(X.ncat, 1));
P = impose_default_value(P, 'beta0', zeros(size(X.c, 2)));

% maximum number of Newton-Raphson iterations to find M-H proposal distribution
P = impose_default_value(P, 'tMH_NR_iter', 100);

% M-H inverse Hessian scaling factor (to tweak M-H proposal dist.)
P = impose_default_value(P, 'beta_nu', 20);
P = impose_default_value(P, 'mutau_ihsf', 1);
P = impose_default_value(P, 'tau_ihsf', 1);
P = impose_default_value(P, 'epsi_ihsf', 1);
P = impose_default_value(P, 'epsi_nu', 5);

%Hyperparameters
%mu, tau ~ normal-gamma(mumu, 1/sqrt(mutau), taua, taub)
P = impose_default_value(P, 'mumu', zeros(X.ncat, 1));
P = impose_default_value(P, 'mutau', ones(X.ncat, 1));

P = impose_default_value(P, 'taua', ones(X.ncat + 1, 1));
P = impose_default_value(P, 'taub', 10*ones(X.ncat + 1, 1));

%beta ~ normal_m(mubeta, covmu) 
P = impose_default_value(P, 'covmu', 2*eye(size(X.c, 2))); X.invcovmu = inv(P.covmu);
P = impose_default_value(P, 'mubeta', zeros(size(X.c, 2), 1));

X.P = P;

%
% Initialization
%

%
%tau/mu
tau = NaN(X.ncat + 1, P.niter); tau(:, 1) = P.t0;
mu = NaN(X.ncat, P.niter); mu(:, 1) = P.m0;

%
%epsilon
fprintf('Initializing epsilon chain: tau0 = %0.2f\n', P.t0);

if ~P.ignore_categories && X.ncat > 1,
  epsi = epsisamp(zeros(X.len, 1), P.beta0, tau(end, 1), mu(X.reidx, 1), zeros(X.len, 1), tau(X.reidx, 1), X);
else
  epsi = zeros(size(X.x));
  epsi_ar = NaN;
end
epsi_h = epsisamp(zeros(X.len, 1), P.beta0, tau(X.reidx, 1), mu(X.reidx, 1), epsi, tau(end, 1), X);

%smooth out initial epsilon samples
for i = 1:5,
  if ~P.ignore_categories && X.ncat > 1,
    [epsi epsi_ar] = epsisamp(epsi, P.beta0, tau(end, 1), mu(X.reidx, 1), epsi_h, tau(X.reidx, 1), X);
  end
  [epsi_h epsi_h_ar] = epsisamp(epsi_h, P.beta0, tau(X.reidx, 1), mu(X.reidx, 1), epsi, tau(end, 1), X);
end

%
%beta
Beta = NaN(size(X.c, 2), P.niter);
if isempty(P.fixed_beta),
  if isempty(P.beta0),
    fprintf('Initializing beta chain ...\n');

    Beta(:, 1) = betasamp(zeros(size(X.c, 2), 1), epsi, P.m0, P.t0, X);
  else
    fprintf('Using user-supplied beta0 ...\n');
    if length(P.beta0) ~= size(X.c, 2), error('Length of P.beta0 must match number of covariates!'); end
    Beta(:, 1) = P.beta0;
  end
else
  if size(P.fixed_beta, 1) ~= size(X.c, 2),
    error('Fixed beta size mismatch!');
  end
  Beta = repmat(P.fixed_beta, 1, P.niter);
end

%whether Metropolis sample was rejected for beta and mu_j,tau_j, tau_0
mrej_b = NaN(P.niter, 1);
mrej_mt = NaN(P.niter, X.ncat);
mrej_t0 = NaN(P.niter, 1);

%full log likelihoods (for temperature diagnostics)
full_lik = NaN(P.niter, 1);

fprintf('\n');

%
% Primary MCMC iterations
%

burnin_flag = 0;
post_abs_ready_flag = 0;
count_abs = 0;
count_beta = 1;
count_PP = 0;
count_PPavg = 0;
last_full_lik_idx = 1;

fprintf('Running full posterior MCMC for %d iterations (burnin %d) ...\n', P.niter, P.burnin)

last_full_lik_idx = 1;

for i = 2:P.niter, 
  %
  % Metropolis/Gibbs draws

  %%draw epsilons via M-H
  if ~P.ignore_categories && X.ncat > 1,
    [epsi epsi_ar] = epsisamp(epsi, Beta(:, i - 1), tau(end, i - 1), mu(X.reidx, i - 1), epsi_h, tau(X.reidx, i - 1), X);
  end
  [epsi_h epsi_h_ar] = epsisamp(epsi_h, Beta(:, i - 1), tau(X.reidx, i - 1), mu(X.reidx, i - 1), epsi, tau(end, i - 1), X);

  %%draw beta via M-H
  if isempty(P.fixed_beta),
    [Beta(:, i), mrej_b(i)] = betasamp(Beta(:, i - 1), epsi, mu(:, i - 1), tau(:, i - 1), epsi_h, X);
  end

  %%draw (mu_j, tau_j) via Gibbs sampling of full conditional
  mt_use_LS = false;
  for j = 1:X.ncat,
    [mu(j, i), tau(j, i), mrej_mt(i, j), use_LS] = mutausamp(mu(:, i - 1), tau(:, i - 1), Beta(:, i), epsi_h, tau(end, i - 1), epsi, j, X);
    mt_use_LS = mt_use_LS || use_LS;
  end

  %%draw tau_0 via Gibbs sampling of full conditional
  if ~P.ignore_categories && X.ncat > 1,
    [tau(end, i) mrej_t0(i)] = tausamp(tau(end, i - 1), mu(:, i), Beta(:, i), epsi, tau(:, i), epsi_h, X);
  else
    tau(end, i) = Inf;
  end

  if ~burnin_flag && i > P.burnin, burnin_flag = 1; end

  if X.ncat <= 4,
    m_string = sprintf('%0.2f ', mu(:, i)); m_string(end) = [];
    t_string = sprintf('%0.2f ', tau(:, i)); t_string(end) = [];
    mrej_t_string = sprintf('%0.2f ', 1 - mean(mrej_mt(max(i - 50, 2):i, :), 1)); mrej_t_string(end) = [];
  else
    m_string = sprintf('<%0.2f>', mean(mu(:, i)));
    t_string = sprintf('<%0.2f>', 1/mean(1./tau(:, i)));
    mrej_t_string = sprintf('<%0.2f>', 1 - mean(mean(mrej_mt(max(i - 50, 2):i, :), 1)));
  end

  if mt_use_LS, warnstr = '[!]'; else warnstr = ''; end

  fprintf('[%d/%d] m = %s, t = %s, e. A.R. = %0.2f, m/t A.R. = %s, b. A.R. = %0.2f, LL = %0.3f %s\n', ...
            i, P.niter, m_string, t_string, epsi_ar, mrej_t_string, 1 - mean(mrej_b(max(i - 50, 2):i)), ...
            full_lik(last_full_lik_idx), warnstr);
end
fprintf('\n');

logmarg_lik = NaN;

end

%
% SUBFUNCTIONS
%

%Generic functions
%%%%%%%%%%%

%log multivariate normal PDF
function p = logmvnpdf(x, mu, icmat) %assumes icmat is already inverted to save compute
  p = -(x - mu)'*icmat*(x - mu)/2;
end

function qrat = mvnqrat(th0, thP, muF, muR, hessF, hessR, hessFinv, hessRinv, ihsf)
  qrat = (-log(det(-hessRinv/ihsf)) + (th0 - muR)'*hessR*ihsf*(th0 - muR) - ...
	 (-log(det(-hessFinv/ihsf)) + (thP - muF)'*hessF*ihsf*(thP - muF)))/2;
end

function qrat = nqrat(th0, thP, muF, muR, sigF, sigR, ihsf)
  qrat = (-2*log(sigR/ihsf) + (th0 - muR)^2/(sigR/ihsf)^2 - ...
	 (-2*log(sigF/ihsf) + (thP - muF)^2/(sigF/ihsf)^2))/2;
end

%returns vector of marginals from multivariate normal proposal ratio, as if we were drawing from
%many univariate distributions
%NOTE: only works in the case of a diagonal covariance matrix
function qrat = margmvnqrat(th0, thP, muF, muR, hessF, hessR, hessFinv, hessRinv, ihsf)
  if size(hessF, 2) ~= 1,
    error('Only diagonal covariance matrices (supplied as column vectors) are supported.');
  end

  qrat = (-log(-hessRinv/ihsf) + (th0 - muR).^2.*hessR*ihsf - ...
         (-log(-hessFinv/ihsf) + (thP - muF).^2.*hessF*ihsf))/2;
end

function qrat = tqrat(th0, thP, muF, muR, sigF, sigR, nu)
  qrat = -log(sigR) - (nu + 1)/2*log(1 + (th0 - muR).^2./(nu*sigR.^2)) - ...
         (-log(sigF) - (nu + 1)/2*log(1 + (thP - muF).^2./(nu*sigF.^2)));
end

function qrat = mvtqrat(th0, thP, muF, muR, hessF, hessR, hessFinv, hessRinv, nu)
  m = length(muF);

  qrat = (-log(det(-hessRinv)) - (nu + m)*log(1 - (th0 - muR)'*hessR*(th0 - muR)/nu) - ...
	 (-log(det(-hessFinv)) - (nu + m)*log(1 - (thP - muF)'*hessF*(thP - muF)/nu)))/2;
end

%Tau_0 M-H functions
%%%%%%%%%%%

function [tau, rej] = tausamp(tau, mu, Beta, epsi, tau_2, epsi_2, X),
  %constant variables
  mu = mu(X.reidx);

  Y = [];
  Y.eps_x = epsi'*X.x;
  Y.eBC = exp(X.c*Beta + X.y); %XXX: does this only work if X.y is a scalar?
  Y.exp_eps_2_tau_2 = exp(epsi_2./sqrt(tau_2(X.reidx)));

  %Newton-Raphson iterations to find proposal density
  %note log transform of tau: tau -> exp(tau) in N-R iterations to allow it to range over all reals.
  [muF, HF, HFinv, is_max] = tauNR(log(tau), mu, Beta, epsi, epsi_2, X, Y);

  %now propose with univariate Gaussian centered at MAP (tau log transformed) with variance from Hessian
  tauP = normrnd(muF, sqrt(-HFinv*X.P.tau_ihsf));

  %if we'd reached a local maximum, then parameterization of forward and reverse jumps identical.
  %if not, then reverse jump will have its own parameterization.
  if ~is_max,
    [muR, HR, HRinv] = tauNR(muF, mu, Beta, epsi, epsi_2, X, Y);
  else
    muR = muF;
    HR = HF; HRinv = HFinv;
  end

  arat = prattau(tau, tauP, mu, Beta, epsi, epsi_2, tau_2, X, Y) + ...
	 nqrat(tau, tauP, muF, muR, sqrt(-HFinv), sqrt(-HRinv), X.P.tau_ihsf);

  if log(rand) < min(0, arat),
    tau = exp(tauP);
    rej = 0;
  else
    rej = 1;
  end
end

%log posterior ratio for tau0 (tau0 log transformed)
function pr = prattau(tau, tauP, mu, Beta, epsi, epsi_2, tau_2, X, Y)
  tauP_e = exp(tauP); tau_e = exp(tau);

  pr = Y.eps_x/sqrt(tauP_e) - (Y.eBC.*Y.exp_eps_2_tau_2)'*(exp(mu).*exp(epsi/sqrt(tauP_e))) + ...
        loggampdf(tauP_e, X.P.taua(end), X.P.taub(end)) + tauP - ...
       (Y.eps_x/sqrt(tau_e) - (Y.eBC.*Y.exp_eps_2_tau_2)'*(exp(mu).*exp(epsi/sqrt(tau_e))) + ...
        loggampdf(tau_e, X.P.taua(end), X.P.taub(end)) + tau);
end

%first derivative of log p(tau0|-)
function gr = gradtau(tau, Beta, epsi, X, Y)
  gr = (-Y.eps_x + Y.d1'*epsi)/(2*tau^(3/2)) + ...
       (X.P.taua(end) - 1)/tau - 1/X.P.taub(end);
end

%second derivative of log p(tau0|-)
function H = hesstau(tau, Beta, epsi, X, Y)
  H = 3/4*tau^(-5/2)*Y.eps_x - ((Y.d1'*epsi)*(3/4)*tau^(-5/2) + ...
                                (Y.d1'*(epsi.^2))*(tau^-3)/4) - ...
      (X.P.taua(end) - 1)/tau^2;
end

%Newton-Raphson iterations for tau
function [tau, H, Hinv, is_max] = tauNR(tau, mu, Beta, epsi, epsi_2, X, Y)
  is_max = 0;
  for i = 1:X.P.tMH_NR_iter,
    tau_e = exp(tau);

    Y.exp_eps_tau = exp(epsi/sqrt(tau_e));
    Y.d1 = exp(mu).*Y.eBC.*Y.exp_eps_tau.*Y.exp_eps_2_tau_2;

    grad = gradtau(tau_e, Beta, epsi, X, Y);
    H = hesstau(tau_e, Beta, epsi, X, Y);

    %change-of-variable chain rule factors
    H = grad*tau_e + tau_e^2*H;
    grad = grad*tau_e;

    %change-of-variable Jacobian factors
    grad = grad + 1;

    %check if we've reached a local maximum (Hessian should be negative definite)
    if norm(grad) <= 1e-6 && H < 0, is_max = 1; break; end

    %N-R update
    tau = tau - grad/H;
  end

  Hinv = 1/H;
end

%Mu_j/Tau_j M-H functions
%%%%%%%%%%%

%a M-H sample from p(mu_j, tau_j|-)
function [mu, tau, rej, use_LS] = mutausamp(mu, tau, Beta, epsi, tau_2, epsi_2, j, X),
  %select jth mu/tau from array
  mu = mu(j);
  tau = tau(j);
  %TODO: select proper hyperparameters in the same fashion

  tidx = logical(X.reidx_bin(:, j));

  epsi = epsi(tidx);
  epsi_2 = epsi_2(tidx);

  %constant variables
  Y = [];
  Y.j = j;
  Y.eps_x = epsi'*X.x(tidx);
  if max(size(X.y)) > 1, y = X.y(tidx); else, y = X.y; end
  Y.eBC = exp(X.c(tidx, :)*Beta + y); %XXX: does this only work if X.y is a scalar?
  
  Y.exp_eps_2_tau_2 = exp(epsi_2/sqrt(tau_2));

  mutau = [mu; log(tau)];

  %Newton-Raphson iterations to find proposal density
  %note log transform of tau: tau -> exp(tau) in N-R iterations to allow it to range over all reals.
  [muF, HF, HFinv, is_max, use_LS] = mutauNR(mu, log(tau), Beta, epsi, epsi_2, tau_2, X, Y);

  %now propose with multivariate Gaussian centered at MAP (tau log transformed) with covariance matrix from Hessian
  mutauP = mvnrnd(muF, -HFinv*X.P.mutau_ihsf)';

  %if we'd reached a local maximum, then parameterization of forward and reverse jumps identical.
  %if not, then reverse jump will have its own parameterization.
  if ~is_max,
    [muR, HR, HRinv] = mutauNR(muF(1), muF(2), Beta, epsi, epsi_2, tau_2, X, Y);
  else
    muR = muF;
    HR = HF; HRinv = HFinv;
  end

  arat = pratmutau(mutau(1), mutau(2), mutauP(1), mutauP(2), Beta, epsi, epsi_2, tau_2, X, Y) + ...
	 mvnqrat(mutau, mutauP, muF, muR, HF, HR, HFinv, HRinv, X.P.mutau_ihsf);

  if log(rand) < min(0, arat),
    mu = mutauP(1);
    tau = exp(mutauP(2));
    rej = 0;
  else
    rej = 1;
  end
end

%full conditional for mu,tau
function fc = fcmutau(mu, tau, Beta, epsi, epsi_2, tau_2, X, Y)
  tau_e = exp(tau);

  fc = mu*X.sums(Y.j) + Y.eps_x/sqrt(tau_e) - exp(mu)*((Y.eBC.*Y.exp_eps_2_tau_2)'*exp(epsi/sqrt(tau_e))) + ...
        normgampdf(mu, tau_e, X.P.taua(Y.j), X.P.taub(Y.j), X.P.mumu(Y.j), X.P.mutau(Y.j), true) + tau;
end

%log posterior ratio for mu,tau (tau log transformed)
function pr = pratmutau(mu, tau, muP, tauP, Beta, epsi, epsi_2, tau_2, X, Y)
  tauP_e = exp(tauP); tau_e = exp(tau);

  pr = fcmutau(muP, tauP, Beta, epsi, epsi_2, tau_2, X, Y) - ...
       fcmutau(mu, tau, Beta, epsi, epsi_2, tau_2, X, Y);
end

%gradient of log p(mu, tau|-)
function gr = gradmutau(mu, tau, Beta, epsi, X, Y)
  gr = [X.sums(Y.j) - exp(mu)*sum(Y.d1) - ...
        X.P.mutau(Y.j)*tau*(mu - X.P.mumu(Y.j)), ...
        (-Y.eps_x + exp(mu)*(Y.d1'*epsi))/(2*tau^(3/2)) + ...
         (2*X.P.taua(Y.j) - 1)/(2*tau) - 1/X.P.taub(Y.j) - X.P.mutau(Y.j)/2*(mu - X.P.mumu(Y.j))^2];
end

%Hessian of log p(mu, tau|-)
function H = hessmutau(mu, tau, Beta, epsi, X, Y)
  H = NaN(2);

  H(1, 1) = -exp(mu)*sum(Y.d1) - ...
            X.P.mutau(Y.j)*tau;
  H(2, 2) = 3/4*tau^(-5/2)*Y.eps_x - exp(mu)*((Y.d1'*epsi)*(3/4)*tau^(-5/2) + ...
                                              (Y.d1'*(epsi.^2))*(tau^-3)/4) - ...
            (2*X.P.taua(Y.j) - 1)/(2*tau^2);
  H(1, 2) = exp(mu)*(Y.d1'*epsi)/(2*tau^(3/2)) - ...
            X.P.mutau(Y.j)*(mu - X.P.mumu(Y.j));
  H(2, 1) = H(1, 2);
end

%Newton-Raphson iterations for mu,tau
function [mutau, H, Hinv, is_max, use_LS] = mutauNR(mu, tau, Beta, epsi, epsi_2, tau_2, X, Y)
  is_max = 0; %if we converged to a maximum
  is_NR_bad = 0; %if regular N-R iterations are insufficent,
                 %so we need to employ line search
  use_LS = 0; %whether we used linesearch (for display purposes only)

  mu_prev = mu;
  tau_prev = tau;

  i = 1;
  while true,
    tau_e = exp(tau);

    Y.exp_eps_tau = exp(epsi/sqrt(tau_e));
    Y.d1 = Y.eBC.*Y.exp_eps_tau.*Y.exp_eps_2_tau_2;

    grad = gradmutau(mu, tau_e, Beta, epsi, X, Y)';
    H = hessmutau(mu, tau_e, Beta, epsi, X, Y);

    %change-of-variable chain rule factors
    H(2, 2) = grad(2)*tau_e + tau_e^2*H(2, 2);
    H(1, 2) = tau_e*H(1, 2); H(2, 1) = H(1, 2);
    grad(2) = grad(2)*tau_e;

    %change-of-variable Jacobian factors
    grad(2) = grad(2) + 1;

    %if Hessian is problematic, rewind an iteration and use line search by default
    if rcond(H) < eps || any(isnan(H(:))) || any(isinf(H(:))),
      %if we're in a problematic region at the first iteration, break and hope for the best
      if i == 1,
	H = eye(2); Hinv = -eye(2);
	disp('WARNING: Full conditional shifted significantly since last iteration!');

	break
      end

      is_NR_bad = 1;
      i = i - 1;

      mu = mu_prev;
      tau = tau_prev;

      continue
    else
      is_NR_bad = 0;
    end

    %check if Hessian is negative definite
    is_ndef = ~is_NR_bad && all(eig(H) < 0);

    %1. if Newton-Raphson will work fine, use it
    Hinv = inv(H);
    step = -Hinv*grad;

    %check if we've reached a local maximum
    if norm(grad) <= 1e-6 && is_ndef, is_max = 1; break; end

    %2. otherwise, employ line search
    fc_step = fcmutau(mu + step(1), tau + step(2), Beta, epsi, epsi_2, tau_2, X, Y);
    if is_NR_bad || isnan(fc_step) || isinf(fc_step) || ...
       fc_step - fcmutau(mu, tau, Beta, epsi, epsi_2, tau_2, X, Y) < -1e-3,
      %indicate that we used linesearch for these iterations
      use_LS = 1;

      %2.1. ensure N-R direction is even ascending.  if not, use direction of gradient
      if ~is_ndef, step = grad; end

      %2.2. regardless of the method, perform line search along direction of step
      s0 = norm(step);
      d_hat = step/s0;

      %bound line search from below at current value
      fc = fcmutau(mu, tau, Beta, epsi, epsi_2, tau_2, X, Y)*[1 1];

      %2.3. do line search
      for l = 0:50,
	s = s0*0.5^l;

	f = fcmutau(mu + s*d_hat(1), tau + s*d_hat(2), Beta, epsi, epsi_2, tau_2, X, Y);
	if fc(mod(l - 1, 2) + 1) > fc(mod(l - 2, 2) + 1) && fc(mod(l - 1, 2) + 1) > f,
	  s = s0*0.5^(l - 1);
	  break;
	else
	  fc(mod(l, 2) + 1) = f;
	end
      end

      step = d_hat*s;
    end

    %update mu,tau
    mu_prev = mu; tau_prev = tau;
    mu = mu + step(1); tau = tau + step(2);

    i = i + 1;
    if i > X.P.tMH_NR_iter,
      if ~is_ndef,
	H = eye(2); Hinv = -eye(2);
	disp('WARNING: Newton-Raphson terminated at non-concave point!');
      end

      break;
    end
  end
  mutau = [mu; tau];
end

%Epsilon M-H functions
%%%%%%%%%%%

%a M-H sample from p(epsi|-)
function [epsi, mrej] = epsisamp(epsi, Beta, tau, mu, epsi_2, tau_2, X),
  %assumes no covariance between epsilons; does not sample as a single block

  %Newton-Raphson iterations to find proposal density
  [muF, HF, HFinv] = epsiNR(epsi, mu, tau, Beta, epsi_2, tau_2, X);

  %now propose with multivariate t centered at epsiMLE with covariance matrix from Hessian
  %note that since Hessian is diagonal, we can just simulate from n univariate t's.
  epsiP = muF + sqrt(-HFinv).*trnd(X.P.epsi_nu, X.len, 1);
  %epsiP = normrnd(muF, -HFinv);

  arat = pratepsi(epsi, epsiP, Beta, tau, mu, epsi_2, tau_2, X) + tqrat(epsi, epsiP, muF, muF, sqrt(-HFinv), sqrt(-HFinv), X.P.epsi_nu);

  ridx = log(rand(X.len, 1)) >= min(0, arat);
  epsi(~ridx) = epsiP(~ridx);
  mrej = mean(~ridx);
end

%log posterior ratios for epsilon vector
function pr = pratepsi(epsi, epsiP, Beta, tau, mu, epsi_2, tau_2, X)
  pr = epsiP.*X.x./sqrt(tau) - exp(mu + X.c*Beta + epsiP./sqrt(tau) + X.y + epsi_2./sqrt(tau_2)) - epsiP.^2/2 - ...
       (epsi.*X.x./sqrt(tau) - exp(mu + X.c*Beta + epsi./sqrt(tau) + X.y + epsi_2./sqrt(tau_2)) - epsi.^2/2);
end

%gradient of log p(epsilon|-)
function gr = gradepsi(epsi, Beta, tau, mu, epsi_2, tau_2, X)
  gr = X.x./sqrt(tau) - exp(mu + X.c*Beta + epsi./sqrt(tau) + X.y + epsi_2./sqrt(tau_2))./sqrt(tau) - epsi;
end

%Hessian of log p(epsilon|-)
function [H, Hinv] = hessepsi(epsi, Beta, tau, mu, epsi_2, tau_2, X)
  %because Hessian is diagonal, we store as a column vector
  H = -exp(mu + X.c*Beta + epsi./sqrt(tau) + X.y + epsi_2./sqrt(tau_2))./tau - 1;
  Hinv = 1./H;
end

%Newton-Raphson iteration for forward and reverse jumps
function [epsi, H, Hinv] = epsiNR(epsi, mu, tau, Beta, epsi_2, tau_2, X),
  for i = 1:100,
    [H, Hinv] = hessepsi(epsi, Beta, tau, mu, epsi_2, tau_2, X);

    %N-R update
    grad = gradepsi(epsi, Beta, tau, mu, epsi_2, tau_2, X);
    epsi = epsi - Hinv.*grad;

    %we've reached a local maximum
    if norm(grad) < 1e-6, break; end
  end
end

%Beta M-H functions
%%%%%%%%%%%

%a M-H sample from p(beta|-)
function [Beta, rej] = betasamp(Beta, epsi, mu, tau, epsi_2, X)
  mu = mu(X.reidx);
  tau_2 = tau(X.reidx);
  tau = tau(end);

  %Newton-Raphson iterations to find proposal density
  [muF, HF, HFinv, is_max] = betaNR(Beta, epsi, mu, tau, epsi_2, tau_2, X);

  %propose from multivariate t centered at muF with covariance matrix from Hessian
  BetaP = mvtrnd_jh(muF, -HFinv, X.P.beta_nu, 1, 0);

  %if we'd reached a local maximum, then parameterization of forward and reverse jumps identical.
  %if not, then reverse jump will have its own parameterization.
  if ~is_max,
    [muR, HR, HRinv] = betaNR(BetaP, epsi, mu, tau, epsi_2, tau_2, X);
  else
    muR = muF;
    HR = HF; HRinv = HFinv;
  end

  arat = pratbeta(Beta, BetaP, epsi, mu, tau, epsi_2, tau_2, X) + mvtqrat(Beta, BetaP, muF, muR, HF, HR, HFinv, HRinv, X.P.beta_nu);

  if log(rand) < min(0, arat),
    Beta = BetaP;
    rej = 0;
  else
    rej = 1;
  end
end

%log posterior ratio for beta vector
function pr = pratbeta(Beta, BetaP, epsi, mu, tau, epsi_2, tau_2, X)
  pr = X.cx*BetaP - sum(exp(mu + X.c*BetaP + epsi/sqrt(tau) + X.y + epsi_2./sqrt(tau_2))) + logmvnpdf(BetaP, X.P.mubeta, X.invcovmu) - ...
      (X.cx*Beta - sum(exp(mu + X.c*Beta + epsi/sqrt(tau) + X.y + epsi_2./sqrt(tau_2))) + logmvnpdf(Beta, X.P.mubeta, X.invcovmu));
end

%gradient of log p(beta|-)
function gr = gradbeta(Beta, epsi, mu, tau, epsi_2, tau_2, X)
  gr = X.cx - exp(mu + X.c*Beta + epsi/sqrt(tau) + X.y + epsi_2./sqrt(tau_2))'*X.c - (X.invcovmu*(Beta - X.P.mubeta))';
end

%Hessian of log p(beta|-)
function [H, Hinv] = hessbeta(Beta, epsi, mu, tau, epsi_2, tau_2, X)
  H = -reshape(X.cc*exp(mu + X.c*Beta + epsi/sqrt(tau) + X.y + epsi_2./sqrt(tau_2)), size(X.c, 2)*[1 1]) - X.invcovmu;

  Hinv = inv(H);
end

%Newton-Raphson iterations for forward and reverse jumps
function [Beta, H, Hinv, is_max] = betaNR(Beta, epsi, mu, tau, epsi_2, tau_2, X),
  is_max = 0;
  for i = 1:X.P.tMH_NR_iter,
    grad = gradbeta(Beta, epsi, mu, tau, epsi_2, tau_2, X)';
    [H, Hinv] = hessbeta(Beta, epsi, mu, tau, epsi_2, tau_2, X); 

    %check if we've reached a local maximum (Hessian should be negative definite)
    if norm(grad) <= 1e-6 && all(eig(H) < 0), is_max = 1; break; end

    %N-R update
    Beta = Beta - Hinv*grad;
  end
end

%Quadrature functions
%%%%%%%%%%%

%return log of Poisson/normal integral for recurrences/covariates
%note: does not include normalizing constants
%      assumes exposure is in eBC term (i.e. eBC = X.c*Beta + X.y)
function g = pois_LN_int(x, eBC, mu, tau, X)
  I = (sqrt(2./tau).*x + mu.*sqrt(2.*tau))*X.h_r' - exp(bsxfun(@plus, eBC, sqrt(2./tau)*X.h_r'));
  m = max(I, [], 2);
  I = I - repmat(m, 1, 100);

  g = m + log(exp(I)*X.q_w);
end

%Likelihood functions
%%%%%%%%%%%

%Full log likelihood
function L = full_loglik(mu, tau, Beta, epsi, X)
  q = X.c*Beta + epsi + X.y;
  L = q'*X.x - sum(exp(q)) + (1/2)*log(tau) - tau/2*sum((epsi - mu).^2);
end

%Marginalized (WRT epsilon) log posterior numerator
function L = marg_postnum(mu, tau, Beta, X)
  eBC = X.c*Beta + X.y;

  L = -X.len/2*log(pi) + ...
      eBC'*X.x - X.len*tau*mu^2/2 - sum(gammaln(X.x + 1)) + ...
      sum(pois_LN_int(X.x, eBC, mu, tau, X)) + ...
      normgampdf(mu, tau, X.P.taua, X.P.taub, X.P.mumu, X.P.mutau, true) + ...
      log(mvnpdf(Beta, X.P.mubeta, X.P.covmu));
end
