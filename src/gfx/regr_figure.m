function regr_figure(M, Mu, C, post_pred, pvalues, mu, tau, Beta, F, P)

figure(F); clf
set(F, 'Position', [1 1 1200 900])

%
%define subplots

%qq
qq_sp = subplot('Position', [0.15 0.6 0.38 0.35]);

%histogram
hist_sp = subplot('Position', [0.6 0.6 0.38 0.35]);

%mu/sigma
musig_sp = subplot('Position', [0.15 0.15 0.38 0.35]);

mumrg_sp = subplot('Position', [0.15 0.05 0.38 0.07]);
sigmrg_sp = subplot('Position', [0.05 0.15 0.07 0.35]);

%regressors
regr_sp = subplot('Position', [0.6 0.15 0.38 0.35]);

%
%Q-Q
set(F, 'currentaxes', qq_sp);
ax = gca;

if P.compute_posterior_predictive,
  if isfield(P, 'pp_recurrence_floor') && P.pp_recurrence_floor > 0,
    error('Not currently supported.')
    %pvs = [ones(nnz(M < P.pp_recurrence_floor), 3); pvalues];
  else
    pvs = pvalues;
  end
  [Mu, si] = sort_struct(Mu, 'p');

  %for pretty QQ plots, p_mid lies randomly between p_min/max
  p_mid = pvalues(:, :, 2) + bsxfun(@times, (pvalues(:, :, 1) - pvalues(:, :, 2)), rand(size(pvalues, 1), 1));

  if isfield(Mu, 'gene'),
    qq_jh(p_mid, Mu.gene(1:10));
  else
    qq_jh(p_mid);
  end 

  xlabel('Uniform theoretical quantiles')
  ylabel('Observed p-value quantiles')
else
  ax.XTick = [];
  ax.YTick = [];
  ax.XAxis.Visible = 'off';
  ax.YAxis.Visible = 'off';
  xlim([0 1])
  ylim([0 1])

  text(0.2, 0.5, 'No posterior predictive computed')
end
ax.Box = 'on';

title('Q-Q Plot')

%
%histogram
set(F, 'currentaxes', hist_sp);
ax = gca;

if P.compute_posterior_predictive,
  hold on

  if isfield(P, 'pp_recurrence_floor') && P.pp_recurrence_floor > 0,
    error('Not currently supported.')
    idx = M >= P.pp_recurrence_floor;
  end

  %observed fraction plot
  h = histc(M, 0:max(M)); h = h./sum(h);
  obs_plot = stem(0:max(M), log10(h), 'Color', 'r');

  %posterior predictives
  regr_plot = scatter(M, log10(post_pred(:, 3)), 'Marker', 's', 'MarkerFaceAlpha', 0.1, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');

  %strawman Poisson regression
  C(isinf(C)) = 0;
  [BetaP, dBeta, stats] = glmfit(C, M, 'poisson');

  pr = poisspdf(M, exp(BetaP(1) + C*BetaP(2:end)));

  prs = accumarray(M + 1, pr, [], @(x) {x});
  pr_mean = cellfun(@sum, cellfun(@mean, prs, 'unif', 0)); pr_mean(pr_mean == 0) = 1e-300;
  pr_max = cellfun(@sum, cellfun(@max, prs, 'unif', 0)); pr_max(pr_max == 0) = 1e-300;
  pr_min = cellfun(@sum, cellfun(@min, prs, 'unif', 0)); pr_min(pr_min == 0) = 1e-300;

  nidx = ~isnan(pr_max);
  pois_plot = scatter(find(nidx) - 1, log10(pr_mean(nidx)), '+', 'MarkerEdgeColor', 'm', 'LineWidth', 1);

  xlim([0 20])
  ylim([min([-20/log(10) log10(min(pr_mean(1:min(20, length(pr_mean))))) - 0.2]) 0])

  ax.XTick = 0:2:20;
  ax.YTickLabel = strsplit(sprintf('10^{%d} ', ax.YTick), ' ');

  xlabel('Sitewise recurrence')
  ylabel('Fraction of sites')

  grid on

  legend([obs_plot regr_plot pois_plot], 'Location', 'SouthWest', 'Observed', 'LN-Pois. regr.', 'Pois. regr.')
else
  ax.XTick = [];
  ax.YTick = [];
  ax.XAxis.Visible = 'off';
  ax.YAxis.Visible = 'off';
  xlim([0 1])
  ylim([0 1])

  text(0.2, 0.5, 'No posterior predictive computed')
end
ax.Box = 'on';

title('Recurrence Histogram')

%
%mu/sigma joint
set(F, 'currentaxes', musig_sp);

hold on
p = plot_vary_color(mu, 1./sqrt(tau), 1:length(mu), 'EdgeAlpha', '0.1');
colormap(jet);
sig = 1./sqrt(tau);
s = scatter(flipud(mu), flipud(sig), 30, 'filled');
alpha(s, 0.01)
scatter(mu(P.burnin), 1./sqrt(tau(P.burnin)), 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'SizeData', 20)
if isfield(P, 'post_abs_mu') && isfield(P, 'post_abs_tau'),
  scatter(P.post_abs_mu, 1./sqrt(P.post_abs_tau), 'o', 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'k', 'SizeData', 20)
end

xl = [min(mu(P.burnin:end)) max(mu(P.burnin:end))];
yl = [min(sig(P.burnin:end)) max(sig(P.burnin:end))];
xlim(xl + diff(xl)*1.3*[-1 1])
ylim(yl + diff(yl)*1.3*[-1 1])

set(gca, 'Box', 'on');

title('\mu,\sigma Joint Posterior')

%
%marginal (mu)
set(F, 'currentaxes', mumrg_sp)
xl = xlim(musig_sp);
mur = xl(1):(diff(xl)/100):xl(2);
muh = histc(mu, mur); muh = muh./sum(muh);
stairs(mur, muh)
xlim(xl)
linkaxes([musig_sp mumrg_sp], 'x')
set(mumrg_sp, 'XTick', [])

set(gca, 'Box', 'on');

xlabel('$\mu$', 'Interpreter', 'latex')

%
%marginals (sigma)
set(F, 'currentaxes', sigmrg_sp)

yl = ylim(musig_sp);
sigr = yl(1):(diff(yl)/100):yl(2);
sigh = histc(sig, sigr); sigh = sigh./sum(sigh);
st = stairs(sigh, sigr);
ylim(yl)
linkaxes([musig_sp sigmrg_sp], 'y')
set(sigmrg_sp, 'YTick', [])

set(gca, 'Box', 'on');

ylabel('$\sigma$', 'Interpreter', 'latex')

%
%regression coefficients
set(F, 'currentaxes', regr_sp);
ax = gca;

boxplot(Beta(:, P.burnin:end)')

ax.XTick = 1:size(Beta, 1);
if isfield(P, 'covar_names'),
  ax.XTickLabel = P.covar_names;
  ax.XTickLabelRotation = 90;
  ax.XAxis.FontSize = 7; 
end
grid on

set(gca, 'Box', 'on');

xlabel('Covariate')
ylabel('Linspace value')

title('Regression Coefficients')
