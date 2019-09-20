function qq_jh(varargin)

usage = 'qq_jh(<p values>, [labels], [P struct])';

pvalues = varargin{1};
labels = [];

if nargin == 2,
  if isstruct(varargin{2}),
    P = varargin{2};
  elseif iscell(varargin{2}),
    labels = varargin{2};
  else
    error(usage);
  end
elseif nargin == 3,
  labels = varargin{2};
  P = varargin{3};
end

if ~exist('P', 'var'), P = []; end

P = impose_default_value(P, 'color', 'b');
P = impose_default_value(P, 'calc_fdr_values', true);
P = impose_default_value(P, 'plot_CI', true);

if size(pvalues, 2) == 1,
  np = length(pvalues);
  [logpvalues, si] = sort(-log10(pvalues));
else
  np = size(pvalues, 1);
  [logpvalues, si] = sortrows(-log10(pvalues), 3);
end
x = as_column(-log10((np:-1:1)/(np + 1)));

hold on

%CI
if P.plot_CI, 
  pts = round([np:(-np/1000):101  101:-1:1]);
  beta95ci = NaN(length(pts), 2);
  for y = [pts; 1:length(pts)],
    i = y(1); j = y(2);
    beta95ci(j, :) = -log10(icdf('beta', [0.025 0.975], i, np - i + 1));
  end

  fill([x(np - pts + 1); flipud(x(np - pts + 1))], [beta95ci(:, 2); flipud(beta95ci(:, 1))], [0.5 0.5 0.5], 'EdgeColor', 'none', 'FaceAlpha', 0.25)
end

%Q-Q plot
cols = repmat([0 0 1], np, 1);
if P.calc_fdr_values,
  qvalues = calc_fdr_value(pvalues);
  qvalues = qvalues(si, :);

  %near significance
  if size(pvalues, 2) == 1, idx = qvalues <= 0.25;
  else idx = qvalues(:, 3) <= 0.25; end
  cols(idx, :) = repmat([0 1 1], nnz(idx), 1);

  %at significance
  if size(pvalues, 2) == 1, idx = qvalues <= 0.1;
  else idx = qvalues(:, 3) <= 0.1; end
  cols(idx, :) = repmat([1 0 0], nnz(idx), 1);
end

[colu, ~, coluj] = unique(cols, 'rows');

if size(pvalues, 2) == 1,
  scatter(x, logpvalues, '.', 'CData', cols, 'MarkerEdgeColor', 'flat')
else
  ci_idx = -diff(logpvalues(:, 1:2), [], 2) > 0.25;
  scatter(x, logpvalues(:, 3), '.', 'CData', cols, 'MarkerEdgeColor', 'flat')

  for i = unique(coluj)', 
    idx = coluj == i & ci_idx;
    line([x(idx) x(idx)]', logpvalues(idx, 1:2)', 'Color', colu(i, :))
  end
end

%1-1 line
xl = xlim;
line(xl, [0 xl(2)], 'LineStyle', '--', 'Color', [0.3 0.3 0.3])

%labels
if ~isempty(labels),
  if max(size(labels)) == length(pvalues),
    labels = labels(si);
    textfit(x((end - 9):end), logpvalues((end - 9):end), labels((end - 9):end), 'fontsize', 8);
  else %assume labels for most significant are preselected
    nl = max(size(labels));
    textfit(x((end - nl + 1):end), logpvalues((end - nl + 1):end), flip(labels), 'fontsize', 8);
  end
end
