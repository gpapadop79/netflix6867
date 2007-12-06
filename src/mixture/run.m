mkdatasets = 0; runem = 0; analyze = 1;

%
%% generate and save the datasets and configs
%

if mkdatasets
  kus = 2.^(1:5); kms = 2.^(1:5);
  maxkm = max(kms); maxku = max(kus);

  % the set of "configurations." A configuration is just a setting of ku and km
  % (the numbers of parameters).
  if 0
    configs = [];
    c = 1;
    for sparsity = [.2 .5 .7]
      for ku = kus
        for km = kms
          configs(c,:) = [sparsity ku km];
          c = c + 1;
        end
      end
    end
  else
    c = 1;
    for sparsity = [.2 .5 .7]
      for k = 1:length(kms)
        configs(c,:) = [sparsity kus(k) kms(k)];
        c = c + 1;
      end
    end
  end
  assert(size(configs,2) == 3);

  load Netflix_subset.mat;
  % rows should be users, cols should be movies
  R = R';
  Rs = [];
  % try different sparsity ratios
  sparsities = unique(configs(:,1));
  for s = 1:length(sparsities)
    [Rs(s).train Rs(s).test] = divideData(R, sparsities(s));
  end
  save datasets Rs configs;
end

%
%% run the EM algorithm
%

if runem
  % number of times to repeat EM for each unique configuration
  tries = 3;

  load datasets Rs configs;
  res = [];
  sparsities = unique(configs(:,1));
  for c = 1:size(configs,1) % configuration counter
    for t = 1:tries % the trial counter
      sparsity = configs(c,1); ku = configs(c,2); km = configs(c,3);
      s = find(sparsities == sparsity);
      fprintf('running with: sparsity = %f, ku = %d, km = %d, try # %d\n', ....
              sparsity, ku, km, t);
      R = Rs(s).train;
      % TODO convert main.m into a function `em`. `em` should return the
      % parameters prum, pui, pmj
      [res(s,c,t).prum, res(s,c,t).pui, res(s,c,t).pmj res(s,c,t).like] = ...
          em(R, ku, km);
    end
  end
  save emresults res;
end

%
%% analyze the results: tables, plots
%

if analyze
  load Netflix_subset.mat R;
  [nu nm] = size(R);
  load datasets Rs configs;

  load emresults res;
  sparsities = unique(configs(:,1));

  % collect results for plotting
  bic = [];
  rmsesml = []; rmseswm = [];

  % print tables of BIC and likelihoods
  fprintf('%14s %5s %5s %14s %9s %9s %14s %14s\n', ...
          'sparsity', 'ku', 'km', 'best LL', 'RMSE-ML', 'RMSE-WM', 'BIC penalty', 'BIC score');
  for c = 1:size(configs,1)
    sparsity = configs(c,1); ku = configs(c,2); km = configs(c,3);
    s = find(sparsities == sparsity);
    R = Rs(s);
    I = find(R.train ~= 0);

    % find the best parameter settings for this configuration
    [like t] = max([res(s,c,:).like]);
    r = res(s,c,t);

    % compute the RMSE
    [rmsesml(c) rmseswm(c)] = rmse(R.test, r.prum, r.pui, r.pmj);

    % compute BIC
    nparams  = numel(r.prum) + numel(r.pui) + numel(r.pmj);
    penalty  = nparams/2 * log(numel(I));
    bic(c)   = like - penalty;

    % print table
    fprintf('%14.5f %5d %5d %14.5f %9.5f %9.5f %14.5f %14.5f\n', ...
            sparsity, ku, km, like, rmsesml(c), rmseswm(c), penalty, bic(c));
  end

  % generate plots
  if 0
    % TODO how to make a 3d plot? define `plotconfigs`
    figure(1); clf; plotconfigs(configs, bic);
    saveas(1, sprintf('em-bic-ku%d-km%d.png', maxku, maxkm));
    saveas(1, sprintf('em-bic-ku%d-km%d.pdf', maxku, maxkm));
    figure(2); clf; plotconfigs(configs, rmses);
    saveas(2, sprintf('em-perf-ku%d-km%d.png', maxku, maxkm));
    saveas(2, sprintf('em-perf-ku%d-km%d.pdf', maxku, maxkm));
  end
end

% vim:et:sw=2:ts=2
