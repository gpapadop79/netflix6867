%
%% dataset initialization
%

load Netflix_subset.mat
% rows should be users, cols should be movies
R = R';
[R Rtest] = divideData(R, 0.5);

% number of times to repeat EM for each unique configuration
tries = 3;

% the set of "configurations." A configuration is just a setting of ku and km
% (the numbers of parameters).
configs = [];
% a 'configuration counter'; each value corresponds to a different
% setting of ku and km
c = 1;
for ku = 2:maxku
  for km = 2:maxkm
    configs(c,:) = [ku km];
    c = c + 1;
  end
end
nc = size(configs, 1);

%% run the EM algorithm
if 0
  [nu nm] = size(R);
  % r is our results; later named emresults, now just r for brevity
  r = [];
  for c = 1:nc % configuration counter
    for t = 1:tries % the trial counter
      ku = configs(c,1); km = configs(c,2);
      fprintf('running with: ku = %d, km = %d, try # %d\n', ku, km, t);
      % TODO convert main.m into a function `em`. `em` should return the
      % parameters prum, pui, pmj
      [r(c,t).prum, r(c,t).pui, r(c,t).pmj r(c,t).like] = ...
          em(R, ku, km, 30, 1e-5);
    end
  end
  emresults = r;
  save emresults;
end

%% analyze the results: tables, plots
if 1
  load emresults;
  rs = emresults;

  [nu nm] = size(R);
  I = find(R ~= 0);

  % collect results for plotting
  bic = [];
  rmses = [];

  % print tables of BIC and likelihoods
  fprintf('%5s %5s %9s %14s %14s %14s\n', 'ku', 'km', 'best LL', 'RMSE', 'BIC penalty', 'BIC score');
  for c = 1:size(configs,1)
%%%%    for t = 1:tries
%%%%      % TODO define `likelihood` function
%%%%      L(c,t) = likelihood(R, prum, pui, pmj);
%%%%    end

    % find the best parameter settings for this configuration
    [like t] = max([rs(c,:).like]);
    r = rs(c,t);

    % compute the RMSE
    rmses(c) = rmse(Rtest, r.prum, r.pui, r.pmj);
%%%%    testlik = likelihood(Rtest, prum, pui, pm);

    % compute BIC
    nparams  = numel(r.prum) + numel(r.pui) + numel(r.pmj);
    penalty  = nparams/2 * log(numel(I));
    bic(c)   = like - penalty;

    % print table
    fprintf('%5d %5d %9.2f %14.5f %14.5f %14.5f\n', ...
        ku, km, like, rmses(c), penalty, bic(c));
  end

  % generate plots
  % TODO how to make a 3d plot? define `plotconfigs`
  figure(1); clf; plotconfigs(configs, bic);
  saveas(1, sprintf('em-bic-ku%d-km%d.png', maxku, maxkm));
  saveas(1, sprintf('em-bic-ku%d-km%d.pdf', maxku, maxkm));
  figure(2); clf; plotconfigs(configs, rmses);
  saveas(2, sprintf('em-perf-ku%d-km%d.png', maxku, maxkm));
  saveas(2, sprintf('em-perf-ku%d-km%d.pdf', maxku, maxkm));
end

% vim:et:sw=2:ts=2
