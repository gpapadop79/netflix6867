%% TODO WARNING UNTESTED CODE, just checking in

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
for ku = 2:5
  for km = 2:5
    configs(c,:) = [ku km];
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
      [r(c,t).prum, r(c,t).pui, r(c,t).pmj] = em(R, ku, km);
    end
  end
  emresults = r;
  save emresults;
end

%% analyze the results: tables, plots
if 0
  load results;
  rs = emresults;

  [nu nm] = size(R);
  I = find(R ~= 0);

  L = []; % likelihoods

  % print tables of BIC and likelihoods
  fprintf('%5s %5s %9s %14s %14s\n', 'ku', 'km', 'best LL', 'BIC penalty', 'BIC score');
  for c = 1:numel(configs)
    for t = 1:tries
      % TODO define `likelihood` function
      L(c,t) = likelihood(R, prum, pui, pmj);
    end

    % find the best parameter settings for this configuration
    [like t] = max(L(c,:));
    r = rs(c,t);

    % compute the test likelihood (TODO use some other test performance
    % metric?)
    testlik = likelihood(Rtest, prum, pui, pm);

    % compute BIC
    nparams = numel(r.prum) + numel(r.pui) + numel(r.pmj);
    penalty = nparams/2 * log(numel(I));
    bic     = like - penalty;

    % print table
    fprintf('%5d %5d %9.2f %14.5f %14.5f 14.5f\n', ku, km, like, testlik, penalty, bic);
  end

  % generate plots
  % TODO how to make a 3d plot? define `plotconfigs`
  figure(1); clf; plotconfigs(configs, bic);
  saveas(fprintf('em-bic-ku%d-km%d.pdf', ku, km));
  figure(2); clf; plotconfigs(configs, perf);
  saveas(fprintf('em-perf-ku%d-km%d.pdf', ku, km));
end

% vim:et:sw=2:ts=2
