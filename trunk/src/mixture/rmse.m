function [err err_mean Rhat Rhatmean] = rmse(R, prum, pui, pmj)

kr = size(prum,1); ku = size(pui,1); km = size(pmj,1);
[nu nm] = size(R);
like = zeros(kr,nu,nm);

% estimate the most likely ratings for each of the test (i,j)
for i = 1:nu
  for j = 1:nm
    if R(i,j) > 0
      % find likelihood of each rating for each user-type/movie-type config.
      % marginal = P(R,U,M|I,J)
      marginal = prum .* shiftdim(repmat(pui(:,i) * pmj(:,j)', [1 1 5]), 2);
      % marginalize over the user and movie types
      like(:,i,j) = sum(sum(marginal, 3), 2);
    end
  end
end

% which r's are most likely? those are our estimates
[l Rhat] = max(like,[],1);
Rhat = reshape(Rhat, nu, nm);


%% use mean as the estimation
Rhatmean = zeros(nu,nm);
for i=1:5
    Rhatmean = Rhatmean + i * reshape(like(i,:,:), nu, nm);
end


% find RMSE between our estimates and the true ratings (ignoring the non-test
% cells)
errs = (Rhat - R).^2 .* (R > 0);
err = sqrt(sum(reshape(errs, numel(errs), 1)) / numel(find(R>0)));

errs_mean = (Rhatmean - R).^2 .* (R > 0);
err_mean = sqrt(sum(reshape(errs_mean, numel(errs_mean), 1)) / numel(find(R>0)));

% vim:et:sw=2:ts=2
