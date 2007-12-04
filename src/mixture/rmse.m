function err = rmse(R, prum, pui, pmj);

kr = size(prum,1); ku = size(pui,1); km = size(pmj,1);
[nu nm] = size(R);
%%%% lprum = log(prum); lpui = log(pui); lpmj = log(pmj);
like = zeros(kr,nu,nm);

% estimate the most likely ratings for each of the test (i,j)
for i = 1:nu
  for j = 1:nm
    if R(i,j) > 0
      % find likelihood of each rating
      for r = 1:kr
        marginal = reshape(prum(r,:,:),ku,km) .* (pui(:,i) * pmj(:,j)');
        like(r,i,j) = like(r,i,j) + sum(reshape(marginal, ku*km, 1));
      end
    end
  end
end

% which r's are most likely? those are our estimates
[l Rhat] = max(like,[],1);
Rhat = reshape(Rhat, nu, nm);

% find RMSE between our estimates and the true ratings (ignoring the non-test
% cells)
errs = (Rhat - R).^2 .* (R > 0);
err = sqrt(sum(reshape(errs, numel(errs), 1)) / numel(find(R>0)));

% vim:et:sw=2:ts=2
