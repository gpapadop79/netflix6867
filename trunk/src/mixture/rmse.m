function err = rmse(R, prum, pui, pmj);

err = 0; kr = size(prum,1); ku = size(pui,1); km = size(pmj,1);
[nu nm] = size(R);

% estimate the most likely ratings for each of the test (i,j)
for i = 1:nu
  for j = 1:nm
    if R(i,j) > 0
      % find likelihood of each rating
      for r = 1:kr
        like(r,i,j) = 0;
        for u = 1:ku
          for m = 1:km
            like(r,i,j) = like(r,i,j) + ...
                log(prum(r,u,m)) + log(pui(u,i)) + log(pmj(m,j));
          end
        end
      end
      % which one is most likely? that's our estimate
      [l r] = max(like(:,i,j));
      % find RMSE between our estimate and the true rating
      err = err + (r - R(i,j))^2;
    end
  end
end

err = sqrt(err / numel(err));

% vim:et:sw=2:ts=2
