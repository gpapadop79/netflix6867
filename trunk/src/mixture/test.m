kr = 5; ku = 2; km = 2;

% prepare our mini dataset
R = [4 5 3; ... % :)
     2 1 1; ... % :(
     5 5 4; ... % :)
     1 1 1];    % :(
R = [5 5 5 0; ... % :)
     2 1 4 5; ... % :(
     5 1 2 5; ... % :)
     1 0 1 5];    % :(

% run EM
[prum pui pmj like pumrij] = oldem2(R, ku, km);

% test the first unknown
i = 1; j = 4;
for r = 1:kr
  prob = 0;
  for u = 1:ku
    for m = 1:km
      prob = prob + prum(r,u,m) * pui(u,i) * pmj(m,j);
    end
  end
  fprintf('P{ R(%d,%d) = %d } = %f\n', i, j, r, prob);
end
% test the second unknown
i = 4; j = 2;
for r = 1:kr
  prob = 0;
  for u = 1:ku
    for m = 1:km
      prob = prob + prum(r,u,m) * pui(u,i) * pmj(m,j);
    end
  end
  fprintf('P{ R(%d,%d) = %d } = %f\n', i, j, r, prob);
end

fprintf('like = %f\n', like);

% vim:et:sw=2:ts=2
