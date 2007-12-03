function [prum pui pmj pumrij] = em(R, ku, km);

%
%% initialization
%

[nu nm] = size(R);
kr = 5; km = 2; ku = 2;

initident = 0;

% the posteriors initialization doesn't matter, since the E-step runs first
pumrij = ones(ku,km,kr,nu,nm);         % P(U,M|R,I,J)
% the priors should be uniformly initialized
pui    = 1/ku * ones(ku,nu);          % P(U|I)
pmj    = 1/km * ones(km,nm);          % P(M|J)
% the likelihood must not be identically initialized; we initialize it randomly
% (but the random distributions must be normalized)
prum   = rand(kr,ku,km);             % P(R|U,M)
for u = ku
  for m = km
    prum(:,u,m) = prum(:,u,m) / sum(prum(:,u,m));
  end
end
if initident, prum = 1/kr * ones(kr,ku,km); end;

for K = 1:100 % TODO when to stop?

  fprintf('iteration %d\n', K);

  %
  %% e-step
  %

  %disp pumrij;
  for i = 1:nu
    for j = 1:nm
      for r = 1:kr
        for m = 1:km
          for u = 1:ku
            pumrij(u,m,r,i,j) = prum(r,u,m) * pui(u,i) * pmj(m,j);
          end
        end
        pumrij(:,:,r,i,j) = pumrij(:,:,r,i,j) / sum(sum(pumrij(:,:,r,i,j)));
      end
    end
  end

  % TODO make sure probs sum to 1; repeat for other mats
  %if find(sum(:,:,:,:,:) != 1), err(); end;
  if find(pumrij < 0 | 1 < pumrij), err(); end;

  %
  %% m-step
  %

  %disp pui
  lastpui = pui;
  pui = zeros(ku,nu);
  for i = 1:nu
    for u = 1:ku
      for j = 1:nm
        if R(i,j) > 0
          pui(u,i) = pui(u,i) + sum(pumrij(u,:,R(i,j),i,j)); % + P(U|R,I,J)
        end
      end
      pui(u,i) = pui(u,i) / numel(find(R(i,:) > 0));
    end
  end

  if find(pui < 0 | 1 < pui), err(); end;

  %disp pmj
  lastpmj = pmj;
  pmj = zeros(km,nm);
  for j = 1:nm
    for m = 1:km
      for i = 1:nu
        if R(i,j) > 0
          pmj(m,j) = pmj(m,j) + sum(pumrij(:,m,R(i,j),i,j)); % + P(M|R,I,J)
        end
      end
      pmj(m,j) = pmj(m,j) / numel(find(R(:,j) > 0));
    end
  end

  if find(pmj < 0 | 1 < pmj), err(); end;

  %disp prum;
  lastprum = prum;
  prum = zeros(kr,ku,km);
  for u = 1:ku
    for m = 1:km
      for i = 1:nu
        for j = 1:nm
          if R(i,j) > 0
            prum(R(i,j),u,m) = prum(R(i,j),u,m) + pumrij(u,m,R(i,j),i,j);
          end
        end
      end
      prum(:,u,m) = prum(:,u,m) / sum(prum(:,u,m));
    end
  end

  if find(prum < 0 | 1 < prum), err(); end;

  if numel(find(pui == lastpui)) > 0 && ...
     numel(find(pmj == lastpmj)) > 0 && ...
     numel(find(prum == lastprum)) > 0
    break;
  end;

end

% DONE added random initialization - this seemed to do the trick
% DONE vectorized some parts of the code
% DONE fixed stupid indexing bug (n,m)
% DONE added support for skipping over unknown values
% TODO make more efficient by not calculating the unneeded values
% TODO make more efficient by vectorizing the for-loops

% vim:et:sw=2:ts=2
