function [prum pui pmj pumrij] = em_hiro_edited(R, ku, km);

%
%% initialization
%

[nu nm] = size(R);
kr = 5; km = 2; ku = 2;

initident = 0;

% the posteriors initialization doesn't matter, since the E-step runs first
pumrij = ones(ku,km,nu,nm);         % P(U,M|R,I,J)
% the priors should be uniformly initialized
pui    = 1/ku * ones(ku,nu);          % P(U|I)
pmj    = 1/km * ones(km,nm);          % P(M|J)
% the likelihood must not be identically initialized; we initialize it randomly
% (but the random distributions must be normalized)
prum   = rand(ku,km,kr);             % P(R|U,M)
for u = ku
  for m = km
    prum(u,m,:) = prum(u,m,:) / sum(prum(u,m,:));
  end
end
if initident, prum = 1/kr * ones(ku,km, kr); end;

for K = 1:100 % TODO when to stop?

  tic;
  fprintf('iteration %d\n', K);

  %
  %% e-step
  %

  %disp pumrij;
  for i = 1:nu
    for j = 1:nm
        if R(i,j) > 0
            pumrij(:,:,i,j) = prum(:,:,R(i,j)) .* (pui(:,i) * pmj(:,j)');
            pumrij(:,:,i,j) = pumrij(:,:,i,j) / sum(sum(pumrij(:,:,i,j)));
        end
    end
  end
  
  toc;

  % TODO make sure probs sum to 1; repeat for other mats
  %if find(sum(:,:,:,:,:) != 1), err(); end;
  if find(pumrij < 0 | 1 < pumrij), err(); end;

  %
  %% m-step
  %

  %disp pui
  lastpui = pui;
  pui = reshape(sum(sum(pumrij,4),2),ku,nu);
  pui = pui./repmat(sum(pui,1),ku,1);
  
  if find(pui < 0 | 1 < pui), err(); end;

  %disp pmj
  lastpmj = pmj;
  pmj = reshape(sum(sum(pumrij,3),1),km,nm);
  pmj = pmj./repmat(sum(pmj,1),km,1);

  if find(pmj < 0 | 1 < pmj), err(); end;

  toc;
  
  %disp prum;
  lastprum = prum;
  prum = zeros(ku,km,kr);
  for u = 1:ku
    for m = 1:km
      for i = 1:nu
        for j = 1:nm
          if R(i,j) > 0
            prum(u,m,R(i,j)) = prum(u,m,R(i,j)) + pumrij(u,m,i,j);
          end
        end
      end
      prum(u,m,:) = prum(u,m,:) / sum(prum(u,m,:));
    end
  end

  if find(prum < 0 | 1 < prum), err(); end;

  if numel(find(pui == lastpui)) > 0 && ...
     numel(find(pmj == lastpmj)) > 0 && ...
     numel(find(prum == lastprum)) > 0
    break;
  end;

  toc;
end

% DONE added random initialization - this seemed to do the trick
% DONE vectorized some parts of the code
% DONE fixed stupid indexing bug (n,m)
% DONE added support for skipping over unknown values
% TODO make more efficient by not calculating the unneeded values
% TODO make more efficient by vectorizing the for-loops

% vim:et:sw=2:ts=2
