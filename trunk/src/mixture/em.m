function [prum pui pmj like pumrij] = em(R, ku, km, maxiter, tol);
% [prum pui pmj pumrij like] = em(R, ku, km, maxiter, tol)
%
% EM algorithm for flexible mixture model (collab filt) MLE
%
% Inputs:
%   R: the dataset
%   ku: the number of user types
%   km: the number of movie types
%   maxiter: cap on the number of iterations (default 30)
%   tol: stop iterating once |L'-L|/L<tol (default 1e-4)
%
% Outputs:
%   prum:   P(R|U,M)
%   pui:    P(U|I) (user type prior)
%   pmj:    P(M|J) (movie type prior)
%   like:   P(R) (likelihood)
%   pumrij: P(U,M|R,I,J) (posterior)

%
%% initialization
%

if nargin < 4, maxiter = 30; end;
if nargin < 5, tol = 1e-4; end; % termination delta ratio tolerance
[nu nm] = size(R);
kr = 5;
like = -inf;
showiters = false;
initident = false;

% the posteriors initialization doesn't matter, since the E-step runs first
pumrij = ones(ku,km,nu,nm);           % P(U,M|R,I,J)
% the priors should be uniformly initialized
pui    = 1/ku * ones(ku,nu);          % P(U|I)
pmj    = 1/km * ones(km,nm);          % P(M|J)
% the likelihood must not be identically initialized; we initialize it randomly
% (but the random distributions must be normalized)
prum   = rand(ku,km,kr);              % P(R|U,M) XXX NOTE prum(u,m,r)
for u = ku
  for m = km
    prum(u,m,:) = prum(u,m,:) / sum(prum(u,m,:));
  end
end
if initident, prum = 1/kr * ones(ku,km,kr); end;

for K = 1:maxiter % TODO when to stop?

  if showiters, fprintf('iteration %d\n', K); end;

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

  % TODO make sure probs sum to 1; repeat for other mats
  assert(0 == numel(find(pumrij < 0 | 1 < pumrij)));

  %
  %% m-step
  %

  %disp pui
  lastpui = pui;
  pui = reshape(sum(sum(pumrij,4),2),ku,nu);
  pui = pui./repmat(sum(pui,1),ku,1);
  
  assert(0 == numel(find(pui < 0 | 1 < pui)));

  %disp pmj
  lastpmj = pmj;
  pmj = reshape(sum(sum(pumrij,3),1),km,nm);
  pmj = pmj./repmat(sum(pmj,1),km,1);

  assert(0 == numel(find(pmj < 0 | 1 < pmj)));
  
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

  
  
  assert(0 == numel(find(prum < 0 | 1 < prum)));

  %
  % calculate likelihood
  %

  lastlike = like;
  like = ones(nu,nm);
  for i = 1:nu
    for j = 1:nm
      if R(i,j) > 0
        like(i,j) = 0;
        for u = 1:ku
          for m = 1:km
            like(i,j) = like(i,j) + prum(u,m,R(i,j)) * pui(u,i) * pmj(m,j);
          end
        end
      end
    end
  end
  like = sum(log(reshape(like,numel(like),1)));
  if showiters, fprintf('%d\n', like); end;
  assert(not(0 < like || like < lastlike));

  if abs((like-lastlike)/lastlike) < tol
    break
  end

  %%%%if numel(find(pui == lastpui)) > 0 && ...
  %%%%   numel(find(pmj == lastpmj)) > 0 && ...
  %%%%   numel(find(prum == lastprum)) > 0
  %%%%  break;
  %%%%end;
end

newprum = zeros(kr,ku,km); % putting the r,u,m in the right spots
for u = 1:ku
  for m = 1:km
    newprum(:,u,m) = prum(u,m,:);
  end
end
prum = newprum;

% DONE added random initialization - this seemed to do the trick
% DONE vectorized some parts of the code
% DONE fixed stupid indexing bug (n,m)
% DONE added support for skipping over unknown values
% TODO make more efficient by not calculating the unneeded values
% TODO make more efficient by vectorizing the for-loops

% vim:et:sw=2:ts=2
