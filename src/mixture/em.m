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
like = -inf;                          % P(all R)
showiters = false;
initident = false;
I = R; I(find(R == 0)) = 1;

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

pumrij = likelihood(R,I,ku,km,nu,nm,prum,pui,pmj);

for K = 1:maxiter % TODO when to stop?

  if showiters, fprintf('iteration %d: ', K); end;

  %
  %% e-step
  %
  % we calculate the likelihood here anyway (the denominator), so combine in
  % this part of the code the e-step and the stopping condition check.
  %
  % actually, this should already be calculated (either from initialization or
  % from the last stage's likelihood calculation), so we have nothing to do
  % here!
  %
  % pseudo-code:
  % for i = 1:nu
  %   for j = 1:nm
  %     if R(i,j) > 0
  %       pumrij(:,:,i,j) = prum(:,:,R(i,j)) .* (pui(:,i) * pmj(:,j)');
  %       pumrij(:,:,i,j) = pumrij(:,:,i,j) / sum(sum(pumrij(:,:,i,j)));
  %     end
  %   end
  % end
  %
  % alternative pseudo-code:
  % like = ones(nu,nm);
  % for i = 1:nu
  %   for j = 1:nm
  %     if R(i,j) > 0
  %       like(i,j) = 0;
  %       for u = 1:ku
  %         for m = 1:km
  %           like(i,j) = like(i,j) + prum(u,m,R(i,j)) * pui(u,i) * pmj(m,j);
  %         end
  %       end
  %     end
  %   end
  % end
  %

  %%%%pumrij = reshape(prum(:,:,I),ku,km,nu,nm) .* ...
  %%%%    (repmat(reshape(pui,ku,1,nu,1), [1 km 1 nm]) .* ...
  %%%%     repmat(reshape(pmj,1,km,1,nm), [ku 1 nu 1]));
  %%%%prij = repmat(sum(sum(pumrij,1),2), [ku km 1 1]);
  %%%%pumrij = pumrij ./ prij;
  %%%%pumrij(find(repmat(reshape(R,1,1,nu,nm),[ku km 1 1]) == 0)) = 1;

  %
  %% m-step
  %

  lastpui = pui;
  pui = reshape(sum(sum(pumrij,4),2),ku,nu);
  pui = pui./repmat(sum(pui,1),ku,1);
  assert(0 == numel(find(pui < 0 | 1 < pui)));

  lastpmj = pmj;
  pmj = reshape(sum(sum(pumrij,3),1),km,nm);
  pmj = pmj./repmat(sum(pmj,1),km,1);
  assert(0 == numel(find(pmj < 0 | 1 < pmj)));

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
  % pseudo-code:
  % like = ones(nu,nm);
  % for i = 1:nu
  %   for j = 1:nm
  %     if R(i,j) > 0
  %       like(i,j) = 0;
  %       for u = 1:ku
  %         for m = 1:km
  %           like(i,j) = like(i,j) + prum(u,m,R(i,j)) * pui(u,i) * pmj(m,j);
  %         end
  %       end
  %     end
  %   end
  % end
  % like = sum(log(reshape(like,numel(like),1)));
  % if showiters, fprintf('%f\n', like); end;
  % assert(not(0 < like || like < lastlike));

  lastlike = like;
  [pumrij prij] = likelihood(R,I,ku,km,nu,nm,prum,pui,pmj);
  like = sum(log(reshape(prij,numel(prij),1)));
  if showiters, fprintf('L = %f\n', like); end;
  assert(like < 0 && lastlike < like);
  if abs((like-lastlike)/lastlike) < tol
    break
  end
end

% newprum = zeros(kr,ku,km); % putting the r,u,m in the right spots
% for u = 1:ku
%   for m = 1:km
%     newprum(:,u,m) = prum(u,m,:);
%   end
% end
% prum = newprum;

prum = shiftdim(prum, 2);

return;

function [pumrij prij pjoint] = likelihood(R,I,ku,km,nu,nm,prum,pui,pmj);
pjoint = reshape(prum(:,:,I),ku,km,nu,nm) .* ...
    (repmat(reshape(pui,ku,1,nu,1), [1 km 1 nm]) .* ...
     repmat(reshape(pmj,1,km,1,nm), [ku 1 nu 1]));
prij = sum(sum(pjoint,1),2);
pumrij = pjoint ./ repmat(prij, [ku km 1 1]);
assert(0 == numel(find(pjoint < 0 | 1 < pjoint)));
assert(0 == numel(find(prij < 0 | 1 < prij)));
assert(0 == numel(find(pumrij < 0 | 1 < pumrij)));
% TODO figure out why this next line makes a difference (eg how do pui/pmj
% depend on the missing ratings?)
pumrij(find(repmat(reshape(R,1,1,nu,nm),[ku km 1 1]) == 0)) = 1;
prij(find(R == 0)) = 1;
return;

% DONE added random initialization - this seemed to do the trick
% DONE vectorized some parts of the code
% DONE fixed stupid indexing bug (n,m)
% DONE added support for skipping over unknown values
% TODO make more efficient by not calculating the unneeded values
% TODO make more efficient by vectorizing the for-loops

% vim:et:sw=2:ts=2
