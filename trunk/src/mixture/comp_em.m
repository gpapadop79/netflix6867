clear;
%load ../dataset/Netflix_subset.mat;

R = [5 5 5 0; ... % :)
     2 1 4 5; ... % :(
     5 1 2 5; ... % :)
     1 0 1 5];    % :(
 
 

ku = 2;
km = 2;
maxiter = 30;
tol = 1e-4;

%% initialization

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



%% Old initialization
% the posteriors initialization doesn't matter, since the E-step runs first
pumrij_old = ones(ku,km,nu,nm);        % P(U,M|R,I,J)
% the priors should be uniformly initialized
pui_old    = 1/ku * ones(ku,nu);          % P(U|I)
pmj_old    = 1/km * ones(km,nm);          % P(M|J)
% the likelihood must not be identically initialized; we initialize it randomly
% (but the random distributions must be normalized)
for r = 1:5
    for u = 1:ku
        for m=1:km
            prum_old(r,u,m) = prum(u,m,r);
        end
    end
end

    
%%





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
 
%% Old E-step
  
    for i = 1:nu
    for j = 1:nm
        if R(i,j)>0
        for m = 1:km
          for u = 1:ku
            pumrij_old(u,m,i,j) = prum_old(R(i,j),u,m) * pui_old(u,i) * pmj_old(m,j);
          end
        end
        pumrij_old(:,:,i,j) = pumrij_old(:,:,i,j) / sum(sum(pumrij_old(:,:,i,j)));
        end
    end
    end

%% Compare E-step result 

tol = 1e-6;

Ind = find(pumrij_old~=pumrij);

if((approx_eq(pumrij, pumrij_old, tol)>0))
    disp('E-step: results are identical.');
else
    disp('E-step: results are different.');
end

    
%% New M step

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

  assert(0 == numel(find(pmj < 0 | 1 < pmj)));
  
%% Old M step

  lastpui = pui_old;
  pui_old = zeros(ku,nu);
  for i = 1:nu
    for u = 1:ku
      for j = 1:nm
        if R(i,j) > 0
          pui_old(u,i) = pui_old(u,i) + sum(pumrij_old(u,:,i,j)); % + P(U|R,I,J)
        end
      end
      pui_old(u,i) = pui_old(u,i) / numel(find(R(i,:) > 0));
    end
  end


  %disp pmj
  lastpmj = pmj_old;
  pmj_old = zeros(km,nm);
  for j = 1:nm
    for m = 1:km
      for i = 1:nu
        if R(i,j) > 0
          pmj_old(m,j) = pmj_old(m,j) + sum(pumrij_old(:,m,i,j)); % + P(M|R,I,J)
        end
      end
      pmj_old(m,j) = pmj_old(m,j) / numel(find(R(:,j) > 0));
    end
  end

  if find(pmj < 0 | 1 < pmj), err(); end;

  
%% Compare

tol = 1e-6;

if(min(min(approx_eq(pui, pui_old, tol)))>0);
    disp('M: pui is identical');
else
    disp('M: pui is different');
end

if(min(min(approx_eq(pmj, pmj_old, tol)))>0);
    disp('M: pmj is identical');
else
    disp('M: pmj is different');
end


%%  new M

 
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

%% old M

  lastprum = prum_old;
  prum_old = zeros(kr,ku,km);
  for u = 1:ku
    for m = 1:km
      for i = 1:nu
        for j = 1:nm
          if R(i,j) > 0
            prum_old(R(i,j),u,m) = prum_old(R(i,j),u,m) + pumrij_old(u,m,i,j);
          end
        end
      end
      prum_old(:,u,m) = prum_old(:,u,m) / sum(prum_old(:,u,m));
    end
  end
  
%% Compare  

if(min(min(approx_eq(shiftdim(prum,2), prum_old, tol)))>0);
    disp('M: prum is identical');
else
    disp('M: prum is different');
end

%%

  
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
  
  input('press any key to proceed');

  
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
