if 1
  ku = 2; nu = 4; km = 3; nm = 5;
  ku = 20; nu = 400; km = 30; nm = 500;
  a = rand(ku,nu); b = rand(km,nm);

  % alt 1
  tic;
  c = zeros(ku,km,nu,nm);
  for i=1:nu, for j=1:nm,
    c(:,:,i,j) = a(:,i) * b(:,j)';
    c(:,:,i,j) = c(:,:,i,j) / sum(sum(c(:,:,i,j)));
  end; end;
  toc;

  % alt 2
  tic;
  d = repmat(reshape(a,ku,1,nu,1), [1 km 1 nm]) .* ...
      repmat(reshape(b,1,km,1,nm), [ku 1 nu 1]);
  d = d ./ repmat(sum(sum(d,1),2), [ku km 1 1]);
  toc;

  assert(numel(find(c ~= d)) == 0);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if 0
  c(x(i)) = a(x(i)) + b(i);
end

% vim;et:sw=2:ts=2
