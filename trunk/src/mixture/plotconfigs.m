%TODO how do i get the right 3d plot?

function plotconfigs(configs, bic);

%plot3(configs(:,1), configs(:,2), bic);
mat = zeros();
for c = 1:size(configs,1)
  mat(configs(c,1), configs(c,2)) = bic(c);
end
mesh(1:max(configs(:,1)), 1:max(configs(:,2)), mat);

% vim:et:sw=2:ts=2
