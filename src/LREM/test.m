clear;

load Netflix_subset;

dim = [1 2 3];
sparsity = [0.3 0.5 0.8];

e = zeros(length(dim),length(sparsity));

for d = dim
    c = 0;
    for i = sparsity
        c = c+1;
        for j = 1:10
            [trainE testE] = main(R,i,d,0);
            e(d,c) = e(d,c) + testE(end);
        end
        disp('hoge');
    end
end

e = e./10;
