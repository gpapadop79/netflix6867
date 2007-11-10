d_kernel = 10;
n_kernel = 30;

sizeA = size(R);

sparse_skirt = 0.2
sparse_other = 0.05

k=1;

%build row skirt
for i=1:d_kernel
    index = randperm(sizeA(2)-n_kernel);
    count=1;
    while count < (sizeA(2)-n_kernel)*sparse_skirt
        data_n(k,1)=i;
        data_n(k,2)=index(count)+n_kernel;
        data_n(k,3)=R(data_n(k,1), data_n(k,2));
        count=count+1;
        k=k+1;
    end
end

sizeD1=size(data_n);
rp = randperm(sizeD1(1));
for i=1:sizeD1(1)
    data(i,:)=data_n(rp(i),:);
end



%build column skirt
for i=1:n_kernel
    index = randperm(sizeA(1)-d_kernel);
    count=1;
    while count < (sizeA(1)-d_kernel)*sparse_skirt
        data_n(k,2)=i;
        data_n(k,1)=index(count)+d_kernel;
        data_n(k,3)=R(data_n(k,1), data_n(k,2));
        count=count+1;
        k=k+1;
    end
end

sizeD2=size(data_n);
rp = randperm(sizeD2(1)-sizeD1(1));
for i=sizeD1(1)+1:sizeD2(1)
    data(i,:)=data_n(rp(i-sizeD1(1))+sizeD1(1),:);
end



%build other part
for i=d_kernel+1:sizeA(1)
    index = randperm(sizeA(2)-n_kernel);
    count=1;
    while count < (sizeA(2)-n_kernel)*sparse_other
        data_n(k,1)=i;
        data_n(k,2)=index(count)+n_kernel;
        data_n(k,3)=R(data_n(k,1), data_n(k,2));
        count=count+1;
        k=k+1;
    end
end

sizeD3=size(data_n);
rp = randperm(sizeD3(1)-sizeD2(1));
for i=sizeD2(1)+1:sizeD3(1)
    data(i,:)=data_n(rp(i-sizeD2(1))+sizeD2(1),:);
end

clear data_n;