function [out mID] = readRatings2(file)

%Format
%  out = | ID1  ID2  ID3  ID4...  |
%        | R1   R2   R3   R4 ...  |
% where IDi is the customer ID of ith customer
% mID is the ID number of the movie
% Ri is ith customer's rating for this movie. 

fid=fopen(file,'rt');
mID=fscanf(fid,'%d');
fscanf(fid,':');
A = fscanf(fid,'%d,%d,%s\n');

temp=size(A);
N = temp(1)/12;
out=zeros(2,N);

for i=1:N
    out(1,i)=A(12*(i-1)+1);
    out(2,i)=A(12*(i-1)+2);
end

fclose(fid);
    