function [An, mIDsn]=findCommon(smID,A,mIDs)

mIDsn = mIDs;
[B mid]=readRatings(smID);

temp = size(find(mIDsn==mid));
pepepe = temp(2);


if pepepe~=0
    display('input video has already been in A.');
    An = A;
else
    mIDsn(end+1)=mid;

    temp = size(A);
    nA = temp(1);

    [An(1,:), IA, IB] = intersect(A(1,:),B(1,:));

    count=1;
    for i=IA
        An(2:nA,count)=A(2:nA,i);
        count=count+1;
    end

    count=1;    
    for i=IB
        An(nA+1,count)=B(2,i);
        count=count+1;
    end
end


    


