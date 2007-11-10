stopCondition=1;

files = ls;
temp=size(files);
N=temp(1);

if exist('A','var')==1 && exist('mid','var')==1
    firstE=1;
    display('A and mid loaded');
else
    firstE=0;
end

for i=4:N
    if strfind(files(i,:),'mv_')==1
        if firstE == 0
            [A mid] = readRatings2(files(i,:));
            firstE=1;
        else
            [A mid] = findCommon2(files(i,:),A,mid);
            if stopCondition==1
                temp=size(mid);
                n=temp(2);
                temp=size(A);
                m=temp(2);
                if n>=m
                    break;
                end
            end
            
        end
    end
end
 
save result;
        
    