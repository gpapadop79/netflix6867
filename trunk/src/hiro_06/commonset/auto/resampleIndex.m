function out = resampleIndex(Wacc, w)

a=1;
b=numel(Wacc)+1;
while b-a~=1
    test = (b+a)/2;
    if Wacc(test)>=w
        b=test;
    else
        a=test;
    end
end

if w<a
    out = a;
else
    out = b;
end


   