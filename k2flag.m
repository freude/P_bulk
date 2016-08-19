function [flag,pth]=k2flag(k1)

if (k1(1)~=0)
    if (k1(1)>0)
        flag='x';
    else
        flag='-x';
    end;
    pth='v0/';
end;
if (k1(2)~=0)
    if (k1(2)>0)
        flag='y';
    else
        flag='-y';
    end;
    pth='v1/';
end;
if (k1(3)~=0)
    if (k1(3)>0)
        flag='z';
    else
        flag='-z';
    end;
    pth='v2/';
end;
