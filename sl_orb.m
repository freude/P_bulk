function f=sl_orb(x,y,z,p0,sigma)

D=sqrt((x-p0(1)).^2+(y-p0(2)).^2+(z-p0(3)).^2);
f=sqrt((sigma^3)/pi)*exp(-D*sigma);