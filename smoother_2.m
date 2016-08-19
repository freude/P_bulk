function V1 = smoother_2(x,G,k1,k2,R_tr)

ij=sqrt(-1);

if (((k1(find(k1))>0)&&(k2(find(k2))<0))||((k2(find(k2))>0)&&(k1(find(k1))<0)))&&(find(k1)==find(k2))
    M1=1;
else
    M1=1;
end;

N=length(x);

flag1=k2flag(k1);
flag2=k2flag(k2);

coord_stps=x(2)-x(1);
kx=linspace(-pi/coord_stps, (pi-2*pi/N)/coord_stps, N);

[kX,kY,kZ]=meshgrid(kx,kx,kx);

bi=br_zone_valley(kX./MyConst.ab, kY./MyConst.ab, kZ./MyConst.ab, flag1, 0);
kk=-(k1-k2);
%kk=-k1;
 
V=zeros(size(kX));
%R_tr=7;

for j=1:50
   qq=sqrt((kX-M1*kk(1)-G(j,1)).^2+(kY-M1*kk(2)-G(j,2)).^2+(kZ-M1*kk(3)-G(j,3)).^2); 
   if flag1==flag2
       VV=4.*pi./(qq.^2).*(1-cos(qq.*R_tr));
       VV(isnan(VV))=4.*pi.*0.5*R_tr^2;       
       V=V-pot_scr(qq).*(G(j,4)+ij*G(j,5)).*VV;      
   else
       V=V-pot_scr(qq).*(G(j,4)+ij*G(j,5)).*4.*pi./(qq.^2);               
   end;
end;

V=bi.*V;
%bi=br_zone_valley(kX./MyConst.ab, kY./MyConst.ab, kZ./MyConst.ab, flag2, 1);
%bi=fftshift(ifftn(ifftshift(bi))).*(abs(kx(2)-kx(1))^3)./((2*pi)^3).*N^3;
%bi=trapz(x,trapz(x,trapz(x,(bi),3),2),1);
V1=fftshift(fftn(ifftshift(V))).*(abs(kx(2)-kx(1))^3)./((2*pi)^3);

