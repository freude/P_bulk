function [V,varargout]=pot_mat(coord,k1,k2)

a=1e-7;

%V=(coord.coord_limits(2)-coord.coord_limits(1)+coord.coord_stps)^3;

x=(coord.coord_limits(1):coord.coord_stps:coord.coord_limits(2));
y=(coord.coord_limits(1):coord.coord_stps:coord.coord_limits(2));
z=(coord.coord_limits(1):coord.coord_stps:coord.coord_limits(2));

[X,Y,Z]=meshgrid(x,y,z);


x1=-3e-9/MyConst.ab;
y1=0;
z1=0;

x2=3e-9/MyConst.ab;
y2=0;
z2=0;

[~,ix1]=min(abs(x-x1));
[~,iy1]=min(abs(y-y1));
[~,iz1]=min(abs(z-z1));

[~,ix2]=min(abs(x-x2));
[~,iy2]=min(abs(y-y2));
[~,iz2]=min(abs(z-z2));

x1=x(ix1)-0.5*coord.coord_stps;
y1=y(iy1)-0.5*coord.coord_stps;
z1=z(iz1)-0.5*coord.coord_stps;

x2=x(ix2)-0.5*coord.coord_stps;
y2=y(iy2)-0.5*coord.coord_stps;
z2=z(iz2)-0.5*coord.coord_stps;

x1=0-0*0.5*coord.coord_stps;
y1=0-0*0.5*coord.coord_stps;
z1=0-0*0.5*coord.coord_stps;

ij=sqrt(-1);
V=real(1/sqrt((X-x1).^2+(Y-y1).^2+(Z-z1).^2+ij*a)).*pot_scr_real(sqrt((X-x1).^2+(Y-y1).^2+(Z-z1).^2));

%V(isinf(V))=0;

V=V.*exp(ij*((k1(1)-k2(1)).*(X-x1)+(k1(2)-k2(2)).*(Y-y1)+(k1(3)-k2(3)).*(Z-z1)));

if nargout>1
    varargaout{1}=X;
    varargaout{2}=Y;
    varargaout{3}=Z;
end;
