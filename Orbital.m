% Plotting hydrogen orbitals

close all;

% Quantum numbers ======================================================

% n=1; 
% l=0;  %  0<= l <n 
 m=0;  % -l<= m <=l

%=======================================================================

a=1; % Bohr radius

% Normalization
N=abs(sign(m)*sqrt(2)+(sign(abs(m))-1)*2);

% Angular part
SphericalYlm=@(l,m,theta,phi) sqrt((2*l+1)/(4*pi)*factorial(l-abs(m))/...
    factorial(l+abs(m)))*AssociatedLegendre(l,m,cos(theta)).*exp(1i*m*phi);
Y=@(l,m,theta,phi) (SphericalYlm(l,m,theta,phi)+SphericalYlm(l,-m,theta,phi))/N;
% Radial part
R=@(n,l,r) sqrt((2/(a*n))^3*factorial(n-l-1)/(2*n*factorial(n+l))).*...
    exp(-r/(a*n)).*(2*r/(a*n)).^l*1/factorial(n-l-1+2*l+1).*...
    AssociatedLaguerre(n-l-1,2*l+1,2*r/(a*n));
% Wave function
psi=@(n,l,m,r,theta,phi) R(n,l,r).*Y(l,m,theta,phi);



% Setting the grid
border=32;
accuracy=100;
[x,y,z]=ndgrid(linspace(-border,border,accuracy),linspace(-border,border,accuracy),linspace(-border,border,accuracy));
% Conversion Cartesian to spherical coordinates
r=sqrt(x.^2+y.^2+z.^2);
theta=acos(z./r);
phi=atan2(y,x);

YY=psi(1,0,0,r,theta,phi);
 
plot(squeeze(YY(:,50,50)))


% Plot orbital,  - and + wave function phase
colors=sign(psi(n,l,m,r,theta,phi));
probability=1E-5;
isosurface(psi(n,l,m,r,theta,phi).^2,probability,colors);
title(['n = ' num2str(n) ', l = ' num2str(l) ', m = ' num2str(m)],'FontName','Times','FontSize',12);
set(gcf,'color',[1 1 1]);
daspect([1 1 1]);
axis off;
view(3);
camlight('left');
camzoom(0.95);
lighting phong;
axis vis3d;
rotate3d on;
brighten(1);

