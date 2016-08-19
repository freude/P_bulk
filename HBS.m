clear all;

% num_cells=4;
% T=80;

%num_cells=60;
%T=8;

N_states = 14;

qn=[1 0  0 -13.6;...  % shell 1
    2 0  0 -3.4;...   % shell 2
    2 1 -1 -3.4;...
    2 1  0 -3.4;...
    2 1  1 -3.4;...
    3 0  0 -1.511;... % shell 3
    3 1 -1 -1.511;...
    3 1  0 -1.511;...
    3 1  1 -1.511;...    
    3 2 -2 -1.511;...
    3 2 -1 -1.511;...
    3 2  0 -1.511;...        
    3 2  1 -1.511 ;...        
    3 2  2 -1.511];

%N_states = 8;

%qn=[1 0  0 -13.6;...
%    2 0  0 -3.4;...
%    3 0  0 -1.511;...
%    4 0  0 -0.85;...
%    5 0  0 -0.544;...
%    6 0  0 -0.3778;...
%    7 0  0 -0.2776;...
%    8 0  0 -0.2125];

num_cells=280;
T=2;

coorsys=CoordSys(num_cells,T,'au');
coorsys.set_origin_cells(num_cells/2+1);
x=coorsys.x();
[X1,Y1,Z1]=meshgrid(x,x,x);

a=1; % Bohr radius

% Normalization

%charge=1.24;
charge=1;

% Angular part
SphericalYlm=@(l,m,theta,phi) sqrt((2*l+1)/(4*pi)*factorial(l-abs(m))/...
    factorial(l+abs(m)))*AssociatedLegendre(l,m,cos(theta)).*exp(1i*m*phi);

Y=@(l,m,theta,phi) (SphericalYlm(l,m,theta,phi)+SphericalYlm(l,-m,theta,phi))/...
    (abs(sign(m)*sqrt(2)+(sign(abs(m))-1)*2));

% Radial part
R=@(n,l,r) sqrt((charge*2/(a*n))^3*factorial(n-l-1)/(2*n*factorial(n+l))).*...
    exp(-charge*r/(a*n)).*(2*charge*r/(a*n)).^l*1/factorial(n-l-1+2*l+1).*...
    AssociatedLaguerre(n-l-1,2*l+1,2*charge*r/(a*n));

% Wave function
psi=@(n,l,m,r,theta,phi) R(n,l,r).*Y(l,m,theta,phi);


% Conversion Cartesian to spherical coordinates
r=sqrt(X1.^2+Y1.^2+Z1.^2);
theta=acos(Z1./r);
theta(281,281,281)=0;    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phi=atan2(Y1,X1);

bas_fun = cell(N_states+1, 1);
bas_fun{1} = x;



for jj = 1:N_states
      
   bas_fun{jj+1} = psi(qn(jj,1),qn(jj,2),qn(jj,3),r,theta,phi);
 
end;

save(strcat(pwd,'/dis_scr/bas_fun.mat'), 'bas_fun', '-v7.3');


a1=qn(1:N_states,4).*40./13.6./2;
Nbands=N_states;
save(strcat(pwd,'/dis_scr/M.mat'), 'Nbands', 'a1')
