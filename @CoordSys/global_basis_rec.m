function [kx,ky,kz]=global_basis_rec(obj)

a0=2*pi/0.5431;
kx=a0.*[-1; 1; 1];
ky=a0.*[1; -1; 1];
kz=a0.*[1; 1; -1];

X=(j1.*kx+j2.*ky+j3.*kz)/obj.arrays_sizes;

kx=X(1)-obj.coord_limits(1)-obj.coord_stps;
ky=X(2)-obj.coord_limits(1)-obj.coord_stps;
kz=X(3)-obj.coord_limits(1)-obj.coord_stps;
