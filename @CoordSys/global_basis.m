function [x,y,z]=global_basis(obj,j1,j2,j3)

a0=0.5431;
x=a0/2.*[0; 1; 1];
y=a0/2.*[1; 0; 1];
z=a0/2.*[1; 1; 0];

X=(j1.*x+j2.*y+j3.*z)/obj.arrays_sizes;

x=X(1)-obj.coord_limits(1)-obj.coord_stps;
y=X(2)-obj.coord_limits(1)-obj.coord_stps;
z=X(3)-obj.coord_limits(1)-obj.coord_stps;
