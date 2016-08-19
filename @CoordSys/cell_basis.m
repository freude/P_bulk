function [x,y,z]=cell_basis(obj,j1,j2,j3)

a0=0.5431/sqrt(2);
x=a0.*[1; 0; 0];
y=a0.*[0; 1; 0];
z=a0.*[0; 0; 1];

X=(j1.*x+j2.*y+j3.*z)/obj.arrays_sizes;

x=X(1)-obj.coord_limits(1)-obj.coord_stps;
y=X(2)-obj.coord_limits(1)-obj.coord_stps;
z=X(3)-obj.coord_limits(1)-obj.coord_stps;