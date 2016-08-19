clear all;

k0=MyConst.k0*MyConst.ab;

kk=k0*[1  0  0;
    -1  0  0;
    0  1  0;
    0 -1  0;
    0  0  1;
    0  0 -1];

pot_for_ff(kk(1,:), kk(1,:),'1');
pot_for_ff(kk(2,:), kk(2,:),'2');
pot_for_ff(kk(3,:), kk(3,:),'3');
