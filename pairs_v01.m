clear all;
bs=load(strcat(pwd,'/dis_scr/bas_fun.mat'));
%bs.bas_fun = bs.bas_fun(1:2);


level_par=0.001;                        % level of reduced noise
ap=350;                                 % number of elements in padding arrays
Rc=3.7;                               % regularization, imaginary number/radius of truncation
save='yes';

[integ,exch] = two_el_int(bs,level_par,ap,Rc,save);