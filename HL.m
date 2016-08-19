clear all;

matObj = matfile(strcat(pwd,'/dis_scr/!bas_fun.mat'),'writable',false);
A=matObj.bas_fun(2,1);
A=A{1};


%bs.bas_fun = bs.bas_fun(1:2);

level_par=0.0001;                        % level of reduced noise
Rc=5.0;                               % regularization, imaginary number/radius of truncation
ap=150;

x = matObj.bas_fun(1,1);
x = x{1};

bas_fun_ms = cell(2, 1);
%ap=500;
bas_fun_ms{1} = fftshift(fftfreq(length(x)+ap, x(3)-x(2)));
    k_tr=200;
    %k_tr=max(bas_fun_ms{1});
    [~,c2]=min(abs(bas_fun_ms{1}-k_tr));
    [~,c1]=min(abs(bas_fun_ms{1}+k_tr));
    bas_fun_ms{1}=bas_fun_ms{1}(c1:c2);

R=7.*(7:17);                                 % number of elements in padding arrays
integ=zeros(length(R),1);

for jjj=1:length(R)
    jjj
    A1=padarray(A,R(jjj),0,'post');
    A2=padarray(A,R(jjj),0,'pre');
    
    Ff=((ifftshift(fftn(fftshift(padarray(A1.*A2,[fix(ap/2) fix(ap/2) fix(ap/2)])))))).*((x(3)-x(2))^3);
    Ff=Ff(c1:c2,c1:c2,c1:c2);
    lev=max(abs(Ff(:)));
    Ff(abs(Ff)<level_par*lev)=0;
    Ff1 = sparse(Ff(:));
    bas_fun_ms{2}= Ff1;
        
    integ(jjj)=el2int_cor_tr1(bas_fun_ms{1},bas_fun_ms{2},bas_fun_ms{2},Rc);    
end;


