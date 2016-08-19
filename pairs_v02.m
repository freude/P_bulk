clear all;
bs=load(strcat(pwd,'/dis_scr/bas_fun.mat'));
%bs.bas_fun = bs.bas_fun(1:2);


level_par=0.001;                        % level of reduced noise
ap=550;                                 % number of elements in padding arrays
Rc=7.7;                               % regularization, imaginary number/radius of truncation
sav_e='yes';

x = bs.bas_fun{1};
N = length(bs.bas_fun) - 1;
%N=5;

bas_fun_ms = cell(N+1, 1);
%ap=500;
bas_fun_ms{1} = fftshift(fftfreq(length(x)+ap, x(3)-x(2)));
ind_len = ind3mat(N,N,N,'sym');

k_tr=120;
k_tr=20;
%k_tr=max(bas_fun_ms{1});
[~,c2]=min(abs(bas_fun_ms{1}-k_tr));
[~,c1]=min(abs(bas_fun_ms{1}+k_tr));
bas_fun_ms{1}=bas_fun_ms{1}(c1:c2);

ind_to_del=1;

for jj=1:ind_len
    jj    
    [i,j]=mat3ind(jj,N);
    if i>ind_to_del
        bs.bas_fun{i}=[];
        ind_to_del=ind_to_del+1;
    end;
    
    Ff=((ifftshift(fftn(fftshift(padarray(bs.bas_fun{i+1}.*bs.bas_fun{j+1},[fix(ap/2) fix(ap/2) fix(ap/2)])))))).*((x(3)-x(2))^3);
    Ff=Ff(c1:c2,c1:c2,c1:c2);
    lev=max(abs(Ff(:)));
    Ff(abs(Ff)<level_par*lev)=0;
    Ff1 = sparse(Ff(:));
    bas_fun_ms{jj+1} = Ff1;
    
end;

clear bs

% [X1,Y1,Z1]=ndgrid(bas_fun_ms{1}, bas_fun_ms{1}, bas_fun_ms{1});
% bas_fun_ms{1} = -k_tr:(2*k_tr/300):k_tr;
% [X2,Y2,Z2]=ndgrid(bas_fun_ms{1}, bas_fun_ms{1}, bas_fun_ms{1});
% 
% for jj=1:ind_len
%     F = griddedInterpolant(X1,Y1,Z1,bas_fun_ms{jj+1}, 'cubic');
%     bas_fun_ms{jj+1}=F(X2,Y2,Z2);    
% end;


EigVec1=load([pwd,'/dis_scr/M.mat'], 'a1');
integ=zeros(ind_len,ind_len);

for j1=1:ind_len
    for j2=1:ind_len
        if j2>=j1
            %integ(j1,j2)=el2int(bas_fun_ms{1},bas_fun_ms{j1+1},bas_fun_ms{j2+1});
            integ(j1,j2)=el2int_cor_tr1(bas_fun_ms{1},bas_fun_ms{j1+1},bas_fun_ms{j2+1},Rc);
        else
            integ(j1,j2)=integ(j2,j1);
        end;
    end;
end;

exch=zeros(N,N);

for j1=1:N
    for j2=1:N
        exch(j1,j2)=0;
        for j3=1:N
            jj1=ind3mat(j1,j3,N,'sym');
            jj2=ind3mat(j3,j2,N,'sym');
            %exch(j1,j2)=exch(j1,j2)-0.5*el2int(bas_fun_ms{1},bas_fun_ms{jj1+1},bas_fun_ms{jj2+1});
            exch(j1,j2)=exch(j1,j2)-0.5*integ(jj1,jj2);
        end;
    end;
end;

% exch(find(eye(N)))=exch(find(eye(N)))+EigVec1.a1(1:N)/40;

if strcmp(sav_e, 'yes')
    save([pwd,'/dis_scr/bas_fun_ms.mat'], 'bas_fun_ms', '-v7.3');
    dlmwrite([pwd,'/dis_scr/E.dat'], EigVec1.a1/40, 'delimiter', ' ');
    dlmwrite([pwd,'/dis_scr/int2e.dat'], real(integ), 'delimiter', ' ');
    dlmwrite([pwd,'/dis_scr/exch.dat'], real(exch), 'delimiter', ' ');
end;
