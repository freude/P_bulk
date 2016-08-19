function integ=el2int_cor_tr1(kk,bas_fun,bas_fun_c,R_tr)

s=length(kk);
dk3=(kk(3)-kk(2))^3;

[ind,~,val] = find(conj(bas_fun).*bas_fun_c);

% for jj=1:length(ind)
%     [j1,j2,j3]=ind2sub([s,s,s],ind(jj));
%     integ=integ+val(jj)/sqrt(kk(j1)^2+kk(j2)^2+kk(j3)^2+1e-10)*dk3;
% end;

[j1,j2,j3]=ind2sub([s,s,s],ind);


if R_tr<0.1
    
    Q=kk(j1).^2+kk(j2).^2+kk(j3).^2+1j*R_tr;        
    V=real(1./Q);
    
else
    
    Q=kk(j1).^2+kk(j2).^2+kk(j3).^2;    
    V=(1-cos(sqrt(Q).*R_tr))./Q;    
    V(Q==0)=0.5*R_tr^2;
    
end;


integ=(2*sum((val.').*V)*dk3/((2*pi)^2));