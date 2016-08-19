function ME=me2(k1,k2,flag)
% if 'me', the function computes the integral 
% over a coupling filtered potential
%
% if 'pot', the function computes the filtered potential itself

if (isequal(k1,k2))&&(strcmp(flag,'mes'))
    
    ME=0; % matrix element is usually needed only for nondiagonal term   
    
else
    
    %-----compute Fourier harmonics of the periodic Bloch fucntion-----
    
%    if (exist(strcat(pwd,'/dis_scr/G.mat'), 'file') == 2)
%        G=load(strcat(pwd,'/dis_scr/G.mat'));
%        G=G.G;
%    else        
        num_cells=1;
        T=70;
        wf=conj(abi_read(num_cells,T,k1)).*(abi_read(num_cells,T,k2));
        G=per_freq(wf);
%        save(strcat(pwd,'dis_scr/G.mat'), 'G');
%    end;
    
    %-----------------------apply filter function---------------------
    
    if isequal(k1,k2)
        num_cells=120;
        coorsys=CoordSys(num_cells,6,'au');
        coorsys.set_origin_cells(num_cells/2+1);
    else
        num_cells=30;
        coorsys=CoordSys(num_cells,16,'au');
        coorsys.set_origin_cells(num_cells/2+1);        
       
    end;

    V = pot_mat(coorsys,k1,k2);
    x = coorsys.x();
         
    %-----------------------apply filter function---------------------
    
   % fprintf('apply filter function...')
    R_tr=7;
    V1sm = smoother_2(x,G,k1,k2,R_tr);
   % fprintf('Done!\n')
    
    if strcmp(flag,'pot')
        
        [~,jj1] = min(abs(x+0.65*R_tr));
        [~,jj2] = min(abs(x-0.65*R_tr));
        
        %     delta=V1sm(jj2,jj2,jj2)-V(jj2,jj2,jj2);
        %     V=V+delta;
        
        V(jj1:jj2,jj1:jj2,jj1:jj2)=-V1sm(jj1:jj2,jj1:jj2,jj1:jj2);
        V1sm=-V;
        
    end;
    
    %----------------------------------------------------------------
    
end;

M2=1;

if (((k1(find(k1))>0)&&(k2(find(k2))<0))||((k2(find(k2))>0)&&(k1(find(k1))<0)))&&(find(k1)==find(k2))
    M1=1;
else
    M1=1;
end;

[X,Y,Z]=meshgrid(x,x,x);
kk=0*k2-0*k1;

if strcmp(flag,'mes')
    ME=2.7*trapz(x,trapz(x,trapz(x,V1sm.*M1.*M2.*exp(1i*(kk(1).*X+kk(2).*Y+kk(3).*Z)),3),2),1);
elseif strcmp(flag,'pot')
    ME{1}=x;
    ME{2}=V1sm;
else
    
end;

