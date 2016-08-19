function pot_for_ff(k1,k2,file_ind)

%-----------------------apply filter function---------------------

pth=pwd

M=me2(k1,k2,'pot');

x=M{1};
V1sm=M{2};

clear M

%----------------------------------------------------------------

fprintf('read the mesh...')
b=gridy();
fprintf('Done!\n')

fprintf('interpolate the potential...')
%V2=-abs(real(V1sm(1:2:end,1:2:end,1:2:end)));

V2i=0;

if isequal(k1,k2)
    [X,Y,Z]=meshgrid(x(1:2:end),x(1:2:end),x(1:2:end));
    V2=(real(V1sm(1:2:end,1:2:end,1:2:end)));
else
    [X,Y,Z]=meshgrid(x,x,x);
    V2=real(V1sm);
    V2i=imag(V1sm);
end;

if length(V2i)==1
    

    if (exist(strcat(pwd,'/dis_scr/F',file_ind,'.mat'), 'file') == 2)
        F=load(strcat(pwd,'/dis_scr/F',file_ind,'.mat'));
        F=F.F;
    else
        P = [2 1 3];
        X = permute(X, P);
        Y = permute(Y, P);
        Z = permute(Z, P);
        V2 = permute(V2, P);
        F = griddedInterpolant(X, Y, Z, V2);  
       save(strcat(pwd,'/dis_scr/F',file_ind,'.mat'), 'F');      
    end;

    fprintf('Done!\n')    
    fprintf('compute the potetnial on the mesh...')
    pot=F(b(:,1),b(:,2),b(:,3));
    pot(isnan(pot))=0;
    fprintf('Done!\n')
    

        fprintf('save the mesh...')
        dlmwrite(strcat(pth, '/dis_scr/pot', file_ind, '.txt'), pot);
        dlmwrite(strcat(pth, '/dis_scr/mesh.dat'), b);
        fprintf('Done!\n')

else

    if ((exist(strcat(pwd,'/dis_scr/F',file_ind,'i.mat'), 'file') == 2)&&(exist(strcat(pwd,'/dis_scr/F',file_ind,'i.mat'), 'file') == 2))
        F=load(strcat(pwd,'/dis_scr/F',file_ind,'r.mat'));
        F=F.F;
        F1=load(strcat(pwd,'/dis_scr/F',file_ind,'i.mat'));
        F1=F1.F1;
    else
        P = [2 1 3];
        X = permute(X, P);
        Y = permute(Y, P);
        Z = permute(Z, P);
        V2 = permute(V2, P);
        V2i = permute(V2i, P);
        F = griddedInterpolant(X, Y, Z, V2);
        F1 = griddedInterpolant(X, Y, Z, V2i);
        save(strcat(pwd,'/dis_scr/F',file_ind,'r.mat'), 'F');     
        save(strcat(pwd,'/dis_scr/F',file_ind,'i.mat'), 'F1');     
    end;


    fprintf('Done!\n')
    

    fprintf('compute the potetnial on the mesh...')
    pot=F(b(:,1),b(:,2),b(:,3));
    pot_i=F1(b(:,1),b(:,2),b(:,3));
    pot(isnan(pot))=0;
    pot_i(isnan(pot_i))=0;
    fprintf('Done!\n')    
    
    ij=sqrt(-1);   
    

        fprintf('save the mesh...')
        dlmwrite(strcat(pth, '/dis_scr/pot_r', file_ind, '.txt'), pot);
        dlmwrite(strcat(pth, '/dis_scr/pot_i', file_ind, '.txt'), pot_i);        
        dlmwrite(strcat(pth, '/dis_scr/mesh.dat'), b);
        fprintf('Done!\n')
  
end;
