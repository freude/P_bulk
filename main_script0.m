clear all;

pth=pwd;
k0=MyConst.k0*MyConst.ab;
kk=k0*[1  0  0;
    -1  0  0;
    0  1  0;
    0 -1  0;
    0  0  1;
    0  0 -1];

bands=[1 3 7 10];
Nbands=length(bands);


%[amp,~]=read_env1();

amp = read_env1(bands,...
    strcat(pth,'/dis_scr/'),...
    kk(1,:),0,[0 0 0]);

amp(abs(amp)<0.01*max(abs(amp(:))))=0;



amp=[amp; amp; amp; amp; amp; amp];

Energy=[-0.87 -0.24];
Energy=[-0.87 -0.24 -0.1 -0.01];
%Energy=[-0.84 -0.29 -0.23];

% ---------------------------------------------------------------------
%
%if (exist(strcat(pwd,'/dis_scr/ME.mat'), 'file') == 2)
%    M=load(strcat(pwd,'/dis_scr/ME.mat'));
%    ME_a=M.M(1);   
%    ME_b=M.M(2);
%else
%    ME_a=me2(kk(1,:),kk(2,:),'mes');
%    ME_b=me2(kk(1,:),kk(3,:),'mes');
%    M=[real(ME_a) real(ME_b)];
%    save(strcat(pwd,'/dis_scr/ME.mat'), 'M');
%end;
%
    ME_a=me2(kk(1,:),kk(2,:),'mes')
    ME_b=me2(kk(1,:),kk(3,:),'mes')
   % ME_b=me2(kk(1,:),kk(4,:),'mes')
   % ME_b=me2(kk(1,:),kk(5,:),'mes')
   % ME_b=me2(kk(1,:),kk(6,:),'mes')

   % ME_a=me2(kk(2,:),kk(1,:),'mes')
   % ME_b=me2(kk(2,:),kk(3,:),'mes')
   % ME_b=me2(kk(2,:),kk(4,:),'mes')
   % ME_b=me2(kk(2,:),kk(5,:),'mes')
   % ME_b=me2(kk(2,:),kk(6,:),'mes')

   % ME_b=me2(kk(3,:),kk(1,:),'mes')
   % ME_b=me2(kk(3,:),kk(2,:),'mes')
   % ME_a=me2(kk(3,:),kk(4,:),'mes')
   % ME_b=me2(kk(3,:),kk(5,:),'mes')
   % ME_b=me2(kk(3,:),kk(6,:),'mes')

   % ME_b=me2(kk(4,:),kk(1,:),'mes')
   % ME_b=me2(kk(4,:),kk(2,:),'mes')
   % ME_a=me2(kk(4,:),kk(3,:),'mes')
   % ME_b=me2(kk(4,:),kk(5,:),'mes')
   % ME_b=me2(kk(4,:),kk(6,:),'mes')

   % ME_b=me2(kk(5,:),kk(1,:),'mes')
   % ME_b=me2(kk(5,:),kk(2,:),'mes')
   % ME_b=me2(kk(5,:),kk(3,:),'mes')
   % ME_b=me2(kk(5,:),kk(4,:),'mes')
   % ME_a=me2(kk(5,:),kk(6,:),'mes')

M=zeros(6*Nbands,6*Nbands);

% ---------------------------------------------------------------------

MEs=zeros(6,6);
for j1=1:6
    for j2=1:6
        if ((j1~=j2)&&(j1<j2))
            if sum(abs(kk(j1,:)+kk(j2,:)))==0
                MEs(j1,j2)=ME_a;
            else
                MEs(j1,j2)=ME_b;
            end;
        else
            MEs(j1,j2)=0;
        end;
    end;
end;

% ---------------------------------------------------------------------

for j1=1:6
    for j2=1:6
        for jj1=1:Nbands
            for jj2=1:Nbands
                if (j1==j2)&&(jj1~=jj2)
                    M(jj1+Nbands*(j1-1),jj2+Nbands*(j2-1))=0;
                elseif (j1==j2)&&(jj1==jj2)
                    M(jj1+Nbands*(j1-1),jj2+Nbands*(j2-1))=Energy(jj1);
                else
                    M(jj1+Nbands*(j1-1),jj2+Nbands*(j2-1))=...
                        MEs(j1,j2)*amp(j1,jj1)*amp(j2,jj2);
                end;
                
            end;
        end;
    end;
end;

for j1=1:6*Nbands
    for j2=1:6*Nbands
        if ((j1~=j2)&&(j1>j2))
            M(j1,j2)=conj(M(j2,j1));
        end;
    end;
end;

% ---------------------------------------------------------------------
% -------------Overlap matrix <k2,i|k1,j> for envelopes----------------
% ---------------------------------------------------------------------

%kin=[Energy(1) 0;...
%    0 Energy(2)];
%
%kin=[Energy(1) 0 0 0;...
%    0 Energy(2) 0 0;...
%    0 0 Energy(3) 0;...
%    0 0 0 Energy(4)];
%% -------------------------------------------
%
%D1=M(1,3);
%D2=M(1,4);
%D3=M(2,3);
%D4=M(2,4);
%
%D_1=[D1 D2; D3 D4];
%
%D1=M(1,5);
%D2=M(1,6);
%D3=M(2,5);
%D4=M(2,6);
%
%D_2=[D1 D2; D3 D4];
%
%D_1=[D1 D2; D3 D4];
%D_2=[D1 D2; D3 D4];
%
%D_1i=conj(D_1)';
%D_2i=conj(D_2)';
%
%Ma=[kin D_1   D_2   D_2   D_2   D_2;...
%    D_1i  kin D_2   D_2   D_2   D_2;...
%    D_2i  D_2i  kin D_1   D_2   D_2;...
%    D_2i  D_2i  D_1i  kin D_2   D_2;...
%    D_2i  D_2i  D_2i  D_2i  kin D_1;...
%    D_2i  D_2i  D_2i  D_2i  D_1i  kin];
%
%% Ma=[Energy(1) D1   D2   D2   D2   D2;...
%%     D1'  Energy(1) D2   D2   D2   D2;...
%%     D2'  D2'  Energy(1) D1   D2   D2;...
%%     D2'  D2'  D1'  Energy(1) D2   D2;...
%%     D2'  D2'  D2'  D2'  Energy(1) D1;...
%%     D2'  D2'  D2'  D2'  D1'  Energy(1)];

% ---------------------------------------------------------------------


[EigVec,a] = eig(M);
a = diag(real(a).*40);
[a1,ind] = sort(a);

% EigVec1=XX*EigVec;
% EigVec1=EigVec1(:,ind);


EigVec1=EigVec(:,ind);


% ---------------------------Save----------------------------

%save(strcat(pth,'/dis_scr/M.mat'), 'EigVec1', 'Nbands', 'kk', 'bands', 'a1')

% dlmwrite('/data/users/mklymenko/abinitio_software/abi_files/tpaw/M.mat', EigVec1);

%construct_wf3
% ---------------------------Add p-orbitals----------------------------
%main_script_p3

% p1 = importdata(['/data/users/mklymenko/abinitio_software/abi_files/tpaw/', 'dis/v0/ff_0.dat']);
% j1=0;
%
% for j=1:(length(squeeze(p1(:,1))))
%     if (p1(j,1)==111)&&(p1(j,2)==111)&&(p1(j,3)==111)
%         j1=j1+1;
%         a(j1)=j;
%     end;
% end;
%
% n1=1;
%
% F1=squeeze(p1(a(n1)+1:a(n1+1)-1,4));
% % ------------------------------------------------------------------
%
% X=squeeze(p1(a(n1)+1:a(n1+1)-1,1));
% Y=squeeze(p1(a(n1)+1:a(n1+1)-1,2));
% Z=squeeze(p1(a(n1)+1:a(n1+1)-1,3));
%
% Fq = TriScatteredInterp(X,Y,Z,F1);
% x=-1.5:0.005:1.5;
% [X,Y,Z]=meshgrid(x,x,x);
% M=Fq(X,Y,Z);

% for j1=1:6
%     for j2=1:6
%         if ((j1~=j2)&&(j1>j2))
%             M(j1,j2)=conj(M(j2,j1));
%         end;
%     end;
% end;


%VO_aniso();

%MM7=ff_data1('/data/users/mklymenko/abinitio_software/abi_files/tpaw/',1,1,[0 0 1],[1 0 0],'_r_xy','_i_xy');
%MM4=ff_data1('/data/users/mklymenko/abinitio_software/abi_files/tpaw/',1,1,[0 1 0],[0 -1 0],potxmx);


