function wf1=read_wf(T,k1)

% associate valley index with abinit wave-function

if k1(1)~=0
    if k1(1)>0    
        indi=5;
    end;
    if k1(1)<0    
        indi=6;
    end;
end;
if k1(2)~=0
    if k1(2)>0    
        indi=2;
    end;
    if k1(2)<0
        indi=4;
    end;
end;
if k1(3)~=0
    if k1(3)>0
        indi=1;
    end;
    if k1(3)<0
        indi=3;
    end;
end;
%if k1(1)~=0
%    if k1(1)>0    
%        indi=6;
%    end;
%    if k1(1)<0    
%        indi=5;
%    end;    
%end;
%if k1(2)~=0
%    if k1(2)>0    
%        indi=4;
%    end;
%    if k1(2)<0
%        indi=2;
%    end;
%end;
%if k1(3)~=0
%    if k1(3)>0
%        indi=3;
%    end;
%    if k1(3)<0
%        indi=1;
%    end;
%end;

ij=sqrt(-1);


% check if the unitcell function is already stored (both real and imaginary parts)

if (exist([pwd, '/dis_scr/wfr_',num2str(indi),'_',num2str(T),'.mat'], 'file') == 2)&&(exist([pwd, '/dis_scr/wfi_',num2str(indi),'_',num2str(T),'.mat'], 'file') == 2)
    
 %   fprintf('(Data is saved)...');
    wfr=load(strcat(pwd,['/dis_scr/wfr_',num2str(indi),'_',num2str(T),'.mat']));
    wfi=load(strcat(pwd,['/dis_scr/wfi_',num2str(indi),'_',num2str(T),'.mat']));
    wf1=(wfr.wfr+ij.*wfi.wfi);
    
else
    
%    fprintf('(Data is not saved)...');
    
    wf=ncread(strcat(pwd,'/dis_scr/tbase3_xo_DS2_AE_WFK-etsf.nc'),...
        'real_space_wavefunctions',[1 1 1 1 1 5 1 1],[Inf Inf Inf Inf 1 1 Inf 1],[1 1 1 1 1 1 1 1]);
    
    sc=1;
    
    wfr=squeeze(wf(1,1:sc:end,1:sc:end,1:sc:end,1,1,indi));
    wfi=squeeze(wf(2,1:sc:end,1:sc:end,1:sc:end,1,1,indi));
    
    wfr=transform_to_uc1(wfr,T);
    wfi=transform_to_uc1(wfi,T);
    
    save(strcat(pwd,['/dis_scr/wfr_',num2str(indi),'_',num2str(T),'.mat']), 'wfr');
    save(strcat(pwd,['/dis_scr/wfi_',num2str(indi),'_',num2str(T),'.mat']), 'wfi');
    
    wf1=wfr+ij.*wfi;
    
end;

x=0:(MyConst.a_Si/MyConst.ab)/T:(MyConst.a_Si/MyConst.ab-(MyConst.a_Si/MyConst.ab)/T);
ME=trapz(x,trapz(x,trapz(x,abs(wf1).^2,3),2),1);
wf1=(1/sqrt(ME))*wf1;
