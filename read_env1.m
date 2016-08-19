function [amp,M]=read_env1(bands,path,k1,indi,coord)

[~,pth]=k2flag(k1);
s = size(coord);
Nimp=s(1);

Nbands=length(bands);

p1 = importdata([path pth 'ff_' num2str(indi) '.dat']);
j1=0;

for j=1:(length(squeeze(p1(:,1))))
    if (p1(j,1)==111)&&(p1(j,2)==111)&&(p1(j,3)==111)
        j1=j1+1;
        a(j1)=j;
    end;
end;

amp=zeros(Nimp,Nbands);

for jjj=1:Nbands
    
    F1=squeeze(p1(a(bands(jjj))+1:a(bands(jjj)+1)-1,4));
    X=squeeze(p1(a(bands(jjj))+1:a(bands(jjj)+1)-1,1));
    Y=squeeze(p1(a(bands(jjj))+1:a(bands(jjj)+1)-1,2));
    Z=squeeze(p1(a(bands(jjj))+1:a(bands(jjj)+1)-1,3));
    
    Fq = TriScatteredInterp(X,Y,Z,F1);
    x=-3:0.02:3;
    
    if nargout > 1        
        [X,Y,Z]=meshgrid(x,x,x);
        M=Fq(X,Y,Z);
        M(isnan(M))=0;        
        ME=trapz(x,trapz(x,trapz(x,abs(M).^2,3),2),1); 
        M=M./sqrt(ME);
        amp(:,jjj)=abs(Fq(coord(:,1),coord(:,2),coord(:,3)))/sqrt(ME);
    else
        amp(:,jjj)=(Fq(coord(:,1),coord(:,2),coord(:,3)));
    end;            
end;
