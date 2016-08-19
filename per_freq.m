function G=per_freq(wf)

    s=size(wf);
    lc=MyConst.a_Si/MyConst.ab;
    dV=(lc/s(1))^3;

    ub=fftshift(fftn(wf)).*dV;    
    
    s=size(ub);    
    x=(1:s(1))-(s(1)/2+1);
    [X,Y,Z]=meshgrid(x,x,x);
    G=zeros((s(1)*s(2)*s(3)),6);
    
    G(:,1)=X(:).*2.*pi./lc;
    G(:,2)=Y(:).*2.*pi./lc;
    G(:,3)=Z(:).*2.*pi./lc;    
    G(:,4)=real(ub(:));
    G(:,5)=imag(ub(:));
    G(:,6)=(abs(ub(:))).^2;
    
    [~,IX] = sort(squeeze(G(:,6)),'descend');
    
    G=G(IX,:);
