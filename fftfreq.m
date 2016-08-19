function f = fftfreq(N,dt)

% f = [0, 1, ...,   n/2-1,     -n/2, ..., -1] / (d*n)   if n is even
% f = [0, 1, ..., (n-1)/2, -(n-1)/2, ..., -1] / (d*n)   if n is odd

% if  mod(N,2)==0
%     a=linspace(0,N/2-1,N/2);
%     b=linspace(1,N/2,N/2);
%     f=[-b(end:-1:1) a]./(dt*N);    
% else
%     a=linspace(0,(N-1)/2,(N-1)/2+1);
%     b=linspace(1,(N-1)/2,(N-1)/2);
%     f=[-b(end:-1:1) a]./(dt*N);    
% end;

if  mod(N,2)==0
    a=linspace(0,N/2-1,N/2);
    b=linspace(1,N/2,N/2);
    f=2*pi*[a -b(end:-1:1)]./(dt*N);    
else
    a=linspace(0,(N-1)/2,(N-1)/2+1);
    b=linspace(1,(N-1)/2,(N-1)/2);
    f=2*pi*[a -b(end:-1:1)]./(dt*N);    
end;