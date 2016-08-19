function [bi]=br_zone_valley(x2, y2, z2, valley, shift)

%-----------------------------------------------------
%-----------------***************---------------------
%---------------*-----------------*-------------------
%-------------*---------------------*-----------------
%-----------*-------------------------*---------------
%---------*-----------------------------*-------------
%---------*-----------------------------*-------------
%---------*-----------------------------*-------------
%---------*--------------0(0,0,0)-*--*--*X(1,0,0)-----
%---------*-------------------*---------*-------------
%---------*----------------*-------*----*-------------
%---------*-------------------*---------*W(1,0.5,0)---
%-----------*--------------------*----*---------------
%-------------*---------------------*K(0.75,0.75,0)---
%---------------*-----------------*-------------------
%-----------------***************---------------------
%-----------------------------------------------------

%%

%     tt=-pi/4;    
% 
%     % Rz=[cos(tt) -sin(tt)  0;
%     %     sin(tt)  cos(tt)  0;
%     %       0       0       1];
%     
%     y1=y2.*cos(tt)-z2.*sin(tt);
%     z1=y2.*sin(tt)+z2.*cos(tt);    
%     x1=x2;   
% 
%     tt=-pi/4; 
%     
%     x2=x1.*cos(tt)-z1.*sin(tt);
%     z2=x1.*sin(tt)+z1.*cos(tt);    
%     y2=y1;   
% 
%     
%     x1=x2.*cos(tt)-y2.*sin(tt);
%     y1=x2.*sin(tt)+y2.*cos(tt);    
%     z1=z2;   
  
x1=x2;
y1=y2;
z1=z2;

%%
if strcmp(valley,'x')
    tt=pi;    
    % Rx=[1     0      0;
    %     0 cos(tt) -sin(tt);
    %     0 sin(tt) cos(tt)];
    %
    % Ry=[cos(tt)  0   -sin(tt);
    %       0      1     0;
    %     sin(tt)  0   cos(tt)];
    %
    % Rz=[cos(tt) -sin(tt)  0;
    %     sin(tt)  cos(tt)  0;
    %       0       0       1];
    
    x=x1.*cos(tt)-y1.*sin(tt);
    y=x1.*sin(tt)+y1.*cos(tt);    
    z=z1;
    x=x-shift*MyConst.k0;       
end;

if strcmp(valley,'y')
    tt=pi/2;    
    % Rx=[1     0      0;
    %     0 cos(tt) -sin(tt);
    %     0 sin(tt) cos(tt)];
    %
    % Ry=[cos(tt)  0   -sin(tt);
    %       0      1     0;
    %     sin(tt)  0   cos(tt)];
    %
    % Rz=[cos(tt) -sin(tt)  0;
    %     sin(tt)  cos(tt)  0;
    %       0       0       1];
    
    x=x1.*cos(tt)-y1.*sin(tt);
    y=x1.*sin(tt)+y1.*cos(tt);  
    z=z1;  
    x=x-shift*MyConst.k0;       
end;

if strcmp(valley,'-y')
    tt=-pi/2;    
    % Rx=[1     0      0;
    %     0 cos(tt) -sin(tt);
    %     0 sin(tt) cos(tt)];
    %
    % Ry=[cos(tt)  0   -sin(tt);
    %       0      1     0;
    %     sin(tt)  0   cos(tt)];
    %
    % Rz=[cos(tt) -sin(tt)  0;
    %     sin(tt)  cos(tt)  0;
    %       0       0       1];
    
    x=x1.*cos(tt)-y1.*sin(tt);
    y=x1.*sin(tt)+y1.*cos(tt);
    z=z1;    
    x=x-shift*MyConst.k0;       
end;

if strcmp(valley,'z')
    tt=pi/2;
    % Rx=[1     0      0;
    %     0 cos(tt) -sin(tt);
    %     0 sin(tt) cos(tt)];
    %
    % Ry=[cos(tt)  0   -sin(tt);
    %       0      1     0;
    %     sin(tt)  0   cos(tt)];
    %
    % Rz=[cos(tt) -sin(tt)  0;
    %     sin(tt)  cos(tt)  0;
    %       0       0       1];
    
    x=x1.*cos(tt)-z1.*sin(tt);
    z=x1.*sin(tt)+z1.*cos(tt);
    y=y1;
    x=x-shift*MyConst.k0;       
end;

if strcmp(valley,'-z')
    tt=-pi/2;
    % Rx=[1     0      0;
    %     0 cos(tt) -sin(tt);
    %     0 sin(tt) cos(tt)];
    %
    % Ry=[cos(tt)  0   -sin(tt);
    %       0      1     0;
    %     sin(tt)  0   cos(tt)];
    %
    % Rz=[cos(tt) -sin(tt)  0;
    %     sin(tt)  cos(tt)  0;
    %       0       0       1];
    
    x=x1.*cos(tt)-z1.*sin(tt);
    z=x1.*sin(tt)+z1.*cos(tt);
    y=y1;    
    x=x-shift*MyConst.k0;       
end;

if strcmp(valley,'-x')        
    x=x1;
    y=y1;
    z=z1;
    x=x-shift*MyConst.k0;    
end;
%%



rho=sqrt(x.^2+y.^2+z.^2);
phi=atan(y./x);
theta=acos(z./rho);

R1=2*pi/0.5431e-9;
R2=0.75*2*pi/0.5431e-9;


s1 = @(phi,theta,rho) (rho.*cos(phi).*sin(theta)-R1)<0;
s2 = @(phi,theta,rho) (rho.*cos(phi).*sin(theta)/(2*R2)+rho.*sin(phi).*sin(theta)/(2*R2)+rho.*cos(theta)/(2*R2)-1)<0;
s3 = @(phi,theta,rho) (rho.*cos(phi).*sin(theta)/(2*R2)-rho.*sin(phi).*sin(theta)/(2*R2)+rho.*cos(theta)/(2*R2)-1)<0;
s4 = @(phi,theta,rho) (rho.*cos(phi).*sin(theta)/(2*R2)+rho.*sin(phi).*sin(theta)/(2*R2)-rho.*cos(theta)/(2*R2)-1)<0;
s5 = @(phi,theta,rho) (rho.*cos(phi).*sin(theta)/(2*R2)-rho.*sin(phi).*sin(theta)/(2*R2)-rho.*cos(theta)/(2*R2)-1)<0;
s6 = @(phi,theta,rho) (cos(phi).*sin(theta)+cos(theta))>0;
s7 = @(phi,theta,rho) (cos(phi).*sin(theta)-cos(theta))>0;
s8 = @(phi,theta,rho) (abs(phi))<pi/4;
s9 = @(phi,theta,rho) (x)<0;

bi=s1(phi,theta,rho).*s2(phi,theta,rho).*s3(phi,theta,rho).*s4(phi,theta,rho).*s5(phi,theta,rho).*s6(phi,theta,rho).*s7(phi,theta,rho).*s8(phi,theta,rho).*s9(phi,theta,rho);

% if (s1(phi,theta,rho)<0)&&...
%         (s2(phi,theta,rho)<0)&&...
%         (s3(phi,theta,rho)<0)&&...
%         (s4(phi,theta,rho)<0)&&...
%         (s5(phi,theta,rho)<0)&&...
%         (s6(phi,theta,rho)>0)&&...
%         (s7(phi,theta,rho)>0)&&...
%         (abs(phi)<pi/4)&&...
%         (x>0)    
%     
%     bi=1;
% else
%     bi=0;
% end;
