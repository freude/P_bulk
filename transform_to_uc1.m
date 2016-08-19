function wf=transform_to_uc1(wf1,L)
% the function transforms a primitive cell wave-function red from ABINIT
% into a unit-cell wave function

a0=0.5431;

s=size(wf1);
T=s(1);

mu=3; % number of primitive cells considered to form unit cell

xx=(0:mu/(mu*T):(mu-mu/(mu*T)))-1;

[x1,y1,z1] = meshgrid(xx,xx,xx);

for j1=1:mu
    for j2=1:mu
        for j3=1:mu
            wf(((j1-1)*T+1):(j1*T),((j2-1)*T+1):(j2*T),((j3-1)*T+1):(j3*T))=wf1;
        end;
    end;
end;

% coordinates transformation

x=(y1+z1)*0.5*a0;
y=(x1+z1)*0.5*a0;
z=(x1+y1)*0.5*a0;

F = TriScatteredInterp(x(:), y(:), z(:), wf(:));
[x1,y1,z1] = meshgrid(0:a0/L:(a0-a0/L),0:a0/L:(a0-a0/L),0:a0/L:(a0-a0/L));

% the origin is placed between nodes
% x1=x1-0.5*a0/L;
% y1=y1-0.5*a0/L;
% z1=z1-0.5*a0/L;

wf = F(x1,y1,z1);




