function wf = abi_read(fac,T,valley)
% the function reads the periodic functions computed by ABINIT
% for fac number of unit cells

wf1 = read_wf(T,valley); % read the wave function for a single unit cells

% if valley(find(valley))<0
%     wf1=-wf1;
% end;

% compose a fac number of cells
wf = zeros(fac*T,fac*T,fac*T);

for j1=1:fac
    for j2=1:fac
        for j3=1:fac
            wf(((j1-1)*T+1):(j1*T),((j2-1)*T+1):(j2*T),((j3-1)*T+1):(j3*T)) = wf1;
        end;
    end;
end;
