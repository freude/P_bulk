function epsilon=pot_scr_real(r)

au=5.29e-11; 
r=r.*MyConst.ab./au;

A=1.175;
alpha=0.757;
betha=0.322;
gamma=2.044;

epsilon = 1+A*MyConst.eps1.*exp(-alpha.*r)+(1-A)*MyConst.eps1.*exp(-betha.*r)-exp(-gamma.*r);