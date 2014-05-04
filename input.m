%# input file
D = 0.;
kcoption = 2;
eta = 0.01;
kappa = 1; %# background dielectric
s = -1; %# band index: either +1 or -1
EFeV = 0.1; %# fermi energy in eV
nus = [-3.5 :0.025: -0.5]'; %# freq. grid for sigma in units of EF
ky = 0.3; %# kpoint in units of kF
