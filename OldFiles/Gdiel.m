more off;
x= [0.01:0.01:4]';

angstrom = 1/0.5291;
sec = 1/2.418884326505e-17; %# wikipedia: atomic units
cm = 10^8*angstrom;
n = 10^13/cm^2;
vF = 10^8*cm/sec;

D0 = sqrt(4*n/pi)/vF;
kF = sqrt(pi*n);
EF = vF*kF;
kappa = 4;
q = x*kF;
Vc = 2*pi./abs(q)/kappa;
vmax = 0.3;
vmin = -15.;
g = 4;
eta = 0.001;

mat = [];
for nu  =  0:0.01:4;

  #nu=0.3;
  w = EF*nu
  
  F  = g/16/pi* vF^2*q.^2 ./sqrt(w.^2 - vF^2*q.^2);
  
  dP = -g*EF/2/pi/vF^2 + F/vF^2 .* ( G((nu+2)./x) -...
	(((2-nu)./x-1)>0).* (G((2-nu)./x) - I*pi) - (((nu-2)./x+1)>0).*G((nu-2)./x) );

  P0 = -I*pi*F/vF^2; 
  Pi  = P0 + dP;
  
  epsilon = 1-Vc .* Pi + I*eta;
  
  v = imag(1./epsilon)/(EF/vF^2);
  v(v>vmax) = vmax; 
  v(v<vmin) = vmin;
  mat = [mat; nu*ones(size(x)) x v];

endfor;
more on;