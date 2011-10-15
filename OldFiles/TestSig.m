more off;
angstrom = 1/0.5291;
sec = 1/2.418884326505e-17; %# wikipedia: atomic units
cm = 10^8*angstrom;
kappa = 4.3747; %# background dielectric constant
kc = 1/angstrom; %# cutoff wave vector
bands = [-1 +1];

n = 10^12/cm^2; %# electron density
s = bands(2); %# band index
ky = 0.001; %# in units of kF
nus = [1.]'; %# freq. for sigma, in units of EF

global vF = 10^8*cm/sec;
global kF = sqrt(pi*n);
ky *= kF;
global EF = vF*kF;
global mu = EF;

%# good convergence parameters:
%# dphi = 0.04 and du=0.001
phis  = [0:0.04:2*pi]';
dphi = phis(2)-phis(1);
us   = [0.001:0.001:1]';
du   = us(2)-us(1);
qs   = kc*us.^2;
Nq   = length(qs);

wcut  = 2.; %# in units of EF
wsf   = [0.0001:0.1:wcut]'*EF;
dwf   = wsf(2)-wsf(1);
wsc   = [wcut: 2. : 300]'*EF;
dwc   = wsc(2)-wsc(1);
ws    = [wsf; wsc];
dw    = [dwf*ones(size(wsf)); dwc*ones(size(wsc))];
Nfreq= length(ws);
global Vc = 2*pi/kappa./abs(qs);

%# get dielectric matrix for line integral
epsC = zeros(Nq,Nfreq);
Ind  = 1;
for wp = ws';

  IepsC(:,Ind) = 1./Cdiel( qs, I*wp ) - 1.;
  Ind += 1;
endfor;

SigV = [];
vec = [];
for nu = nus';

  nu
  w = nu*EF;
  k = [0 ky];

  Ind = 1;
  for n = 1:Nfreq;
    
    wp=ws(n)
    SigL = 0;
    for sp = bands;
      for phi = phis';
	
	kpq = qs*[cos(phi) sin(phi)] + ones(Nq,1)*k;
	xi  = sp*vF*sqrt(sum(kpq.^2,2)) - mu;
	uFss= us .* F(s,sp,k,kpq);
	SigL += sum( 2* uFss.*(xi-w)./( (xi-w).^2 + wp.^2) .* IepsC(:,Ind) );
      endfor;     
    endfor;

    Ind += 1;
    vec = [vec; SigL];
  endfor;
  
  SigL *= -kc/kappa*dphi*2*du/(2*pi)^2;
  SigV = [SigV; SigL];
endfor;
more on;