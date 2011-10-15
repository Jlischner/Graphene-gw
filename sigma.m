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
nus = [0.2:0.2:3.]'; %# freq. for sigma, in units of EF

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
qs   = kc*us.^2;
Nq   = length(qs);
global Vc = 2*pi/kappa./abs(qs);

wcut  = 2.; %# in units of EF                           
wsf   = [0.0001:0.1:wcut]'*EF;
wsc   = [wcut: 2. : 300]'*EF;
ws    = [wsf; wsc];
Nfreq= length(ws);

%# get dielectric matrix for line integral
epsC = zeros(Nq,Nfreq);
for n = 1:Nfreq;
  IepsC(:,n) = 1./Cdiel( qs, I*ws(n) ) - 1.;
endfor;

SigXw = [];
SigPw = [];
SigLw = [];

for nu = nus';

  nu
  w = nu*EF;
  k = [0 ky];

  SigX = 0.;
  SigP = 0.;
  SigL = 0.;

  for sp = bands;
    for phi = phis';
    
      kpq = qs*[cos(phi) sin(phi)] + ones(Nq,1)*k;
      xi  = sp*vF*sqrt(sum(kpq.^2,2)) - mu;
      uFss= us .* F(s,sp,k,kpq);
      SigX += trapz(us, (xi<0) .* uFss ) ;
      
      Ieps = 1./diel( qs, xi-w ) - 1.;
      SigP += trapz(us, ( ((w-xi)>0) - (xi<0) ) .* uFss .* Ieps );

      SigLv = zeros(Nfreq,1);
      for n = 1:Nfreq;
	SigLv(n) += trapz(us, uFss.*(xi-w)./( (xi-w).^2 + ws(n).^2) .* IepsC(:,n) );
      endfor;
      SigL += trapz(ws, SigLv);
      
    endfor;%# end phis loop
  endfor;%# end bands loop

  SigX *= -kc/kappa*dphi*2/(2*pi);
  SigP *=  kc/kappa*dphi*2/(2*pi);
  SigL *= -kc/kappa*dphi*4/(2*pi)^2;
  
  SigXw = [SigXw; SigX];
  SigPw = [SigPw; SigP];
  SigLw = [SigLw; SigL];

endfor;%# end nus loop
more on;