more off;
angstrom = 1/0.5291;
sec = 1/2.418884326505e-17; %# wikipedia: atomic units
cm = 10^8*angstrom;
kc = 1/(2.46*angstrom); %# cutoff wave vector
bands = [-1 +1];

%# read input
input;

global e0 = 1.0
global vF = 10^8*cm/sec;
global EF = EFeV/27.21;
global kF = EF/vF;
global mu = EF;
k  = [0 ky*kF];
n   = kF^2/pi*cm^2; %# electron density 

%# good convergence parameters:
%# dphi = 0.04 and du=0.001
phis  = [0:0.01:2*pi]';
dphi = phis(2)-phis(1);
us   = [0.001:0.001:1]';
qs   = kc*us.^2;
Nq   = length(qs);
global Vc = 2*pi/kappa./abs(qs);

wts = [0.0:0.2:15]'*sqrt(EF);
ws = wts.^2;
Nfreq= length(wts);

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
#  k = [0 ky];

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
      SigL += trapz(wts,2*wts.* SigLv);
      
    endfor;%# end phis loop
  endfor;%# end bands loop

  SigX *= -kc/kappa*dphi*2/(2*pi);
  SigP *=  kc/kappa*dphi*2/(2*pi);
  SigL *= -kc/kappa*dphi*4/(2*pi)^2;
  
  SigXw = [SigXw; SigX];
  SigPw = [SigPw; SigP];
  SigLw = [SigLw; SigL];

endfor;%# end nus loop

Sig = SigXw + SigPw + SigLw;
xi  = s*vF*norm(k)-EF;
ws  = nus*EF;
save output Sig ws xi EF kF SigXw SigPw SigLw
more on;
