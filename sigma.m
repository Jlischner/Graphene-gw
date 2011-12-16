more off;
angstrom = 1/0.5291;
sec = 1/2.418884326505e-17; %# wikipedia: atomic units
cm = 10^8*angstrom;
bands = [-1 +1];

%# read input
input;

global e0 = 1.0
global vF = 10^8*cm/sec;
global EF = EFeV/27.21;
global kF = EF/vF;
global mu = EF;

if( kcoption == 1); 
  kc = 1/(2.46*angstrom); %# cutoff wave vector
elseif( kcoption == 2);
  kc = 100*kF;
endif;

k  = [0 ky*kF];
n   = kF^2/pi*cm^2; %# electron density 

%# good convergence parameters:
%# dphi = 0.01 and du=0.001
phis  = [0:0.01:2*pi]';
dphi = phis(2)-phis(1);
us   = [0.001:0.001:1.2]';
qs   = kc*us.^2;
Nq   = length(qs);
#global Vc = 2*pi/kappa./abs(qs);
epsB = 2*kappa-1;
global Vc = 2*pi./abs(qs) .*(1-(epsB-1)/(epsB+1).*exp(-2*D*abs(qs)) ); 
qVc = qs .* Vc;

%# was using wts = [0.0:0.2:15]'*sqrt(EF)
wts = [0.0:0.2:25]'*sqrt(EF);
ws = wts.^2;
Nfreq= length(wts);

%# get dielectric matrix for line integral
epsC = zeros(Nq,Nfreq);
for n = 1:Nfreq;
  IepsC(:,n) = 1./Cdiel( qs, ws(n) ) -1.;
endfor;

%# first compute interaction contribution to mu
SigL = 0.;
SigX = 0.;
w = 0;
kapp = k;
sapp = s;
k = [0 kF];
s = +1;

for sp = bands;
  for phi = phis';
    
    kpq = qs*[cos(phi) sin(phi)] + ones(Nq,1)*k;
    abskpq = sqrt(sum(kpq.^2,2));
    xi  = sp*vF*abskpq - mu;
    uFss= us .* F(s,sp,k,kpq) .* (kc>abskpq);
    SigX += trapz(us, (xi<0) .* uFss .* qVc );   
    
    SigLv = zeros(Nfreq,1);
    for n = 1:Nfreq;
      SigLv(n) += trapz(us, uFss.*qVc .*(xi-w)./( (xi-w).^2 + ws(n).^2) .* IepsC(:,n) );
    endfor;
    SigL += trapz(wts,2*wts.* SigLv);
    
  endfor;
endfor;

SigX *= -kc*dphi*2/(2*pi)^2;                                                                                          
SigL *=  kc*dphi*4/(2*pi)^3; 
dmu = SigX + SigL;

%# now back to kpoint of interest
k = kapp;
s = sapp;
SigPw = [];
SigLw = [];
SigXw = [];

for nu = nus';

  printf("%f \n",nu);
  w = nu*EF;

  SigX = 0.;
  SigP = 0.;
  SigL = 0.;

  for sp = bands;
    for phi = phis';
    
      kpq = qs*[cos(phi) sin(phi)] + ones(Nq,1)*k;
      abskpq = sqrt(sum(kpq.^2,2));
      xi  = sp*vF*abskpq - mu;
      uFss= us .* F(s,sp,k,kpq) .* (kc>abskpq);
      SigX += trapz(us, (xi<0) .* uFss .* qVc );

      Ieps = 1./diel( qs, xi-w ) - 1.;
      SigP += trapz(us, ( ((w-xi)>0) - (xi<0) ) .* uFss .* qVc .* Ieps );

      SigLv = zeros(Nfreq,1);
      for n = 1:Nfreq;
	SigLv(n) += trapz(us, uFss.*qVc .*(xi-w)./( (xi-w).^2 + ws(n).^2) .* IepsC(:,n) );
      endfor;
      SigL += trapz(wts,2*wts.* SigLv);
      
    endfor;%# end phis loop
  endfor;%# end bands loop
 
  SigX *= -kc*dphi*2/(2*pi)^2;                                                                                          
  SigP *=  kc*dphi*2/(2*pi)^2;                                                                                          
  SigL *=  kc*dphi*4/(2*pi)^3; 

  SigXw = [SigXw; SigX];
  SigPw = [SigPw; SigP];
  SigLw = [SigLw; SigL];

endfor;%# end nus loop

Sig = SigPw + SigLw + SigXw;
xi  = s*vF*norm(k)-EF;
ws  = nus*EF;
ef  = EF;
A   = 1/pi * abs(imag(Sig)) ./ ( (ws-xi-real(Sig-dmu)).^2 + imag(Sig).^2 );
save output ws xi ef kF SigPw SigLw Sig SigXw dmu A;
more on;
