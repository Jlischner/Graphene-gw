more off;
addpath("/auto/jlischner/2Sarma/DasSarma/Graphene-gw");

angstrom = 1/0.5291;
sec = 1/2.418884326505e-17; %# wikipedia: atomic units
cm = 10^8*angstrom;
bands = [-1 +1];

%# read input
input;

global vF = 0.85*10^8*cm/sec;
global EF = EFeV/27.21;
global kF = EF/vF;

if( kcoption == 1);
  printf("Das Sarma cutoff \n");
  kc = 1/(2.46*angstrom); %# Das Sarma cutoff wave vector
elseif( kcoption == 2);
  printf("V_BZ conserving cutoff \n");
  kc = 2* sqrt( 2*pi/(3*sqrt(3)) ) / (1.42*angstrom); %# V_BZ conserving cutoff
elseif( kcoption == 3);
  printf("kc=30*kF cutoff \n");
  kc = 30*kF;
endif;
printf("kc/kF = %f \n",kc/kF);

Nk = length(kys);
Nf = length(nus);
freqs = nus*EF;
n   = kF^2/pi*cm^2; %# electron density in cm^-2

%# good convergence parameters:
%# dphi = 0.01 and du=0.001
phis  = [0:0.01:2*pi]';
us   = [0.01:0.001:1.2]';

dphi = phis(2)-phis(1);
qs   = kc*us.^2;
Nq   = length(qs);

%# parameter for sigma -> pi screening model for eps^{-1}:
if( wehlingflag == 1);
  eps1 = 2.4;
  d1 = 2.8*angstrom;
  epsW = 1/eps1*(eps1+1+(eps1-1)*exp(-d1*abs(qs)))./(eps1+1-(eps1-1)*exp(-d1*abs(qs)));
  vchiW = 1 - 1./epsW; 
else;
  vchiW = 0;
endif;

%# screening from substrate:
epsB = 2*kappa-1;
epsS = 1-(epsB-1)/(epsB+1).*exp(-2*D*abs(qs));
vchiS = 1 - 1./epsS;
epsI = 1./(1 - vchiS - vchiW); %# includes contributions from substrate and sigma -> pi stuff
global Vc = 2*pi./abs(qs).*epsI;
qVc = qs .* Vc;

%# calculate self energies
ImSigs  = zeros(Nf,Nk);
SigXs   = zeros(Nk,1);

for ik = 1:Nk;
  
  k = [0 kys(ik)*kF];

  for iw = 1:Nf;

    printf("freq: %d out of %d, kpoint: %d of ouf %d \n",iw,Nf,ik,Nk);
    w = freqs(iw);

    SigX = 0.;
    SigP = 0.;

    for sp = bands;
      for phi = phis';
    
	kpq = qs*[cos(phi) sin(phi)] + ones(Nq,1)*k;
	abskpq = sqrt(sum(kpq.^2,2));
	xi  = sp*vF*abskpq - EF;
	uFss= us .* F(s,sp,k,kpq) .* (kc>abskpq);

	if(iw == 1);
	  SigX += trapz(us, (xi<0) .* uFss .* qVc );
	endif;

	Ieps = 1./diel( qs, xi-w , eta) - 1.;
	SigP += trapz(us, ( ((w-xi)>0) - (xi<0) ) .* uFss .* qVc .* Ieps );
	
      endfor;%# end phis loop      
    endfor;%# end bands loop
    
    SigX *= -kc*dphi*2/(2*pi)^2; 
    SigP *=  kc*dphi*2/(2*pi)^2;             
    ImSigs(iw,ik)  = imag(SigP);

    if(iw == 1);
      SigXs(ik) = SigX;
    endif;

  endfor;%# end nus loop
endfor; %# end kys loop

xis = s*vF*abs(kys)*kF-EF;
ws  = freqs;
ef  = EF;
vf  = vF 

save output ws xis ef kF ImSigs SigXs kappa D vF kc s vf epsI;
more on;
