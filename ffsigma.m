%# this code does not calculate dmu any more:
%# need postprocessing to calculate spectral functions
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
n   = kF^2/pi*cm^2; %# electron density 

%# good convergence parameters:
%# dphi = 0.01 and du=0.001
phis  = [0:0.01:2*pi]';
us   = [0.01:0.001:1.2]';
%# for kappa=4, use dwts=0.2
wts = [0.1:0.2:25]'*sqrt(EF);

dphi = phis(2)-phis(1);
qs   = kc*us.^2;
Nq   = length(qs);

%# parameter for sigma -> pi screening model:
eps1 = 2.4;
d1 = 2.8*angstrom;
epsI = 1/eps1*(eps1+1+(eps1-1)*exp(-d1*abs(qs)))./(eps1+1-(eps1-1)*exp(-d1*abs(qs)));

%# screening from substrate:
epsB = 2*kappa-1;
#global Vc = 2*pi./abs(qs) .*(1-(epsB-1)/(epsB+1).*exp(-2*D*abs(qs)) ); 
global Vc = 2*pi./abs(qs).*epsI;
qVc = qs .* Vc;
ws = wts.^2;
Nfreq= length(wts);

%# get dielectric matrix for line integral
epsC = zeros(Nq,Nfreq);
for n = 1:Nfreq;
  IepsC(:,n) = 1./Cdiel( qs, ws(n) ) -1.;
endfor;

%# calculate self energies
SigPw = zeros(Nf,Nk);
SigLw = zeros(Nf,Nk);
SigXw = zeros(Nf,Nk);
Sigs  = zeros(Nf,Nk);

for ik = 1:Nk;
  
  k = [0 kys(ik)*kF];

  for iw = 1:Nf;

    printf("freq: %d out of %d, kpoint: %d of ouf %d \n",iw,Nf,ik,Nk);
    w = freqs(iw);

    SigX = 0.;
    SigP = 0.;
    SigL = 0.;

    for sp = bands;
      for phi = phis';
    
	kpq = qs*[cos(phi) sin(phi)] + ones(Nq,1)*k;
	abskpq = sqrt(sum(kpq.^2,2));
	xi  = sp*vF*abskpq - EF;
	uFss= us .* F(s,sp,k,kpq) .* (kc>abskpq);
	SigX += trapz(us, (xi<0) .* uFss .* qVc );
	
	Ieps = 1./diel( qs, xi-w , eta) - 1.;
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
    
    SigXw(iw,ik) = SigX;
    SigPw(iw,ik) = SigP;
    SigLw(iw,ik) = SigL;
    Sigs(iw,ik)  = SigX + SigP + SigL;

  endfor;%# end nus loop
endfor; %# end kys loop

xis = s*vF*abs(kys)*kF-EF;
ws  = freqs;
ef  = EF;

save output ws xis ef kF SigPw SigLw SigXw Sigs kappa D vF kc s kF vF;
more on;
