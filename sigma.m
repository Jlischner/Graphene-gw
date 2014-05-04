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

#Nk = length(kys);
k  = [0 ky*kF];
n   = kF^2/pi*cm^2; %# electron density 

%# good convergence parameters:
%# dphi = 0.01 and du=0.001
phis  = [0:0.01:2*pi]';
us   = [0.01:0.001:1.2]';
%# for kappa=4, use dwts=0.2
wts = [0.1:0.2:25]'*sqrt(EF);
#eta = 0.02;


%# print info
#printf("doing k=%f *kF \n",ky);
printf("using Efermi = %f eV \n", EFeV);
printf("using vF = %f * 10^8 cm/sec \n",vF*sec/cm/10^8);
printf("using kc = %f bohr \n",kc);
printf("using kappa %f \n",kappa);
printf("using D=%f bohr \n",D);
printf("convergence parameter: \n");
printf("dphi = %f \n",phis(2)-phis(1));
printf("du = %f \n",us(2)-us(1));
printf("u_max = %f \n", us(length(us)));
printf("dwts = %f \n", ( wts(2)-wts(1) ) /sqrt(EF));
printf("wts_max = %f \n",wts(length(wts))/sqrt(EF));
printf("eta = %f \n",eta);

dphi = phis(2)-phis(1);
qs   = kc*us.^2;
Nq   = length(qs);
epsB = 2*kappa-1;
%# parameter for sigma -> pi screening model
eps1 = 2.4;
d1 = 2.8*angstrom;
epsI = 1/eps1*(eps1+1+(eps1-1)*exp(-d1*abs(qs)))./(eps1+1-(eps1-1)*exp(-d1*abs(qs)));
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

%# first compute interaction contribution to dmu
SigL = 0.;
SigX = 0.;
SigP = 0.;
w = 0;
kapp = k;
sapp = s;
k = [0 kF];
s = +1;

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
    
  endfor;
endfor;

SigX *= -kc*dphi*2/(2*pi)^2;
printf("SigX = %f \n",SigX); 
SigL *=  kc*dphi*4/(2*pi)^3;
printf("SigL = %f \n",SigL); 
SigP *=  kc*dphi*2/(2*pi)^2;
dmu = SigX + SigL + SigP;

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

  SigXw = [SigXw; SigX];
  SigPw = [SigPw; SigP];
  SigLw = [SigLw; SigL];

endfor;%# end nus loop

Sig = SigPw + SigLw + SigXw;
xi  = s*vF*norm(k)-EF;
ws  = nus*EF;
ef  = EF;
A   = 1/pi * abs(imag(Sig)) ./ ( (ws-xi-real(Sig-dmu)).^2 + imag(Sig).^2 );

save output ws xi ef kF SigPw SigLw Sig SigXw dmu A kappa D vF kc s kF vF;
more on;
