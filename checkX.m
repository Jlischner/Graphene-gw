input;

angstrom = 1/0.5291;
cm = 10^8*angstrom;
sec = 1/2.418884326505e-17; %# wikipedia: atomic units

EF = EFeV/27.21;
vF = 10^8*cm/sec;
kF = EF/vF;


if( kcoption == 1);
  kc = 1/(2.46*angstrom); %# cutoff wave vector
elseif( kcoption == 2);
  kc = 100*kF;
endif;

SigXiv = [];
SigXev = [];

ks = ky*kF;

for k = ks';
  
  xi = k/kc;
  xe = k/kF
  
  [ki,ei] = ellipke( xi^2 );
  [ke,ee] = ellipke( xe^2 );
  fi = ei;
  fe = ee;
  
  hi = xi*(pi/4*log(4/xi) -pi/8 );
  he = xe*(pi/4*log(4/xe) -pi/8 );
  
  yt = [0.002:0.00001:1]';
  ys = xi * yt;
  [Ks,Es] = ellipke( ys.^2 );
  Ii  = trapz(yt, (Ks-Es-pi/4*ys.^2)./ys.^3 );
  hi -= Ii*xi^2;
  
  ys = xe * yt;
  [Ks,Es] = ellipke( ys.^2 );
  Ie  = trapz(yt, (Ks-Es-pi/4*ys.^2)./ys.^3 );
  he -= Ie*xe^2;
  
  
  SigXi = -fi + s*hi;
  SigXe = -fe - s*he;
  
  SigXiv = [SigXiv; SigXi];
  SigXev = [SigXev; SigXe];
  
endfor;

SigXiv *= kc/pi/kappa;
SigXev *= kF/pi/kappa;
SigX = SigXiv + SigXev;

#Sigi = 1/kappa*(-kc/2 + s*ks/4.*(log(4*kc./ks) - 0.5 ));
#Sige = 1/kappa*(-kF/2 - s*ks/4.*(log(4*kF./ks) - 0.5 ));

save xout SigX
