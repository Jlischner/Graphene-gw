more off;
addpath("/auto/jlischner/2Sarma/DasSarma/Graphene-gw");

eta = 0.0001;
%# get dmu from upper band at kF
load Upper/output;

nus = ws/ef;
wss = [min(nus)+0.01 : 0.01 : max(nus)-0.01]'*ef;
Nfreq = length(wss);
Nk = length(xis);
Ass = zeros(Nfreq,Nk);
Sigss = zeros(Nfreq,Nk);

%# get dmu:
[a,b] = min(abs(xis));
ikf = b;
printf("smallest xi=%f eV\n",a*27.21);

ImSigF = interp1(ws,ImSigs(:,ikf),wss);
ReSigF = getReSig(wss,ImSigF,eta);
ReSigF += SigXs(ikf);

[a,b] = min(abs(wss/ef));
if( abs(a) < 0.001 );
  dmu = ReSigF(b);
  printf("dmu = %f eV \n",dmu*27.21);
else;
  printf("w=0 not in frequency grid! \n");
  exit;
endif;

%# problem in GW+C at k close to kF
Nkmax = Nk;
load Lower/output;
ImSigMax = max(abs(ImSigF));
shift = ImSigMax/100;

for ik = 1:Nkmax;
  
  printf("doing ik=%d of %d states \n",ik,Nkmax);
  xi = xis(ik);
  ImSig = interp1(ws,ImSigs(:,ik),wss);
  ReSig = getReSig(wss,ImSig,eta);
  ReSig += SigXs(ik);

#  Sig = Sigs(:,ik);
#  ReSig = interp1(ws,real(Sig),wss);

  ImSig += shift;

  Sig = ReSig + I*ImSig;
  A   = 1/pi * abs(imag(Sig)) ./ ( (wss-xi-real(Sig-dmu)).^2 + imag(Sig).^2 );
  
  Ass(:,ik) = A;
  Sigss(:,ik) = Sig;

end;

more on;