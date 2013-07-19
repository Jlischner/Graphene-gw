%# Graphene-gw takes EFeV and kappa as input params
%# convert from interaction strength alpha and density n

alpha = 0.5;
ncm = 10^(12); 

angstrom = 1/0.5291;
sec = 1/2.418884326505e-17; %# wikipedia: atomic units
cm = 10^8*angstrom;
vF = 0.85*10^8*cm/sec; %# LDA value

kappa = 1/alpha/vF;
EF = sqrt( pi*ncm * vF^2/cm^2 );
EFeV = EF*27.21;

printf("alpha = %f and n=%f x 10^12 1/cm^2 \n",alpha,ncm*10^-12);
printf("result in (using LDA value for vF): \n")
printf("kappa = %f and EF = %f eV \n",kappa,EFeV);