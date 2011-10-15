more off;
x= [0.01:0.05:3]';

angstrom = 1/0.5291;
sec = 1/2.418884326505e-17; %# wikipedia: atomic units
cm = 10^8*angstrom;
n = 10^13/cm^2;
vF = 10^8*cm/sec;

D0 = sqrt(4*n/pi)/vF;
kF = sqrt(pi*n);
EF = vF*kF;
kappa = 1;
q = x*kF;
Vc = 2*pi./abs(q)/kappa;
vmax = 0.5;

mat = [];
#for nu  =  0:0.05:3;

nu=1.3;
mnu = -nu;

f1 = (2+nu).*sqrt((2+nu).^2 -x.^2) - ...
    x.^2.*log( (sqrt( (2+nu).^2 -x.^2) + (2+nu))./abs(sqrt(nu.^2-x.^2)+nu));

f2 = x.^2 .* log( (nu-sqrt(nu.^2 - x.^2))./x);

f3 = (2+nu).*sqrt(x.^2-(2+nu).^2) + x.^2.*asin((2+nu)./x);

f4 = (2+nu).*sqrt( (2+nu).^2 - x.^2) - ...
    x.^2.* log( ( sqrt((2+nu).^2 -x.^2) + (2+nu))./x );

f1m = (2+mnu).*sqrt((2+mnu).^2 -x.^2) - ...
    x.^2.*log( (sqrt( (2+mnu).^2 -x.^2) + (2+mnu))./abs(sqrt(mnu.^2-x.^2)+mnu));

f3m = (2+mnu).*sqrt(x.^2-(2+mnu).^2) + x.^2.*asin((2+mnu)./x);

f4m = (2+mnu).*sqrt( (2+mnu).^2 - x.^2) - ...
    x.^2.* log( ( sqrt((2+mnu).^2 -x.^2) + (2+mnu))./x );

RePi1p = 1.-1./8./sqrt(nu.^2-x.^2).*(f1 .* ((abs(2+nu)-x)>0)+...
    sign(nu-2+x).*f1m.*((abs(2-nu)-x)>0)+f2.*( ((x+2-nu)>0) + ((2-x-nu)>0)));

RePi2p = 1.-1./8./sqrt(x.^2-nu.^2).*(f3 .* ((x-abs(nu+2))>0)+...
    f3m.*((x-abs(nu-2))>0)+pi*x.^2/2.*( ((abs(nu+2)-x)>0) + ((abs(nu-2)-x)>0)));

ImPi1p = -1./8./sqrt(nu.^2-x.^2) .*( f3m.*((x-abs(nu-2))>0)+... 
    pi*x.^2/2.* ( ((x+2-nu)>0) +((2-x-nu)>0) ));

ImPi2p = ((nu-x+2)>0)./8./sqrt(x.^2-nu.^2) .* ( f4 - f4m.*((2-x-nu)>0));

Pip = (RePi1p + I*ImPi1p).*((nu-x)>0) + ...
    (RePi2p + I*ImPi2p).*((x-nu)>0);

Pim = pi*x.^2.*((x-nu)>0)./8./sqrt(x.^2-nu.^2) + ...
    I*pi*x.^2.*((nu-x)>0)./8./sqrt(nu.^2-x.^2);

Pi = D0*( Pip + Pim );
epsilon = 1-Vc .* Pi;

v = real(1./epsilon);
v(v>vmax) = vmax; 
mat = [mat; nu*ones(size(x)) x v];
#endfor;
more on;