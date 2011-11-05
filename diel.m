function epsilon = diel(q,w)

  g = 4;
  global vF;
  global EF;
  global kF;
  global Vc;
  global e0;

  eta = 0.03;
  x = q/kF;

  aw = abs(w);
  nu = aw/EF;

  F  = g/16/pi* vF^2*q.^2 ./sqrt(aw.^2 - vF^2*q.^2);
  
  dP = -g*EF/2/pi/vF^2 + F/vF^2 .* ( G((nu+2)./x) -...
	(((2-nu)./x-1.)>0).* (G((2-nu)./x) - I*pi) - (((nu-2)./x+1.)>0).*G((nu-2)./x) );

  P0 = -I*pi*F/vF^2; 
  Pi  = P0 + dP;

  epsilon = e0 -Vc .* Pi + I*eta;
  epsilon(w<0) = epsilon(w<0)';
  #epsilon = Pi;
endfunction;
