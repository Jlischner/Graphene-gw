function epsilon = Cdiel(qs, w);
  global vF;
  global EF;
  global Vc;

  Omega = w;

  chiMDF = -qs.^2./16./sqrt( Omega.^2 + vF^2*qs.^2 ) - ...
      EF/2/pi/vF^2 + qs.^2./8/pi./sqrt( Omega.^2 + vF^2*qs.^2 ) .*...
      real( asin( (2*EF+I*Omega)/vF./qs ) + (2*EF+I*Omega)/vF./qs .*...
	   sqrt(1-(( 2*EF + I*Omega )/vF./qs ).^2 ) );
  
  Pi0 = -qs.^2/16./sqrt(Omega.^2+vF^2*qs.^2);
  Pi = 4*chiMDF;
  epsilon = 1 - Vc .* Pi;
#  epsilon = Pi;

endfunction;
