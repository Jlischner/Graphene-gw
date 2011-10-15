function out = G(x);

  out = x.*sqrt(x.^2-1) - log(x+sqrt(x.^2-1));
endfunction;