function ReSig = getReSig(ws,ImSig,eta);

  Nf = length(ws);
  ReSig = zeros(size(ImSig));
  for ii = 1:Nf;
    ReSig(ii) += trapz(ws, ImSig./( ws-ws(ii) + I*eta));
  endfor;
  ReSig = -real(ReSig)/pi;
  
endfunction;