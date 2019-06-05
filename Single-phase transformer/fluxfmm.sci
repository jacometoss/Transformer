function  [Flux,Fmm]=fluxfmm(Voc,Ioc,Np,l,A)
    fhz=60;
    B=Voc/(sqrt(2)*%pi*fhz*Np*A);
    Flux=B*A;
    Fmm=Np*Ioc;
endfunction
