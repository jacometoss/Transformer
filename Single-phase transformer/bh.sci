function  [B,H]=bh(Voc,Ioc,Np,l,A)
        fhz=60;
        B=Voc/(sqrt(2)*%pi*fhz*Np*A);
        H=(Ioc*Np)/l;
endfunction
