function [Peq7,Peq6,Ldp]=pdtotal(u,Np)
     
     //----------------------------------------------------------------//
    uo=4*%pi*1e-7;
    A1=(1.9/100)*(4.6/100);
    A2=(3.8/100)*(4.6/100);
    l1=4.75/100;
    l2=7.6/100;    
      
    PTadyizq=6.385*1e-9;
    PTadyder=6.385*1e-9;
    PTvizq=7.515*1e-9;
    PTvder=7.515*1e-9;
    
    //PTadyizq=0.0;
    //PTadyder=0.0;
    //PTvizq=0.0;
    //PTvder=0.0;
    
    Pyoke=(l1/(u*A1))^-1;
    Pleg=(l2/(u*A1))^-1;
    Pmleg=(l2/(u*A2))^-1;
    Peq1=(Pyoke*Pleg)/(2*Pleg+Pyoke);
    Peq2=(Pyoke*Pleg)/(2*Pleg+Pyoke);
    Peq3=PTadyizq+PTadyder;
    Peq4=PTvizq+PTvder;
    Peq5=Peq1+Peq2;
    Peq6=(Peq5*Pmleg)/(Pmleg+Peq5);
    Peq7=Peq6+Peq3+Peq4;
    //****Inductancia de Dispersion****//
    Ldp=(Peq7)*Np^2;

endfunction
