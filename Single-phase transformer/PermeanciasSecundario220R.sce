//****************************************************************// 
// .Rutina: Permeancia Devanado Primario 127V                  ...//
// .Autor: Marco Polo Jacome Toss                              ...//
// .Version : 1.0                                              ...//
// .Plataforma : Scilab (https://www.scilab.org) 6.0           ...//
// .Fecha : 2018.06.10                                         ...//
// .Nota  : Apendice A - Dimensiones en Devanados              ...//
// ...      Capitulo 2 - Permeancia de Fuga                    ...//
// .... ..... ..... .... .... .... .... .... .... .... .... ... ..// 
// .PW1CorteAA : Permeancia interior (Sólo una ventana)        ...//
// .PW1CorteBB : Permeancia Adyacente al nucleo (Sólo un lado) ...//
// .PW1CorteCC : Permeancia En rel Radio (Un radio)            ...//
// ...... ..... .... .... .... .... .... .... .... .... .... .....//
// .PW1_Total_CorteAA : Permeancia Total en las dos ventanas   ...//
// .PW1_Total_CorteBB : Permeancia Total adyacentes (izq. y der.).//
// .... ..... ..... .... .... .... .... .... .... .... .... ......// 
uo=4*%pi*1e-7;
d1=5.3/100; 
d2=0.2/100; 
d3=0.05/100; //Distancia influye mucho, unidades metros.
d4=0.4/100;
d5=0.425/100;
//-----------------------------------------//
la=0.638/100;
lb=5.5/100;
lc=5/100; 
ld=4.2/100;
lrm=2.278/100;
//----------------------------------------//
k1=abs(d1-d4);
k2=min(d1,d4);
//Permeancnia del Devanado de 127V 
PW1CorteAA=((uo*lc)/(128*(d4^2)*(d1^2)))*((4*k2^4)+(8*k1*k2^3)+(2*k1^2*k2^2)-(2*k1^3*k2)+(k1^4*log((2*k2/k1+1))));
PW1CorteBB=((uo*ld)/(128*(d4^2)*(d1^2)))*((4*k2^4)+(8*k1*k2^3)+(2*k1^2*k2^2)-(2*k1^3*k2)+(k1^4*log((2*k2/k1+1))));
PW1CorteCC=((uo*lrm)/(128*(d4^2)*(d1^2)))*((4*k2^4)+(8*k1*k2^3)+(2*k1^2*k2^2)-(2*k1^3*k2)+(k1^4*log((2*k2/k1+1))));
//disp (PWinding1);

//Permeancia Horizontal En Aire Devanado 127V
Pw1ahCorteAA=((uo*lc)/((2*d4+d3+d5)/(2*d2))) 
Pw1ahCorteBB=((uo*ld)/((2*d4+d3+d5)/(2*d2)))
Pw1ahCorteCC=((uo*lrm)/((2*d4+d3+d5)/(2*d2))) 
//disp (Pw1ah);
//Permeancia VErtical En Aire Devanado 127V
Pw1avCorteAA=((uo*lc)/((d1+d2)/(d1)))
Pw1avCorteBB=((uo*ld)/((d1+d2)/(d1)))
Pw1avCorteCC=((uo*lrm)/((d1+d2)/(d1))) 
//disp (Pw1av);
//Permeancia Total Devanado Primario 127 V  Aire
PW1AirCorteAA=(2*((1/Pw1ahCorteAA)+(1/Pw1avCorteAA)))^-1;
PW1AirCorteBB=(2*((1/Pw1ahCorteBB)+(1/Pw1avCorteBB)))^-1;
PW1AirCorteCC=(2*((1/Pw1ahCorteCC)+(1/Pw1avCorteCC)))^-1;
//Permeancia Total (Aire y Devanado)
PW1CorteAA=PW1AirCorteAA+PW1CorteAA;
PW1CorteBB=PW1AirCorteBB+PW1CorteBB;
PW1CorteCC=PW1AirCorteCC+PW1CorteCC;
disp('****** Permeancia  Devanado Primario 127V Corte A-A Lado Izquierdo********')
disp (PW1CorteAA);
disp('****** Permeancia  Devanado Primario 127V Exterior Corte B-B Lado Izquierdo********')
disp (PW1CorteBB);
disp('****** Permeancia  Devanado Primario 127V Exterior Corte C-C Un Radio********')
disp (PW1CorteCC);
disp('****** Permeancia Total Vista Corte A-A********')
PW1_Total_CorteAA=2*PW1CorteAA+2*PW1CorteCC;;
disp(PW1_Total_CorteAA)
disp('****** Permeancia Total Vista Corte B-B********')
PW1_Total_CorteBB=2*PW1CorteBB+2*PW1CorteCC;
disp(PW1_Total_CorteBB)
