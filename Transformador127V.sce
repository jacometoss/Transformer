//**********************************************************************// 
// ... Rutina: Transformador Monofásico 127V SIN CARGA               ...//
// ... Autor: Marco Polo Jacome Toss                                 ...//
// ... Version : 0.01                                                ...//
// ... Plataforma : Scilab (https://www.scilab.org) 6.0              ...//
// ... Fecha : 2018.06.10                                            ...//
// .... ..... ..... .... .... .... .... .... .... .... .... .... .... ..//
// ... MAGCURVE_127_TA : Curva de Saturacion [V-I]                   ...//
// ... function [B,H]=BH(Voc,Ioc,Np,l,A) : Funcion Curva BH          ...//
// ... function [Flux,Fmm]=FluxFMM(Voc,Ioc,Np,l,A) Funcion Curva FmmFlux
// ...           >Flux: Flujo Magnetico (Wb)                         ...//
// ...           >Fmm: Fuerza Magnetomotriz (A-V)                    ...//
// ... function [Peq7,Peq6,Ldp]=PermeanciaTotal(u,Np) Permeancia Total .//
// ...           >Peq7: Permeancia Total                             ...//
// ...           >Peq6: Permeancia                                   ...//
// ...           >Ldp : Inductancia de dispersion                    ...//
// ... function Iaeq                                                 ...//
// ...          >    Model: (1) para modelo T (2) para Pi            ...//
// ...          >  fluxmag: Flujo Magnertico (Wb)                    ...//
// ...          >         u: Permeabilidad                           ...//
// ...          >        Ia: Corriente (A) Inductancia Saturable     ...//
// ...          >   dfluxdt: Derivada del flujo magnetico            ...//
// ...          >      didt: Derivada de la corriente malla primaria ...//
// ...          >         t: Tiempo de simulacion                    ...//
// ...          >      Bvar: Densidad de flujo magnetico (T)         ...//
// ...          >       Irn: Corriente en R magnetizacion            ...//
// ...          >      Eind: Voltaje Inducido                        ...//
// ...          >     Iload: Corriente de carga                      ...//
// ...          >  SeccionI: Seccion de columna central en nucleo    ...//
// ...          > SeccionII: Seccion de columnas externas en nucleo  ...//
// ...          >    Np_127: Vueltas en lado primario (127V)         ...//
// ...          >    Ns_127: Vueltas en lado secundario (127V)       ...//
// ...          >    lmFlux: Longitud media circuito magnetico (m)   ...//
// ...          >        Rn: Resistencia del nucleo  Prueba C.O.     ...//
// ...          >        Xm: Inductancia del nucleo  Prueba  C.O     ...//
// ...          >       Req: Resistencia del conductor (Ohm)         ...//
// ...          >       Xeq: Resistencia del conductor (Ohm)         ...//
// ...          >     Rload: Resistencia de carga      (Ohm)         ...//
// ...          >     Xload: Inductancia de carga      (Ohm)         ...//
// ...          >    Vtotal: Voltaje de suministro comprobacion      ...//
// ... function PermeanciaTotal                                      ...//
// ...          >  Peq1: Permeancia Eq1 Columna y yugos Izquierdos   ...//
// ...          >  Peq2: Permeancia Eq2 Columna y yugos Derechos     ...//
// ...          >  Peq3: Permeancias Externas (Devanado adyacente)   ...//
// ...          >  Peq4: Permeancia Internas (Devanado interna )     ...//
// ...          >  Peq5: Permeancia Eq5 yugos y piernas externas     ...//
// ...          >  Peq6: Permeancia Eq6 (Pierna central y Eq5)       ...//
// ...          >  Peq7: Permeancia Eq7 (Total)                      ...//
// ...          >  PTadyizq: Permeancia Adyacente al nucleo Izquierda...//
// ...          >  PTadyder: Permeancia Adyacente al nucleo Derecha  ...//
// ...          >  PTvizq:   Permeancia Ventana Izquierda            ...//
// ...          >  PTvder:   Permeancia Ventana Derecha              ...//
// ...          >  Nota: PTadyizq,PTadyder,PTvizq,PTvder             ...//
// ...                   Se calculan a partir de otros archivos      ...//
//**********************************************************************//
global MAGCURVE_127_TA;
MAGCURVE_127_TA=[
-270,-10;
-237,-7.00;
-233,-6.7;
-220,-5.5;
-215,-5.10;
-208,-4.50;
-203,-4.00;
-199,-3.6;
-195,-3.2;
-190,-2.70; 
-183,-2.05;
-175,-1.58;
-168,-1.25;
-128,-0.7;
-80.3,-0.19;
-55.00,-0.15;
55.00,0.15;
80.3,0.19;
128,0.7;
168,1.25;
175,1.58;
183,2.05;
190,2.70; 
195,3.2;
199,3.6;
203,4.00;
208,4.50;
215,5.10;
220,5.5;
233,6.7;
237,7.00;
270,10];

funcprot(0)
    function  [B,H]=BH(Voc,Ioc,Np,l,A)
        fhz=60;
        B=Voc/(sqrt(2)*%pi*fhz*Np*A);
        H=(Ioc*Np)/l;
endfunction
    
funcprot(0)
function  [Flux,Fmm]=FluxFMM(Voc,Ioc,Np,l,A)
    fhz=60;
    B=Voc/(sqrt(2)*%pi*fhz*Np*A);
    Flux=B*A;
    Fmm=Np*Ioc;
endfunction
       
//************************ATRIBUTOS*****************************//
	   
    SeccionI=(3.8/100)*(4.6/100); 
    SeccionII=(1.9/100)*(4.6/100);
    Np_127=262;
    Ns_127=249;
    lmFlux=0.248;
//*************************PRIMARIO*****************************//
    Rnp=560.56;
    Xmp=279.32; 
    Reqp=3.29;
    Xeqp=0.001908303;        
//*************************SECUNDARIO****************************//
    Rns=560.56;
    Xms=279.32; 
    Reqs=3.29;
    Xeqs=0.001987752;        
//-----------------------------------------------------//
    Rload=0;
    Xload=0;

function [Peq7,Peq6,Ldp]=PermeanciaTotal(u,Np)
     
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
    /*
    PTadyizq=0.0;
    PTadyder=0.0;
    PTvizq=0.0;
    PTvder=0.0;
    */
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

function [LMAG]=LM(u,Np)
     
    uo=4*%pi*1e-7;
    A1=(1.9/100)*(4.6/100);
    A2=(3.8/100)*(4.6/100);
    l1=4.75/100;
    l2=7.6/100;    
      
    PTadyizq=0.0;
    PTadyder=0.0;
    PTvizq=0.0;
    PTvder=0.0;
   
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
    LMAG=(Peq7)*Np^2;
  
endfunction

function [diferencial]=diff(M,h)
    
    diferencial=M;
    diferencial(1)=((M(2))-M(1))/h;
for i=2:max(size(M))-1
    diferencial(i)=((M(i+1))-(M(i-1)))/(2*h);
end
    diferencial(max(size(M)))=(M(max(size(M)))-(M(max(size(M-1)))))/(h);
    
endfunction 


function [fluxmag,Peq7,u,Ia,Ldp,dfluxdt,didt,t,E,Bvar,Irn,Eind,I,Vsum]=Iaeq(Np,Ns,Req,Xeq,fluxmag,SeccionI,SeccionII,Theta,B,H,h,t0,N,Vm,W,Rn,Model)  
  
        t(1) = t0;

for i = 1:N  

select Model
    
case 1 then
    
       Bvar(i)=(fluxmag(i))/SeccionI;
       Bvar(i)=abs(Bvar(i))
       Hvar(i)=interp1(B,H,Bvar(i));
       u(i)=(Bvar(i)/Hvar(i)); 
       [Peq7(i),Peq6(i),Ldp(i)]=PermeanciaTotal(u(i),Np)    
       disp([(t(i)),Vm,Ldp(i),Peq7(i)])
       Ia(i)=(Np*fluxmag(i))/Ldp(i);
       Lambda(i)=fluxmag(i)*Np;
    
        k1 = h*((1/Np)*(Vm*cos(W*t(i)+Theta)-((Req*Np*fluxmag(i))/(Ldp(i)))));
        k2 = h*(((1/Np)*(Vm*cos(W*(t(i)+h/2)+Theta)-((Req*Np*fluxmag(i))/(Ldp(i))))+0.5*k1));
        k3 = h*(((1/Np)*(Vm*cos(W*(t(i)+h/2)+Theta)-((Req*Np*fluxmag(i))/(Ldp(i))))+0.5*k2));                    
        k4 = h*(((1/Np)*(Vm*cos(W*(t(i)+h)+Theta)-((Req*Np*fluxmag(i))/(Ldp(i))))+k3));
        fluxmag(i+1) = fluxmag(i) + (k1 + 2*k2 + 2*k3 + k4)/6;
        
       Ia(i+1)=(Np*fluxmag(i))/(abs(Peq7(i))*Np^2);
       Peq7(i+1)=Peq7(i);
       Ldp(i+1)=Ldp(i);
       t(i+1) = t0 + i*h;
       u(i+1)=u(i);
       Lambda(i+1)=fluxmag(i)*Np;;        

case 2 then
    
       Bvar(i)=(fluxmag(i))/SeccionI;
       Bvar(i)=abs(Bvar(i))
       Hvar(i)=interp1(B,H,Bvar(i));
       u(i)=(Bvar(i)/Hvar(i)); 
       [LMAG(i)]=LM(u(i),Np)
       [Peq7(i),Peq6(i),Ldp(i)]=PermeanciaTotal(u(i),Np) 
       disp([(t(i)),Vm,LMAG(i),Peq7(i)])
       Ia(i)=(Np*fluxmag(i))/LMAG(i);
       Lambda(i)=fluxmag(i)*Np;

        k1 = h*((1/(Np*((Req+Rn)/Rn))*(Vm*cos(W*t(i)+Theta)-((Req*Np*fluxmag(i))/(LMAG(i))))));
        k2 = h*(((1/(Np*((Req+Rn)/Rn))*(Vm*cos(W*(t(i)+h/2)+Theta)-((Req*Np*fluxmag(i))/(LMAG(i))))+0.5*k1)));
        k3 = h*(((1/(Np*((Req+Rn)/Rn))*(Vm*cos(W*(t(i)+h/2)+Theta)-((Req*Np*fluxmag(i))/(LMAG(i))))+0.5*k2)));                    
        k4 = h*(((1/(Np*((Req+Rn)/Rn))*(Vm*cos(W*(t(i)+h)+Theta)-((Req*Np*fluxmag(i))/(LMAG(i))))+k3)));
        fluxmag(i+1) = fluxmag(i) + (k1 + 2*k2 + 2*k3 + k4)/6;
        
       Ia(i+1)=Ia(i);
       Peq7(i+1)=Peq7(i);
       Ldp(i+1)=Ldp(i);
       t(i+1) = t0 + i*h;
       u(i+1)=u(i);
       LMAG(i+1)=LMAG(i);
       Lambda(i+1)=fluxmag(i)*Np;
 
end
end 
[diferencial]=diff(fluxmag,h)
dfluxdt=diferencial;

[diferencial]=diff(Lambda,h)
dfLambda=diferencial;
    
 
select Model
    case 1 then
        Eind=Np*dfluxdt;
        E=Eind;
        I=Ia;
        Irn=Eind/Rn;
        for i=1:N
            Irn(i)=Irn(i)*0;
            Irn(i+1)=Irn(i);
        end  
        [diferencial]=diff(Ia,h);
        didtXm=diferencial; 
        didt=didtXm;
        Eind=(Vm*cos(W*t+Theta))-(Req*I+Xeq*didt);
        Vsum=(Vm*cos(W*t+Theta));       
        
    case 2 then
        Eind=Np*dfluxdt;
        E=Eind;
        Irn=Eind/Rn;
        I=Ia+Irn;
        [diferencial]=diff(Irn,h)
        didtRn=diferencial;
        [diferencial]=diff(Ia,h)
        didtXm=diferencial; 
        didt=didtRn+didtXm;
        Eind=(Vm*cos(W*t+Theta))-(Req*I);
        Vsum=(Vm*cos(W*t+Theta));
end   

endfunction 


function [X,armonicas,amplitud,FFT_50_HARM,FFT_50_HARM_SEC_ZERO,FluxT,FluxTSC,HARM]=ARMONICAS(N,h,Muestra,tk,fhz,Np,lmFlux,SeccionI,Senal,nHarm)

w=2*%pi*fhz;
//      ... ... ... ... ...  funciones internas ... ... ... ... ...  // 
funcprot(0)
function  [Flux,Fmm]=FluxFMM(Voc,Ioc,Np,l,A)
    B=Voc/(sqrt(2)*%pi*60*Np*A);
    Flux=B*A;
    Fmm=Np*Ioc;
endfunction   
//................CURVA DE SATURACION A 127V..........................//  
    
tm=(1:N);
Muestra=Muestra';
A=Muestra(1:N);
X=fft(A',-1);
fs=1/(h);
f=round(fs*(0:(N/2))/N); //:Asociar Frecuencia a un Vector
//*********************************************************************//
//... .... .... .... .... .... .... .... .... .... ... .... .... .... .//
//....Proporciona modulo y argumento en cada espectro de frecuencia....//
//... .... .... .... .... .... .... .... .... .... ... .... .... .... .//
//*********************************************************************//
for i=1:1:202
   ang(i)=atan(imag(X(i)),real(X(i))) 
   //ang(i)=(ang(i)*180)/%pi
   //disp([i,round(f(i)),abs(X(i)),ang(i)])
        
   if i==1 then 
     // disp([i,abs(X(i))/N,ang(i)])
      X(i)=(X(i))/N     
   else
     // disp([i,abs(X(i))*2/N,ang(i)])
      X(i)=(X(i))*2/N     
   end 
   
end

X=X';
ejeX=nHarm;
for i=1:1:ejeX
    Eje_f(i)=f(i);
end

armonicas=(0:nHarm-1);
amplitud=abs(X(:,1:max(size(Eje_f))));

HARM=zeros(nHarm,4);

select Senal
case 1 then
    ConsPi=%pi/2;
    kPi=1;
case 2 then 
   kPi=-1;   
   ConsPi=-%pi/2; 
end

for n=1:nHarm

if  n==1 then

     disp([n-1,round(Eje_f(n,1)),abs(X(n)),atan(imag(X(n)),real(X(n)))])
     HARM(n,1)=(n-1);
     HARM(n,2)=round(Eje_f(n,1));
     HARM(n,3)=abs(X(1,n)');
     HARM(n,4)=kPi*atan(imag(X(n)),real(X(n)))
     
     elseif (modulo(n,2)==0) then
    
     disp([n-1,round(Eje_f(n,1)),abs(X(n)),atan(imag(X(n)),real(X(n)))-ConsPi])
     HARM(n,1)=(n-1);
     HARM(n,2)=Eje_f(n,1);
     HARM(n,3)=abs(X(1,n)');
     HARM(n,4)=kPi*atan(imag(X(n)),real(X(n)))-ConsPi;
    
     
     elseif  (modulo(n,2)==1) then
  
     disp([n-1,round(Eje_f(n,1)),abs(X(n)),atan(imag(X(n)),real(X(n)))+ConsPi])
     HARM(n,1)=(n-1);
     HARM(n,2)=round(Eje_f(n));
     HARM(n,3)=abs(X(1,n)');
     HARM(n,4)=kPi*atan(imag(X(n)),real(X(n)))+ConsPi;
     
end
     
end
//*** ***** ***** ****** ****** ***** ***** ****** ****** ****** *****//
//... .... .... .... ARMÓNICAS PARES E IMPARES ... .... ..... .... ..//
//*** ***** ***** ****** ****** ***** ***** ****** ****** ****** *****//
for i=1:1:nHarm
     for j=1:1:max(size(tk))
     FFT_HARM(j,i)=HARM(i,3)*sin(w*tk(j)*HARM(i,1)+HARM(i,4));
end
end
FFT_50_HARM=zeros(max(size(tk)),1)
for i=1:1:nHarm 
      FFT_50_HARM=FFT_50_HARM+FFT_HARM(:,i)
end
for i=4:3:nHarm
     for j=1:1:max(size(tk))
     HARM_SEC_CERO(j,i)=HARM(i,3)*sin(w*tk(j)*HARM(i,1)+HARM(i,4));
     end
end
FFT_50_HARM_SEC_ZERO=zeros(max(size(tk)),1)
for i=1:1:nHarm
      FFT_50_HARM_SEC_ZERO=FFT_50_HARM_SEC_ZERO+HARM_SEC_CERO(:,i)
end
 
[Flux,Fmm]=FluxFMM(MAGCURVE_127_TA(:,1),MAGCURVE_127_TA(:,2),Np,lmFlux,SeccionI)
FluxT=interp1(Fmm,Flux,FFT_50_HARM'*Np)
FluxTSC=interp1(Fmm,Flux,(FFT_50_HARM'-FFT_50_HARM_SEC_ZERO')*Np)
endfunction

function [HARM_1,HARM_3,HARM_5,HARM_7,HARM_9,HARM_11,HARM_13]=HARM_GRAPH(HARM,w,tk)
    
HARM_1=HARM(2,3)*sin(w*tk*HARM(2,1)+HARM(2,4));
HARM_3=HARM(4,3)*sin(w*tk*HARM(4,1)+HARM(4,4));
HARM_5=HARM(6,3)*sin(w*tk*HARM(6,1)+HARM(6,4));
HARM_7=HARM(8,3)*sin(w*tk*HARM(8,1)+HARM(8,4));
HARM_9=HARM(10,3)*sin(w*tk*HARM(10,1)+HARM(10,4));
HARM_11=HARM(12,3)*sin(w*tk*HARM(12,1)+HARM(12,4));
HARM_13=HARM(14,3)*sin(w*tk*HARM(14,1)+HARM(14,4));

endfunction


t0=-1;
tf=3/60;
Ciclos=3;
N=25200; // Numero de muestras
Vrms=127; //Voltaje de fuente 
fhz=60; // Frecuencia 
Model=2; //Modelo T 1  %pi 2
flux_almacenado=1e-16; //Flujo Almacenado
fluxmag(1)=flux_almacenado; // Estado Inicial    
h=(tf-t0)/N; //Paso de integracion
W=2*%pi*fhz; 
Vm=Vrms*sqrt(2); //Voltaje Maximo  

Theta=0; // Angulo [Fase A]
[B,H]=BH(MAGCURVE_127_TA(:,1),MAGCURVE_127_TA(:,2),Np_127,lmFlux,SeccionI);
BTA=B;
HTA=H;
[fluxmag,Peq7,u,Ia,Ldp,dfluxdt,didt,t,E,Bvar,Irn,Eind,I,Vsum]=Iaeq(Np_127,Ns_127,Reqp,Xeqp,fluxmag,SeccionI,SeccionII,Theta,B,H,h,t0,N,Vm,W,Rnp,Model);
  FluxTAp=fluxmag;
  EindTAp=Eind;    
  IaTAp=I;
  VAp=Vsum;
  IaTApRn=Irn;
  IaTApXm=Ia;

// Muestra t=0 a tf=1/60
 xx=24001
 yy=25200

 scf(1)
 plot(t(xx:yy),FluxTAp(xx:yy),'r');
 title("Flujo Magnético en Unidad Monofásica")
 xlabel("Tiempo, Seg")
 ylabel("Flujo, Wb")
 legend("Flujo Magnético TA ")
 xgrid(35);  
 
 scf(2)
 plot(t(xx:yy),IaTAp(xx:yy),'r');
 title("Corriente de Vacio en Unidad Monofásica ")
 xlabel("Tiempo, Seg")
 ylabel("Io, A")
 legend("Io")
 xgrid(35)
 
 scf(3)
 plot(t(xx:yy),[VAp(xx:yy) EindTAp(xx:yy)]);
 title("Voltaje Primario ")
 xlabel("Tiempo, Seg")
 ylabel("Voltaje Aplicado, V")
 legend("Vp","Eind")
 xgrid(35)
 
 
Nm=400; // Muestras x Ciclo
tm=t(xx:yy);
tk=0:h:3/60;
fhz=60;
w=2*%pi*fhz;
nHarm=16; // 0 1 2 3 ... 15

Muestra=FluxTAp(xx:yy);
Senal=1  //(1 Normal) (2 Desfasada)
[X,armonicas,amplitud,FFT_50_HARM,FFT_50_HARM_SEC_ZERO,FluxT,FluxTSC,HARM]=ARMONICAS(Nm,h,Muestra,tk,fhz,Np_127,lmFlux,SeccionI,Senal,nHarm);
data_1a=X;
data_1b=armonicas;
data_1c=amplitud;
data_1d=FFT_50_HARM';
data_1e=FluxTSC';//Desconexion de Tierra
data_1f=HARM; 

Muestra=IaTAp(xx:yy);
Senal=2 //(1 Normal) (2 Desfasada)
[X,armonicas,amplitud,FFT_50_HARM,FFT_50_HARM_SEC_ZERO,FluxT,FluxTSC,HARM]=ARMONICAS(Nm,h,Muestra,tk,fhz,Np_127,lmFlux,SeccionI,Senal,nHarm);
data_2a=X;
data_2b=armonicas;
data_2c=amplitud;
data_2d=FFT_50_HARM';
data_2e=FluxTSC';//Desconexion de Tierra
data_2f=HARM;   


scf(4)
bar(data_1b,data_1c,'red')
title('Espectro de Frecuencias, Flujo Magnético Unidad Monofasica TA  a Vn')
xlabel('Armónicas')
ylabel('Módulo WB')
legend("Espectro de Frecuencias, Flujo")
xgrid

scf(5)
bar(data_2b,data_2c,'red')
title('Espectro de Frecuencias, Corriente de Vacio en TA ')
xlabel('Armónicas')
ylabel('Módulo A')
legend("Espectro de Frecuencias, Io")
xgrid

scf(6)
plot(t(xx:yy),IaTAp(xx:yy),'b');
plot(t(xx:yy),data_2d(1:Nm*Ciclos),'r')
title("Comparacion de Onda Original y Recuperada ")
xlabel("Tiempo, Seg")
ylabel("Io, A")
legend("Io","Io FTT")
xgrid(35)


 
 
  
  
