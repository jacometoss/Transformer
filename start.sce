getd .;
t0=-30/60;
tf=3/60;
Ciclos=3;
N=13200; // Numero de muestras
Vrms=127; //Voltaje de fuente 
fhz=60; // Frecuencia 
Model=2; //Modelo T 1  %pi 2
flux_almacenado=1e-16; //Flujo Almacenado
fluxmag(1)=flux_almacenado; // Estado Inicial    
h=(tf-t0)/N; //Paso de integracion
W=2*%pi*fhz; 
Vm=Vrms*sqrt(2); //Voltaje Maximo  

Theta=0; // Angulo [Fase A]
[B,H]=bh(MAGCURVE_127_TA(:,1),MAGCURVE_127_TA(:,2),Np_127,lmFlux,SeccionI);
BTA=B;
HTA=H;
[fluxmag,Peq7,u,Ia,dfluxdt,didt,t,E,Bvar,Irn,Eind,I,Vsum]=io(Np_127,Ns_127,Reqp,Xeqp,fluxmag,SeccionI,SeccionII,Theta,B,H,h,t0,N,Vm,W,Rnp,Model);
  FluxTAp=fluxmag;
  EindTAp=Eind;    
  IaTAp=I;
  VAp=Vsum;
  IaTApRn=Irn;
  IaTApXm=Ia;

// Muestra t=0 a tf=1/60
  xx=12001
  yy=13200
  
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
[X,armonicas,amplitud,FFT_50_HARM,FFT_50_HARM_SEC_ZERO,FluxT,FluxTSC,HARM]=harmonic(Nm,h,Muestra,tk,fhz,Np_127,lmFlux,SeccionI,Senal,nHarm);
data_1a=X;
data_1b=armonicas;
data_1c=amplitud;
data_1d=FFT_50_HARM';
data_1e=FluxTSC';//Desconexion de Tierra
data_1f=HARM; 

Muestra=IaTAp(xx:yy);
Senal=2 //(1 Normal) (2 Desfasada)
[X,armonicas,amplitud,FFT_50_HARM,FFT_50_HARM_SEC_ZERO,FluxT,FluxTSC,HARM]=harmonic(Nm,h,Muestra,tk,fhz,Np_127,lmFlux,SeccionI,Senal,nHarm);
data_2a=X;
data_2b=armonicas;
data_2c=amplitud;
data_2d=FFT_50_HARM';
data_2e=FluxTSC';//Desconexion de Tierra
data_2f=HARM;   

Muestra=EindTAp(xx:yy);
Senal=2 //(1 Normal) (2 Desfasada)
[X,armonicas,amplitud,FFT_50_HARM,FFT_50_HARM_SEC_ZERO,FluxT,FluxTSC,HARM]=harmonic(Nm,h,Muestra,tk,fhz,Np_127,lmFlux,SeccionI,Senal,nHarm);
data_3a=X;
data_3b=armonicas;
data_3c=amplitud;
data_3d=FFT_50_HARM';
data_3e=FluxTSC';//Desconexion de Tierra
data_3f=HARM;  


scf(4)
bar(data_1b,data_1c,'red')
title('Espectro de Frecuencias, Flujo Magnético Unidad Monofasica TA  a Vn')
xlabel('Armónicas')
ylabel('Módulo WB')
xgrid(35)

scf(5)
bar(data_2b,data_2c,'red')
title('Espectro de Frecuencias, Corriente de Vacio en TA ')
xlabel('Armónicas')
ylabel('Módulo A')
legend("Espectro de Frecuencias, Io")
xgrid(35)
scf(6)
plot(t(xx:yy),IaTAp(xx:yy),'b');
plot(t(xx:yy),data_2d(1:Nm*Ciclos)','r')
title("Comparacion de vacío ( Original y Recuperada ) ")
xlabel("Tiempo, Seg")
ylabel("Io, A")
legend("Io","Io FTT")
xgrid(35)
d1=t(xx:yy);
d2=FluxTAp(xx:yy);
d3=IaTAp(xx:yy);
d4=EindTAp(xx:yy);
d5=VAp(xx:yy);
ruta=strcat([pwd(),'\result.csv'])
write_csv([d1,d2,d3,d4,d5],ruta)
ruta2=strcat([pwd(),'\harm.csv'])
write_csv([data_2f data_3f],ruta2)

