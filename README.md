# Transformador
Simulación Transformador Monofásico 127 V
Rutina: Transformador Monofásico 127V SIN CARGA 
Autor: Marco Polo Jacome Toss 
Version : 0.01                                           
Plataforma : Scilab 6.0 (https://www.scilab.org)      
Fecha : 2018.06.10          

## Funcion : harmonic
El archivo require el paso, numero de muestras, tiempo de simulación, frecuencia, Vueltas del bobinado, longitud media de la seccion del circuito magnético, sección de las piernas, seccion de la pierna central.

Las variables de salida 
* X descomposicion de armonicas pares e impares.
* Armónicas 
* Amplitud de armónicas


```scilab

function [X,armonicas,amplitud,FFT_50_HARM,FFT_50_HARM_SEC_ZERO,FluxT,FluxTSC,HARM]=harmonic(N,h,Muestra,tk,fhz,Np,lmFlux,SeccionI,Senal,nHarm)

w=2*%pi*fhz;
    
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
 
[Flux,Fmm]=fluxfmm(MAGCURVE_127_TA(:,1),MAGCURVE_127_TA(:,2),Np,lmFlux,SeccionI)
FluxT=interp1(Fmm,Flux,FFT_50_HARM'*Np)
FluxTSC=interp1(Fmm,Flux,(FFT_50_HARM'-FFT_50_HARM_SEC_ZERO')*Np)
endfunction

```
