function [fluxmag,Peq7,u,Ia,dfluxdt,didt,t,E,Bvar,Irn,Eind,I,Vsum]=io(Np,Ns,Req,Xeq,fluxmag,SeccionI,SeccionII,Theta,B,H,h,t0,N,Vm,W,Rn,Model)  
  
        t(1) = t0;

for i = 1:N  

select Model
    
case 1 then
    
       Bvar(i)=(fluxmag(i))/SeccionI;
       Bvar(i)=abs(Bvar(i))
       Hvar(i)=interp1(B,H,Bvar(i));
       u(i)=(Bvar(i)/Hvar(i)); 
       [Peq7(i),Peq6(i),Ldp(i)]=pdtotal(u(i),Np)    
       disp([(t(i)),Vm,Ldp(i),Peq7(i)])
       Ia(i)=(Np*fluxmag(i))/Ldp(i);
       Lambda(i)=fluxmag(i)*Np;
    
        k1 = h*((1/Np)*(Vm*cos(W*t(i)+Theta)-((Req*Np*fluxmag(i))/(Ldp(i)))));
        k2 = h*(((1/Np)*(Vm*cos(W*(t(i)+h/2)+Theta)-((Req*Np*fluxmag(i))/(Ldp(i))))+0.5*k1));
        k3 = h*(((1/Np)*(Vm*cos(W*(t(i)+h/2)+Theta)-((Req*Np*fluxmag(i))/(Ldp(i))))+0.5*k2));                    
        k4 = h*(((1/Np)*(Vm*cos(W*(t(i)+h)+Theta)-((Req*Np*fluxmag(i))/(Ldp(i))))+k3));
        fluxmag(i+1) = fluxmag(i) + (k1 + 2*k2 + 2*k3 + k4)/6;
       
       t(i) = t0 + i*h;
   

case 2 then
    
       Bvar(i)=(fluxmag(i))/SeccionI;
       Bvar(i)=abs(Bvar(i))
       Hvar(i)=interp1(B,H,Bvar(i));
       u(i)=(Bvar(i)/Hvar(i)); 
       [LMAG(i),Peq7(i)]=lm(u(i),Np)
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
       LMAG(i+1)=LMAG(i);
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
