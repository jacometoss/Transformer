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
