***--------------------------------------------------------------------------***
***                           DEFINITION of VARIABLES                        ***
***--------------------------------------------------------------------------***

  Variables
* objective values
  welfareSpot            "welfare in spot market"
  costRedispatch         "cost at redispatch level"
  welfareRedispatch      "welfare at redispatch level"

* Spot Market
  f_sp(L,T)            "trade flow in spot market"
* Redispatch
  f_rd(L,T)            "transmission flows redispatch"
  angle(N,T)           "phase angle in redispatch model"
  ;

  Positive Variables
* Spot Market
  d_sp(D,T)            "demand spot market"
  g_sp(G,T)            "generation amount spot market"
  ig_sp(G)             "installed capacity of generators in spot market"

* Redispatch
  d_rd(D,T)            "demand redispatcht"
  g_rd(G,T)            "generation amount redispatch"
  gb_rd(B,T)           "generation backup capacity redispatch"
  ib_rd(B)             "investment in backup capacity redispatch"
  g_n_rd(G,T)          "negative generation redispatch"
  g_p_rd(G,T)          "positive generation redispatch"
  ls_rd(D,T)           "load shedding redispatch"
  ;

***--------------------------------------------------------------------------***
***                          SPOT MARKET MODEL                               ***
***--------------------------------------------------------------------------***

*** Objective function
  Equation welfSpot;
  welfSpot..         welfareSpot =e= (sum((D,T), periodScale(T)*( consObjA(D,T) * d_sp(D,T)
                                         - 0.5 * consObjB(D,T) * d_sp(D,T) * d_sp(D,T) ))
                                         - sum((G,T), genVarInv(G) * g_sp(G,T) * periodScale(T) ) ) * Year
                                         - sum(G, genFixInv(G) * ig_sp(G) ) ;

*** Zonal First Kirchhoffs Law

  Equation ZFKL;
  ZFKL(Z,T)..

         sum(D$(ConsInZone(D) = Z.val), d_sp(D,T)) =e=
                     sum(G$(sum(N$(genAtNode(G) = N.val), NodeInZone(N)) = Z.val), g_sp(G,T) )
                   - sum(L$(sum(N$(lineStart(L) = N.val), NodeInZone(N)) = Z.val and lineInter(L) = 1), f_sp(L,T))
                   + sum(L$(sum(N$(lineEnd(L) = N.val),   NodeInZone(N)) = Z.val and lineInter(L) = 1), f_sp(L,T)) ;

*** Market Coupling Flow Restrictions

  Equation MCF1;
  MCF1(L,T)$(lineInter(L) = 1 and lineIsNew(L) = 0)..   f_sp(L,T) =l= lineUB(L);
  Equation MCF2;
  MCF2(L,T)$(lineInter(L) = 1 and lineIsNew(L) = 0).. - lineUB(L)=l= f_sp(L,T);
  Equation MCF3;
  MCF3(L,T)$(lineInter(L) = 1 and lineIsNew(L) = 1)..   f_sp(L,T) =l= lineB(L) * lineUB(L);
  Equation MCF4;
  MCF4(L,T)$(lineInter(L) = 1 and lineIsNew(L) = 1).. - lineB(L) * lineUB(L) =l= f_sp(L,T);

***Generation Capacity Limits

  Equation GCLSpot ;
  GCLSpot(G,T)..  g_sp(G,T) =l= avail(T,G) * ig_sp(G) ;


***--------------------------------------------------------------------------***
***                     NETWORK- and REDISPATCH LEVEL                        ***
***--------------------------------------------------------------------------***

  Equation costRed ;
  costRed..         costRedispatch =e=  sum((G,T), GenVarInv(G) * ( g_rd(G,T) - SP_GEN_G(G,T) ) * periodScale(T) ) * YEAR
                                           + sum((B,T), buVarInv * gb_rd(B,T) * periodScale(T) ) * YEAR
                                           + sum(B, buFixInv * ib_rd(B) )
                                           + sum(L$(lineIsNew(L) = 1), lineFixInv(L) * lineB(L)) ;

  Equation welfRed;
  welfRed..         welfareRedispatch =e= (sum((D,T), (consObjA(D,T) * SP_dem(D,T) - 0.5 * consObjB(D,T) * SP_dem(D,T) * SP_dem(D,T)) * periodScale(T))
                                            - sum((G,T), GenVarInv(G) * g_rd(G,T) * periodScale(T))
                                            - sum((B,T), buVarInv * gb_rd(B,T) * periodScale(T))) * Year
                                            - sum(G,genFixInv(G) * SP_CAP_G(G))
                                            - sum(B,buFixInv * ib_rd(B) )
                                            - sum(L$(lineIsNew(L) = 1), lineFixInv(L) * lineB(L));

***First Kirchhoffs Law

  Equation FKL;
  FKL(N,T)..   sum(D$(consAtNode(D) = N.val), SP_DEM(D,T)) =e=
                 sum(G$(genAtNode(G)   = N.val), g_rd(G,T))
                 + sum(B$(buAtNode(B)  = N.val), gb_rd(B,T))
                 + sum(L$(lineEnd(L)   = N.val), f_rd(L,T))
                 - sum(L$(lineStart(L) = N.val), f_rd(L,T)) ;

***Second Kirchhoffs Law

  Equation SKL1;
  SKL1(L,T)$(lineIsNew(L) = 0).. f_rd(L,T) + lineGamma(L) * (sum(N$(lineStart(L) = N.val), angle(N,T)) - sum(N$(lineEnd(L) = N.val), angle(N,T))) =e= 0;
  Equation SKL2;
  SKL2(L,T)$(lineIsNew(L) = 1).. - M * (1 - lineB(L)) =l=  f_rd(L,T) + lineGamma(L) * (sum(N$(lineStart(L) = N.val), angle(N,T)) - sum(N$(lineEnd(L) = N.val), angle(N,T)));
  Equation SKL3;
  SKL3(L,T)$(lineIsNew(L) = 1).. f_rd(L,T) + lineGamma(L) * (sum(N$(N.val = lineStart(L)), angle(N,T)) - sum(N$(N.val = lineEnd(L)), angle(N,T))) =l= M * (1 - lineB(L));

***Voltage Phase Angle

  Equation VPA;
  VPA(N,T)$(N.val = 1).. angle(N,T) =e= 0;

***Trasmission Flow Limits

  Equation TFL1;
  TFL1(L,T)$(lineIsNew(L) = 0)..   f_rd(L,T) =l= lineUB(L);
  Equation TFL2;
  TFL2(L,T)$(lineIsNew(L) = 0).. - lineUB(L) =l= f_rd(L,T);
  Equation TFL3;
  TFL3(L,T)$(lineIsNew(L) = 1)..   f_rd(L,T) =l= lineB(L) * lineUB(L);
  Equation TFL4;
  TFL4(L,T)$(lineIsNew(L) = 1).. - lineB(L) * lineUB(L) =l= f_rd(L,T);

***Generation Capacity Limits (Redispatch Level)

  Equation GCLRed;
  GCLRed(G,T).. g_rd(G,T) =l= avail(T,G) * SP_CAP_G(G) ;

  Equation GCLBu;
  GCLBu(B,T)..  gb_rd(B,T) =l= ib_rd(B) ;

*** Fix Spot Market and Redispatch Quantities

  Equation fixDem;
  fixDem(D,T).. d_rd(D,T) =e= SP_DEM(D,T)
  ;

  Equation fixGen;
  fixGen(G,T).. g_rd(G,T) =e= SP_GEN_G(G,T) + g_p_rd(G,T) - g_n_rd(G,T);

***--------------------------------------------------------------------------***
***                           DEFINITION MODELS                              ***
***--------------------------------------------------------------------------***

  Model Spotmarket
  / welfspot,
    ZFKL,
    MCF1,
    MCF2,
    MCF3,
    MCF4,
    GCLSpot /;

  Model Redispatch
  / welfRed,
*   costRed,
    FKL,
    SKL1,
    SKL2,
    SKL3,
    VPA,
    TFL1,
    TFL2,
    TFL3,
    TFL4,
    GCLRed,
    GCLBu,
    fixDem,
    fixGen /;

***--------------------------------------------------------------------------***
***                           END MODEL SECTION                              ***
***--------------------------------------------------------------------------***
