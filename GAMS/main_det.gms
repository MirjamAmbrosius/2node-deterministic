***--------------------------------------------------------------------------***
***                             GENERAL OPTIONS                              ***
***--------------------------------------------------------------------------***

option  optcr = 0.0001
        limrow = 0,
*equations listed per block */
        limcol = 0
*variables listed per block */
        solprint = off,
*solver's solution output printed */
        sysout = off,
*define standard solver
        QCP = Gurobi,
        LP = Gurobi
;

***--------------------------------------------------------------------------***
***            OPTIONS FOR DIFFERENT SCENARIOS & LINE INVESTMENT             ***
***--------------------------------------------------------------------------***

Sets
         L "indices for power lines"     / 1 * 5 /
         LineInvest                      / 1 * 5 /
         Loop_Probability                / 1  /
         results                         / "CS", "PS", "CR", "ALL"
                                           "RD_G", "RD_B","RD_L", "SP_G", "SP_P"
                                           "C_L", "C_B" , "C_G",
                                           "p_sp", "p_rd", "p_cr" /  ;

***--------------------------------------------------------------------------***
***             LOAD DATA AND SETUP FOR LOOP WITH PROBABILITIES              ***
***--------------------------------------------------------------------------***

$include input_det.gms
$include parameters_det.gms
$include model_det.gms

*** read gurobi.opt
*  Spotmarket.OptFile = 1 ;
*  RedispatchWelfare.OptFile = 1 ;

*** time after whcih the solver terminates:
 Spotmarket.reslim = 10000;
 Redispatch.reslim = 36000;

 Alias(LineInvest,LineInvest2) ;

***--------------------------------------------------------------------------***
***           START MODEL LOOP FOR PROBABILITIES AND LINE INVEST             ***
***--------------------------------------------------------------------------***


     Loop(LineInvest,

       lineB(L) = 1$(LineInvest.val=L.val);

***--------------------------------------------------------------------------***
***                        SOLVE SPOT MARKET MODEL                           ***
***--------------------------------------------------------------------------***

  option clear= welfareSpot ;
  option clear= d_sp        ;
  option clear= g_sp        ;
  option clear= ig_sp       ;
  option clear= f_sp        ;

  SOLVE Spotmarket USING QCP MAXIMIZE welfareSpot ;

  SP_DEM(D,T)   = d_sp.l(D,T)          ;
  SP_GEN_G(G,T) = g_sp.l(G,T)          ;
  SP_CAP_G(G)     = ig_sp.l(G)             ;
  SP_FLOW(L,T)  = f_sp.l(L,T)          ;
  SP_CAP_L(L)     = lineB(L) * lineUB(L)   ;


$ontext
  option clear= welfareSpot ;
  option clear= d_sp        ;
  option clear= g_sp        ;
  option clear= ig_sp       ;
  option clear= f_sp        ;
$offtext

***--------------------------------------------------------------------------***
***                         SOLVE REDISPATCH MODEL                           ***
***--------------------------------------------------------------------------***


  option clear = costRedispatch ;
  option clear = f_rd   ;
  option clear = angle  ;
  option clear = d_rd   ;
  option clear = g_rd   ;
  option clear = gb_rd  ;
  option clear = ib_rd  ;
  option clear = g_n_rd ;
  option clear = g_p_rd ;

*  SOLVE Redispatch USING LP MINIMIZE costRedispatch ;
   SOLVE Redispatch USING LP MAXIMIZE welfareRedispatch;

  RD_GEN_G(G,T) = g_rd.l(G,T)  ;
  RD_CAP_B(B)     = ib_rd.l(B)     ;
  RD_GEN_B(B,T) = gb_rd.l(B,T) ;


$ontext
  option clear = costRedispatch ;
  option clear = f_rd   ;
  option clear = angle  ;
  option clear = d_rd   ;
  option clear = g_rd   ;
  option clear = gb_rd  ;
  option clear = ib_rd  ;
  option clear = g_n_rd ;
  option clear = g_p_rd ;
$offtext

***--------------------------------------------------------------------------***
***                        CALCULATION OF RESULTS                            ***
***--------------------------------------------------------------------------***


*Welfare after redispatch
  wf_all    = ( sum((D,T), ( consObjA(D,T) * SP_dem(D,T) - 0.5 * consObjB(D,T) * SP_dem(D,T) * SP_dem(D,T) ) * periodScale(T) )
                                   - sum((G,T), GenVarInv(G) * RD_GEN_G(G,T) * periodScale(T) )
                                   - sum((B,T), buVarInv * RD_GEN_B(B,T) * periodScale(T) )) * Year
                                     - sum(G,genFixInv(G) * SP_CAP_G(G) )
                                     - sum(B,buFixInv * RD_CAP_B(B) )
                                     - sum(L$SP_CAP_L(L), lineFixInv(L) )
;
  total_generation = sum((G,T), RD_GEN_G(G,T));
  total_bu_generation = sum((B,T),RD_GEN_B(B,T));
  total_spot_generation = sum((G,T),SP_GEN_G(G,T));



***--------------------------------------------------------------------------***
***                       RESULTS to LOOP-PARAMETER                          ***
***--------------------------------------------------------------------------***

  Loop_welfare_all(LineInvest)         = wf_all        ;

  Loop_genInv(LineInvest, G)           = SP_CAP_G(G)   ;
  Loop_lineInv(LineInvest)              = sum(l, SP_CAP_L(L) ) ;

***--------------------------------------------------------------------------***
***                     CLEAR PARAMETERs OF MODEL RUN                        ***
***--------------------------------------------------------------------------***

* Clear Spot Resuls
  option clear= SP_DEM           ;
  option clear= SP_GEN_G         ;
  option clear= SP_CAP_G         ;
  option clear= SP_FLOW          ;
  option clear= SP_CAP_L         ;
*  option clear= wf_all           ;

  );

***--------------------------------------------------------------------------***
***                            END OF MODEL LOOP                             ***
***--------------------------------------------------------------------------***

***--------------------------------------------------------------------------***
***                 OUTPUT WITH RESULTS FOR BEST LINE INVEST                 ***
***--------------------------------------------------------------------------***

  Parameter
  Results_genInv(G)
  Results_lineInv

  maxWelfare
  Results_welfare_all
  ;


  Loop(LineInvest,

    maxWelfare$(Loop_welfare_all(LineInvest)=smax(LineInvest2, Loop_welfare_all(LineInvest2) )) = LineInvest.val             ;

  );

  Results_genInv(G)                     = sum(LineInvest$(ord(LineInvest)=maxWelfare), Loop_genInv( LineInvest, G) )         ;
  Results_lineInv                      = sum(LineInvest$(ord(LineInvest)=maxWelfare), Loop_lineInv( LineInvest) )           ;

  Results_welfare_all                  = sum(LineInvest$(ord(LineInvest)=maxWelfare), Loop_welfare_all(LineInvest) )        ;

*$include OutputWriter.gms
