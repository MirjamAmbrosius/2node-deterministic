  Parameters

*welfare
  wf_sp                  "spot market welfare"
  wf_all                 "final welfare"

*generation and demand
  SP_DEM(D,T)          "demand Spot"
  SP_GEN_G(G,T)        "generation amount Spot"
  RD_GEN_G(G,T)        "generation after Redispatch"
  RD_GEN_B(B,T)        "generation of backup capacity b in scenario s and time t"
  lineB(L)
  SP_FLOW(L,T)
  SP_CAP_L(L)
  SP_CAP_G(G)
  RD_CAP_B(B)

  Loop_welfare_all(LineInvest)          "total welfare"

  Loop_genInv(LineInvest, G)           "generation investment"
  Loop_lineInv(LineInvest)              "cost of line investement"
  total_generation               "total generation by private firms"
  total_bu_generation            "total bu generation"
  total_spot_generation          "total spot generation"
;
