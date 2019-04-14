*** General Sets ***

  Sets
  T    "indices for times"       / 1 *  400 /
  N    "indices for nodes"       / 1 *   2 /
  G    "indices for generators"  / 1 *   8 /
  D(N) "indices for consumers"   / 1 *   2 /
  Z    "indices for zones"       / 1 *   2 /
  S    "indices for scenarios"   /low_co2, high_co2/
  B    "indices for backup"      / 1 *   2 /
  ;

  Alias (n,nn) ;

  Scalars
  M              bigM                                  / 10000    /
  epsilon        elacticity of demand                  /    -0.10 /
  Year           Hours per year                        /  8760    /
  buFixInv       Annuity per 1 MW backup capacity      / 32000    /
  DSM            Load Shedding costs                   /  3000    /
  buVarInv       Variable cost per MWh for backup      /    79    /
  L_step         Capacity steps for lines              /   0.2 /
  L_cost         Cost for 0.01 line capacity           /   250    /
  ;

*** Line Parameters ***

  Parameters
  lineIsNew(L)   "candidate line"
  lineGamma(L)   "susceptance"
  lineUB(L)      "thermal capacity"
  lineStart(L)   "start node"
  lineEnd(L)     "end node"
  lineFixInv(L)  "line investment cost (candidate lines)"
  lineInter(L) "line is interzone link"

*** Generator Parameters ***

  genIsRES(G)    "renewable generator"                   / 7 1, 8 1 /
  genAtNode(G)   "location (node)"                       / 1 1, 2 1, 3 1, 4 2, 5 2, 6 2, 7 1, 8 2 / 

  genVarInv(G)   "variable cost"               /1 51, 2 50, 3 79, 4 54, 5 50, 6 79, 7 0, 8 0/

  genFixInv(G)   "investment cost"                       / 1 93000, 2 58000, 3 32000, 4 93000, 5 58000, 6 32000, 7 78000, 8 93000 /
  buAtNode(B)    "location (node) of backup"             / 1 1, 2 2 /
  probability(S) "probability for scenario S"
  avail(T,G)     "availability of generators"

*** Demand Parameters ***
  consAtNode(D)   "location (node)"                      / 1 1, 2 2 /
  qPeak(D)        "peak consumption at consumer D"       / 1 0.3, 2 0.7 /
  consObjB(D,T)   "slope demand function"
  consObjA(D,T)   "intercept demand function"
  periodScale(T)  "occurence of scenarios"
  dRef(T)         "reference demand per season"
  pRef(T)         "reference price"
  ConsInZone(D) "consumer in zone" /1 1, 2 1/
  GenInZone(G) "generation in zone"    / 1 1, 2 1, 3 1, 4 1, 5 1, 6 1, 7 1, 8 1 /
  NodeInZone(N) "node in zone" /1 1, 2 1/
  ;

  lineIsNew(L)     = 1 ;
  lineGamma(L)     = 1 ;
  lineUB(L)        = ( L.Val - 1 ) * L_step ;
  lineStart(L)     = 1 ;
  lineEnd(L)       = 2 ;
  lineFixInv(L)    = L_cost * L_step / 0.01 * ( L.Val - 1 ) ;
  lineInter(L) = 0

  ;


*** Read.csv Input Data

$call csv2gdx Data/InputRES_wind.txt id=avail Index=1 Value='(2..9)' UseHeader=Y StoreZero=Y FieldSep=Tab Output=input.gdx
$gdxin input.gdx
$load avail
$gdxin

$call csv2gdx Data/Input_HourlyValues_wind.txt id=periodScale Index=1 Value=2 UseHeader=Y StoreZero=Y FieldSep=Tab Output=input.gdx
$gdxin input.gdx
$load periodScale
$gdxin

$call csv2gdx Data/Input_HourlyValues_wind.txt id=dRef Index=1 Value=3 UseHeader=Y StoreZero=Y FieldSep=Tab Output=input.gdx
$gdxin input.gdx
$load dRef
$gdxin

$call csv2gdx Data/Input_HourlyValues_wind.txt id=pRef Index=1 Value=4 UseHeader=Y StoreZero=Y FieldSep=Tab Output=input.gdx
$gdxin input.gdx
$load pRef
$gdxin

*** Demand Curves ***

  consObjB(D,T) = (-1) * pRef(T) / ( dRef(T) * qPeak(D) * epsilon )      ;
  consObjA(D,T) =  pRef(T) + consObjB(D,T) * dRef(T) * qPeak(D)          ;
