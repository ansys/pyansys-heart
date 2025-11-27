*KEYWORD
*SET_PART_LIST_TITLE
$#                                                                         title
myocardium_all                                                                  
$#     sid       da1       da2       da3       da4    solver
         1       0.0       0.0       0.0       0.0MECH      
$#   parts     parts
         1         2
*SET_PART_LIST_TITLE
$#                                                                         title
septum                                                                          
$#     sid       da1       da2       da3       da4    solver
         2       0.0       0.0       0.0       0.0MECH      
$#   parts
         3
*SET_NODE_ADD_TITLE
$#                                                                         title
all_endocardium_segments                                                        
$#     sid       da1       da2       da3       da4    solver
        10       0.0       0.0       0.0       0.0MECH      
$#   nodes     nodes     nodes
         5         7         9
*SET_NODE_ADD_TITLE
$#                                                                         title
all_epicardium_segments                                                         
$#     sid       da1       da2       da3       da4    solver
        11       0.0       0.0       0.0       0.0MECH      
$#   nodes     nodes
         6         8
*SET_NODE_ADD_TITLE
$#                                                                         title
all_but_left_endocardium                                                        
$#     sid       da1       da2       da3       da4    solver
        12       0.0       0.0       0.0       0.0MECH      
$#   nodes     nodes     nodes
         6         8         9
*SET_NODE_LIST_TITLE
$#                                                                         title
base nodes                                                                      
$#     sid       da1       da2       da3       da4    solver       its    unused
        13       0.0       0.0       0.0       0.0MECH                          
$#   nodes     nodes     nodes     nodes     nodes     nodes     nodes     nodes
     13855     14269     12396     12144     11748     11636     18488     17678
     13775     13226     19356     19557     13046     12537     19508     19643
     14612     14301     17262     17518     14540     15115     15585     15958
     19179     18798     18415     18775     19681     19752     15121     14917
     19116     19428     16463     17020     14240     13888     12302     11226
     10942     15707     15328     18824     12016     19658     13130     10931
     11118      9657      9749     11945     12470     10418      9983      9479
      9361      9776      9428      9303      9581     11335     12585     10156
     16645     16318     10287     10440     10576     10644     18013     17009
     13841     13242      6813      7495     14383     14218      7868      8075
     14171     14161     10762     11176      8595      9104      5986      6337
      5660     14406      9744     10387     14159     13968     11730     12810
     13114     10222      9402      8331      7888      7359      6658     10686
      5402      5084      5143      4995      5708      5574     12090     11712
     12907     12525      6143      5879      6213      5245     11371      8884
      8568     13635     12229     13394     13209      1812      1752      1915
      4392      4554      4726      4905       905       804      1015      5050
      5131      4575      5382      5258      2994      2886      3724      3515
      4039      3784      3342      3146      5134      5431      4375      4146
      4456      4355      5558      5501      3481      3437      5125      5443
      4524       586       624      2681      1581       239       250       430
       470       286       321       389       352      1246      1128       195
       176      2334      2160      1517      1454       205       202       564
       530       213       236      2064       334       240       637       468
       276       233      2499      4295      4001      4005      1391      1489
      3611      3755       832       898      1684      1785      1951      2159
      3198      3404      2469      2643      2827      3003      1186       742
      4211      2316      1037       704     18137     18602     17830     15985
     16352     19279     19762     12950     13493     20407     20813     13824
     14499     19521     19299     14122     13554      5740      5768      6634
      7142      9494      9924     11650     12133      8846      8298      7620
      7126      6324      5953     20604     20573     20637     20626     12221
     11364      6017      6222     10035      9457     10254     10661     20327
     20060     16719     16364     19036     18539     12633     20000      5853
     17670     15351      7604      8157      8425      8914      6623     11094
     10607      5801     12970     20839     17341     18653      5786     15656
     17096     14701     20210
*SET_NODE_LIST_TITLE
$#                                                                         title
apex node                                                                       
$#     sid       da1       da2       da3       da4    solver       its    unused
        14       0.0       0.0       0.0       0.0MECH                          
$#   nodes
     40504
*EM_EP_FIBERINITIAL
$#      id    partid     stype     ssid1     ssid2
         1         1         2        13        14
*EM_EP_FIBERINITIAL
$#      id    partid     stype     ssid1     ssid2
         2         1         2        11        10
*EM_EP_FIBERINITIAL
$#      id    partid     stype     ssid1     ssid2
         3         2         2        12         5
*EM_EP_CREATEFIBERORIENTATION
$# partsid   solvid1   solvid2     alpha      beta     wfile    prerun
         1         1         2      -101      -102         1         1
*EM_EP_CREATEFIBERORIENTATION
$# partsid   solvid1   solvid2     alpha      beta     wfile    prerun
         2         1         3      -101      -103         1         1
*DEFINE_FUNCTION
$#     fid                                                               heading
       101                                                                      
$#                                                                      function
float alpha(
            float x_ele, float y_ele, float z_ele,
            float phi_len, float phi_thi)
{ 
  float alpha1;
  float pi;
  float alpha_endo;
  float alpha_epi;
  pi=3.14159265359;
  alpha_endo=-60.00*pi/180;
  alpha_epi=60.00*pi/180;
  alpha1=alpha_endo*(1-phi_thi)+alpha_epi*phi_thi;
  return alpha1;
}
*DEFINE_FUNCTION
$#     fid                                                               heading
       102                                                                      
$#                                                                      function
    float beta(
            float x_ele, float y_ele, float z_ele,
            float phi_len, float phi_thi)
{  
  float beta1;
  float pi;
  float beta_endo;
  float beta_epi;
  pi=3.14159265359;
  beta_endo=-65.00*pi/180;
  beta_epi=25.00*pi/180;
  beta1=beta_endo*(1-phi_thi)+beta_epi*phi_thi;
  return beta1;
}
*DEFINE_FUNCTION
$#     fid                                                               heading
       103                                                                      
$#                                                                      function
    float betaW(
            float x_ele, float y_ele, float z_ele,
            float phi_len, float phi_thi)
{  
  float beta1;
  float pi;
  float beta_endo;
  float beta_epi;
  pi=3.14159265359;
  beta_endo=-65.00*pi/180;
  beta_epi=25.00*pi/180;
  beta1=beta_endo*(1-phi_thi)+beta_endo*phi_thi;
  return beta1;
}
*END