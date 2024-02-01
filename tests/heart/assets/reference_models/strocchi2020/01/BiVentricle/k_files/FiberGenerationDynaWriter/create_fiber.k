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
         8       0.0       0.0       0.0       0.0MECH      
$#   nodes     nodes     nodes
         3         5         7
*SET_NODE_ADD_TITLE
$#                                                                         title
all_epicardium_segments                                                         
$#     sid       da1       da2       da3       da4    solver
         9       0.0       0.0       0.0       0.0MECH      
$#   nodes     nodes
         4         6
*SET_NODE_ADD_TITLE
$#                                                                         title
all_but_left_endocardium                                                        
$#     sid       da1       da2       da3       da4    solver
        10       0.0       0.0       0.0       0.0MECH      
$#   nodes     nodes     nodes
         4         6         7
*SET_NODE_LIST_TITLE
$#                                                                         title
base nodes                                                                      
$#     sid       da1       da2       da3       da4    solver
        11       0.0       0.0       0.0       0.0MECH      
$#      n1        n2        n3        n4        n5        n6        n7        n8
      2420      2514      2515      2526      2527      2496      2497      2529
      2530      2486      2487      2634      2517      2518      2511      2512
      2535      2536      2502      2503      2505      2506      2468      2469
      2474      2475      2480      2481      2607      2508      2509      2490
      2491      2547      2548      2553      2554      2632      2419        27
       151       152       148       149       165        17      1439      1317
      2424      2538      2539      2635      2541      2542      2477      2478
      2520      2521      2617      2464      2465      2461      2462      2493
      2494      2601      2499      2500      2556      2557      2471      2472
      2550      2551      2627      2544      2545      2456      2457      2422
      1329      1407      1408      1485      1426      1427      1436      1437
      1404      1405      1414      1415      1423      1424      1410      1411
      1417      1418      1480      1401      1402      1429      1430      1398
      1399      1482      1386      1387     10521     10520     10522     10409
     10410     10328     10329     10355     10356     10481     10367     10368
     10385     10386     10489     10370     10371     10406     10407     10376
     10377     10373     10374     10391     10392     10397     10398     10400
     10401     10388     10389     10394     10395     10508     10304     10305
     10319     10320     10473     10343     10344     10358     10359     10346
     10347     10423     10424     10454     10455     10379     10380     10307
     10308     10458     10334     10335     10364     10365     10259     10260
     10258     11934     11894     11895     11909     11910     11879     11880
     11876     11877     11981     11900     11901     11863     11864     11983
     11882     11883     11994     11873     11874     11970     11903     11904
     11897     11898     11995     11885     11886     11860     11861     11848
     11849     11906     11907     11737     13940     13976     13977     13987
     13988     14038     14039     14095     14032     14033     14035     14036
     14041     14042     14025     14026     14022     14023     13984     13985
     13990     13991     13994     13995     14015     14016     14000     14001
     14018     14019     14003     14004     14085     14006     14007     14029
     14030     14076     14009     14010     13980     13981     14012     14013
     13997     13998     14046     11730     11731                              
*SET_NODE_LIST_TITLE
$#                                                                         title
apex node                                                                       
$#     sid       da1       da2       da3       da4    solver
        12       0.0       0.0       0.0       0.0MECH      
$#      n1        n2        n3        n4        n5        n6        n7        n8
      6230                                                                      
*EM_EP_FIBERINITIAL
$#      id    partid     stype     ssid1     ssid2
         1         1         2        11        12
*EM_EP_FIBERINITIAL
$#      id    partid     stype     ssid1     ssid2
         2         1         2         9         8
*EM_EP_FIBERINITIAL
$#      id    partid     stype     ssid1     ssid2
         3         2         2        10         3
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