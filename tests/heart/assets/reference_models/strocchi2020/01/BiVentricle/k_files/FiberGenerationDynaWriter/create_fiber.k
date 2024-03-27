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
       165        20      2398      2494      2495      2506      2507      2479
      2480      2509      2510      2466      2467      2615      2497      2498
      2491      2492      2515      2516      2482      2483      2485      2486
      2447      2448      2454      2455      2460      2461      2588      2488
      2489      2470      2471      2527      2528      2533      2534      2579
      2397        27       151       152       148       149      1407      1312
      2401      2518      2519      2599      2524      2525      2457      2458
      2500      2501      2601      2444      2445      2441      2442      2473
      2474      2584      2476      2477      2536      2537      2450      2451
      2530      2531      2611      2521      2522      2433      2434      2399
      1323      1409      1410      1484      1424      1425      1431      1432
      1403      1404      1415      1416      1421      1422      1412      1413
      1418      1419      1481      1400      1401      1427      1428      1387
      1388      1383      1384      1406     10287     10288     10323     10324
     10329     10330     10443     10314     10315     10495     10492     10427
     10368     10369     10332     10333     10398     10311     10312     10299
     10300     10353     10354     10473     10493     10494     10411     10412
     10383     10384     10296     10297     10326     10327     10455     10338
     10339     10356     10357     10482     10341     10342     10380     10381
     10347     10348     10344     10345     10362     10363     10371     10372
     10374     10375     10359     10360     10365     10366     10457     10275
     10276     11929     11888     11889     11903     11904     11873     11874
     11870     11871     11976     11894     11895     11857     11858     11979
     11876     11877     11990     11867     11868     11965     11897     11898
     11891     11892     11967     11879     11880     11854     11855     11842
     11843     11900     11901     11731     13938     13974     13975     13985
     13986     14036     14037     14093     14030     14031     14033     14034
     14039     14040     14023     14024     14020     14021     13982     13983
     13988     13989     13992     13993     14013     14014     13998     13999
     14016     14017     14001     14002     14083     14004     14005     14027
     14028     14074     14007     14008     13978     13979     14010     14011
     13995     13996     14044     11724     11725                              
*SET_NODE_LIST_TITLE
$#                                                                         title
apex node                                                                       
$#     sid       da1       da2       da3       da4    solver
        12       0.0       0.0       0.0       0.0MECH      
$#      n1        n2        n3        n4        n5        n6        n7        n8
      6619                                                                      
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