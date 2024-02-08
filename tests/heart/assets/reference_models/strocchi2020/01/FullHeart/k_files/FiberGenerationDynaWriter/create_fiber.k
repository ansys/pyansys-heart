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
      7862      7863     11607     11712     11621     11622     11486     11487
     11488     11466     11467     11550     11409     11410     11396     11397
     11518      9615      9616      9780      9698      9699      9713      9714
      9796      9704      9705      9775      9701      9702      9716      9717
      9695      9696      9719      9720      9684      9685      9722      9723
      9618     11364     11510     11511     11406     11407     11592     11436
     11437     11497     11498     11504     11505     11500     11501     11413
     11414     11481     11482     11416     11417     11403     11404     11432
     11433     11439     11440     11401      7865     11453     11454     11442
     11443     11493     11494     11513     11514     11507     11508     11484
     11485     11571     11475     11476     11478     11479     11490     11491
     11459     11460     11575     11472     11473     11462     11463     11546
     11426     11427     11585     11448     11449     11469     11470     11566
     11445     11446     11429     11430     11554     11456     11457     11423
     11424     11419     11420     28453     17920     28432     28433     28488
     28489     28539     28458     28459     28536     28549     28537     28450
     28451     28321     30018     30101     30157     30152     30153     30148
     30111     30085     30086     30107     30054     30055     30073     30074
     30168     30149     30113     30092     30093     30103     30104     30105
     30151     30150     30089     30090     30070     30071     30079     30080
     30136     30154     30095     30096     30023     28315     28316     28438
     28439     28532     28533     28467     28407     28408     28498     28413
     28414     28428     28429     28398     28399     28535     17137     28517
     28939     29046     29047     29049     29050     29116     29043     29044
     29020     29021     29099     29017     29018     29008     29009     29011
     29012     29092     29027     29028     29014     29015     29030     29031
     29033     29034     29115     29024     29025     29036     29037     29040
     29041     28999     29000     29084     29002     29003     29102     28936
     28937     31113     31184     31185     31181     31182     31260     31178
     31179     31175     31176     31256     31187     31188     31190     31191
     31202     31203     31167     31168     31196     31197     31275     31206
     31207     31108     28337     28435     28436     28425     28426     28416
     28417     28395     28396     28500     28404     28405     28410     28411
     28447     28448     28327     28328                                        
*SET_NODE_LIST_TITLE
$#                                                                         title
apex node                                                                       
$#     sid       da1       da2       da3       da4    solver
        12       0.0       0.0       0.0       0.0MECH      
$#      n1        n2        n3        n4        n5        n6        n7        n8
     13458                                                                      
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