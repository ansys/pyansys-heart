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
     11614     11739     11626     11490     11491     11475     11476     11436
     11437     11403     11404     11414     11372      9614      9615      9775
      9698      9699      9713      9714      9791      9704      9705      9770
      9701      9702      9716      9717      9695      9696      9719      9720
      9684      9685      9722      9723      9617     11371     11514     11515
     11411     11412     11599     11443     11444     11500     11501     11508
     11509     11503     11504     11417     11418     11472     11473     11420
     11421     11547     11408     11409     11579     11439     11440     11558
      7864      7862     11457     11458     11446     11447     11496     11497
     11517     11518     11511     11512     11487     11488     11577     11481
     11482     11484     11485     11493     11494     11463     11464     11583
     11478     11479     11466     11467     11552     11430     11431     11564
     11453     11454     11469     11470     11572     11449     11450     11433
     11434     11560     11460     11461     11427     11428     11423     11424
     28455     17945     28434     28435     28490     28491     28541     28460
     28461     28538     28551     28539     28452     28453     28323     30020
     30103     30159     30154     30155     30150     30113     30087     30088
     30109     30056     30057     30075     30076     30170     30151     30115
     30094     30095     30105     30106     30107     30153     30152     30091
     30092     30072     30073     30081     30082     30138     30156     30097
     30098     30025     28317     28318     28440     28441     28534     28535
     28469     28409     28410     28500     28415     28416     28430     28431
     28400     28401     28537     17162     28519     28941     29048     29049
     29051     29052     29118     29045     29046     29022     29023     29101
     29019     29020     29010     29011     29013     29014     29094     29029
     29030     29016     29017     29032     29033     29035     29036     29117
     29026     29027     29038     29039     29042     29043     29001     29002
     29086     29004     29005     29104     28938     28939     31115     31186
     31187     31183     31184     31262     31180     31181     31177     31178
     31258     31189     31190     31192     31193     31204     31205     31169
     31170     31198     31199     31277     31208     31209     31110     28339
     28437     28438     28427     28428     28418     28419     28397     28398
     28502     28406     28407     28412     28413     28449     28450     28329
     28330                                                                      
*SET_NODE_LIST_TITLE
$#                                                                         title
apex node                                                                       
$#     sid       da1       da2       da3       da4    solver
        12       0.0       0.0       0.0       0.0MECH      
$#      n1        n2        n3        n4        n5        n6        n7        n8
     13460                                                                      
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