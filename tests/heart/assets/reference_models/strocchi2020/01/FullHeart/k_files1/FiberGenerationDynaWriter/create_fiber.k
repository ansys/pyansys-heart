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
         7
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
     30205     30604     38655     38644     32439     32752     30831     31230
     38015     38366     33841     34349     38533     38672     37517     37786
     35009     35646     31745     33276     38753     36978     36091     36480
     35721     35513     30156     30162     37079     36294     31819     31211
     28555     28859     33363     33209     29595     29120     27594     28023
     26759     26454     32307     32118     34899     34438     27934     27686
     27141     26745     30817     30561     28890     28374     38342     38002
     26521     26313     37757     37740     33113     38311     28995     26863
     27195     33785     26427     29577     32552     32123     26627     26621
     26795     29898     21031     20372     21728     28365     27786     22870
     22151     32591     32104     20217     20319     24362     24068     22829
     23792     27600     28105     32421     32361     22027     22346     25273
     25895     20071     20429     29768     29467     31017     30439     24230
     24732     32554     32648     30618     31090     22285     28894     28927
     29323     27226     31519     32094     26746     24149     24985     31461
     23928     26728     21419     20799     31783     21765     22082     25906
     29943     31801     29874     10359     10415     10025     14405     14189
     20167     20513     20156     19794     19572     19107     15729     16324
     10914     11243     19022     18729     11794     12102     17339     17786
     16458     16951     14252     14852     20654     20629     16005     15485
     18385     18073     18910     19255     16910     16648     20638     20597
     18243     18499     12644     13219     17277     17129     20847     19962
     15473     17248     18883     17598     14939     13774     19622     19878
      9746      7681      7964      7317      7355      8651      8448      9479
      9149      8208      7929     11312     11036      7251      7117      8466
      8804      7141      7127     10146      9674      7230      7380     13612
     13089      7803      7670      7467      7568      7337     12162     11847
      9042      9203      9070      7126     12683      8313      7387      8219
     10576      9082     23421     24115     21733     21569     21335     21416
     21500     21830     22102     22779     22224     40243     39564     39581
     32959     32124     36430     36458     36999     37324     37660     39018
     39635     39665     39968     27389     28024     39674     39845     38306
     38829     29892     30467     26105     26575     38583     38321     24863
     25314     24931     24150     34253     33620     35487     35173     30851
     29965     37932     37431     33993     35268     30774     31569     31346
     32038     23371     23645     36562     36997     24242     39038     25591
     28566     29284     35833     22792     29051     28103     27561     35997
     26471     33204     21455
*SET_NODE_LIST_TITLE
$#                                                                         title
apex node                                                                       
$#     sid       da1       da2       da3       da4    solver       its    unused
        14       0.0       0.0       0.0       0.0MECH                          
$#   nodes
     59705
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