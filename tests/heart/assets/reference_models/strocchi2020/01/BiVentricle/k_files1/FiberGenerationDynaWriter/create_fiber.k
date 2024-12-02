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
     23442      8944     26972     26965     26485     26650     24236     24164
     26953     26888     19313     19314     26916     26889     24201     26930
     13707     24521     22158      1699     21998     21515     22050     19611
     20934     26601     19296     13273     19668      1198     23930     24117
     19740     20842     16670     21575     26999     26954     26464     27004
     26994     13871     23933     24234     26490     13776     23441     26736
     26767     26788     26772     26782     26781     26797     16748     14853
     26764      4033     26725     21491     26766     26742     17218     25990
     25952      1262     26792     19383     26796     26768     22238     23806
     26616     26004     12972     18938     26552     24802     14511     24629
     25148     25147     15123     15122     25995     24014     13581     20105
     25448     25398     24508     23481     24651     25951     26049     24102
     23988     21169     18889     25796     13092      4738      4737     18790
     21983      1501     13105     18662     20433     26421     20279     25184
     26901     20281     11283     11284     18170     26424     26804     22355
     19253     23989     13710     23942     24269     39980     38933     39923
     30130     36288     36261     35100     34865     30326     30994     35102
     31152     31151     35343     34836     36392     36050     35037     35455
     30887     31853     38655     36321     35252     35128     34105     27964
     31140     35731     34652     35381     32323     40952     34476     35358
     32437     37589     36434     32144     38932     33927     31501     28538
     33846     31131     32310     38267     35643     28239     30835     36182
     31751     39614     41996     41477     30643     33225     39087     32879
     32936     39596     29167     33683     33763     31245     32876     30961
     32271     28296     37659     32500     31316     29790     32483     27179
     29656     36260     31187     30185     27180     34196     35287     31543
     28125     33137     34202     35333     34848     33149     31849     31180
     27902     34203     31142     31845     36387     39392     30315     32448
     32332     38919     33447     39010     28787     32380     33458     30221
     40991     38347     35804     33167     33166     33961     36561     40145
     39116     31481     36257     40959     38634     31171     31172     31510
     38031     41309     28198     34080     41866     41833     33904     33153
     28950     30336     39417     32902     28275     29150     40533     33985
     35958     37541     37050     37876     30529     30017     39743     36625
     30530     36485     34822     38752     37098     32539     38382     37658
     32901     32495     36877     33814     32333     40873     41909     36875
     32816     36957     33184
*SET_NODE_LIST_TITLE
$#                                                                         title
apex node                                                                       
$#     sid       da1       da2       da3       da4    solver       its    unused
        14       0.0       0.0       0.0       0.0MECH                          
$#   nodes
     18916
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