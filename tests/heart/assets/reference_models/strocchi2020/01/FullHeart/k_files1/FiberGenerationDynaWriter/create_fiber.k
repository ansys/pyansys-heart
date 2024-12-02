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
     22768     11469     12516     20572      7042      7044     10519      4913
      8843     21382     18057      6798     24026     19144     20651      6084
     22872     24139      9837     11703     12515     25315     20519     15781
     20776      8572     14838     23751     13128     23024     22310     22309
     17387     19340     22962     10365     22521     23349     23452     21816
     26404     23155     22603     22401       118       201     26441     26430
     23807     23805     22475     22399     23616     26419      5342      5343
     25891      2240       256      5312     14701      1769      6995     22795
     22797     22562     26079     23908     23355     22314     25836     21927
     21946     22520      6496      4521      1917     26012     25479      8362
     25541      5916     20745      8143     26129     26018     26019     25885
     22325      7709     23363     23244     23168     19830     17203     22953
      6118     25255      5412     26442     25354     22331      8226      8569
     22954     23021     18499     21921     12081      8170     25246     25916
     23365     24201     16517      8737      6117      4446     25990     23002
      4448     24200      8174     25544     23003      8243      9209     23777
      8240      5494     23792     30051     34497     36006     33357     33355
     27254     27273     35044     32187     33174     33119     31644     29885
     35273     38593     27692     27690     40995     40667     37179     37180
     29023     28575     34841     32966     41111     38646     32634     30588
     33586     40043     37985     37768     34181     32633     33600     33695
     36724     37249     40668     37883     38301     34182     40195     34704
     33377     40135     36062     34139     32790     33891     38286     38303
     30072     35951     34683     32201     33456     33753     27864     38464
     36144     33354     33027     34971     34998     38234     32830     28042
     38789     36014     34116     36022     38463     31685     31687     38612
     38609     29405     30402     34647     34008     30362     38552     38528
     39700     36335     33754     33151     35244     27873     33572     32833
     34950     36334     34270     34263     33560     34579     29846     29842
     36142     41175     34352     27731     35038     30116     31245     30047
     32304     31935     32417     27425     27426     28450     35179     33833
     27611     27558     30317     36285     34422     32320     33450     34506
     34505     34490     31648     35531     34226     26752     30101     36263
     38648     35385     38742     30048     32302     30268     29367     31911
     34511     30102     30453     27984     27547     28377     34545     34379
     36561     38494     39300     36560     32124     34641     26469     35384
     34484     34486     30267     38493     34479     34476     32617     32418
     32618     34856     36454
*SET_NODE_LIST_TITLE
$#                                                                         title
apex node                                                                       
$#     sid       da1       da2       da3       da4    solver       its    unused
        14       0.0       0.0       0.0       0.0MECH                          
$#   nodes
     25299
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