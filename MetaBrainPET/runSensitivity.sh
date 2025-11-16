#MODIFIED TO PERFORM SIMULATIONS FOR DIFFERENT CRYSTAL TYPES
#./runnew.sh
#!/bin/bash
#set -x

#Gate -a [crystal,LYSO-Proteus][axial_pos,-5] mainMacro-Sensitivity.mac&
# source ./runSensitivity.sh > outputSensitivity.txt &

for Material in  BGO_EJ232_20 #LYSO-Proteus_BaF2_30 BGO_BaF2_20 BGO_BaF2_25 BGO_EJ232_20 BGO_EJ232_25 LYSO-Proteus_EJ232_25 LYSO-Proteus_EJ232_30 LYSO-Proteus_BaF2_25 LYSO-Proteus_BaF2_30 LYSO-Proteus_BaF2_35 BGO_BaF2_15 LYSO-Proteus BGO
do
 echo "Crystal material" is $Material
 for z in -75 -60 -45 -30 -15 0 15 30 45 60 75
 do
  echo "axial_pos: $z"
  mkdir -p output/Sensitivity/$Material/$z
  Gate -a [crystal,$Material][axial_pos,$z] mainMacro-Sensitivity.mac&
 done
done
