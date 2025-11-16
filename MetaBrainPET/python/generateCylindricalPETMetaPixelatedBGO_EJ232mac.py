# Python script to generate a GATE macro to describe a PET scanner based on monolithic crystal
import math

## User section
##Parameters to be edited by the user
## The unit of length is always mm

#distance from the center of the field of view to the crystal surface
CFOV = 145.00
#SiPM array dimensions
SiPMPixSizeZ = 2.90
SiPMPixSizeY = 3.00
SiPMPitchSizeZ = 3.60
SiPMPitchSizeY = 3.60
nSiPMZ = 8
nSiPMY = 8
#crystal dimensions
CSizeZ = SiPMPitchSizeZ*nSiPMZ-(SiPMPitchSizeZ-SiPMPixSizeZ)
CSizeY = SiPMPitchSizeY*nSiPMY-(SiPMPitchSizeY-SiPMPixSizeY)
CSizeX = 22.7 #crystal thickness
#metapixel dimensions
MetaPixelZBGO = 0.3
MetaPixelZEJ232 = 0.1
NMetaPixelZBGO = 7
NMetaPixelZEJ232 = 8
#module array
NModZ = 6 #number of rings
ModPitchSizeZ = CSizeZ + 1.2
NModY = 1
ModPitchSizeY = CSizeY + 1.2
#ring or rsector
#NrsectorPerRing = int(math.floor(math.pi/math.atan(CSizeY/(2*CFOV)))) #number of heads per ring (automatic method)
NrsectorPerRing = 30 #number of heads per ring (manual method)
#optical grease thickness
OptCouplingThick = 0.01
#silicon thickness
SiThick = 0.40
# SiPM epoxy cover
SiPMCoverThick = 0.26
# Module Cover Thickness
ModCoverThick = 2.00
##Do not edit the parameters beyond this line!!!
##End of the user section

##GATE parameters section
##Do not edit the following parameters!!!

#cylindricalPET
cylindricalPETRmax = CFOV + CSizeX + 20.00
cylindricalPETRmin = CFOV - ModCoverThick
cylindricalPETHeight = NModZ*ModPitchSizeZ
#rsector
Xrsector = CSizeX + 2*ModCoverThick
Yrsector = NModY*ModPitchSizeY
Zrsector = NModZ*ModPitchSizeZ
TranslationXrsector = cylindricalPETRmin + Xrsector/2
#module
Xmodule = CSizeX + 2*ModCoverThick
Ymodule = ModPitchSizeY
Zmodule = ModPitchSizeZ
#crystal size
Xcrystal = CSizeX
Ycrystal = SiPMPixSizeY
Zcrystal = SiPMPixSizeZ
#optical coupling
Xoptcoup = OptCouplingThick
Yoptcoup = CSizeY
Zoptcoup = CSizeZ
TranslationXoptcoup = CSizeX/2+Xoptcoup/2
# SiPM epoxy cover
XSiPMCover = SiPMCoverThick
YSiPMCover = SiPMPixSizeY
ZSiPMCover = SiPMPixSizeZ
TranslationXSiPMCover = CSizeX/2 + XSiPMCover/2 + Xoptcoup
# SiPM array
XSiPM = SiThick
YSiPM = SiPMPixSizeY
ZSiPM = SiPMPixSizeZ
TranslationXSiPM = CSizeX/2 + XSiPM/2 + XSiPMCover + Xoptcoup

f = open('cylindricalPETMetapy.mac','w')
f.write('# This macro describes a PET scanner based on pixelated crystal detector blocks')
f.write('\n' + '# each block consists of an array of ' + str(nSiPMZ) + ' x ' + str(nSiPMY) + ' elements with ' + str(Ycrystal) + 'x' + str(Zcrystal) + ' x ' + str(Xcrystal) + ' mm3 crystal coupled to a SiPM array')
f.write('\n' + '# of ' + str(nSiPMZ) + ' x ' + str(nSiPMY) + ' elements with '  + str(SiPMPixSizeZ) + ' x ' + str(SiPMPixSizeZ) + ' mm2 each one.' )
f.write('\n' + '# The distance from the center of the field of view to the crystal surface is ' + str(CFOV) + ' mm.')
f.write('\n')
f.write('\n' + '#### Using cylindricalPET system ####')
f.write('\n' + '/gate/world/daughters/name cylindricalPET')
f.write('\n' + '/gate/world/daughters/insert cylinder')
f.write('\n' + '/gate/cylindricalPET/setMaterial G4_AIR')
f.write('\n' + '/gate/cylindricalPET/placement/setTranslation  0.0 0.0 0.0 mm')
f.write('\n' + '/gate/cylindricalPET/geometry/setRmax ' + str(cylindricalPETRmax) + ' mm')
f.write('\n' + '/gate/cylindricalPET/geometry/setRmin ' + str(cylindricalPETRmin) + ' mm')
f.write('\n' + '/gate/cylindricalPET/geometry/setHeight ' + str(cylindricalPETHeight) + ' mm')
f.write('\n' + '/gate/cylindricalPET/vis/forceWireframe 1')
f.write('\n' + '/gate/cylindricalPET/vis/setVisible 0')
f.write('\n')
f.write('\n' + '/gate/cylindricalPET/daughters/name rsector')
f.write('\n' + '/gate/cylindricalPET/daughters/insert box')
f.write('\n' + '/gate/rsector/placement/setTranslation ' + str(TranslationXrsector) + ' 0. 0. mm')
f.write('\n' + '/gate/rsector/geometry/setXLength ' + str(Xrsector) + ' mm')
f.write('\n' + '/gate/rsector/geometry/setYLength ' + str(Yrsector) + ' mm')
f.write('\n' + '/gate/rsector/geometry/setZLength ' + str(Zrsector) + ' mm')
f.write('\n' + '/gate/rsector/setMaterial G4_AIR')
f.write('\n' + '/gate/rsector/vis/forceWireframe 1')
f.write('\n' + '/gate/rsector/vis/setVisible 0')
f.write('\n')
f.write('\n' + '/gate/rsector/daughters/name module')
f.write('\n' + '/gate/rsector/daughters/insert box')
f.write('\n' + '/gate/module/geometry/setXLength ' + str(Xmodule) + ' mm')
f.write('\n' + '/gate/module/geometry/setYLength ' + str(Ymodule) + ' mm')
f.write('\n' + '/gate/module/geometry/setZLength ' + str(Zmodule) + ' mm')
f.write('\n' + '/gate/module/setMaterial PTFE')
f.write('\n' + '/gate/module/vis/forceWireframe')
f.write('\n' + '/gate/module/vis/setVisible 0')
f.write('\n')
f.write('\n' + '/gate/module/daughters/name crystal')
f.write('\n' + '/gate/module/daughters/insert box')
f.write('\n' + '/gate/crystal/geometry/setXLength ' + str(Xcrystal) + ' mm')
f.write('\n' + '/gate/crystal/geometry/setYLength ' + str(Ycrystal) + ' mm')
f.write('\n' + '/gate/crystal/geometry/setZLength ' + str(Zcrystal) + ' mm')
f.write('\n' + '/gate/crystal/setMaterial G4_AIR')
f.write('\n' + '/gate/crystal/placement/setTranslation    0. 0. 0. mm')
f.write('\n' + '/gate/crystal/vis/setColor blue')
f.write('\n' + '/gate/crystal/vis/setVisible 0')
f.write('\n')
f.write('\n' + '# BGO layer {Setting the dimensions of each crystal layer}')
f.write('\n' + '/gate/crystal/daughters/name BGO')
f.write('\n' + '/gate/crystal/daughters/insert box')
f.write('\n' + '/gate/BGO/placement/setTranslation 0.0 0.0 0.0 mm')
f.write('\n' + '/gate/BGO/geometry/setXLength ' + str(Xcrystal) + ' mm')
f.write('\n' + '/gate/BGO/geometry/setYLength ' + str(Ycrystal) + ' mm')
f.write('\n' + '/gate/BGO/geometry/setZLength ' + str(MetaPixelZBGO) + ' mm')
f.write('\n' + '/gate/BGO/setMaterial BGO')
f.write('\n' + '/gate/BGO/vis/forceWireframe')
f.write('\n' + '/gate/BGO/vis/setColor green')
f.write('\n' + '/gate/BGO/repeaters/insert cubicArray')
f.write('\n' + '/gate/BGO/cubicArray/setRepeatNumberX 1')
f.write('\n' + '/gate/BGO/cubicArray/setRepeatNumberY 1')
f.write('\n' + '/gate/BGO/cubicArray/setRepeatNumberZ ' + str(NMetaPixelZBGO))
f.write('\n' + '/gate/BGO/cubicArray/setRepeatVector 0. 0. ' + str(MetaPixelZEJ232+MetaPixelZBGO) + ' mm')
f.write('\n' + '/gate/BGO/vis/setVisible 1')
f.write('\n')
f.write('\n' + '# EJ232 layer {Setting the dimensions of each crystal layer}')
f.write('\n' + '/gate/crystal/daughters/name EJ232')
f.write('\n' + '/gate/crystal/daughters/insert box')
f.write('\n' + '/gate/EJ232/placement/setTranslation 0.0 0.0 0.0 mm')
f.write('\n' + '/gate/EJ232/geometry/setXLength ' + str(Xcrystal) + ' mm')
f.write('\n' + '/gate/EJ232/geometry/setYLength ' + str(Ycrystal) + ' mm')
f.write('\n' + '/gate/EJ232/geometry/setZLength ' + str(MetaPixelZEJ232) + ' mm')
f.write('\n' + '/gate/EJ232/setMaterial EJ232')
f.write('\n' + '/gate/EJ232/vis/forceWireframe')
f.write('\n' + '/gate/EJ232/vis/setColor red')
f.write('\n' + '/gate/EJ232/repeaters/insert cubicArray')
f.write('\n' + '/gate/EJ232/cubicArray/setRepeatNumberX 1')
f.write('\n' + '/gate/EJ232/cubicArray/setRepeatNumberY 1')
f.write('\n' + '/gate/EJ232/cubicArray/setRepeatNumberZ ' + str(NMetaPixelZEJ232))
f.write('\n' + '/gate/EJ232/cubicArray/setRepeatVector 0. 0. ' + str(MetaPixelZEJ232+MetaPixelZBGO) + ' mm')
f.write('\n' + '/gate/EJ232/vis/setVisible 1')
f.write('\n')
f.write('\n' + '# === Optical coupling ===')
f.write('\n' + '/gate/module/daughters/name OpticalCoupling')
f.write('\n' + '/gate/module/daughters/insert box')
f.write('\n' + '/gate/OpticalCoupling/geometry/setXLength ' + str(Xoptcoup) + ' mm')
f.write('\n' + '/gate/OpticalCoupling/geometry/setYLength ' + str(Yoptcoup) + ' mm')
f.write('\n' + '/gate/OpticalCoupling/geometry/setZLength ' + str(Zoptcoup) + ' mm')
f.write('\n' + '/gate/OpticalCoupling/setMaterial BC-630')
f.write('\n' + '/gate/OpticalCoupling/vis/setColor yellow')
f.write('\n' + '/gate/OpticalCoupling/placement/setTranslation ' + str(TranslationXoptcoup) + ' 0. 0.  mm')
f.write('\n' + '/gate/OpticalCoupling/vis/setVisible 0')
f.write('\n')
f.write('\n' + '# === SiPM cover window ===')
f.write('\n' + '/gate/module/daughters/name SiPMCover')
f.write('\n' + '/gate/module/daughters/insert box')
f.write('\n' + '/gate/SiPMCover/geometry/setXLength ' + str(XSiPMCover) + ' mm')
f.write('\n' + '/gate/SiPMCover/geometry/setYLength ' + str(YSiPMCover) + ' mm')
f.write('\n' + '/gate/SiPMCover/geometry/setZLength ' + str(ZSiPMCover) + ' mm')
f.write('\n' + '/gate/SiPMCover/setMaterial Epoxy')
f.write('\n' + '/gate/SiPMCover/vis/setColor blue')
f.write('\n' + '/gate/SiPMCover/placement/setTranslation ' + str(TranslationXSiPMCover) + ' 0. 0.  mm')
f.write('\n' + '/gate/SiPMCover/repeaters/insert              cubicArray')
f.write('\n' + '/gate/SiPMCover/cubicArray/setRepeatNumberX   1')
f.write('\n' + '/gate/SiPMCover/cubicArray/setRepeatNumberY ' + str(nSiPMY))
f.write('\n' + '/gate/SiPMCover/cubicArray/setRepeatNumberZ ' + str(nSiPMZ))
f.write('\n' + '/gate/SiPMCover/cubicArray/setRepeatVector 0 ' + str(SiPMPitchSizeY) + ' ' + str(SiPMPitchSizeZ) + ' mm')
f.write('\n' + '/gate/SiPMCover/vis/setVisible 0')
f.write('\n')
f.write('\n' + '# === SiPM ===')
f.write('\n' + '/gate/module/daughters/name             SiPM')
f.write('\n' + '/gate/module/daughters/insert           box')
f.write('\n' + '/gate/SiPM/geometry/setXLength ' + str(XSiPM) + ' mm')
f.write('\n' + '/gate/SiPM/geometry/setYLength ' + str(YSiPM) + ' mm')
f.write('\n' + '/gate/SiPM/geometry/setZLength ' + str(ZSiPM) + ' mm')
f.write('\n' + '/gate/SiPM/setMaterial                   Silicon')
f.write('\n' + '/gate/SiPM/vis/setColor                  red')
f.write('\n' + '/gate/SiPM/placement/setTranslation ' + str(TranslationXSiPM) + ' 0. 0.  mm')
f.write('\n' + '/gate/SiPM/repeaters/insert              cubicArray')
f.write('\n' + '/gate/SiPM/cubicArray/setRepeatNumberX   1')
f.write('\n' + '/gate/SiPM/cubicArray/setRepeatNumberY ' + str(nSiPMY))
f.write('\n' + '/gate/SiPM/cubicArray/setRepeatNumberZ ' + str(nSiPMZ))
f.write('\n' + '/gate/SiPM/cubicArray/setRepeatVector 0. ' + str(SiPMPitchSizeY) + ' ' + str(SiPMPitchSizeZ) + ' mm')
f.write('\n' + '/gate/SiPM/vis/setVisible 0')
f.write('\n')
f.write('\n' + '/gate/crystal/repeaters/insert              cubicArray')
f.write('\n' + '/gate/crystal/cubicArray/setRepeatNumberX   1')
f.write('\n' + '/gate/crystal/cubicArray/setRepeatNumberY ' + str(nSiPMY))
f.write('\n' + '/gate/crystal/cubicArray/setRepeatNumberZ ' + str(nSiPMZ))
f.write('\n' + '/gate/crystal/cubicArray/setRepeatVector 0 ' + str(SiPMPitchSizeY) + ' ' + str(SiPMPitchSizeZ) + ' mm')
f.write('\n')
f.write('\n' + '/gate/module/repeaters/insert cubicArray')
f.write('\n' + '/gate/module/cubicArray/setRepeatNumberX 1')
f.write('\n' + '/gate/module/cubicArray/setRepeatNumberY ' + str(NModY))
f.write('\n' + '/gate/module/cubicArray/setRepeatNumberZ ' + str(NModZ))
f.write('\n' + '/gate/module/cubicArray/setRepeatVector 0. ' + str(ModPitchSizeY) + ' ' + str(ModPitchSizeZ) + ' mm')
f.write('\n')
f.write('\n' + '/gate/rsector/repeaters/insert ring')
f.write('\n' + '/gate/rsector/ring/setRepeatNumber ' + str(NrsectorPerRing))
f.write('\n')
f.write('\n' + '/gate/systems/cylindricalPET/rsector/attach rsector')
f.write('\n' + '/gate/systems/cylindricalPET/module/attach module')
f.write('\n' + '#/gate/systems/cylindricalPET/module/attach SiPMCover')
f.write('\n' + '/gate/systems/cylindricalPET/crystal/attach crystal')
f.write('\n' + '/gate/systems/cylindricalPET/layer0/attach BGO')
f.write('\n' + '/gate/systems/cylindricalPET/layer1/attach EJ232')
f.write('\n')
f.write('\n' + '# Attach CrystalSD')
f.write('\n' + '/gate/BGO/attachCrystalSD')
f.write('\n' + '/gate/EJ232/attachCrystalSD')
f.write('\n' + '#/gate/crystal/attachCrystalSD')
f.write('\n' + '#/gate/SiPMCover/attachCrystalSD')
f.write('\n')
f.write('\n' + '/gate/systems/cylindricalPET/describe')
f.write('\n')
f.write('\n')

f.close()

f = open('opticalsurfacespy.mac','w')
f.write('\n' + '# === surfaces ===')
f.write('\n' + '/gate/crystal/surfaces/name              crystal_module')
f.write('\n' + '/gate/crystal/surfaces/insert            module')
f.write('\n' + '#CrystalSurface (PTFE_wrapped or Air or smooth)')
f.write('\n' + '/gate/crystal/surfaces/crystal_module/setSurface PTFE_wrapped')
f.write('\n')
f.write('\n' + '/gate/crystal/surfaces/name              optcoupsurf')
f.write('\n' + '/gate/crystal/surfaces/insert            OpticalCoupling')
f.write('\n' + '/gate/crystal/surfaces/optcoupsurf/setSurface smooth')
f.write('\n')
f.write('\n' + '/gate/OpticalCoupling/surfaces/name              crystalbacksurf')
f.write('\n' + '/gate/OpticalCoupling/surfaces/insert            crystal')
f.write('\n' + '/gate/OpticalCoupling/surfaces/crystalbacksurf/setSurface smooth')
f.write('\n')
f.write('\n' + '/gate/OpticalCoupling/surfaces/name              SiPMCoversurf')
f.write('\n' + '/gate/OpticalCoupling/surfaces/insert            SiPMCover')
f.write('\n' + '/gate/OpticalCoupling/surfaces/SiPMCoversurf/setSurface smooth')
f.write('\n')
f.write('\n' + '/gate/OpticalCoupling/surfaces/name              lateralsurf')
f.write('\n' + '/gate/OpticalCoupling/surfaces/insert            module')
f.write('\n' + '/gate/OpticalCoupling/surfaces/lateralsurf/setSurface smooth')
f.write('\n')
f.write('\n' + '/gate/SiPMCover/surfaces/name              SiPMCoverbacksurf')
f.write('\n' + '/gate/SiPMCover/surfaces/insert            OpticalCoupling')
f.write('\n' + '/gate/SiPMCover/surfaces/SiPMCoverbacksurf/setSurface smooth')
f.write('\n')
f.write('\n' + '/gate/SiPMCover/surfaces/name              lateralsurf')
f.write('\n' + '/gate/SiPMCover/surfaces/insert            module')
f.write('\n' + '/gate/SiPMCover/surfaces/lateralsurf/setSurface smooth')
f.write('\n')
f.write('\n' + '/gate/SiPMCover/surfaces/name              SiPMsurf')
f.write('\n' + '/gate/SiPMCover/surfaces/insert            SiPM')
f.write('\n' + '/gate/SiPMCover/surfaces/SiPMsurf/setSurface SiPM')

f.close()

