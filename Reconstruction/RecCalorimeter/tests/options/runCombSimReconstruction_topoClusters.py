import os
#import numpy as np
from GaudiKernel.SystemOfUnits import MeV,GeV

#loads array of random seeds from file                                                                                                                                                                                                   
#seed_array = np.loadtxt('/afs/cern.ch/user/c/cneubuse/FCCSW/condor/seeds.txt',dtype='int',delimiter=',')                                                                                                                                   
#set these in the .sh script                                                                                                                                                                                                                

energy=10*GeV
num_events=1
bfield=0
i=1
particle=1
eta=0.

particleType = "pi-"
if particle==0:
    particleType = "e-"
if particle==2:
    particleType = "mu-"
if particle==3:
    particleType = "pi0"
print particleType

from Gaudi.Configuration import *

from Configurables import ApplicationMgr, FCCDataSvc, PodioOutput

podioevent = FCCDataSvc("EventDataSvc")
# input="/eos/experiment/fcc/users/c/cneubuse/FccHcal/fullModelSim/combCalo/deltaEta.001/output_combCalo_"+str(particleType)+str(energy)+"GeV_bfield"+str(bfield)+"_eta"+str(eta)+"_part"+str(i)+".root")

from Configurables import SimG4Svc, GeoSvc
geoservice = GeoSvc("GeoSvc", detectors=[  'file:Detector/DetFCChhBaseline1/compact/FCChh_DectEmptyMaster.xml',
                                           'file:Detector/DetFCChhECalInclined/compact/FCChh_ECalBarrel_withCryostat.xml',
                                           'file:Detector/DetFCChhHCalTile/compact/FCChh_HCalBarrel_TileCal.xml'],
                    OutputLevel = INFO)

geantservice = SimG4Svc("SimG4Svc", detector='SimG4DD4hepDetector', physicslist="SimG4FtfpBert", actions="SimG4FullSimActions")
geantservice.g4PostInitCommands += ["/run/setCut 0.1 mm"]

# Magnetic field                                                                                                                                                                                                                           
from Configurables import SimG4ConstantMagneticFieldTool
field = SimG4ConstantMagneticFieldTool("SimG4ConstantMagneticFieldTool",FieldOn=False)

from Configurables import SimG4Alg, SimG4SaveCalHits, InspectHitsCollectionsTool
saveecaltool = SimG4SaveCalHits("saveECalHits", readoutNames = ["ECalBarrelEta"],
                                positionedCaloHits = "ECalPositionedHits",
                                caloHits = "ECalHits")
savehcaltool = SimG4SaveCalHits("saveHCalHits", readoutNames = ["BarHCal_Readout"],
                                positionedCaloHits = "HCalPositionedHits",
                                caloHits = "HCalHits")

from Configurables import SimG4SingleParticleGeneratorTool
pgun = SimG4SingleParticleGeneratorTool("SimG4SingleParticleGeneratorTool",saveEdm=True,
                                        particleName=particleType,energyMin=energy,energyMax=energy,etaMin=-0.5,etaMax=0.5,
                                        OutputLevel =DEBUG)

geantsim = SimG4Alg("SimG4Alg",
                    outputs= ["SimG4SaveCalHits/saveECalHits", "SimG4SaveCalHits/saveHCalHits"],
                    eventProvider=pgun,
                    OutputLevel=DEBUG)

# common CAL specific information
# readout name
ecalReadoutName = "ECalBarrelPhiEta"
# active material identifier name
ecalIdentifierName = ["layer"]
# active material volume name
ecalVolumeName = ["layer"]
ecalNumberOfLayers = [130]
# ECAL bitfield names & values system:4,cryo:1,type:3,subtype:3,layer:8,eta:9,phi:10
ecalFieldNames = ["system"]
ecalFieldValues = [5]
# readout name
hcalReadoutName = "BarHCal_Readout_phieta"
# active material identifier name
hcalIdentifierName = ["row","layer"]
# active material volume name
hcalVolumeName = ["wedgeVolume","layerVolume"]
hcalNumberOfLayers = [510,10]
## HCAL bitfield names & values system:4,row:9,layer:5,eta:-9,phi:-10
hcalFieldNames = ["system"]
hcalFieldValues = [8]

#Configure tools for calo reconstruction
from Configurables import CalibrateInLayersTool
calibEcells = CalibrateInLayersTool("Calibrate",
                                    # sampling fraction obtained using SamplingFractionInLayers from DetStudies package                                                                                                                  
                                    samplingFraction = [0.12125] * 4 + [0.14283] * 18 + [0.16354] * 18 + [0.17662] * 18 + [0.18867] * 18 + [0.19890] * 18 + [0.20637] * 18 + [0.20802] * 18,
                                    readoutName = ecalReadoutName,
                                    layerFieldName = "layer")
#Configure tools for calo reconstruction
from Configurables import CalibrateCaloHitsTool
calibHcells = CalibrateCaloHitsTool("CalibrateHCal", invSamplingFraction="34.5 ")

from Configurables import CreateCaloCells
createEcells = CreateCaloCells("CreateECaloCells",
                               doCellCalibration=True,
                               calibTool=calibEcells,
                               addCellNoise=False, filterCellNoise=False,
                               OutputLevel=DEBUG,
                               hits="ECalHits",
                               cells="ECalCells")
createHcells = CreateCaloCells("CreateHCaloCells",
                               doCellCalibration=True,
                               calibTool=calibHcells,
                               addCellNoise=False, filterCellNoise=False,
                               OutputLevel=INFO,
                               hits="HCalHits",
                               cells="HCalCells")

from Configurables import CreateVolumeCaloPositions
# Ecal cell positions
positionsEcal = CreateVolumeCaloPositions("positionsEcal", OutputLevel = INFO)
positionsEcal.hits.Path = "ECalCells"
positionsEcal.positionedHits.Path = "ECalPositions"
positionsHcal = CreateVolumeCaloPositions("positionsHcal", OutputLevel = INFO)
positionsHcal.hits.Path = "HCalCells"
positionsHcal.positionedHits.Path = "HCalPositions"

from Configurables import RedoSegmentation
resegmentHcal = RedoSegmentation("ReSegmentationHcal",
                                 # old bitfield (readout)
                                 oldReadoutName = "BarHCal_Readout",
                                 # specify which fields are going to be altered (deleted/rewritten)
                                 oldSegmentationIds = ["module","tile"],
                                 # new bitfield (readout), with new segmentation
                                 newReadoutName = hcalReadoutName,
                                 debugPrint = 10,
                                 OutputLevel = DEBUG,
                                 inhits = "HCalPositions",
                                 outhits = "newHCalCells")
resegmentEcal = RedoSegmentation("ReSegmentationEcal",
                                 # old bitfield (readout)
                                 oldReadoutName = "ECalBarrelEta",
                                 # specify which fields are going to be altered (deleted/rewritten)
                                 oldSegmentationIds = ["module"],
                                 # new bitfield (readout), with new segmentation
                                 newReadoutName = ecalReadoutName,
                                 debugPrint = 10,
                                 OutputLevel = DEBUG,
                                 inhits = "ECalPositions",
                                 outhits = "newECalCells")

from Configurables import CreateVolumeCaloPositions
# Ecal cell positions
positionsEcal2 = CreateVolumeCaloPositions("newPositionsEcal", OutputLevel = INFO)
positionsEcal2.hits.Path = "newECalCells"
positionsEcal2.positionedHits.Path = "newECalPositions"
positionsHcal2 = CreateVolumeCaloPositions("newPositionsHcal", OutputLevel = INFO)
positionsHcal2.hits.Path = "newHCalCells"
positionsHcal2.positionedHits.Path = "newHCalPositions"

#Create topo clusters
from Configurables import TubeLayerPhiEtaCaloTool, CombinedCaloTopoCluster
ecalgeo = TubeLayerPhiEtaCaloTool("EcalGeo",
                                  readoutName = ecalReadoutName,
                                  activeVolumeNames = ecalVolumeName,
                                  activeFieldNames = ecalIdentifierName,
                                  fieldNames = ecalFieldNames,
                                  fieldValues = ecalFieldValues,
                                  # to make it working with MergeLayers algorithm
                                  activeVolumesNumber = ecalNumberOfLayers,
                                  OutputLevel = INFO)
hcalgeo = TubeLayerPhiEtaCaloTool("HcalGeo",
                                  readoutName = hcalReadoutName,
                                  activeVolumeNames = hcalVolumeName,
                                  activeFieldNames = hcalIdentifierName,
                                  fieldNames = hcalFieldNames,
                                  fieldValues = hcalFieldValues,
                                  # to make it working with MergeLayers algorithm
                                  activeVolumesNumber = hcalNumberOfLayers,
                                  OutputLevel = INFO)
createTopoClusters = CombinedCaloTopoCluster("CreateTopoClusters",
                                             ecalCells = "newECalCells",
                                             hcalCells = "newHCalCells",
                                             geometryTool = ecalgeo,
                                             ecalReadoutName = ecalReadoutName,
                                             ecalFieldNames = ecalFieldNames,
                                             ecalFieldValues = ecalFieldValues,
                                             hcalReadoutName = hcalReadoutName,
                                             hcalFieldNames = hcalFieldNames,
                                             hcalFieldValues = hcalFieldValues,
                                             seedThresholdEcal =  7.5, #in MeV
                                             neighbourThresholdEcal = 3, # in MeV
                                             OutputLevel = DEBUG)
createTopoClusters.clusters.Path = "topoClusters"

out = PodioOutput("out", filename = "~/FCCSW/condor/output_ecal_reconstructionTopoClusters_"+str(particleType)+str(energy)+"GeV_bfield"+str(bfield)+"_part"+str(i)+".root",
                  OutputLevel=DEBUG)
out.outputCommands = ["keep *"]#,"drop ECalHits"]

#CPU information
from Configurables import AuditorSvc, ChronoAuditor
chra = ChronoAuditor()
audsvc = AuditorSvc()
audsvc.Auditors = [chra]
geantsim.AuditExecute = True
createEcells.AuditExecute = True
createHcells.AuditExecute = True
positionsEcal.AuditExecute = True
positionsHcal.AuditExecute = True
resegmentEcal.AuditExecute = True
resegmentHcal.AuditExecute = True
positionsEcal2.AuditExecute = True
positionsHcal2.AuditExecute = True
createTopoClusters.AuditExecute = True
out.AuditExecute = True

ApplicationMgr(
    TopAlg = [geantsim,
              createEcells,createHcells,
              positionsEcal,positionsHcal,
              resegmentEcal,resegmentHcal,
              positionsEcal2,positionsHcal2,
              createTopoClusters,
              out
              ],
    EvtSel = 'NONE',
    EvtMax   = int(num_events),
    ExtSvc = [geoservice, podioevent, audsvc],
)
