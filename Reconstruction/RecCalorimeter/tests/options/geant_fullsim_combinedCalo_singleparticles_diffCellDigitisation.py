
import os
import numpy as np

#loads array of random seeds from file                                                                                                                                                                          
#seed_array = np.loadtxt('/afs/cern.ch/user/c/cneubuse/FCCSW/condor/seeds.txt',dtype='int',delimiter=',')

#set these in the .sh script                                                                                                                                                                                                                
energy=10000
num_events=1
magnetic_field=0
i=1
particle=1
etaMin=.5
etaMax=.5

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

podioevent   = FCCDataSvc("EventDataSvc")

from Configurables import GeoSvc
geoservice = GeoSvc("GeoSvc", detectors=[  'file:Detector/DetFCChhBaseline1/compact/FCChh_DectEmptyMaster.xml',
                                           'file:Detector/DetFCChhECalInclined/compact/FCChh_ECalBarrel_withCryostat.xml',
                                           'file:Detector/DetFCChhHCalTile/compact/FCChh_HCalBarrel_TileCal.xml',
                                           'file:Detector/DetFCChhHCalTile/compact/FCChh_HCalExtendedBarrel_TileCal.xml',
                                           'file:Detector/DetFCChhCalDiscs/compact/Endcaps_coneCryo.xml' 
                                           ],
                    OutputLevel = INFO)
# Geant4 service                                                                                                                                                                         
# Configures the Geant simulation: geometry, physics list and user actions
from Configurables import SimG4Svc
geantservice = SimG4Svc("SimG4Svc", detector='SimG4DD4hepDetector', physicslist="SimG4FtfpBert", actions="SimG4FullSimActions")

#Setting random seeds for Geant4 simulations
#geantservice.g4PreInitCommands  += ["/random/setSeeds "+str(x)+" 0"] #where x is the number you want

# range cut
geantservice.g4PostInitCommands += ["/run/setCut 0.1 mm"]

# Magnetic field
from Configurables import SimG4ConstantMagneticFieldTool
field = SimG4ConstantMagneticFieldTool("SimG4ConstantMagneticFieldTool",FieldOn=False)

# Geant4 algorithm                                                                                                                                                                                                       
# Translates EDM to G4Event, passes the event to G4, writes out outputs via tools                                                                                                                                          
# and a tool that saves the calorimeter hits                                                                                                                                                                                         
from Configurables import SimG4Alg, SimG4SaveCalHits, InspectHitsCollectionsTool
saveecaltool = SimG4SaveCalHits("saveECalHits", readoutNames = ["ECalBarrelEta"],
                                positionedCaloHits = "ECalPositionedHits",
                                caloHits = "ECalHits")
savehcaltool = SimG4SaveCalHits("saveHCalHits",readoutNames = ["BarHCal_Readout"],
                                positionedCaloHits="HCalPositionedHits",
                                caloHits="HCalHits")                                
saveexthcaltool = SimG4SaveCalHits("saveExtHCalHits",readoutNames = ["ExtBarHCal_Readout"],
                                   positionedCaloHits="extHCalPositionedHits",
                                   caloHits="ExtHCalHits")
savecalendcaptool = SimG4SaveCalHits("saveECalEndcapHits", readoutNames = ["EMECPhiEta"])
savecalendcaptool.positionedCaloHits.Path = "ECalEndcapPositionedHits"
savecalendcaptool.caloHits.Path = "ECalEndcapHits"

# next, create the G4 algorithm, giving the list of names of tools ("XX/YY")                                                                                                                                                     
from Configurables import SimG4SingleParticleGeneratorTool
pgun = SimG4SingleParticleGeneratorTool("SimG4SingleParticleGeneratorTool",saveEdm=True,
                particleName=particleType,energyMin=energy,energyMax=energy,etaMin=etaMin,etaMax=etaMax,
                OutputLevel =DEBUG)

geantsim = SimG4Alg("SimG4Alg",
                       outputs= ["SimG4SaveCalHits/saveECalHits", "SimG4SaveCalHits/saveHCalHits", "SimG4SaveCalHits/saveExtHCalHits", "SimG4SaveCalHits/saveECalEndcapHits"],
                       eventProvider=pgun,
                       OutputLevel=DEBUG)

# common ECAL specific information
# readout name
ecalReadoutName = "ECalBarrelPhiEta"
# common HCAL specific information
# hcal readout name
hcalReadoutName = "BarHCal_Readout_phieta"
# extHcal readout name
extHcalReadoutName = "ExtBarHCal_Readout_phieta"

# Configure tools for calo reconstruction                                                                                                                                                                    
from Configurables import CalibrateInLayersTool
calibEcells = CalibrateInLayersTool("Calibrate",
                                    # sampling fraction obtained using SamplingFractionInLayers from DetStudies package
                                    samplingFraction = [0.12125] + [0.14283] + [0.16354] + [0.17662] + [0.18867] + [0.19890] + [0.20637] + [0.20802],
                                    readoutName = "ECalBarrelEta",
                                    layerFieldName = "layer")

calibcellsEndcap = CalibrateInLayersTool("CalibrateEndcap",
                                         # sampling fraction obtained using SamplingFractionInLayers from DetStudies package
                                         samplingFraction = [0.15] * 118,
                                         readoutName = "EMECPhiEta",
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
                               addCellNoise = False, filterCellNoise = False,
                               OutputLevel = DEBUG,
                               hits="HCalHits",
                               cells="HCalCells")

createExtHcells = CreateCaloCells("CreateExtHCaloCells",
                                  doCellCalibration=True,
                                  calibTool=calibHcells,
                                  addCellNoise = False, filterCellNoise = False,
                                  OutputLevel = DEBUG,
                                  hits="ExtHCalHits",
                                  cells="ExtHCalCells")

createcellsEndcap = CreateCaloCells("CreateCaloCellsEndcap",
                                    doCellCalibration=True,
                                    calibTool=calibcellsEndcap,
                                    addCellNoise=False, filterCellNoise=False,
                                    OutputLevel=DEBUG,
                                    hits="ECalEndcapHits",
                                    cells="ECalEndcapCells")

# additionally for HCal
from Configurables import CreateVolumeCaloPositions, CreateCellCaloPositions
volPositionsEcal = CreateVolumeCaloPositions("volPositionsEcal", OutputLevel = INFO)
volPositionsEcal.hits.Path = "ECalCells"
volPositionsEcal.positionedHits.Path = "ECalPositions"
volPositionsEmec = CreateVolumeCaloPositions("volPositionsEmec", OutputLevel = DEBUG)
volPositionsEmec.hits.Path = "ECalEndcapCells"
volPositionsEmec.positionedHits.Path = "EmecPositions"
volPositionsHcal = CreateVolumeCaloPositions("volPositionsHcal", OutputLevel = INFO)
volPositionsHcal.hits.Path = "HCalCells"
volPositionsHcal.positionedHits.Path = "HCalPositions"
# HCAL cell positions retrievable for original segmentation (w/o additional eta segmentation)
positionsHcal = CreateCellCaloPositions("cellPositionsHcal", readoutName="BarHCal_Readout", OutputLevel = INFO)
positionsHcal.hits.Path = "HCalCells"
positionsHcal.positionedHits.Path = "cellHCalPositions"
positionsExtHcal = CreateCellCaloPositions("cellPositionsExtHcal", readoutName="ExtBarHCal_Readout", OutputLevel = INFO)
positionsExtHcal.hits.Path = "ExtHCalCells"
positionsExtHcal.positionedHits.Path = "cellExtHCalPositions"
positionsEmec = CreateCellCaloPositions("cellPositionsEmec", readoutName="EMECPhiEta", OutputLevel = DEBUG)
positionsEmec.hits.Path = "ECalEndcapCells"
positionsEmec.positionedHits.Path = "cellEmecPositions"

from Configurables import RedoSegmentation
resegmentEcal = RedoSegmentation("ReSegmentationEcal",
                                 # old bitfield (readout)
                                 oldReadoutName = "ECalBarrelEta",
                                 # specify which fields are going to be altered (deleted/rewritten)
                                 # oldSegmentationIds = [],
                                 # new bitfield (readout), with new segmentation
                                 newReadoutName = ecalReadoutName,
                                 debugPrint = 10,
                                 OutputLevel = DEBUG,
                                 inhits = "ECalPositions",
                                 outhits = "newECalCells")
resegmentHcal = RedoSegmentation("ReSegmentationHcal",
                                 # old bitfield (readout)
                                 oldReadoutName = "BarHCal_Readout",
                                 # specify which fields are going to be altered (deleted/rewritten)
                                 oldSegmentationIds = ["eta"],
                                 # new bitfield (readout), with new segmentation
                                 newReadoutName = hcalReadoutName,
                                 debugPrint = 10,
                                 OutputLevel = DEBUG,
                                 inhits = "HCalPositions",
                                 outhits = "newHCalCells")

from Configurables import CreateCellCaloPositions
# Ecal cell positions
positionsEcal2 = CreateCellCaloPositions("cellPositionsEcal", readoutName=ecalReadoutName, OutputLevel = INFO)
positionsEcal2.hits.Path = "newECalCells"
positionsEcal2.positionedHits.Path = "cellECalPositions"
positionsHcal2 = CreateCellCaloPositions("cellPositionsHcalEtaPhi", readoutName=hcalReadoutName, OutputLevel = INFO)
positionsHcal2.hits.Path = "newHCalCells"
positionsHcal2.positionedHits.Path = "cellHCalPositionsEtaPhi"

out = PodioOutput("out", 
                  OutputLevel=DEBUG)
out.outputCommands = ["keep *", "drop ECalHits", "drop HCalHits", "drop ExtHCalHits", "drop ECalEndcapHits"]
#out.filename = "testHCALSegmentations_eta05.root"
out.filename = "output_combCalo_"+str(particleType)+str(int(energy/1e3))+"GeV_etaMin"+str(etaMin)+"_etaMax"+str(etaMax)+"_part"+str(i)+".root"
print "output_combCalo_"+str(particleType)+str(int(energy/1e3))+"GeV_etaMin"+str(etaMin)+"_etaMax"+str(etaMax)+"_part"+str(i)+".root"

#CPU information
from Configurables import AuditorSvc, ChronoAuditor
chra = ChronoAuditor()
audsvc = AuditorSvc()
audsvc.Auditors = [chra]
geantsim.AuditExecute = True
createEcells.AuditExecute = True
createHcells.AuditExecute = True
createExtHcells.AuditExecute = True
#createcellsEndcap.AuditExecute = True
volPositionsEcal.AuditExecute = True
volPositionsHcal.AuditExecute = True
#volPositionsEmec.AuditExecute = True
positionsHcal.AuditExecute = True
positionsExtHcal.AuditExecute = True
#positionsEmec.AuditExecute = True
resegmentEcal.AuditExecute = True
resegmentHcal.AuditExecute = True
positionsEcal2.AuditExecute = True
positionsHcal2.AuditExecute = True
out.AuditExecute = True

ApplicationMgr(
    TopAlg = [geantsim,
              createEcells,createHcells,createExtHcells,
              volPositionsEcal,volPositionsHcal,
              positionsHcal,positionsExtHcal,
              resegmentEcal,resegmentHcal,
              positionsEcal2,positionsHcal2,
              out
              ],
    EvtSel = 'NONE',
    EvtMax   = int(num_events),
    ExtSvc = [podioevent, geoservice, geantservice, audsvc],
 )

