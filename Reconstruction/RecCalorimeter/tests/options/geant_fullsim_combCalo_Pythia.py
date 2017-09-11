import os
import numpy as np

#loads array of random seeds from file                                                                                                                                                                          
#seed_array = np.loadtxt('/afs/cern.ch/user/c/cneubuse/FCCSW/condor/seeds.txt',dtype='int',delimiter=',')

#set these in the .sh script                                                                                                                                                                                                                
energy=100000
num_events=500
magnetic_field=0
i=1
particle=1
eta=1.0

particleType = "pi-"
if particle==0:
    particleType = "e-"
if particle==2:
    particleType = "mu-"
if particle==3:
    particleType = "pi0"
if particle==4:
    particleType = "gamma"
print particleType

from Gaudi.Configuration import *
messageLevelPythia =INFO 

from Configurables import ApplicationMgr, FCCDataSvc, PodioOutput

podioevent   = FCCDataSvc("EventDataSvc")

from Configurables import GeoSvc
geoservice = GeoSvc("GeoSvc", detectors=[  'file:/eos/experiment/fcc/hh/simulation/FCCSW/Detector/DetFCChhBaseline1/compact/FCChh_DectEmptyMaster.xml',
                                           'file:/eos/experiment/fcc/hh/simulation/FCCSW/Detector/DetFCChhECalInclined/compact/FCChh_ECalBarrel_withCryostat.xml',
                                           'file:/eos/experiment/fcc/hh/simulation/FCCSW/Detector/DetFCChhHCalTile/compact/FCChh_HCalBarrel_TileCal.xml',
                                           'file:/eos/experiment/fcc/hh/simulation/FCCSW/Detector/DetFCChhHCalTile/compact/FCChh_HCalExtendedBarrel_TileCal.xml',
                                           'file:/eos/experiment/fcc/hh/simulation/FCCSW/Detector/DetFCChhCalDiscs/compact/Endcaps_coneCryo.xml',
                                           'file:/eos/experiment/fcc/hh/simulation/FCCSW/Detector/DetFCChhCalDiscs/compact/Forward_coneCryo.xml'],
                    OutputLevel = INFO)
# Geant4 service                                                                                                                                                                         
# Configures the Geant simulation: geometry, physics list and user actions
from Configurables import SimG4Svc
geantservice = SimG4Svc("SimG4Svc", detector='SimG4DD4hepDetector', physicslist="SimG4FtfpBert", actions="SimG4FullSimActions")

# range cut                                                                                                                                                                                                                        
geantservice.g4PostInitCommands += ["/run/setCut 0.1 mm"]

# Magnetic field                                                                                                                                                                                                                           
from Configurables import SimG4ConstantMagneticFieldTool
field = SimG4ConstantMagneticFieldTool("SimG4ConstantMagneticFieldTool",FieldOn=False)

# Geant4 algorithm                                                                                                                                                                                                       
# Translates EDM to G4Event, passes the event to G4, writes out outputs via tools                                                                                                                                          
# and a tool that saves the calorimeter hits                                                                                                                                                                                         
from Configurables import SimG4Alg, SimG4SaveCalHits, InspectHitsCollectionsTool
saveecaltool = SimG4SaveCalHits("saveECalBarrelHits", readoutNames = ["ECalBarrelPhiEta"],
                                positionedCaloHits = "ECalBarrelPositionedHits",
                                caloHits = "ECalBarrelHits")
savecalendcaptool = SimG4SaveCalHits("saveECalEndcapHits", readoutNames = ["EMECPhiEta"])
savecalendcaptool.positionedCaloHits.Path = "ECalEndcapPositionedHits"
savecalendcaptool.caloHits.Path = "ECalEndcapHits"
savecalfwdtool = SimG4SaveCalHits("saveECalFwdHits", readoutNames = ["EMFwdPhiEta"])
savecalfwdtool.positionedCaloHits.Path = "ECalFwdPositionedHits"
savecalfwdtool.caloHits.Path = "ECalFwdHits"
savehcaltool = SimG4SaveCalHits("saveHCalHits",readoutNames = ["BarHCal_Readout_phieta"],
                                positionedCaloHits="HCalPositionedHits",
                                caloHits="HCalHits")                                
saveexthcaltool = SimG4SaveCalHits("saveExtHCalHits",readoutNames = ["ExtBarHCal_Readout_phieta"],
                                   positionedCaloHits="extHCalPositionedHits",
                                   caloHits="extHCalHits")


pythiaConfFile="Pythia_pp_Higgs_100TeV.cmd"
#pythiaConfFile="Generation/data/Pythia_minbias_pp_100TeV.cmd"
#pythiaConfFile="/eos/experiment/fcc/hh/utils/pythiacards/pythia_pp_jj_M_500_1000.cmd"
#pythiaConfFile="/eos/experiment/fcc/hh/utils/pythiacards/pythia_pp_Zprime_40TeV_ttbar.cmd"
from Configurables import PythiaInterface, GenAlg
### PYTHIA algorithm
pythia8gentool = PythiaInterface("Pythia8Interface", Filename=pythiaConfFile, OutputLevel=messageLevelPythia)
pythia8gen = GenAlg("Pythia8", SignalProvider=pythia8gentool)
## Write the HepMC::GenEvent to the data service
pythia8gen.hepmc.Path = "hepmc"

# reads an HepMC::GenEvent from the data service and writes a collection of EDM Particles
from Configurables import HepMCToEDMConverter
hepmc_converter = HepMCToEDMConverter("Converter")
hepmc_converter.hepmc.Path="hepmc"
hepmc_converter.genparticles.Path="allGenParticles"
hepmc_converter.genvertices.Path="allGenVertices"

from Configurables import SimG4PrimariesFromEdmTool
particle_converter = SimG4PrimariesFromEdmTool("EdmConverter")
particle_converter.genParticles.Path = "allGenParticles"

from Configurables import SimG4SingleParticleGeneratorTool
pgun = SimG4SingleParticleGeneratorTool("SimG4SingleParticleGeneratorTool",saveEdm=True,
                particleName=particleType,energyMin=energy,energyMax=energy, #etaMin=eta,etaMax=eta,
                OutputLevel =DEBUG)

geantsim = SimG4Alg("SimG4Alg",
                    outputs= ["SimG4SaveCalHits/saveECalBarrelHits", "SimG4SaveCalHits/saveECalEndcapHits", "SimG4SaveCalHits/saveECalFwdHits", "SimG4SaveCalHits/saveHCalHits", "SimG4SaveCalHits/saveExtHCalHits"],
                    eventProvider=particle_converter,
#                    eventProvider=pgun,
                    OutputLevel=DEBUG)

# Configure tools for calo reconstruction                                                                                                                                                                    
from Configurables import CalibrateInLayersTool
calibcellsBarrel = CalibrateInLayersTool("CalibrateBarrel",
                                   # sampling fraction obtained using SamplingFractionInLayers from DetStudies package
                                   samplingFraction = [0.12125] * 4 + [0.14283] * 18 + [0.16354] * 18 + [0.17662] * 18 + [0.18867] * 18 + [0.19890] * 18 + [0.20637] * 18 + [0.20802] * 18,
                                   readoutName = "ECalBarrelEta",
                                   layerFieldName = "layer")
calibcellsEndcap = CalibrateInLayersTool("CalibrateEndcap",
                                         # sampling fraction obtained using SamplingFractionInLayers from DetStudies package
                                    samplingFraction = [0.15] * 118,
                                    readoutName = "EMECPhiEta",
                                    layerFieldName = "layer")
calibcellsFwd = CalibrateInLayersTool("CalibrateFwd",
                                         # sampling fraction obtained using SamplingFractionInLayers from DetStudies package
                                    samplingFraction = [0.00056] * 22,
                                    readoutName = "EMFwdPhiEta",
                                    layerFieldName = "layer")

#Configure tools for calo reconstruction
from Configurables import CalibrateCaloHitsTool
calibHcells = CalibrateCaloHitsTool("CalibrateHCal", invSamplingFraction="34.5 ")

from Configurables import CreateCaloCells
createEcellsBarrel = CreateCaloCells("CreateECaloCellsBarrel",
                                    doCellCalibration=True,
                                    calibTool=calibcellsBarrel,
                                    addCellNoise=False, filterCellNoise=False,
                                    OutputLevel=DEBUG)
createEcellsBarrel.hits.Path="ECalBarrelHits"
createEcellsBarrel.cells.Path="ECalBarrelCells"
createEcellsEndcap = CreateCaloCells("CreateECaloCellsEndcap",
                                    doCellCalibration=True,
                                    calibTool=calibcellsEndcap,
                                    addCellNoise=False, filterCellNoise=False,
                                    OutputLevel=DEBUG)
createEcellsEndcap.hits.Path="ECalEndcapHits"
createEcellsEndcap.cells.Path="ECalEndcapCells"
createEcellsFwd = CreateCaloCells("CreateECaloCellsFwd",
                                    doCellCalibration=True,
                                    calibTool=calibcellsFwd,
                                    addCellNoise=False, filterCellNoise=False,
                                    OutputLevel=DEBUG)
createEcellsFwd.hits.Path="ECalFwdHits"
createEcellsFwd.cells.Path="ECalFwdCells"

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
                               hits="extHCalHits",
                               cells="extHCalCells")

# additionally for HCal
from Configurables import CreateVolumeCaloPositions
positions = CreateVolumeCaloPositions("positions", OutputLevel = INFO)
positions.hits.Path = "HCalCells"
positions.positionedHits.Path = "HCalPositions"

# Ecal cell positions
positionsEcal = CreateVolumeCaloPositions("positionsEcal", OutputLevel = INFO)
positionsEcal.hits.Path = "ECalBarrelCells"
positionsEcal.positionedHits.Path = "ECalBarrelPositions"

# extHcal cell positions                                                                                                                                                                                       
positionsExtHcal = CreateVolumeCaloPositions("positionsExtHcal", OutputLevel = INFO)
positionsExtHcal.hits.Path = "extHCalCells"
positionsExtHcal.positionedHits.Path = "extHCalPositions"

out = PodioOutput("out", 
                  OutputLevel=DEBUG)
out.outputCommands = ["keep *"]
out.filename = "output_combCalo_"+str(particleType)+str(int(energy/1e3))+"GeV_eta"+str(eta)+"_part"+str(i)+".root"
out.filename = "output_higgs.root"

#CPU information
from Configurables import AuditorSvc, ChronoAuditor
chra = ChronoAuditor()
audsvc = AuditorSvc()
audsvc.Auditors = [chra]
#pythia8gen.AuditExecute = True
#hepmc_converter.AuditExecute = True
geantsim.AuditExecute = True
createEcellsBarrel.AuditExecute = True
createEcellsEndcap.AuditExecute = True
createEcellsFwd.AuditExecute = True
createHcells.AuditExecute = True
createExtHcells.AuditExecute = True
#positions.AuditExecute = True
#positionsExtHcal.AuditExecute = True
out.AuditExecute = True

ApplicationMgr(
    TopAlg = [pythia8gen, 
              hepmc_converter,
              geantsim,
              createEcellsBarrel,
              createEcellsEndcap,
              createEcellsFwd,
              createHcells,
              createExtHcells,
#              positions,
#              positionsEcal,
#              positionsExtHcal,
              out
              ],
    EvtSel = 'NONE',
    EvtMax   = int(num_events),
    ExtSvc = [podioevent, geoservice, geantservice, audsvc],
 )

