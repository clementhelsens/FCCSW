
import os
import numpy as np
from FWCore.joboptions import parse_standard_job_options
args = parse_standard_job_options()

#loads array of random seeds from file                                                                                                                                                                          
#seed_array = np.loadtxt('/afs/cern.ch/user/c/cneubuse/FCCSW/condor/seeds.txt',dtype='int',delimiter=',')

#set these in the .sh script                                                                                                                                                                                                                
num_events=100
if args.nevents is not None:
    num_events = args.nevents

from Gaudi.Configuration import *
messageLevelPythia =INFO 

from Configurables import ApplicationMgr, FCCDataSvc, PodioOutput

podioevent   = FCCDataSvc("EventDataSvc")

from Configurables import GeoSvc
geoservice = GeoSvc("GeoSvc", detectors=[  'file:Detector/DetFCChhBaseline1/compact/FCChh_DectEmptyMaster.xml',
                                           'file:Detector/DetFCChhECalInclined/compact/FCChh_ECalBarrel_withCryostat.xml',
                                           'file:Detector/DetFCChhHCalTile/compact/FCChh_HCalBarrel_TileCal.xml',
                                           'file:Detector/DetFCChhHCalTile/compact/FCChh_HCalExtendedBarrel_TileCal.xml',
                                           'file:Detector/DetFCChhCalDiscs/compact/Endcaps_coneCryo.xml',
                                           'file:Detector/DetFCChhTrackerTkLayout/compact/Tracker.xml'
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
field = SimG4ConstantMagneticFieldTool("SimG4ConstantMagneticFieldTool",FieldOn=True)

# Geant4 algorithm                                                                                                                                                                                                       
# Translates EDM to G4Event, passes the event to G4, writes out outputs via tools                                                                                                                                          
# and a tool that saves the calorimeter hits                                                                                                                                                                                         
from Configurables import SimG4Alg, SimG4SaveCalHits, InspectHitsCollectionsTool,SimG4SaveTrackerHits
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

savetrackertool = SimG4SaveTrackerHits("saveTrackerHits", readoutNames = ["TrackerBarrelReadout","TrackerEndcapReadout"])
savetrackertool.positionedTrackHits.Path = "positionedHits"
savetrackertool.trackHits.Path = "hits"



pythiaConfFile="/afs/cern.ch/user/h/helsens/FCCsoft/Calo/FCCSW/pythia_pp_DrellYann.cmd"
if args.inputfile != '':
    pythiaConfFile = args.inputfile

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


from Configurables import GenParticleFilter
### Filters generated particles
# accept is a list of particle statuses that should be accepted
genfilter = GenParticleFilter("StableParticles", accept=[1], OutputLevel=INFO)
genfilter.allGenParticles.Path = "allGenParticles"
genfilter.filteredGenParticles.Path = "GenParticles"


from Configurables import SimG4PrimariesFromEdmTool
particle_converter = SimG4PrimariesFromEdmTool("EdmConverter")
particle_converter.genParticles.Path = "GenParticles"


# next, create the G4 algorithm, giving the list of names of tools ("XX/YY")                                                                                                                                                     
#from Configurables import SimG4SingleParticleGeneratorTool
#gun = SimG4SingleParticleGeneratorTool("SimG4SingleParticleGeneratorTool",saveEdm=True,
#                particleName=12,energyMin=100,energyMax=100,etaMin=1,etaMax=1,
#                OutputLevel=DEBUG)


geantsim = SimG4Alg("SimG4Alg",
                    outputs= ["SimG4SaveCalHits/saveECalHits", 
                              "SimG4SaveCalHits/saveHCalHits", 
                              "SimG4SaveCalHits/saveExtHCalHits",
                              "SimG4SaveCalHits/saveECalEndcapHits", 
                              "SimG4SaveTrackerHits/saveTrackerHits"],
                    eventProvider=particle_converter,
                    OutputLevel=INFO)


################################################################
################################################################
#pileupFilenames = ["/eos/experiment/fcc/hh/simulation/samples/tracker_calobarrel/minbias/output_helsens_2017-09-14_15:01:14.071.root"]
import glob
All_pileupfilenames = glob.glob("/eos/experiment/fcc/hh/simulation/samples/tracker_calobarrel/minbias/output_helsens*.root")
pileupFilenames = [f for f in All_pileupfilenames if 'converted' not in f]
import random
random.shuffle(pileupFilenames)
# the collections to be read from the signal event
signalCollections = ["allGenVertices", "allGenParticles", "hits", "positionedHits", "ECalPositionedHits", "HCalPositionedHits", "extHCalPositionedHits"]


from Configurables import PileupParticlesMergeTool
particlemergetool = PileupParticlesMergeTool("MyPileupParticlesMergeTool")
# branchnames for the pileup
particlemergetool.genParticlesBranch = "allGenParticles"
particlemergetool.genVerticesBranch = "allGenVertices"
# branchnames for the signal
particlemergetool.signalGenVertices.Path = "allGenVertices"
particlemergetool.signalGenParticles.Path = "allGenParticles"
# branchnames for the output
particlemergetool.mergedGenParticles.Path = "overlaidGenParticles"
particlemergetool.mergedGenVertices.Path = "overlaidGenVertices"

# edm data from simulation: hits and positioned hits for the tracker
from Configurables import PileupTrackHitMergeTool
trackhitsmergetool = PileupTrackHitMergeTool("MyTrackHitMergeTool")
# branchnames for the pileup
trackhitsmergetool.pileupHitsBranch = "hits"
trackhitsmergetool.pileupPositionedHitsBranch = "positionedHits"
# branchnames for the signal
trackhitsmergetool.signalHits = "hits"
trackhitsmergetool.signalPositionedHits = "positionedHits"
# branchnames for the output
trackhitsmergetool.mergedHits = "overlaidTrackHits"
trackhitsmergetool.mergedPositionedHits = "overlaidPositionedTrackHits"


# edm data from simulation: hits and positioned hits for the EM Calo
from Configurables import PileupCaloHitMergeTool
ECalhitsmergetool = PileupCaloHitMergeTool("MyECalHitMergeTool")
# branchnames for the pileup
ECalhitsmergetool.pileupHitsBranch = "ECalHits"
ECalhitsmergetool.pileupPositionedHitsBranch = "ECalPositionedHits"
# branchnames for the signal
ECalhitsmergetool.signalHits = "ECalHits"
ECalhitsmergetool.signalPositionedHits = "ECalPositionedHits"
# branchnames for the output
ECalhitsmergetool.mergedHits = "overlaidECalHits"
ECalhitsmergetool.mergedPositionedHits = "overlaidECalPositionedHits"

# edm data from simulation: hits and positioned hits for the Had Calo
from Configurables import PileupCaloHitMergeTool
HCalhitsmergetool = PileupCaloHitMergeTool("MyHCalHitMergeTool")
# branchnames for the pileup
HCalhitsmergetool.pileupHitsBranch = "HCalHits"
HCalhitsmergetool.pileupPositionedHitsBranch = "HCalPositionedHits"
# branchnames for the signal
HCalhitsmergetool.signalHits = "HCalHits"
HCalhitsmergetool.signalPositionedHits = "HCalPositionedHits"
# branchnames for the output
HCalhitsmergetool.mergedHits = "overlaidHCalHits"
HCalhitsmergetool.mergedPositionedHits = "overlaidHCalPositionedHits"


# use the pileuptool to specify the number of pileup
from Configurables import ConstPileUp
pileuptool = ConstPileUp("MyPileupTool", numPileUpEvents=200)

# algorithm for the overlay
from Configurables import PileupOverlayAlg
overlay = PileupOverlayAlg()
overlay.pileupFilenames = pileupFilenames
overlay.randomizePileup = False
overlay.mergeTools = ["PileupParticlesMergeTool/MyPileupParticlesMergeTool",
                      "PileupTrackHitMergeTool/MyTrackHitMergeTool",
                      "PileupCaloHitMergeTool/MyECalHitMergeTool",
                      "PileupCaloHitMergeTool/MyHCalHitMergeTool"]
overlay.PileUpTool = pileuptool


################################################################
################################################################





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
                               OutputLevel=INFO,
                               #hits="ECalHits",
                               hits="overlaidECalHits",
                               cells="ECalCells")

createHcells = CreateCaloCells("CreateHCaloCells",
                               doCellCalibration=True,
                               calibTool=calibHcells,
                               addCellNoise = False, filterCellNoise = False,
                               OutputLevel = INFO,
                               hits="overlaidHCalHits",
                               cells="HCalCells")

createExtHcells = CreateCaloCells("CreateExtHCaloCells",
                                  doCellCalibration=True,
                                  calibTool=calibHcells,
                                  addCellNoise = False, filterCellNoise = False,
                                  OutputLevel = INFO,
                                  hits="ExtHCalHits",
                                  cells="ExtHCalCells")

createcellsEndcap = CreateCaloCells("CreateCaloCellsEndcap",
                                    doCellCalibration=True,
                                    calibTool=calibcellsEndcap,
                                    addCellNoise=False, filterCellNoise=False,
                                    OutputLevel=INFO,
                                    hits="ECalEndcapHits",
                                    cells="ECalEndcapCells")

# additionally for HCal
from Configurables import CreateVolumeCaloPositions, CreateCellCaloPositions
volPositionsEcal = CreateVolumeCaloPositions("volPositionsEcal", OutputLevel = INFO)
volPositionsEcal.hits.Path = "ECalCells"
volPositionsEcal.positionedHits.Path = "ECalPositions"

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

from Configurables import RedoSegmentation
resegmentEcal = RedoSegmentation("ReSegmentationEcal",
                                 # old bitfield (readout)
                                 oldReadoutName = "ECalBarrelEta",
                                 # specify which fields are going to be altered (deleted/rewritten)
                                 # oldSegmentationIds = [],
                                 # new bitfield (readout), with new segmentation
                                 newReadoutName = ecalReadoutName,
                                 debugPrint = 10,
                                 OutputLevel = INFO,
                                 inhits = "ECalPositions",
                                 outhits = "newECalCells")

from Configurables import CreateCellCaloPositions
# Ecal cell positions
positionsEcal2 = CreateCellCaloPositions("cellPositionsEcal", readoutName=ecalReadoutName, OutputLevel = INFO)
positionsEcal2.hits.Path = "newECalCells"
positionsEcal2.positionedHits.Path = "cellECalPositions"


out = PodioOutput("out", 
                  OutputLevel=INFO)
out.outputCommands = ["keep *", "drop ECalHits", "drop HCalHits", "drop ExtHCalHits", "drop ECalEndcapHits", 
                      "drop ECalPositionedHits", "drop HCalPositionedHits", "drop extHCalPositionedHits", "drop ECalEndcapPositionedHits", 
                      "drop ECalCells", "drop HCalCells", "drop ExtHCalCells", "drop ECalPositions", "drop HCalPositions",
                      "drop overlaidECalHits", "drop overlaidHCalHits", "drop overlaidextHCalHits",
                      "drop overlaidECalPositionedHits",  "drop overlaidHCalPositionedHits"]
out.filename = "output.root"
if args.outputfile != '':
    out.filename = args.outputfile


#CPU information
from Configurables import AuditorSvc, ChronoAuditor
chra = ChronoAuditor()
audsvc = AuditorSvc()
audsvc.Auditors = [chra]
pythia8gen.AuditExecute = True
hepmc_converter.AuditExecute = True
geantsim.AuditExecute = True
createEcells.AuditExecute = True
createHcells.AuditExecute = True
createExtHcells.AuditExecute = True
#createcellsEndcap.AuditExecute = True
volPositionsEcal.AuditExecute = True
volPositionsHcal.AuditExecute = True
positionsHcal.AuditExecute = True
positionsExtHcal.AuditExecute = True
positionsEcal2.AuditExecute = True
out.AuditExecute = True

ApplicationMgr(
    TopAlg = [pythia8gen, 
              hepmc_converter,
              genfilter,
              geantsim,
              overlay,
              createEcells,createHcells,createExtHcells,
              volPositionsEcal,volPositionsHcal,
              positionsHcal,positionsExtHcal,
              resegmentEcal,
              positionsEcal2,
              out
              ],
    EvtSel = 'NONE',
    EvtMax   = int(num_events),
    ExtSvc = [podioevent, geoservice, geantservice, audsvc],
 )

