
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

savetrackertool = SimG4SaveTrackerHits("saveTrackerHits", readoutNames = ["TrackerBarrelReadout", "TrackerEndcapReadout"])
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
genfilter = GenParticleFilter("StableParticles", accept=[1], OutputLevel=DEBUG)
genfilter.allGenParticles.Path = "allGenParticles"
genfilter.filteredGenParticles.Path = "GenParticles"


from Configurables import SimG4PrimariesFromEdmTool
particle_converter = SimG4PrimariesFromEdmTool("EdmConverter")
particle_converter.genParticles.Path = "GenParticles"


# next, create the G4 algorithm, giving the list of names of tools ("XX/YY")                                                                                                                                                   
geantsim = SimG4Alg("SimG4Alg",
                       outputs= ["SimG4SaveCalHits/saveECalHits", 
                                 "SimG4SaveCalHits/saveHCalHits", 
                                 "SimG4SaveCalHits/saveExtHCalHits",
                                 "SimG4SaveCalHits/saveECalEndcapHits", 
                                 "SimG4SaveTrackerHits/saveTrackerHits"],
                       eventProvider=particle_converter,
                       OutputLevel=DEBUG)




out = PodioOutput("out", 
                  OutputLevel=DEBUG)

#keeping the hits for MinBias only (to mix pileup OTF)
out.outputCommands = ["keep *",  "drop ECalEndcapHits",
                      "drop ECalEndcapPositionedHits"]

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
out.AuditExecute = True

ApplicationMgr(
    TopAlg = [pythia8gen, 
              hepmc_converter,
              genfilter,
              geantsim,
              out
              ],
    EvtSel = 'NONE',
    EvtMax   = int(num_events),
    ExtSvc = [podioevent, geoservice, geantservice, audsvc],
 )

