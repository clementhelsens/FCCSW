
### \file
### \ingroup SimulationExamples
### | **input (alg)**               | other algorithms                 |                                                         |                        |                                  | **output (alg)**                              |
### |-------------------------------|----------------------------------|---------------------------------------------------------|------------------------|----------------------------------|-----------------------------------------------|
### | generating single particle events from a given list of types, with momentum, phi and theta from a given range, saving to HepMC | convert `HepMC::GenEvent` to EDM | geometry parsed from XML (TestHCal.xml) by DD4hep using GeoSvc | FTFP_BERT physics list | saving HCal hits | write the EDM output to ROOT file using PODIO |



from Gaudi.Configuration import *

from Configurables import FCCDataSvc
## Data service
podioevent = FCCDataSvc("EventDataSvc")

from Configurables import ParticleGunAlg
## Particle Gun using MomentumRangeParticleGun tool and FlatSmearVertex
# MomentumRangeParticleGun generates particles of given type(s) within given momentum, phi and theta range
# FlatSmearVertex smears the vertex with uniform distribution
gen = ParticleGunAlg("ParticleGun", ParticleGunTool="MomentumRangeParticleGun", VertexSmearingToolPGun="FlatSmearVertex")
gen.DataOutputs.hepmc.Path = "hepmc"

from Configurables import Gaudi__ParticlePropertySvc
## Particle service
# list of possible particles is defined in ParticlePropertiesFile
ppservice = Gaudi__ParticlePropertySvc("ParticlePropertySvc", ParticlePropertiesFile="Generation/data/ParticleTable.txt")

from Configurables import HepMCConverter
## Reads an HepMC::GenEvent from the data service and writes a collection of EDM Particles
hepmc_converter = HepMCConverter("Converter")
hepmc_converter.DataInputs.hepmc.Path="hepmc"
hepmc_converter.DataOutputs.genparticles.Path="allGenParticles"
hepmc_converter.DataOutputs.genvertices.Path="allGenVertices"

from Configurables import GeoSvc
## DD4hep geometry service
# Parses the given xml file
geoservice = GeoSvc("GeoSvc", detectors=['file:Detector/DetFCChhTrackerSimple/compact/FCChh_TrackerStandalone.xml'],
                    OutputLevel = DEBUG)

from Configurables import G4SimSvc
## Geant4 service
# Configures the Geant simulation: geometry, physics list and user actions
geantservice = G4SimSvc("G4SimSvc", detector='G4DD4hepDetector', physicslist="G4FtfpBert",
                        actions="G4FullSimActions")

from Configurables import G4SimAlg, G4SaveTrackerHits, G4SaveCalHits
## Geant4 algorithm
# Translates EDM to G4Event, passes the event to G4, writes out outputs via tools
# first, create a tool that saves the tracker hits
# Name of that tool in GAUDI is "XX/YY" where XX is the tool class name ("G4SaveTrackerHits")
# and YY is the given name ("saveTrackerHits")
savetrackertool = G4SaveTrackerHits("saveTrackerHits")
savetrackertool.DataOutputs.trackClusters.Path = "clusters"
savetrackertool.DataOutputs.trackHits.Path = "hits"
savetrackertool.DataOutputs.trackHitsClusters.Path = "hitClusterAssociation"
# and a tool that saves the calorimeter hits with a name "G4SaveCalHits/saveHCalHits"
savehcaltool = G4SaveCalHits("saveHCalHits", caloType = "HCal")
savehcaltool.DataOutputs.caloClusters.Path = "caloClusters"
savehcaltool.DataOutputs.caloHits.Path = "caloHits"
# next, create the G4 algorithm, giving the list of names of tools ("XX/YY")
geantsim = G4SimAlg("G4SimAlg",
                        outputs= ["G4SaveTrackerHits/saveTrackerHits", "G4SaveCalHits/saveHCalHits"])
geantsim.DataInputs.genParticles.Path="allGenParticles"

from Configurables import PodioOutput
out = PodioOutput("out",
                   OutputLevel=DEBUG)
out.outputCommands = ["keep *"]

from Configurables import ApplicationMgr
ApplicationMgr( TopAlg=[gen, hepmc_converter, geantsim, out],
                EvtSel='NONE',
                EvtMax=1,
                ## order is important, as GeoSvc is needed by G4SimSvc
                ExtSvc=[podioevent, geoservice, geantservice, ppservice],
                OutputLevel=DEBUG
 )
