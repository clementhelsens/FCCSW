#include "CombinedCaloTopoCluster.h"

// FCCSW
#include "DetInterface/IGeoSvc.h"

// datamodel
#include "datamodel/PositionedCaloHitCollection.h"
#include "datamodel/CaloCluster.h"
#include "datamodel/CaloHitCollection.h"
#include "datamodel/CaloHit.h"

// DD4hep
#include "DD4hep/LCDD.h"
#include "DD4hep/Readout.h"

#include "DetCommon/DetUtils.h"

#include <algorithm>

DECLARE_ALGORITHM_FACTORY(CombinedCaloTopoCluster)

CombinedCaloTopoCluster::CombinedCaloTopoCluster(const std::string& name, ISvcLocator* svcLoc)
: GaudiAlgorithm(name, svcLoc) {
  declareProperty("ecalCells", m_ecalCells, "calo/ecalCells (input)");
  declareProperty("hcalCells", m_hcalCells, "calo/hcalCells (input)");
  declareProperty("ECalPositions", m_ecalPositions, "ecalCell positions (input)");
  declareProperty("HCalPositions", m_hcalPositions, "hcalCell positions (input)");
  declareProperty("ecalReadoutName", m_ecalReadoutName);
  declareProperty("hcalReadoutName", m_hcalReadoutName);
  // the default value to calculate the position of clusters
  
}

StatusCode CombinedCaloTopoCluster::initialize() {
  if (GaudiAlgorithm::initialize().isFailure()) {
    return StatusCode::FAILURE;
  }
  m_geoSvc = service("GeoSvc");
  if (!m_geoSvc) {
    error() << "Unable to locate Geometry Service. "
            << "Make sure you have GeoSvc and SimSvc in the right order in the configuration." << endmsg;
    return StatusCode::FAILURE;
  }
  // check if readouts exist
  if (m_geoSvc->lcdd()->readouts().find(m_ecalReadoutName) == m_geoSvc->lcdd()->readouts().end()) {
    error() << "Readout <<" << m_ecalReadoutName << ">> does not exist." << endmsg;
    return StatusCode::FAILURE;
  }
  if (m_geoSvc->lcdd()->readouts().find(m_hcalReadoutName) == m_geoSvc->lcdd()->readouts().end()) {
    error() << "Readout <<" << m_hcalReadoutName << ">> does not exist." << endmsg;
    return StatusCode::FAILURE;
  }
  // retrieve PhiEta segmentation
  m_ecalSegmentation = dynamic_cast<DD4hep::DDSegmentation::GridPhiEta*>(
      m_geoSvc->lcdd()->readout(m_ecalReadoutName).segmentation().segmentation());
  if (m_ecalSegmentation == nullptr) {
    error() << "There is no phi-eta segmentation in the electromagnetic calorimeter." << endmsg;
    return StatusCode::FAILURE;
  }
  m_hcalSegmentation = dynamic_cast<DD4hep::DDSegmentation::GridPhiEta*>(
      m_geoSvc->lcdd()->readout(m_hcalReadoutName).segmentation().segmentation());
  if (m_hcalSegmentation == nullptr) {
    error() << "There is no segmentation in the hadronic calorimeter." << endmsg;
    return StatusCode::FAILURE;
  }

  return StatusCode::SUCCESS;
}

bool myFunction (fcc::CaloCluster hit1, fcc::CaloCluster hit2){
  return hit1.core().energy < hit2.core().energy;
}

StatusCode CombinedCaloTopoCluster::execute() {
  const fcc::PositionedCaloHitCollection* cellsEcal = m_ecalPositions.get();
  const fcc::PositionedCaloHitCollection* cellsHcal = m_hcalPositions.get();

  // Finds seeds and fills the list of allCells
  CombinedCaloTopoCluster::findingSeeds(cellsEcal, m_seedThr_ecal, m_firstSeedsEcal, m_allCellsEcal);
  CombinedCaloTopoCluster::findingSeeds(cellsHcal, m_seedThr_hcal, m_firstSeedsHcal, m_allCellsHcal);
  
  std::cout << "Number of seeds found in ECAL = " << m_firstSeedsEcal.size() << std::endl;
  std::cout << "Number of seeds found in HCAL = " << m_firstSeedsHcal.size() << std::endl;
  
  //decending order of seeds
  std::sort (m_firstSeedsEcal.begin(), m_firstSeedsEcal.end(), myFunction);
  std::sort (m_firstSeedsHcal.begin(), m_firstSeedsHcal.end(), myFunction);

  //  CombinedCaloTopoCluster::buildingProtoCluster(m_firstSeedsEcal, m_allCellsEcal, cellsEcal);

  return StatusCode::SUCCESS;
}

void CombinedCaloTopoCluster::findingSeeds(const fcc::PositionedCaloHitCollection* cells, double threshold, std::vector<fcc::CaloCluster>& seeds, std::vector<fcc::CaloHit>& allCells){
  std::cout << "seed threshold  = " << threshold << std::endl;
  for (const auto& cell : *cells) {
    allCells.push_back(cell);
    if (cell.core().energy / dd4hep::MeV > threshold)
      fcc::CaloCluster cluster;
    cluster.addHits(cell);
      seeds.push_back(cluster);
  }
} 


void CombinedCaloTopoCluster::buildingProtoCluster(std::vector<fcc::CaloCluster> seeds, std::vector<fcc::CaloHit>& allCells, fcc::CaloClusterCollection* clusterCollection){

  // Take readout bitfield decoder from GeoSvc                                                             
  auto decoderEcal = m_geoSvc->lcdd()->readout(m_ecalReadoutName).idSpec().decoder();
  auto decoderHcal = m_geoSvc->lcdd()->readout(m_hcalReadoutName).idSpec().decoder();
 
  const std::vector<std::pair<int, int>> m_fieldExtremesEcal = det::utils::bitfieldExtremes((*decoderEcal), m_fieldNamesEcal);
  const std::vector<std::pair<int, int>> m_fieldExtremesHcal = det::utils::bitfieldExtremes((*decoderHcal), m_fieldNamesHcal);
 
  // Loop over every seed in Ecal to get neighbouring cells
  std::vector<fcc::CaloCluster>::iterator itSeed = seeds.begin();
  uint iSeeds = 0;
  while(itSeed != seeds.end(), ++itSeed, iSeeds++){
    auto seed = *itSeed;
    // the seed is first entry in collection to form proto-clusters
    // global position of the cell
    uint64_t id = seed.core().cellId;
    auto position = seed.core().position;
    auto edmPos = fcc::Point();
    edmPos.x = position.x() * 10.;
    edmPos.y = position.y() * 10.;
    edmPos.z = position.z() * 10.;

    auto positionedHit = clusterCollection->create(edmPos, seed.core());
   
    // remove seed from cell list
    allCells.erase(std::remove(allCells.begin(), allCells.end(), id), allCells.end());
    
    // retrieve the neighbours of the seed
    std::vector<uint64_t> Neighbours = det::utils::neighbours((*decoderEcal), m_fieldNamesEcal, m_fieldExtremesEcal, seed.core().cellId);
    std::vector<uint64_t>::iterator itNeighbour = Neighbours.begin();
    uint iSeeds = 0;

    // Loop over all neighbours
    while(itNeighbour != Neighbours.end(), ++itNeighbour, iNeighbours++){
      // Find the neighbours of the seeds in the ECal cell collection 
      bool foundNeighbour = false; 
      for (const auto& cell : *cells) {
	foundNeighbour = std::find(cells.core().cellID == ;
 std::find (clusterCollection.begin(), clusterCollection.end(), cell ) != clusterCollection.end();

      if (foundNeighbour){
//	// Check if Neighbour is already assigned to another proto-cluster
//	bool foundHit = false;
//	std::vector<fcc::CaloClusterCollection>::iterator itProtoClusters;
//	while (itProtoClusters != m_seedsEcalCollection.end()){
//	  auto clusterCollection = *itProtoClusters;
//	  foundHit = std::find (clusterCollection.begin(), clusterCollection.end(), cell ) != clusterCollection.end();
//	}
	// If the neighouring cell possesses more than threshold energy, the cell is added to proto-cluster collection
	if (neighbourHit.core().energy / dd4hep::MeV > m_neighbourThr_ecal ){
	// 	neighbourHit = m_seedsEcalCollection.at(iSeeds).create();   
	// }
	//      else if (neighbourHit.core().energy / dd4hep::MeV > 0){ // TODO: find the seed that is closest and assign!!! 
	if (foundHit == false){
	  //neighbourHit = seeds.at(iSeeds).create();   
	  // remove the neighbour from cell list
	  allCellIDs.erase(std::remove(allCellIDs.begin(), allCellIDs.end(), neighbourHit.core().cellId), allCellIDs.end());
	}
      }
      else {
	std::cout << "neighbouring cell could not be found in cell collection!! Soemthing is wrong" << std::endl;
	break;
      }
    }  
  }
}

StatusCode CombinedCaloTopoCluster::finalize() { 
  return GaudiAlgorithm::finalize(); 
}


