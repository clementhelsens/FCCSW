#include "CombinedCaloTopoCluster.h"

// FCCSW
#include "DetCommon/DetUtils.h"
#include "DetInterface/IGeoSvc.h"

// datamodel
#include "datamodel/CaloCluster.h"
#include "datamodel/CaloClusterCollection.h"
#include "datamodel/CaloHit.h"
#include "datamodel/CaloHitCollection.h"
#include "datamodel/PositionedCaloHit.h"
#include "datamodel/PositionedCaloHitCollection.h"

// DD4hep
#include "DD4hep/LCDD.h"
#include "DD4hep/Readout.h"

#include <algorithm>
#include <map>

DECLARE_ALGORITHM_FACTORY(CombinedCaloTopoCluster)

CombinedCaloTopoCluster::CombinedCaloTopoCluster(const std::string& name, ISvcLocator* svcLoc)
    : GaudiAlgorithm(name, svcLoc) {
  declareProperty("ecalCells", m_ecalCells, "calo/ecalCells (input)");
  declareProperty("hcalCells", m_hcalCells, "calo/hcalCells (input)");
  declareProperty("ecalReadoutName", m_ecalReadoutName);
  declareProperty("hcalReadoutName", m_hcalReadoutName);
  declareProperty("ecalFieldNames", m_ecalFieldNames);
  declareProperty("hcalFieldNames", m_hcalFieldNames);
  declareProperty("ecalFieldValues", m_ecalFieldValues);
  declareProperty("hcalFieldValues", m_hcalFieldValues);
  declareProperty("clusters", m_preClusterCollection, "Handle for calo clusters (output collection)");
  declareProperty("geometryToolEcal", m_geoToolEcal, "Handle for the geometry tool of the Ecal");
  declareProperty("geometryToolHcal", m_geoToolHcal, "Handle for the geometry tool of the Hcal");
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

  // Take readout bitfield decoder from GeoSvc
  auto decoderEcal = m_geoSvc->lcdd()->readout(m_ecalReadoutName).idSpec().decoder();
  auto decoderHcal = m_geoSvc->lcdd()->readout(m_hcalReadoutName).idSpec().decoder();
  
  const std::vector<std::pair<int, int>> m_fieldExtremesEcal =
    det::utils::bitfieldExtremes((*decoderEcal), m_fieldNamesEcal);
  const std::vector<std::pair<int, int>> m_fieldExtremesHcal =
    det::utils::bitfieldExtremes((*decoderHcal), m_fieldNamesHcal);
  
  // Geometry settings
  if (!m_geoToolEcal.retrieve()) {
    error() << "Unable to retrieve the ECAL geometry tool!!!" << endmsg;
    return StatusCode::FAILURE;
  }
  if (!m_geoToolHcal.retrieve()) {
    error() << "Unable to retrieve the HCAL geometry tool!!!" << endmsg;
    return StatusCode::FAILURE;
  }
   std::unordered_map<uint64_t, double> ecalDetCells;
  std::unordered_map<uint64_t, double> hcalDetCells;
  
  // Prepare map of all existing cells in calorimeter
  info() << "Initialising ECAL CellID-Neighbours map, this can take some time!" << endmsg;
  StatusCode sc_prepareEcalCells = m_geoToolEcal->prepareEmptyCells(ecalDetCells);
  if (sc_prepareEcalCells.isFailure()) {
    error() << "Unable to create empty cells!" << endmsg;
    return StatusCode::FAILURE;
  }
  // Filling of map for each Cell to its neighbours
  for (std::unordered_map<uint64_t, double>::iterator itCellID = ecalDetCells.begin(); itCellID != ecalDetCells.end(); itCellID++) {
    uint64_t cellID =  itCellID->first;
    NeighboursMapEcal[cellID] = det::utils::neighbours((*decoderEcal), m_fieldNamesEcal, m_fieldExtremesEcal, cellID);
  }
  
  // Prepare map of all existing cells in calorimeter
  info() << "Initialising HCAL CellID-Neighbours map, this can take some time!" << endmsg;
  StatusCode sc_prepareHcalCells = m_geoToolHcal->prepareEmptyCells(hcalDetCells);
  if (sc_prepareHcalCells.isFailure()) {
    error() << "Unable to create empty cells!" << endmsg;
    return StatusCode::FAILURE;
  }
  info() << "Number of cells in HCAL: " << hcalDetCells.size() << endmsg;
   // Filling of map for each Cell to its neighbours
  for (std::unordered_map<uint64_t, double>::iterator itCellID = hcalDetCells.begin(); itCellID != hcalDetCells.end(); itCellID++) {
    uint64_t cellID =  itCellID->first;
    NeighboursMapHcal[cellID] = det::utils::neighbours((*decoderHcal), m_fieldNamesHcal, m_fieldExtremesHcal, cellID);
  }

  return StatusCode::SUCCESS;
}

bool myFunction(fcc::PositionedCaloHit hit1, fcc::PositionedCaloHit hit2) { return hit1.core().energy < hit2.core().energy; }

StatusCode CombinedCaloTopoCluster::execute() {
  const fcc::PositionedCaloHitCollection* ecalCells = m_ecalCells.get();
  const fcc::PositionedCaloHitCollection* hcalCells = m_hcalCells.get();

  // Finds seeds and fills the list of allCells
  CombinedCaloTopoCluster::findingSeeds(ecalCells, m_seedThr_ecal, firstSeedsEcal, allCellsEcal);
  CombinedCaloTopoCluster::findingSeeds(hcalCells, m_seedThr_hcal, firstSeedsHcal, allCellsHcal);

  info() << "Number of seeds found in ECAL = " << firstSeedsEcal.size() << endmsg;
  info() << "Number of seeds found in HCAL = " << firstSeedsHcal.size() << endmsg;
  info() << "All Cells in ECAL             = " << allCellsEcal.size() << endmsg;
  info() << "All Cells in HCAL             = " << allCellsHcal.size() << endmsg;
  
  // decending order of seeds
  std::sort(firstSeedsEcal.begin(), firstSeedsEcal.end(), myFunction);
  std::sort(firstSeedsHcal.begin(), firstSeedsHcal.end(), myFunction);

  fcc::CaloClusterCollection* edmClusters = m_preClusterCollection.createAndPut();
  //while (allCellsEcal.size() > 0) {
  CombinedCaloTopoCluster::buildingProtoCluster(NeighboursMapEcal, firstSeedsEcal, allCellsEcal, edmClusters);
  //}
  // while (firstSeedsHcal.size() > 0) {
  CombinedCaloTopoCluster::buildingProtoCluster(NeighboursMapHcal, firstSeedsHcal, allCellsHcal, edmClusters);
    //}
  info() << "number of reconstructed clusters: "<< edmClusters->size() << endmsg;
  return StatusCode::SUCCESS;
}

void CombinedCaloTopoCluster::findingSeeds(const fcc::PositionedCaloHitCollection* cells,
					   double threshold,
                                           std::vector<fcc::PositionedCaloHit>& seeds,
                                           std::map<uint64_t, fcc::PositionedCaloHit>& allCells) {
  //info() << "cells : " << cells->size() << endmsg;
  info() << "seed threshold  = " << threshold << "MeV " << endmsg;
  for (const auto& cell : *cells) {
    allCells[cell.cellId()] = cell;
    if (cell.core().energy / dd4hep::MeV > threshold) {
      seeds.push_back(cell);
    }
  }
}

void CombinedCaloTopoCluster::buildingProtoCluster(const std::map<uint64_t, std::vector<uint64_t> >& neighboursMap,
						   std::vector<fcc::PositionedCaloHit>& seeds,
                                                   std::map<uint64_t, fcc::PositionedCaloHit>& allCells,
                                                   fcc::CaloClusterCollection* preClusterCollection) {
  std::map<unsigned int, fcc::CaloCluster> clusterOfCell;
  // New seed list for neighbours of original seeds
  std::vector<fcc::PositionedCaloHit> newSeeds;
  std::vector<fcc::PositionedCaloHit> cellsToFormCluster;

  // Loop over every seed in Cal to get neighbouring cells
  auto itSeed = seeds.begin();
  uint iSeeds = 0;
  debug() << "seeds to loop over : " << seeds.size() << endmsg;
  while (itSeed != seeds.end(), ++itSeed, iSeeds++) {
    auto seedCell = *itSeed;
    uint64_t id = seedCell.cellId();
    cellsToFormCluster.push_back(seedCell);
    
    // remove cell added to cluster from (free) cell list
    allCells.erase(id);

    // retrieve the neighbours of the seed
    std::vector<uint64_t> Neighbours = neighboursMap.find(id)->second;
    std::vector<uint64_t>::iterator itNeighbour = Neighbours.begin();
    debug() << "Number of found neigbours: " << Neighbours.size() << endmsg;
    
    int iNeighbours = 0;
    // Loop over all neighbour cellID to check if it was hit
    while (itNeighbour != Neighbours.end(), ++itNeighbour, iNeighbours++) {

      // Find the neighbours of the seeds in the Cal cell collection
      bool foundNeighbour = false;
      uint64_t neighbourID = *itNeighbour;
      auto itAllCells = allCells.find(neighbourID);
      if (itAllCells != allCells.end()) {
        debug() << "Found neighbour with CellID: " << itAllCells->first << endmsg;
        foundNeighbour = true;
        auto neighbouringCell = allCells[neighbourID];

        // Check if Neighbour is already assigned to another proto-cluster
        bool foundCellInAnotherCluster = false;
        auto it = clusterOfCell.find(neighbourID);
        if (it != clusterOfCell.end()) {
          foundCellInAnotherCluster = true;
          debug() << "neighbour found in another cluster, move on" << endmsg;
          continue;
        } else {
          // check if the neighouring cell possesses more than threshold energy
          // if yes, the cell is added to proto-cluster collection
          // and erased from cell list
          if (neighbouringCell.core().energy / dd4hep::MeV > m_neighbourThr_ecal) {
            cellsToFormCluster.push_back(neighbouringCell);
            newSeeds.push_back(neighbouringCell);
            allCells.erase(neighbourID);
          } else if (neighbouringCell.core().energy / dd4hep::MeV > 0) {
            // TODO: find the seed that is closest and assign!!!
          }
        }
      } else {
        debug() << "neighbour not found in cell list" << endmsg;
      }     
    }
  }
    
//  auto cluster = fcc::CaloCluster();
//  double posX = 0.;
//  double posY = 0.;
//  double posZ = 0.;
//  double energy = 0.;
//
//  for (int i=0; i<cellsToFormCluster.size(); i++) {
//    posX += cellsToFormCluster[i].position().x * cellsToFormCluster[i].energy();
//    posY += cellsToFormCluster[i].position().y * cellsToFormCluster[i].energy();
//    posZ += cellsToFormCluster[i].position().z * cellsToFormCluster[i].energy();
//    energy += cellsToFormCluster[i].core().energy;
//    cluster.addhits(cellsToFormCluster[i]);
//  }
//  cluster.core().energy(energy);
//  cluster.core().position.x(posX/energy);
//  cluster.core().position.y(posY/energy);
//  cluster.core().position.z(posZ/energy);
//
//  auto edmCluster = preClusterCollection->create();
//  auto& edmClusterCore = edmCluster.core();
//  edmClusterCore.position.x = cluster.core().position.x;
//  edmClusterCore.position.y = cluster.core().position.y;
//  edmClusterCore.position.z = cluster.core().position.z;
//  edmClusterCore.energy = cluster.core().energy;
//  debug() << "Cluster energy:     " << cluster.core().energy << endmsg;
//  debug() << "Cluster position x: " << cluster.core().position.x << endmsg;
//  debug() << "Cluster position y: " << cluster.core().position.y << endmsg;
//  debug() << "Cluster position z: " << cluster.core().position.z << endmsg;
  
  // Entries in seeds cleared and replace by the list of new seeds (neighbouring cells)
  seeds.clear();
  seeds = newSeeds;
  debug() << "Number of seeds for next iteration: " << newSeeds.size() << endmsg;
}

StatusCode CombinedCaloTopoCluster::finalize() { return GaudiAlgorithm::finalize(); }
