#ifndef RECCALORIMETER_COMBINEDCALOTOPOCLUSTER_H
#define RECCALORIMETER_COMBINEDCALOTOPOCLUSTER_H

// from Gaudi
#include "GaudiAlg/GaudiAlgorithm.h"
#include "GaudiKernel/ToolHandle.h"

// FCCSW
#include "FWCore/DataHandle.h"
#include "DetSegmentation/GridPhiEta.h"

class IGeoSvc;

// datamodel
namespace fcc {
  class CaloHitCollection;
  class CaloHit;
  class CaloClusterCollection;
  class CaloCluster;
  class PositionedCaloHitCollection;
}

namespace DD4hep {
namespace DDSegmentation {
    class Segmentation;
}
}

/** @class CombinedCaloTopoClusterAlgorithm Reconstruction/RecCalorimeter/src/components/CombinedCaloTopoCluster.h CombinedCaloTopoCluster.h
 *
 *  Algorithm building the topological clusters for the energy reconstruction and as input for ATLAS PFA.
 *
 *  @author Coralie Neubueser
 */

class CombinedCaloTopoCluster : public GaudiAlgorithm {
 public:
  CombinedCaloTopoCluster(const std::string& name, ISvcLocator* svcLoc);

  StatusCode initialize();

  StatusCode execute();

  StatusCode finalize();
  /**  Find cells with a signal to noise ratio > 6 for ECal and > 4 for HCal, following ATLAS note ATL-LARG-PUB-2008-002.
   *   For simulation without electronic and pile-up noise, the average noise levels are taken as reference for seeding (1.5 and 3.5MeV/cell for E and HCAL, electronic noise only), (2.5 and 100MeV/cell for E and HCAL, added pile-up).
   *   @return list of seed cells ("proto-clusters").
   */
  virtual void findingSeeds(const fcc::CaloHitCollection* cells, double threshold, std::vector<fcc::CaloCluster>& seeds, std::map<uint64_t, fcc::CaloHit>& allCells);
  /**Building proto-clusters
  */
  virtual void buildingProtoCluster(const std::vector<fcc::CaloCluster>& seeds, std::map<uint64_t, fcc::CaloHit>& allCells, fcc::CaloClusterCollection* preClusterCollection);

 private:
  /// Handle for electromagnetic calorimeter cells (input collection)
  DataHandle<fcc::CaloHitCollection> m_ecalCells{"ecalCells", Gaudi::DataHandle::Reader, this};
  /// Handle for hadronic calorimeter cells (input collection)
  DataHandle<fcc::CaloHitCollection> m_hcalCells{"hcalCells", Gaudi::DataHandle::Reader, this};
  /// Handle for electromagnetic calorimeter cells (input collection)
  DataHandle<fcc::PositionedCaloHitCollection> m_ecalPositions{"ECalPositions", Gaudi::DataHandle::Reader, this};
  /// Handle for hadronic calorimeter cells (input collection)
  DataHandle<fcc::PositionedCaloHitCollection> m_hcalPositions{"HCalPositions", Gaudi::DataHandle::Reader, this};
  // Pre-cluster collections
  DataHandle<fcc::CaloClusterCollection> m_preClusterCollection{"calo/clusters", Gaudi::DataHandle::Writer, this};
  /// Pointer to the geometry service
  SmartIF<IGeoSvc> m_geoSvc;
  /// Name of the electromagnetic calorimeter readout
  Gaudi::Property<std::string> m_ecalReadoutName{this, "ecalReadoutName", "name of the ecal readout"};
  /// Name of thehadronic calorimeter readout
  Gaudi::Property<std::string> m_hcalReadoutName{this, "hcalReadoutName", "name of the hcal readout"};

  /// PhiEta segmentation of the electromagnetic detector (owned by DD4hep)
  DD4hep::DDSegmentation::GridPhiEta* m_ecalSegmentation;
  /// PhiEta segmentation of the hadronic detector (owned by DD4hep)
  DD4hep::DDSegmentation::GridPhiEta* m_hcalSegmentation;

  /// Seed threshold Ecal
  Gaudi::Property<double> m_seedThr_ecal{this, "seedThresholdEcal", 7.5, "seed threshold estimate [MeV]"};
  /// Seed threshold hcal
  Gaudi::Property<double> m_seedThr_hcal{this, "seedThresholdHcal", 11.5, "seed threshold estimate [MeV]"};
  /// Seed threshold Ecal
  Gaudi::Property<double> m_neighbourThr_ecal{this, "neighbourThresholdEcal", 3, "neighbour threshold estimate [MeV]"};
  /// Seed threshold hcal
  Gaudi::Property<double> m_neighbourThr_hcal{this, "neighbourThresholdHcal", 3.5, "neighbour threshold estimate [MeV]"};

  /// allCells
  std::map<uint64_t, fcc::CaloHit> allCellsEcal;
  std::map<uint64_t, fcc::CaloHit> allCellsHcal;
  /// First list of CaloCells above seeding threshold 
  std::vector<fcc::CaloCluster> firstSeedsEcal;
  std::vector<fcc::CaloCluster> firstSeedsHcal;

  /// Collection of CaloCells above neighbouring threshold associated to seeds (used for proto-clustering)
  std::vector<fcc::CaloHitCollection> seedsEcalCollection;
  std::vector<fcc::CaloHitCollection> seedsHcalCollection;

  /// Name of the bit-fields (in the readout) describing the volume ECAL                      
  const std::vector<std::string> m_fieldNamesEcal{"active","cell"};
  std::string m_activeVolumeNameEcal = "LAr_sensitive";
   /// Name of the bit-fields (in the readout) describing the volume HCAL                      
  const std::vector<std::string> m_fieldNamesHcal{"layer","row"};
  std::string m_activeVolumeNameHcal = "Polystyrene_tile";
  //std::vector<std::pair<int, int>> m_fieldExtremesHcal;
  //hcalVolumeName = ["moduleVolume", "wedgeVolume", "layerVolume", "modCompVolume"]
  //hcalIdentifierName = ["module", "row", "layer", "tile"]
  //hcalFieldValues=[8]

};

#endif /* RECCALORIMETER_COMBINEDCALOTOPOCLUSTER_H */
