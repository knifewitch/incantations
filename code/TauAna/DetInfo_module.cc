// Tau analysis module by J. Hewes <jhewes15@fnal.gov>

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"

#include "larcore/Geometry/Geometry.h"

namespace dune {

  class DetInfo : public art::EDAnalyzer {

  public:

    explicit DetInfo(fhicl::ParameterSet const& pset);
    //virtual ~DetInfo();

    void beginJob();
    void analyze(const art::Event& evt);

  private:

  }; // class dune::DetInfo
} // namespace dune

dune::DetInfo::DetInfo(fhicl::ParameterSet const& pset) :
  EDAnalyzer(pset)
{

} // dune::DetInfo::DetInfo()

void dune::DetInfo::beginJob()
{

  art::ServiceHandle<geo::Geometry> geom;

  double min_x = 999999;
  double max_x = -999999;
  double min_y = 999999;
  double max_y = -999999;
  double min_z = 999999;
  double max_z = -999999;

  int it_tpc = 0;

  for (geo::TPCGeo const& tpc : geom->IterateTPCs()) {

    geo::BoxBoundedGeo bbox = tpc.ActiveBoundingBox();

    double t_min_x = bbox.MinX();
    double t_max_x = bbox.MaxX();

    double t_min_y = bbox.MinY();
    double t_max_y = bbox.MaxY();

    double t_min_z = bbox.MinZ();
    double t_max_z = bbox.MaxZ();

    ++it_tpc;

    std::cout << "-------------------------" << std::endl;
    std::cout << "TPC " << it_tpc << " dimensions:" << std::endl;
    std::cout << "X: " << t_min_x << " -> " << t_max_x << std::endl;
    std::cout << "Y: " << t_min_y << " -> " << t_max_y << std::endl;
    std::cout << "Z: " << t_min_z << " -> " << t_max_z << std::endl;
    std::cout << "-------------------------" << std::endl;

    if (min_x > t_min_x) min_x = t_min_x;
    if (max_x < t_max_x) max_x = t_max_x;

    if (min_y > t_min_y) min_y = t_min_y;
    if (max_y < t_max_y) max_y = t_max_y;

    if (min_z > t_min_z) min_z = t_min_z;
    if (max_z < t_max_z) max_z = t_max_z;
  }
  std::cout << "-------------------------" << std::endl;
  std::cout << "10kt dimensions:" << std::endl;
  std::cout << "X: " << min_x << " -> " << max_x << std::endl;
  std::cout << "Y: " << min_y << " -> " << max_y << std::endl;
  std::cout << "Z: " << min_z << " -> " << max_z << std::endl;
  std::cout << "-------------------------" << std::endl;

} // dune::DetInfo::beginJob

void dune::DetInfo::analyze(const art::Event& evt)
{

} // dune::DetInfo::analyze()

namespace dune {

  DEFINE_ART_MODULE(DetInfo)

}
