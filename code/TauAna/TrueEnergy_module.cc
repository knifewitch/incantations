// Tau analysis module by J. Hewes <jhewes15@fnal.gov>

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Services/Optional/TFileService.h"

#include "nusimdata/SimulationBase/MCTruth.h"

#include "TH1.h"

namespace dune {

  class TrueEnergy : public art::EDAnalyzer {

  public:

    explicit TrueEnergy(fhicl::ParameterSet const& pset);
    //virtual ~TrueEnergy();

    void analyze(const art::Event& evt);

  private:

    TH1* fTrueEnergy;

  }; // class dune::TrueEnergy
} // namespace dune

dune::TrueEnergy::TrueEnergy(fhicl::ParameterSet const& pset) :
  EDAnalyzer(pset)
{
  // Get a handle to the TFile service

  art::ServiceHandle<art::TFileService> tfs;

  // Create true energy hist
  
  fTrueEnergy = tfs->make<TH1D>("true_energy",
    ";True #nu energy (GeV);Events", 50, 0, 16); 

} // dune::TrueEnergy::TrueEnergy()

void dune::TrueEnergy::analyze(const art::Event& evt)
{
  // First, get an art handle to MC truth information

  art::Handle<std::vector<simb::MCTruth>> mct;
  std::vector<art::Ptr<simb::MCTruth>> truthlist;
  if (evt.getByLabel("generator", mct))
    art::fill_ptr_vector(truthlist, mct);
  simb::MCTruth truth(*(truthlist[0]));

  // Then fill true neutrino energy

  const simb::MCNeutrino & nu = truthlist[0]->GetNeutrino();
  double e_true = nu.Nu().E();
  fTrueEnergy->Fill(e_true);

} // dune::TrueEnergy::analyze()

namespace dune {

  DEFINE_ART_MODULE(TrueEnergy)

}
