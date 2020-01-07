// Tau analysis module by J. Hewes <jhewes15@fnal.gov>

// Framework includes

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Event.h"
#include "larcore/Geometry/Geometry.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

// Data product includes

#include "nusimdata/SimulationBase/MCTruth.h"
#include "dune/FDSensOpt/FDSensOptData/EnergyRecoOutput.h"

// ROOT includes

#include "TTree.h"
#include "TH1.h"
#include "TH2.h"

// c++ includes

#include <vector>

namespace dune {

  class LeptonInfo : public art::EDAnalyzer {

  public:

    explicit LeptonInfo(fhicl::ParameterSet const& pset);
    //virtual ~LeptonInfo();

    void analyze(const art::Event& evt);

  private:

    void ProcessChildren(const simb::MCParticle* p);

    std::set<int> fLeptonChild;
    double fNuEnergy;
    double fLepEnergy;
    bool fIsCC;
    int fInteractionType;

    TTree* fTree;

  }; // class dune::LeptonInfo
} // namespace dune

dune::LeptonInfo::LeptonInfo(fhicl::ParameterSet const& pset) :
  EDAnalyzer(pset),
  fTree(nullptr)
{
  // Get a handle to the TFile service

  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("lep_tree", "lep tree");

  fTree->Branch("lepton_child", &fLeptonChild);
  fTree->Branch("nu_energy", &fNuEnergy);
  fTree->Branch("lep_energy", &fLepEnergy);
  fTree->Branch("is_cc", &fIsCC);
  fTree->Branch("interaction_type", &fInteractionType);

} // dune::LeptonInfo::LeptonInfo()

void dune::LeptonInfo::analyze(const art::Event& evt)
{
  // First, get an art handle to MC truth information
  art::Handle<std::vector<simb::MCTruth>> mct;
  std::vector<art::Ptr<simb::MCTruth>> truthlist;
  if (evt.getByLabel("generator", mct))
    art::fill_ptr_vector(truthlist, mct);
  simb::MCTruth truth(*(truthlist[0]));

  // Fill event level information
  simb::MCParticle lep = truth.GetNeutrino().Lepton();
  fNuEnergy = truth.GetNeutrino().Nu().E();
  fLepEnergy = lep.E();
  fIsCC = (truth.GetNeutrino().CCNC() == simb::kCC);
  int lepTrackID = -1;
  std::vector<std::pair<int, double>> lepKids;
  if (fIsCC) {
    // Find CC tau in GENIE truth
    for (int iP = 0; iP < truth.NParticles(); ++iP) {
      simb::MCParticle p = truth.GetParticle(iP);
      if (p.StatusCode() == 3 && p.PdgCode() == 15) {
        lepTrackID = p.TrackId();
        break;
      }
    }
    // Find children of tau
    for (int iP = 0; iP < truth.NParticles(); ++iP) {
      simb::MCParticle p = truth.GetParticle(iP);
      if (p.StatusCode() == 1 && p.Mother() == lepTrackID) {
        // std::cout << "Lepton child here with PDG " << p.PdgCode() << ", energy " << p.E() << " and track ID " << p.TrackId() << std::endl;
        lepKids.push_back(std::make_pair(p.PdgCode(), p.E()));
      }
    }
  }

  // for (int iP = 0; iP < truth.NParticles(); ++iP) {
    // simb::MCParticle p = truth.GetParticle(iP);
    // std::cout << "True particle has PDG " << p.PdgCode() << ", status code " << p.StatusCode() << ", track ID " << p.TrackId() << " and parent " << p.Mother() << std::endl;
  // }

  // Get G4 truth information
  art::Handle<std::vector<simb::MCParticle>> h_mcpart;
  if (!evt.getByLabel("largeant", h_mcpart))
    throw std::runtime_error("Could not get MC particle information!");
  std::vector<simb::MCParticle> const& mcpart(*h_mcpart);

  // Find all G4 primaries
  // MCParticle* lepPrim = nullptr;
  for (const simb::MCParticle &p : mcpart) {
    if (p.Mother() == 0) {
      // std::cout << "Primary particle found with PDG " << p.PdgCode() << ", energy " << p.E() << " and track ID " << p.TrackId() << std::endl;
      for (auto kid : lepKids) {
        if (p.PdgCode() == kid.first && (p.E() - kid.second < 1e-5)) {
          // std::cout << "Identified this one as a child, calling recursive loop." << std::endl;
          ProcessChildren(&p);
          break;
        }
      }
    }
  }

  // Is this a neutral current event?
  if (!fIsCC) fInteractionType = 0;

  else {

    // Loop over all MC particles
    bool foundLepton = false;
    for (int i = 0; i < truthlist[0]->NParticles(); ++i) {
      simb::MCParticle p = truthlist[0]->GetParticle(i);

      // We only care about final state particles here
      if (p.StatusCode() == 1) {

        // Is there an electron?
        if (abs(p.PdgCode()) == 11) {
          fInteractionType = 1;
          foundLepton = true;
        }

        // Is there a muon?
        if (abs(p.PdgCode()) == 13) {
          fInteractionType = 2;
          foundLepton = true;
        }
      } // if p.StatusCode() == 1
    } // for i in truthlist

    if (!foundLepton) fInteractionType = 4;

  } // if cc

  // Fill the tree and then reset values
  fTree->Fill();
  fLeptonChild.clear();

} // dune::LeptonInfo::analyze()

void dune::LeptonInfo::ProcessChildren(const simb::MCParticle* p) {

  art::ServiceHandle<cheat::ParticleInventoryService> pi;

  // Add particle track ID to children
  std::cout << "Adding particle " << p->TrackId() << " to lepton children." << std::endl;
  fLeptonChild.insert(p->TrackId());
  // If there are no children, quit out of the recursive loop
  if (p->NumberDaughters() == 0) return;
  for (int i = 0; i < p->NumberDaughters(); ++i) {
    int childID = p->Daughter(i);
    // Recursive loop over each child
    ProcessChildren(pi->TrackIdToParticle_P(childID));
  }

}

namespace dune {

  DEFINE_ART_MODULE(LeptonInfo)

}
