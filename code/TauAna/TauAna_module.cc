// Tau analysis module by J. Hewes <jhewes15@fnal.gov>

// Framework includes

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Event.h"
#include "larcore/Geometry/Geometry.h"

// Data product includes

#include "nusimdata/SimulationBase/MCTruth.h"
// #include "nusimdata/SimulationBase/GTruth.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "dune/FDSensOpt/FDSensOptData/EnergyRecoOutput.h"

// ROOT includes

#include "TTree.h"
#include "TH1.h"
#include "TH2.h"

// c++ includes

#include <vector>

namespace dune {

  enum InteractionType {
    kNue,
    kNumu,
    kNueLike,
    kNumuLike,
    kHadronic,
    kNC
  };

  enum EnergyType {
    kNueEnergy,
    kNumuEnergy,
    kCaloEnergy,
    kTrueEnergy
  };

  class TauAna : public art::EDAnalyzer {

  public:

    explicit TauAna(fhicl::ParameterSet const& pset);
    //virtual ~TauAna();

    void analyze(const art::Event& evt);
    void beginJob();

  private:

    std::vector<double> energies;

    TTree* fTree;
    std::vector<TH2*> energy_hists;

    TH1* fTrueEnergy;

    // Energy migration matrices

    TH2* fContainedNCEnergy;
    TH2* fContainedCCEnergy;
    TH2* fContainedEnergy;

    TH2* fAllNCEnergy;
    TH2* fAllCCEnergy;
    TH2* fAllEnergy;

    TH2* fNDContainedNCEnergy;
    TH2* fNDContainedCCEnergy;
    TH2* fNDContainedEnergy;

    bool fContained;
    bool fNDContained;

    int fNuPDG;
    bool fIsNC;
    int fInteractionType;
    int fGenieMode;
    int fGenieInteractionType;

    bool fVertexContained;

    double fVertexX;
    double fVertexY;
    double fVertexZ;

    // Atmospheric information!

    bool fVertexContainedAtmospheric;
    bool fEventContainedAtmospheric;

    // Directionality information

    double fCosThetaVis;

    geo::BoxBoundedGeo fFullVolume;
    geo::BoxBoundedGeo fNDVolume;
    geo::BoxBoundedGeo fAtmosphericVolume;

    const std::vector<std::string> interaction_types = { "nuelike", "numulike", "hadronic", "nc" };
    const std::vector<std::string> reco_types = { "nue", "numu", "calo", "true" };

    void ResetEnergies();

  }; // class dune::TauAna
} // namespace dune

dune::TauAna::TauAna(fhicl::ParameterSet const& pset) :
  EDAnalyzer(pset),
  fTree(nullptr),
  fFullVolume(),
  fNDVolume(),
  fAtmosphericVolume()
{
  // Get a handle to the TFile service

  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("tau_tree", "tau tree");

  fTree->Branch("contained", &fContained);
  fTree->Branch("nu_pdg", &fNuPDG);
  fTree->Branch("is_nc", &fIsNC);
  fTree->Branch("interaction_type", &fInteractionType);
  fTree->Branch("genie_mode", &fGenieMode);
  fTree->Branch("genie_interaction_type", &fGenieInteractionType);

  fTree->Branch("vertex_contained", &fVertexContained);

  fTree->Branch("vertex_x", &fVertexX);
  fTree->Branch("vertex_y", &fVertexY);
  fTree->Branch("vertex_z", &fVertexZ);

  fTree->Branch("vertex_contained_atmospheric", &fVertexContainedAtmospheric);
  fTree->Branch("event_contained_atmospheric",  &fEventContainedAtmospheric);

  fTree->Branch("cos_theta_vis", &fCosThetaVis);

  // Set size of energy vector

  energies.resize(reco_types.size());

  // Set size of energy histogram

  energy_hists.resize(reco_types.size() - 1);

  for (unsigned int i = 0; i < reco_types.size(); ++i) {
    std::string branch_name = "e_nu_" + reco_types[i];
    fTree->Branch(branch_name.c_str(), &(energies[i]));
    if (reco_types[i] != "true") {
      std::string hist_name = "e_true_vs_" + reco_types[i];
      energy_hists[i] = tfs->make<TH2D>(hist_name.c_str(),
        ";True #nu energy (GeV);Reco #nu energy (GeV)", 50, 0, 100, 50, 0, 100);
    }
  }

  // Create true energy hist
  
  fTrueEnergy = tfs->make<TH1D>("true_energy",
    ";True #nu energy (GeV);Events", 50, 0, 16); 

  // Contained energy migration matrices

  fContainedNCEnergy = tfs->make<TH2D>("contained_nc_true_vs_calo_energy",
    ";True #nu energy (GeV);Reco #nu energy (GeV)", 50, 0, 100, 50, 0, 100);
  fContainedCCEnergy = tfs->make<TH2D>("contained_cc_true_vs_calo_energy",
    ";True #nu energy (GeV);Reco #nu energy (GeV)", 50, 0, 100, 50, 0, 100);
  fContainedEnergy = tfs->make<TH2D>("contained_true_vs_calo_energy",
    ";True #nu energy (GeV);Reco #nu energy (GeV)", 50, 0, 100, 50, 0, 100);

  // All energy migration matrices

  fAllNCEnergy = tfs->make<TH2D>("all_nc_true_vs_calo_energy",
    ";True #nu energy (GeV);Reco #nu energy (GeV)", 50, 0, 100, 50, 0, 100);
  fAllCCEnergy = tfs->make<TH2D>("all_cc_true_vs_calo_energy",
    ";True #nu energy (GeV);Reco #nu energy (GeV)", 50, 0, 100, 50, 0, 100);
  fAllEnergy = tfs->make<TH2D>("all_true_vs_calo_energy",
    ";True #nu energy (GeV);Reco #nu energy (GeV)", 50, 0, 100, 50, 0, 100);

  // ND contained energy migration matrices
  fNDContainedNCEnergy = tfs->make<TH2D>("nd_contained_nc_true_vs_calo_energy",
    ";True #nu energy (GeV);Reco #nu energy (GeV)", 50, 0, 100, 50, 0, 100);
  fNDContainedCCEnergy = tfs->make<TH2D>("nd_contained_cc_true_vs_calo_energy",
    ";True #nu energy (GeV);Reco #nu energy (GeV)", 50, 0, 100, 50, 0, 100);
  fNDContainedEnergy = tfs->make<TH2D>("nd_contained_true_vs_calo_energy",
    ";True #nu energy (GeV);Reco #nu energy (GeV)", 50, 0, 100, 50, 0, 100);

} // dune::TauAna::TauAna()

void dune::TauAna::beginJob()
{
  art::ServiceHandle<geo::Geometry> geom;

  double min_x = 999999;
  double max_x = -999999;
  double min_y = 999999;
  double max_y = -999999;
  double min_z = 999999;
  double max_z = -999999;

  // Loop over the TPCs and get the full boundaries

  for (geo::TPCGeo const& tpc : geom->IterateTPCs()) {
    geo::BoxBoundedGeo box = tpc.ActiveBoundingBox();
    if (min_x > box.MinX()) min_x = box.MinX();
    if (max_x < box.MaxX()) max_x = box.MaxX();
    if (min_y > box.MinY()) min_y = box.MinY();
    if (max_y < box.MaxY()) max_y = box.MaxY();
    if (min_z > box.MinZ()) min_z = box.MinZ();
    if (max_z < box.MaxZ()) max_z = box.MaxZ();
  }
  std::cout << "-------------------------" << std::endl;
  std::cout << "10kt dimensions:" << std::endl;
  std::cout << "X: " << min_x << " -> " << max_x << std::endl;
  std::cout << "Y: " << min_y << " -> " << max_y << std::endl;
  std::cout << "Z: " << min_z << " -> " << max_z << std::endl;
  std::cout << "-------------------------" << std::endl;

  // Set bounding box dimensions

  fFullVolume.SetBoundaries(min_x, max_x, min_y, max_y, min_z, max_z);

  fNDVolume.SetBoundaries(-150., 150., -350., 350., min_z, max_z);

  fAtmosphericVolume.SetBoundaries(-600., 600., -600., 600., 0., 1200.);

} // function TauAna::beginJob

void dune::TauAna::analyze(const art::Event& evt)
{
  // First, get an art handle to MC truth information

  art::Handle<std::vector<simb::MCTruth>> mct;
  std::vector<art::Ptr<simb::MCTruth>> truthlist;
  if (evt.getByLabel("generator", mct))
    art::fill_ptr_vector(truthlist, mct);
  simb::MCTruth truth(*(truthlist[0]));

  // Next, get a handle on MC reco information

  art::Handle<std::vector<sim::MCTrack>> h_mctrack;
  if (!evt.getByLabel("mcreco", h_mctrack))
    throw std::runtime_error("Could not get MCTrack information!");
  std::vector<sim::MCTrack> const & mctrack(*h_mctrack);

  art::Handle<std::vector<sim::MCShower>> h_mcshower;
  if (!evt.getByLabel("mcreco", h_mcshower))
    throw std::runtime_error("Could not get MCShower information!");
  std::vector<sim::MCShower> const& mcshower(*h_mcshower);

  // Check containment

  fContained = true;
  fNDContained = true;
  fEventContainedAtmospheric = true;
  art::ServiceHandle<geo::Geometry> geom;

  // Get the true neutrino vertex and see if it's contained in
  // various geometry configurations

  TVector3 nu_vtx(truth.GetNeutrino().Nu().Position().Vect());
  fVertexContained = fFullVolume.ContainsPosition(nu_vtx);
  fVertexContainedAtmospheric = fAtmosphericVolume.ContainsPosition(nu_vtx);

  // Save the vertex position

  fVertexX = nu_vtx.X();
  fVertexY = nu_vtx.Y();
  fVertexZ = nu_vtx.Z();

  TVector3 nu_direction(truth.GetNeutrino().Nu().Momentum().Vect().Unit());

  // Get G4 truth information

  art::Handle<std::vector<simb::MCParticle>> h_mcpart;
  if (!evt.getByLabel("largeant", h_mcpart))
    throw std::runtime_error("Could not get MC particle information!");
  std::vector<simb::MCParticle> const& mcpart(*h_mcpart);

  // Make a map of TrackID to MCParticle

  std::map<unsigned int, simb::MCParticle> part_map;
  for (simb::MCParticle p : mcpart)
    part_map.insert(std::make_pair(p.TrackId(), p));

  // Loop over MCTracks

  TVector3 visible_direction(0, 0, 0);

  for (sim::MCTrack mct : mctrack) {

    // Get associated MCParticle

    simb::MCParticle p = part_map[mct.TrackID()];

    // Check if MCTrack is contained

    TVector3 start = mct.Start().Position().Vect();
    TVector3 end   = mct.End().Position().Vect();

    // MCTrack visible momentum

    TVector3 start_momentum = mct.Start().Momentum().Vect();
    TVector3 end_momentum = mct.End().Momentum().Vect();

    visible_direction += start_momentum - end_momentum;

    // Is event contained in full volume?

    if (fContained) {
      bool start_contained = fFullVolume.ContainsPosition(start);
      bool end_contained   = fFullVolume.ContainsPosition(end);
      if (!start_contained || !end_contained) fContained = false;
    }

    // Is event contained in ND-equivalent volume?

    if (fNDContained) {
      bool start_contained_nd = fNDVolume.ContainsPosition(start);
      bool end_contained_nd = fNDVolume.ContainsPosition(end);
      if (!start_contained_nd || !end_contained_nd) fNDContained = false;
    }

    // Is event contained in atmospheric cube volume?

    if (fEventContainedAtmospheric) {
      bool start_contained_atmos = fAtmosphericVolume.ContainsPosition(start);
      bool end_contained_atmos = fAtmosphericVolume.ContainsPosition(end);
      if (!start_contained_atmos || !end_contained_atmos) fEventContainedAtmospheric = false;
    }

    if (!fContained && !fNDContained && !fEventContainedAtmospheric) break;

  } // for each MCTrack

  // Loop over MCShowers

  for (sim::MCShower mcs : mcshower) {

    // Get associated MCParticle

    simb::MCParticle p = part_map[mcs.TrackID()];

    // Check if MCShower is contained

    TVector3 start = mcs.Start().Position().Vect();
    TVector3 end   = mcs.End().  Position().Vect();

    // MCShower deposited momentum

    TVector3 start_momentum = mcs.Start().Momentum().Vect();
    TVector3 end_momentum = mcs.End().Momentum().Vect();
    // std::cout << "MCShower has start momentum (" << start_momentum.X()
    //   << "," << start_momentum.Y() << "," << start_momentum.Z()
    //   << ") and end momentum (" << end_momentum.X() << ","
    //   << end_momentum.Y() << "," << end_momentum.Z() << ")" << std::endl;
    visible_direction += start_momentum - end_momentum;

    // Is event contained in full volume?

    if (fContained) {
      bool start_contained = fFullVolume.ContainsPosition(start);
      bool end_contained   = fFullVolume.ContainsPosition(end);
      if (!start_contained || !end_contained) fContained = false;
    }

    // Is event contained in ND-equivalent volume?

    if (fNDContained) {
      bool start_contained_nd = fNDVolume.ContainsPosition(start);
      bool end_contained_nd = fNDVolume.ContainsPosition(end);
      if (!start_contained_nd || !end_contained_nd) fNDContained = false;
    }

    // Is event contained in atmospheric cube volume?

    if (fEventContainedAtmospheric) {
      bool start_contained_atmos = fAtmosphericVolume.ContainsPosition(start);
      bool end_contained_atmos = fAtmosphericVolume.ContainsPosition(end);
      if (!start_contained_atmos || !end_contained_atmos) fEventContainedAtmospheric = false;
    }

    // If we've determined the event is uncontained by all metrics, may as well quit now

    if (!fContained && !fNDContained && !fEventContainedAtmospheric) break;

  } // for each mcshower

  // Now get costheta of visible energy
  fCosThetaVis = nu_direction.Dot(visible_direction.Unit());

  // We want to know what kind of tau neutrino interaction
  // this is, so here we check.

  const simb::MCNeutrino & nu = truthlist[0]->GetNeutrino();
  double e_true = nu.Nu().E();
  fTrueEnergy->Fill(e_true);
  fIsNC = nu.CCNC();
  fNuPDG = nu.Nu().PdgCode();
  InteractionType interaction_type;

  // Is this a neutral current event?

  if (fIsNC) interaction_type = kNC;

  // Otherwise it's a charged current event.

  else if (fabs(fNuPDG) == 12) interaction_type = kNue;
  else if (fabs(fNuPDG) == 14) interaction_type = kNumu;

  else {

    // Loop over all MC particles

    bool found_lepton = false;
    for (int i = 0; i < truthlist[0]->NParticles(); ++i) {
      simb::MCParticle p = truthlist[0]->GetParticle(i);

      // We only care about final state particles here

      if (p.StatusCode() == 1) {

        // Is there an electron?

        if (abs(p.PdgCode()) == 11) {
          interaction_type = kNueLike;
          found_lepton = true;
        }

        // Is there a muon?

        if (abs(p.PdgCode()) == 13) {
          interaction_type = kNumuLike;
          found_lepton = true;
        }
      } // if p.StatusCode() == 1
    } // for i in truthlist

    if (!found_lepton) interaction_type = kHadronic;

  } // if cc

  // Set interaction type

  fInteractionType = interaction_type;

  fGenieMode = nu.Mode();
  fGenieInteractionType = nu.InteractionType();

  // Set true energy

  energies[kTrueEnergy] = e_true;

  // Get nue reconstructed energy

  art::Handle<dune::EnergyRecoOutput> handle_energy_reco_nue;
  evt.getByLabel("energynue", handle_energy_reco_nue);
  if (handle_energy_reco_nue->recoMethodUsed == 1) {
    double e_reco = handle_energy_reco_nue->fNuLorentzVector.E();
    energies[kNueEnergy] = e_reco;
    energy_hists[kNueEnergy]->Fill(e_true, e_reco);
  }

  // Get numu reconstructed energy

  art::Handle<dune::EnergyRecoOutput> handle_energy_reco_numu;
  evt.getByLabel("energynumu", handle_energy_reco_numu);
  if (handle_energy_reco_numu->recoMethodUsed == 2) {
    double e_reco = handle_energy_reco_numu->fNuLorentzVector.E();
    energies[kNumuEnergy] = e_reco;
    energy_hists[kNumuEnergy]->Fill(e_true, e_reco);
  }

  // Get calorimetric reconstructed energy

  art::Handle<dune::EnergyRecoOutput> handle_energy_reco_calo;
  evt.getByLabel("energynutau", handle_energy_reco_calo);
  double e_reco = handle_energy_reco_calo->fNuLorentzVector.E();
  energies[kCaloEnergy] = e_reco;
  energy_hists[kCaloEnergy]->Fill(e_true, e_reco);

  // Fill calo reco energy migration matrices

  fAllEnergy->Fill(e_true, e_reco);
  if (fIsNC) fAllNCEnergy->Fill(e_true, e_reco);
  else fAllCCEnergy->Fill(e_true, e_reco);

  if (fContained) {
    fContainedEnergy->Fill(e_true, e_reco);
    if (fIsNC) fContainedNCEnergy->Fill(e_true, e_reco);
    else fContainedCCEnergy->Fill(e_true, e_reco);
  }

  if (fNDContained) {
    fNDContainedEnergy->Fill(e_true, e_reco);
    if (fIsNC) fNDContainedNCEnergy->Fill(e_true, e_reco);
    else fNDContainedCCEnergy->Fill(e_true, e_reco);
  }

  // Fill the tree and then reset values

  fTree->Fill();
  ResetEnergies();

} // dune::TauAna::analyze()

void dune::TauAna::ResetEnergies() {

  for (double & val : energies)
    val = -1;

} // dune::TauAna::ResetEnergies()

namespace dune {

  DEFINE_ART_MODULE(TauAna)

}
