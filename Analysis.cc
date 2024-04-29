#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include <sstream>
#include <iomanip> 


using namespace std;


//object definition
class Analysis : public edm::one::EDAnalyzer<> {
public:

  //constructor, function is called when new object is created
  explicit Analysis(const edm::ParameterSet& conf);

  //destructor, function is called when object is destroyed
  ~Analysis();

  //edm filter plugin specific functions
  virtual void beginJob();

  double muonEnergy( const pat::Muon& );
  double invariantMass ( const pat::Muon& , const pat::Muon& );
  double invariantMassTest( const pat::Muon& , const pat::Muon& );
  void muonParams ( const pat::Muon& );

  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();

private:

  edm::ParameterSet theConfig;
  unsigned int theEventCount;
  TH1D *histo;
  TH1D *hID ;
  TH1D *hDisp ;
  TH1D *hCan ;
  TH1D *hVz ;

  edm::EDGetTokenT< vector<pat::Muon> >            theMuonToken;
  edm::EDGetTokenT< vector<pat::Muon> >            theDisplacedMuonToken;
  edm::EDGetTokenT< vector<pat::PackedCandidate> > theCandidateToken;

};


Analysis::Analysis(const edm::ParameterSet& conf)
  : theConfig(conf), theEventCount(0)
{
  cout <<" CTORXX" << endl;
  theMuonToken          = consumes< vector<pat::Muon> >            ( theConfig.getParameter<edm::InputTag>("muonSrc"));
  theDisplacedMuonToken = consumes< vector<pat::Muon> >            ( theConfig.getParameter<edm::InputTag>("displacedSrc"));
  theCandidateToken     = consumes< vector<pat::PackedCandidate> > ( theConfig.getParameter<edm::InputTag>("candidateSrc"));
}

Analysis::~Analysis()
{
  cout <<" DTOR" << endl;
}

void Analysis::beginJob()
{
  //create a histogram
  histo = new TH1D("histo","Invariant~mass~of~\\mu^+\\mu^-~pairs; \\sqrt{s}~[GeV]; Counts", 2400 , 0. , 12. );

  hID = new TH1D("hID","Particle ID; ID; Counts", 225 , -9.5, 215.5 );

  hDisp = new TH1D("hDisp","Invariant~mass~of~\\mu^+\\mu^-~pairs~(SlimmedDisplacedMuons); \\sqrt{s}~[GeV]; Counts", 2400 , 0. , 12. );

  hCan = new TH1D("hCan","Invariant~mass~of~\\mu^+\\mu^-~pairs~(Packed Candidate); \\sqrt{s}~[GeV]; Counts", 2400 , 0. , 12. );

  hVz = new TH1D("hVz","The~distance~along~the~Z-axis~between~vertices~for ~\\mu^+\\mu^-~pairs; Distance~[1/GeV]; Counts", 2400 , 0. , 12. );

  cout << "HERE Analysis::beginJob()" << endl;
}

void Analysis::endJob()
{
  //make a new Root file
  TFile myRootFile( theConfig.getParameter<std::string>("outHist").c_str(), "RECREATE");

  //write histogram data
  histo -> Write();
  hID   -> Write();
  hDisp -> Write();
  hCan  -> Write();
  hVz   -> Write();

  myRootFile.Close();

  delete histo;
  delete hID;
  delete hDisp;
  delete hVz;
  delete hCan;

  cout << "HERE Cwiczenie::endJob()" << endl;
}

double Analysis::muonEnergy(const pat::Muon& muon ){
  
  double mMu = 0.1056583755; // GeV

  double energy = sqrt( pow(muon.p() , 2.) + mMu * mMu);

  return energy ;
}

double Analysis::invariantMass(const pat::Muon& muon1, const pat::Muon& muon2 ){

  TLorentzVector p4_1(muon1.px(), muon1.py(), muon1.pz(), muonEnergy( muon1 ));
  TLorentzVector p4_2(muon2.px(), muon2.py(), muon2.pz(), muonEnergy( muon2 ));

  TLorentzVector sum = p4_1 + p4_2;

  return sum.M() ;
}

double Analysis::invariantMassTest(const pat::Muon& muon1, const pat::Muon& muon2 ){
  
  double inv = sqrt( pow( muonEnergy(muon1) + muonEnergy(muon2), 2.) - pow( muon1.pz() + muon2.pz() , 2.) - pow( muon1.px() + muon2.px() , 2.) - pow( muon1.py() + muon2.py() , 2.));
  
  return inv ;
}

void Analysis::muonParams ( const pat::Muon& muon ){

    std::cout << std::setw(10) << std::left << "q/m/px/py/pz/vx/vy/vz: " << muon.charge() << " " << muon.mass() << " " 
              << muon.px() << " " << muon.py() << " " << muon.pz() << " " 
              << muon.vx() << " " << muon.vy() << " " << muon.vz() << " \n" ;

}

void Analysis::analyze(
    const edm::Event& ev, const edm::EventSetup& es)
{
  std::cout << " -------------------------------- HERE Cwiczenie::analyze "<< std::endl;
  bool debug = true;

  const vector<pat::Muon> & muons = ev.get(theMuonToken);
  const vector<pat::Muon> & dispMuons = ev.get(theDisplacedMuonToken);
  
  if (debug) std::cout <<" number of muons (slimmedMuons): " << muons.size() <<std::endl;
  if (debug) std::cout <<" number of muons (slimmedDisplacedMuons): " << dispMuons.size() <<std::endl;

  /*
  for (const auto & muon : muons) {
    if (debug) std::cout <<" reco muon pt: "<<muon.pt() << ",   muon q: " << muon.charge() <<std::endl;
    //histo->Fill(muon.eta());
  }
  */

  if ( theEventCount == 1 ){


    const vector<pat::PackedCandidate> & candidates = ev.get(theCandidateToken);
    unsigned int numberOfParticles = 0;

    vector<pat::PackedCandidate> packedMuons;

    //std::cout <<  "vx / vy / vz / mass: \n";

    for (const auto & can : candidates) {

      
      
      // find muons in candidates
      if ( std::abs( can.pdgId() ) == 13 ){

        packedMuons.push_back( can );
      }

      hID -> Fill( abs( can.pdgId() ) );

      //if (debug) std::cout <<" Particle id: "<<can.pdgId() << std::endl;
      numberOfParticles++;
    }

      
    cout << "Number of particles in thi event: " << numberOfParticles << endl;
    cout << "Number of muons in PackedCandidate / slimmed Muons / slimmed Displaced Muons: " 
         << packedMuons.size() << "/" << muons.size() << "/" << dispMuons.size() << endl; 


  }

  // SlimmedMuons - histogram
  for (unsigned int i = 0; i < muons.size(); ++i) {
        
        if( theEventCount == 1 ){

        //cout << "----------------------SLIMMED MUONS:\n";
        //muonParams( muons[i] );

        }

        // Invariant mass of muon pairs + vz difference

        for (unsigned int j = i + 1; j < muons.size(); ++j) {

  
            // invariant mass
          if (muons[i].charge() * muons[j].charge() < 0) {

              // vz
              if ( muons[i].pt() >= 2 && muons[j].pt() >= 2 ) {

                  hVz -> Fill( abs( muons[i].vz() - muons[j].vz() ) );

              }
              
              histo -> Fill( invariantMass( muons[i], muons[j] ) );
          }
        }
  }


  // SlimmedDisplacedMuons - histogram
  for (unsigned int i = 0; i < dispMuons.size(); ++i) {

        /*if( theEventCount == 1 ){

          //cout << "------------SLIMMED DISPLACED MUONS:\n";
          //muonParams( dispMuons[i] );

        }*/

        for (unsigned int j = i + 1; j < dispMuons.size(); ++j) {
            
            if (dispMuons[i].charge() * dispMuons[j].charge() < 0) {
                
                hDisp -> Fill(invariantMass( dispMuons[i] , dispMuons[j] ) );
            }
        }
  }


  //write std io
  cout <<"*** Analyze event: " << ev.id()<<" analysed event count:"<<++theEventCount << endl;
}

DEFINE_FWK_MODULE(Analysis);

