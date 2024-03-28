#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/PatCandidates/interface/Muon.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include <sstream>


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
  double InvariantMass( const pat::Muon& , const pat::Muon&);
  double InvariantMassStupid( const pat::Muon& , const pat::Muon&);
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();

private:

  edm::ParameterSet theConfig;
  unsigned int theEventCount;
  TH1D *histo;

  edm::EDGetTokenT< vector<pat::Muon> > theMuonToken;

};


Analysis::Analysis(const edm::ParameterSet& conf)
  : theConfig(conf), theEventCount(0)
{
  cout <<" CTORXX" << endl;
  theMuonToken = consumes< vector<pat::Muon> >( theConfig.getParameter<edm::InputTag>("muonSrc"));
}

Analysis::~Analysis()
{
  cout <<" DTOR" << endl;
}

void Analysis::beginJob()
{
  //create a histogram
  histo =new TH1D("histo","test; X; #events",150, 0., 200.);
  cout << "HERE Analysis::beginJob()" << endl;
}

void Analysis::endJob()
{
  //make a new Root file
  TFile myRootFile( theConfig.getParameter<std::string>("outHist").c_str(), "RECREATE");
  //write histogram data
  histo->Write();
  myRootFile.Close();
  delete histo;
  cout << "HERE Cwiczenie::endJob()" << endl;
}

double Analysis::InvariantMass(const pat::Muon& muon1, const pat::Muon& muon2 ){
  
  TLorentzVector p4_1(muon1.px(), muon1.py(), muon1.pz(), muon1.energy());
  TLorentzVector p4_2(muon2.px(), muon2.py(), muon2.pz(), muon2.energy());
  //cout << "px: " << muon1.px() << " " << muon2.px() <<endl;
  //cout << "py: " << muon1.py() << " " << muon2.py() <<endl;
  //cout << "pz: " << muon1.pz() << " " << muon2.pz() <<endl;
  //cout << "energy: " << muon1.energy() << " " << muon2.energy() <<endl;

  TLorentzVector sum = p4_1 + p4_2;

  return sum.M() ;
}

double Analysis::InvariantMassStupid(const pat::Muon& muon1, const pat::Muon& muon2 ){
  
  double inv = sqrt( pow( muon1.energy() + muon2.energy() , 2.) - pow( muon1.pz() + muon2.pz() , 2.) - pow( muon1.px() + muon2.px() , 2.) - pow( muon1.py() + muon2.py() , 2.));
  
  return inv ;
}

void Analysis::analyze(
    const edm::Event& ev, const edm::EventSetup& es)
{
  std::cout << " -------------------------------- HERE Cwiczenie::analyze "<< std::endl;
  bool debug = true;
  const vector<pat::Muon> & muons = ev.get(theMuonToken);
  if (debug) std::cout <<" number of muons: " << muons.size() <<std::endl;
  
  for (const auto & muon : muons) {
    if (debug) std::cout <<" reco muon pt: "<<muon.pt() << ",   muon q: " << muon.charge() <<std::endl;
    //histo->Fill(muon.eta());
  }


  for (unsigned int i = 0; i < muons.size(); ++i) {
        for (unsigned int j = i + 1; j < muons.size(); ++j) {
            
            if (muons[i].charge() * muons[j].charge() < 0) {
                //cout <<" invariant mass: " <<  InvariantMass(muons[i], muons[j]) <<std::endl;
                
                
                histo -> Fill(InvariantMass(muons[i], muons[j]));
            }
        }
    }


  //write std io
  cout <<"*** Analyze event: " << ev.id()<<" analysed event count:"<<++theEventCount << endl;
}

DEFINE_FWK_MODULE(Analysis);

