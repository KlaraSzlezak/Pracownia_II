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

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
// #include "TrackingTools/IPTools/interface/IPTools.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include <sstream>
#include <iomanip> 
#include <utility>


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

  double muonEnergy     ( const pat::Muon& );
  double invariantMass  ( const pat::Muon& , const pat::Muon& );
  double invariantMass  ( const pat::Muon& , const pat::Muon& , const pat::Muon& , const pat::Muon& );
  double invariantMass  ( const pat::PackedCandidate& , const pat::PackedCandidate& , double );
  double pTJPsi         ( const pat::Muon& , const pat::Muon& );
  double invariantMass3 ( const pat::Muon& , const pat::Muon&, const pat::PackedCandidate& , double& ); //only kaon kaon mu mu
  double invariantMass4 ( const pat::Muon& , const pat::Muon&, const pat::PackedCandidate& , const pat::PackedCandidate& , double& );


  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();

private:

//PDG:
  double mMu = 0.1056583755; // GeV
  double mK = 0.493677; //GeV
  double mPi = 0.13957039; //GeV
  double zAv = 0.038164028445384024; // cm
  double mPhi =1.019461; //GeV

  edm::ParameterSet theConfig;
  unsigned int theEventCount;

  TH1D *histo;
  TH1D *hInv;
  TH1D *hInvG;
  TH1D *hID ;
  TH1D *hDisp ;
  TH1D *hCan ;
  TH1D *hVz ;
  TH1D *hVzPC ;
  TH1D *hProb ;
  TH1D *hJpsiProb ;
  TH1D *hVz01;
  TH1D *hVz08;
  //TH1D *hInvK;
  TH1D *hkkProb ;
  TH1D *hMuMuKK ;
  TH1D *hInvMuMuK ;
  TH1D *hKK ;
  TH1D *hTest ;

  edm::EDGetTokenT< vector<pat::Muon> >            theMuonToken;
  edm::EDGetTokenT< vector<pat::Muon> >            theDisplacedMuonToken;
  edm::EDGetTokenT< vector<pat::PackedCandidate> > theCandidateToken;

  // common vertex
  edm::EDGetTokenT< vector<reco::Track> > theDsplTrkToken; 
  edm::EDGetTokenT< vector<reco::Track> > theDsplMuToken;
  edm::EDGetTokenT< vector<reco::Vertex> > theVertexToken, theVertexWithBSToken; //?
  edm::EDGetTokenT< vector<reco::VertexCompositePtrCandidate> > theVertexCPCToken; 
  edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> theTrackBuilderToken;

};


Analysis::Analysis(const edm::ParameterSet& conf)
  : theConfig(conf), theEventCount(0)
{
  cout <<" CTORXX" << endl;
  theMuonToken          = consumes< vector<pat::Muon> >            ( theConfig.getParameter<edm::InputTag>("muonSrc"));
  theDisplacedMuonToken = consumes< vector<pat::Muon> >            ( theConfig.getParameter<edm::InputTag>("displacedSrc"));
  theCandidateToken     = consumes< vector<pat::PackedCandidate> > ( theConfig.getParameter<edm::InputTag>("candidateSrc"));

  //theMuonDsplToken     = consumes< vector<pat::Muon>    >( edm::InputTag("slimmedDisplacedMuons"));
  theDsplMuToken       = consumes< vector<reco::Track>  >( edm::InputTag("displacedGlobalMuons"));
  theDsplTrkToken      = consumes< vector<reco::Track>  >( edm::InputTag("displacedTracks"));
  theVertexToken       = consumes< vector<reco::Vertex> >( edm::InputTag("offlineSlimmedPrimaryVertices"));
  theVertexWithBSToken = consumes< vector<reco::Vertex> >( edm::InputTag("offlineSlimmedPrimaryVerticesWithBS"));

  theVertexCPCToken    = consumes< vector<reco::VertexCompositePtrCandidate> >(edm::InputTag("slimmedSecondaryVertices"));

  theTrackBuilderToken = esConsumes(edm::ESInputTag("", "TransientTrackBuilder"));

}

Analysis::~Analysis()
{
  cout <<" DTOR" << endl;
}

void Analysis::beginJob()
{
  //create a histogram
  histo = new TH1D("histo","test; X; #events",10, 0., 10.);

  hInv  = new TH1D("hInv","Invariant~mass~of~\\mu^+\\mu^-~pairs; M_{inv} [GeV]; Counts", 4800 , 0. , 12. );
  hInvG = new TH1D("hInvG","Invariant~mass~of~\\mu^+\\mu^-~pairs~(global~muons~only); M_{inv} [GeV]; Counts", 4800 , 0. , 12. );
  hID   = new TH1D("hID","Particle ID; ID; Counts", 225 , -9.5, 215.5 );
  hDisp = new TH1D("hDisp","Invariant~mass~of~\\mu^+\\mu^-~pairs~(SlimmedDisplacedMuons); M_{inv} [GeV]; Counts", 4800 , 0. , 12. );
  hCan  = new TH1D("hCan","Invariant~mass~of~\\mu^+\\mu^-~pairs~(Packed Candidate); M_{inv} [GeV]; Counts", 4800 , 0. , 12. );
  hVz   = new TH1D("hVz","The~distance~along~the~Z-axis~between~vertices~for~\\mu^+\\mu^-~pairs; Distance [cm]; Counts", 2400 , 0. , 6. );
  hVzPC = new TH1D("hVzPC","The~distance~along~the~Z-axis~between~vertices~for ~\\mu^+\\mu^-~pairs~(PC); Distance [cm]; Counts", 2400 , 0. , 6. );
  hProb = new TH1D("hProb","The~probability~of~a~common~vertex~for~\\mu^+\\mu^-~pairs; Probability; Counts", 200 , 0. , 1. );
  hJpsiProb = new TH1D("hJpsiProb","The~probability~of~a~common~vertex~for~\\mu^+\\mu^-~pairs,~M_{inv}~=~m(J/\\psi); Probability; Counts", 200 , 0. , 1. );
  hVz01 = new TH1D("hVz01","The~distance~along~the~Z-axis~between~vertices~for~\\mu^+\\mu^-~pairs,~M_{inv}~=~m(J/\\psi); Distance [cm]; Counts", 200000 , 0. , 20. );
  hVz08 = new TH1D("hVz08","The~distance~along~the~Z-axis~between~vertices~for~\\mu^+\\mu^-~pairs,~M_{inv}~=~m(J/\\psi); Distance [cm]; Counts", 200000 , 0. , 20. );
  //hInvK  = new TH1D("hInvK","Invariant~mass~of~K^+K^-~pairs; M_{inv}~[GeV]; Counts", 2400 , 0. , 12. );
  hkkProb = new TH1D("hkkProb","The~probability~of~a~common~vertex~for~K^+K^-~pairs; Probability; Counts", 200 , 0. , 1. );
  hInvMuMuK = new TH1D("hInvMuMuK","Invariant~mass~of~K^+J/\\psi; M_{inv} [GeV]; Counts", 2400 , 0. , 12. );
  hKK = new TH1D("hKK","Invariant~mass~of~\\Phi ; M_{inv} [GeV]; Counts", 2400 , 0.8 , 2. );
  hMuMuKK  = new TH1D("hMuMuKK ","Invariant~mass~of~K^+K^-\\mu^+\\mu^-; M_{inv} [GeV]; Counts", 1200 , 4. , 7. );
  hTest = new TH1D("hTest","Invariant~mass~of~\\pi^+J/\\psi; M_{inv} [GeV]; Counts", 2400 , 0. , 12. );

  cout << "HERE Analysis::beginJob()" << endl;
}

void Analysis::endJob()
{
  //make a new Root file
  TFile myRootFile( theConfig.getParameter<std::string>("outHist").c_str(), "RECREATE");

  //write histogram data
  histo -> Write();
  hInv  -> Write();
  hInvG  -> Write();
  hID   -> Write();
  hDisp -> Write();
  hCan  -> Write();
  hVz   -> Write();
  hVzPC -> Write();
  hProb -> Write();
  hJpsiProb -> Write();
  hVz01 -> Write();
  hVz08 -> Write();
 // hInvK  -> Write();
  hkkProb -> Write();
  hInvMuMuK -> Write();
  hKK -> Write();
  hMuMuKK -> Write();
  hTest -> Write();

  myRootFile.Close();

  delete histo;
  delete hInv;
  delete hInvG;
  delete hID;
  delete hDisp;
  delete hCan;
  delete hVz;
  delete hVzPC;
  delete hProb;
  delete hJpsiProb ;
  delete hVz01;
  delete hVz08;
  //delete hInvK;
  delete hkkProb;
  delete hInvMuMuK;
  delete hKK;
  delete hMuMuKK;
  delete hTest;

  cout << "HERE Cwiczenie::endJob()" << endl;
}


double Analysis::muonEnergy(const pat::Muon& muon ){

  double energy = sqrt( pow(muon.p() , 2.) + mMu * mMu);

  return energy ;
}

double Analysis::invariantMass(const pat::Muon& muon1, const pat::Muon& muon2 ){

  TLorentzVector p4_1(muon1.px(), muon1.py(), muon1.pz(), muonEnergy( muon1 ));
  TLorentzVector p4_2(muon2.px(), muon2.py(), muon2.pz(), muonEnergy( muon2 ));

  TLorentzVector sum = p4_1 + p4_2;

  return sum.M() ;
}

double Analysis::invariantMass(const pat::Muon& muon1, const pat::Muon& muon2, const pat::Muon& muon3, const pat::Muon& muon4 ){

  TLorentzVector p4_1(muon1.px(), muon1.py(), muon1.pz(), muonEnergy( muon1 ));
  TLorentzVector p4_2(muon2.px(), muon2.py(), muon2.pz(), muonEnergy( muon2 ));

  TLorentzVector p4_3(muon3.px(), muon3.py(), muon3.pz(), muonEnergy( muon3 ));
  TLorentzVector p4_4(muon4.px(), muon4.py(), muon4.pz(), muonEnergy( muon4 ));

  TLorentzVector sum = p4_1 + p4_2+ p4_3 + p4_4;

  return sum.M() ;
}

double Analysis::invariantMass(const pat::PackedCandidate& par1, const pat::PackedCandidate& par2, double m ){

  double energy1 = sqrt( pow(par1.p() , 2.) + m * m);
  double energy2 = sqrt( pow(par2.p() , 2.) + m * m);
  
  TLorentzVector p4_1(par1.px(), par1.py(), par1.pz(), energy1);
  TLorentzVector p4_2(par2.px(), par2.py(), par2.pz(), energy2);

  TLorentzVector sum = p4_1 + p4_2;

  return sum.M() ;
}

double Analysis::pTJPsi( const pat::Muon& muon1 , const pat::Muon& muon2 ){

  TLorentzVector p4_1(muon1.px(), muon1.py(), muon1.pz(), muonEnergy( muon1 ));
  TLorentzVector p4_2(muon2.px(), muon2.py(), muon2.pz(), muonEnergy( muon2 ));

  TLorentzVector sum = p4_1 + p4_2;

  return sum.Pt() ;
}

double Analysis::invariantMass3(const pat::Muon& muon2, const pat::Muon& muon1,const pat::PackedCandidate& par1, double& m){

  double energy1 = sqrt( pow(par1.p() , 2.) + m * m);
  //double energy2 = sqrt( pow(par2.p() , 2.) + mK * mK);
  
  TLorentzVector p4_1(muon1.px(), muon1.py(), muon1.pz(), muonEnergy( muon1 ));
  TLorentzVector p4_2(muon2.px(), muon2.py(), muon2.pz(), muonEnergy( muon2 ));
  TLorentzVector p4_3(par1.px(), par1.py(), par1.pz(), energy1);
  //TLorentzVector p4_4(par2.px(), par2.py(), par2.pz(), energy2);

  TLorentzVector sum = p4_1 + p4_2 + p4_3 ;

  return sum.M() ;
}

double Analysis::invariantMass4(const pat::Muon& muon2, const pat::Muon& muon1,const pat::PackedCandidate& par1,const pat::PackedCandidate& par2, double& mK){

  double energy1 = sqrt( pow(par1.p() , 2.) + mK * mK);
  double energy2 = sqrt( pow(par2.p() , 2.) + mK * mK);
  
  TLorentzVector p4_1(muon1.px(), muon1.py(), muon1.pz(), muonEnergy( muon1 ));
  TLorentzVector p4_2(muon2.px(), muon2.py(), muon2.pz(), muonEnergy( muon2 ));
  TLorentzVector p4_3(par1.px(), par1.py(), par1.pz(), energy1);
  TLorentzVector p4_4(par2.px(), par2.py(), par2.pz(), energy2);

  TLorentzVector sum = p4_1 + p4_2 + p4_3 + p4_4;

  return sum.M() ;
}


void Analysis::analyze(
    const edm::Event& ev, const edm::EventSetup& es)
{
  std::cout << " -------------------------------- HERE Cwiczenie::analyze "<< std::endl;
  //bool debug = true;
  
  const vector<pat::Muon> & muons = ev.get(theMuonToken);
  //const vector<pat::Muon> & dispMuons = ev.get(theDisplacedMuonToken);
  const vector<pat::PackedCandidate> & candidates = ev.get(theCandidateToken); 


  // for muon pairs with a common vertex

  std::vector< std::pair< pat::Muon , pat::Muon > > commonPairsJPsi;
  std::vector< std::pair< pat::Muon , pat::Muon > > commonPairsPhi;


  const auto & trackBuilder = es.getData(theTrackBuilderToken); 


  for (std::vector<pat::Muon>::const_iterator im1 = muons.begin(); im1 < muons.end(); im1++) {
    //const pat::Muon & muon = *im1; 

    if( !im1->isGlobalMuon() ) continue; // skip if mu1 is not GLB (PASS IF LooseMuId)

    //if (debug) {}
    //reco::TrackRef mu1Ref = im1->get<reco::TrackRef>();
    reco::TrackRef mu1Ref = im1->track();
    if (!mu1Ref) continue;


    for (std::vector<pat::Muon>::const_iterator im2 = im1+1; im2 < muons.end(); im2++) {
      
      if ( im2-> charge() * im1->charge() >= 0) continue;
      if (!im2->isGlobalMuon() ) continue; // skip if mu1 is not GLB  (PASS IF LooseMuId)
      //reco::TrackRef mu2Ref = im1->get<reco::TrackRef>();
       
      
          
      reco::TrackRef mu2Ref = im2->track();
      if (!mu2Ref)continue;
      //cout << "after getting track of muon 2: " << !mu2Ref << endl;
      
      std::vector<reco::TransientTrack> muonTTs;
      //cout << "before putting muon 1 and muon 2 into muonTTs : " << muonTTs.size() << endl;
      
      muonTTs.push_back(trackBuilder.build(mu1Ref));
      muonTTs.push_back(trackBuilder.build(mu2Ref));
      //cout << "after putting muon 1 and muon 2 into muonTTs : " << muonTTs.size() << endl;


      KalmanVertexFitter kvf(true); //true?
      reco::Vertex vjPsi(TransientVertex(kvf.vertex(muonTTs)));

      double  prob = TMath::Prob(vjPsi.chi2(),vjPsi.ndof());

      hProb -> Fill( prob ) ;

      double mumuInv = invariantMass( *im1 , *im2 );

      double delta_vz = abs( im1 -> vz() - im2 -> vz() );


      if ( prob > 0.2 && delta_vz < 0.3 ){
        hInvG -> Fill( mumuInv ); //only global muons
      }

      if( mumuInv > 3.03 && mumuInv < 3.15 ){ // J/Psi - > a loop over 221
          hJpsiProb -> Fill( prob );

          if( prob > 0.1 && pTJPsi(*im1, *im2) > 8.){


              for( std:: vector<pat::PackedCandidate>::const_iterator ic1 = candidates.begin(); ic1 < candidates.end(); ic1++){
                  
                  if( ic1->pdgId() != 211 || !ic1->hasTrackDetails() || ic1->pt() < 2. ) continue; //K+
                  const reco::Track & rtrk1 = ic1 -> pseudoTrack();
                  if( fabs(vjPsi.position().z() - rtrk1.vz()) < 0.3 ) continue;
                  std::vector<reco::TransientTrack> tracks3 = muonTTs;
                  //cout<< "muonTTs size: " << muonTTs.size() << endl;
                  tracks3.push_back(trackBuilder.build(&rtrk1));
                  //cout<< "muonTTs size: " << muonTTs.size() <<  endl;
                  //KalmanVertexFitter kvf(true); //true?
                  reco::Vertex tv3(TransientVertex(kvf.vertex(tracks3)));
                  //muonTTs.pop_back();
                  //cout<< "muonTTs size: " << muonTTs.size() << endl;
                  double prob = TMath::Prob(tv3.chi2(),tv3.ndof());
                  //double vzMuons = ( im1->vz() + im2->vz() ) / 2;
                  if(prob <0.1 ) continue;
                  hInvMuMuK -> Fill( invariantMass3(*im1, *im2, *ic1, mK) );
                  //cout << invariantMass3(*im1, *im2, *ic1, mK) << endl;
                  hTest -> Fill( invariantMass3(*im1, *im2, *ic1, mPi) );
                  //cout << invariantMass3(*im1, *im2, *ic1, mPi) << endl;
                  for ( std:: vector<pat::PackedCandidate>::const_iterator ic2 = ic1+1; ic2 < candidates.end(); ic2++){
                    if( ic2->pdgId() != 211 || !ic2->hasTrackDetails() || ic2->pt() < 2. ) continue; //K+
                    const reco::Track & rtrk2 = ic2 -> pseudoTrack();
                    if( rtrk1.charge()*rtrk2.charge() != -1 ) continue;
                    if( fabs(vjPsi.position().z() - rtrk2.vz()) < 0.3) continue;
                    std::vector<reco::TransientTrack> tracks4 = tracks3;
                    //cout<< "muonTTs size: " << muonTTs.size() << endl;
                    tracks4.push_back(trackBuilder.build(&rtrk2));
                    //cout<< "muonTTs size: " << muonTTs.size() <<  endl;
                    //KalmanVertexFitter kvf(true); //true?
                    reco::Vertex tv4(TransientVertex(kvf.vertex(tracks4)));
                    //muonTTs.pop_back();
                    //cout<< "muonTTs size: " << muonTTs.size() << endl;
                    double prob = TMath::Prob(tv4.chi2(),tv4.ndof());
                    if(prob < 0.1 ) continue;
                    double massPhi = invariantMass(*ic1, *ic2, mK);
                    hKK -> Fill(massPhi);
                    if( fabs(massPhi - mPhi ) > 0.02) continue;
                    double massBs = invariantMass4(*im1, *im2, *ic1, *ic2, mK);
                    hMuMuKK -> Fill(massBs);

                  }

              }

              hVz01 -> Fill( delta_vz );

              commonPairsJPsi.push_back( std::make_pair( *im1 , *im2) );

          }

          if( prob > 0.8){
              hVz08 -> Fill( delta_vz );

          } 

      }else if( mumuInv > 1. && mumuInv < 1.04 && prob > 0.1 && delta_vz < 0.3 ){ // Phi
  

              commonPairsPhi.push_back( std::make_pair( *im1 , *im2) );

          

      }

    }
  }

/*
for (unsigned int i = 0; i < commonPairsPhi.size(); ++i) {
    for (unsigned int j = i+1; j < commonPairsJPsi.size(); ++j) {
        hExp->Fill(invariantMass(commonPairsPhi[i].first, commonPairsPhi[i].second, commonPairsJPsi[j].first, commonPairsJPsi[j].second));
    }
}
*/



  //cout << "Muon pairs, after muon loop " << commonPairs.size() << endl;

/*
  for (std::vector<pat::PackedCandidate>::const_iterator ic1 = candidates.begin(); ic1 < candidates.end(); ic1++) {
    
    if( abs(ic1->pdgId()) != 211 || !ic1->hasTrackDetails() || ic1->pt() < 1 ) continue;



    for (std::vector<pat::PackedCandidate>::const_iterator ic2 = ic1+1; ic2 < candidates.end(); ic2++) {

      if( abs(ic2->pdgId()) != 211 || (ic1->charge())*(ic2->charge()) >= 0 || !ic2->hasTrackDetails() || ic2->pt() < 1) continue;

      //reco::TrackRef mu1Ref = im1->get<reco::TrackRef>();

      const reco::Track & rtrk1 = ic1 -> pseudoTrack();
      const reco::Track & rtrk2 = ic2 -> pseudoTrack();

      //const reco::Track &  rtrk2 = packedCand.pseudoTrack();
      //TransientTrack PiMTT2 = trackBuilder.build(&rtrk2);
      
      std::vector<reco::TransientTrack> patTTs;
     
      patTTs.push_back(trackBuilder.build(&rtrk1));
      patTTs.push_back(trackBuilder.build(&rtrk2));

      KalmanVertexFitter kvf(true); //true?
      reco::Vertex tv(TransientVertex(kvf.vertex(patTTs)));

      double kkInv = invariantMass( *ic1, *ic2 , mK );

      
      double prob = TMath::Prob(tv.chi2(),tv.ndof());
      double vz = abs(ic2 -> vz() - ic1 -> vz()) ; 
      //cout << "distance between vertices/ prob : " << vz << "/" << prob << endl;

      
      hkkProb -> Fill ( prob );
      //cout << "hkkProb filled " <<  endl;
      //cout << "hMuMuKK, hInvK filled if: " << vz << "<= 0.038... AND " << prob << ">0.1" <<  endl;


      if( vz <= zAv && prob > 0.1 ){
          //cout << " hInvK filled  " << kkInv << endl;
          hInvK -> Fill( kkInv );
          
          for (const auto& pair : commonPairs) {
            const pat::Muon& muon1 = pair.first;
            const pat::Muon& muon2 = pair.second;


            double vzMuons = ( muon1.vz() + muon2.vz() ) / 2;
            double vzKaons = ( ic1 -> vz() + ic2 -> vz() ) / 2;

            if ( abs( vzMuons - vzKaons ) <= zAv ) hMuMuKK -> Fill(invariantMass4( muon1 , muon2 , *ic1, *ic2));
            cout << "hMuMuKK filled: " << invariantMass4(commonV1[i], commonV2[i] , *ic1, *ic2) << endl;
          }
      }
    }
  }  



  cout << "Number of muons in PackedCandidate / slimmed Muons / slimmed Displaced Muons: " 
       << packedMuons.size() << "/" << muons.size() << "/" << dispMuons.size() << endl; 


  unsigned int numberOfParticles = 0;

  vector<pat::PackedCandidate> packedMuons;


  for (const auto & can : candidates) {
    
    numberOfParticles++;   
      
    // find muons in candidates
    if ( std::abs( can.pdgId() ) == 13 ){
         
        packedMuons.push_back( can );
    }

    if ( theEventCount == 1 ){
      hID -> Fill( abs( can.pdgId() ) );
    }

  }

  cout << "Number of particles in thi event: " << numberOfParticles << endl;


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

              hInv -> Fill( invariantMass( muons[i], muons[j] ) );
          }
        }
  }

  // Packed Candidates Muons - histogram
  for (unsigned int i = 0; i < packedMuons.size(); ++i) {

        // Invariant mass of muon pairs + vz difference

        for (unsigned int j = i + 1; j < packedMuons.size(); ++j) {

  
            // invariant mass
          if (packedMuons[i].charge() * packedMuons[j].charge() < 0) {

              // vz

              if ( packedMuons[i].pt() >= 2 && packedMuons[j].pt() >= 2 ) {
                  cout << "VzPC: " << abs( packedMuons[i].vz() - packedMuons[j].vz() ) << endl;
                  hVzPC -> Fill( abs( packedMuons[i].vz() - packedMuons[j].vz() ) );

              }

              hCan -> Fill( invariantMass( packedMuons[i], packedMuons[j], mMu ) );
          }
        }
  }

  // SlimmedDisplacedMuons - histogram
  for (unsigned int i = 0; i < dispMuons.size(); ++i) {


        for (unsigned int j = i + 1; j < dispMuons.size(); ++j) {
            
            if (dispMuons[i].charge() * dispMuons[j].charge() < 0) {
                
                hDisp -> Fill(invariantMass( dispMuons[i] , dispMuons[j] ) );
            }
        }
  }
  
*/
  //write std io
  cout <<"*** Analyze event: " << ev.id()<<" analysed event count:"<<++theEventCount << endl;
}

DEFINE_FWK_MODULE(Analysis);

