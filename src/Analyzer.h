#ifndef Analyzer_h
#define Analyzer_h



// system include files
#include <memory>

// user include files
#include <Math/VectorUtil.h>
#include <fstream>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TTree.h>
#include <string>
#include <vector>
#include <unordered_map>
#include <map>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <TRandom3.h>
#include <TMath.h>
#include <iostream>
#include <iomanip>
#include <utility>
#include <TROOT.h>
#include <TBranch.h>
#include <TApplication.h>
#include <TChain.h>
#include <TDirectory.h>
#include <TLorentzVector.h>
#include <TEnv.h>
#include "Particle.h"
#include "Histo.h"
#include "./svfit/SVfitStandaloneAlgorithm.h"

//#include <ctime>

//#define const
using namespace std;

static const int nTrigReq = 2;

class Analyzer {

 public:
  Analyzer(string, string);
  ~Analyzer();
  void clear_values();
  void preprocess(int);
  void fillCuts();
  void printCuts();
  void writeout();
  int nentries;
  void fill_histogram();

 private:
  void fill_Folder(string, int);

  void getInputs();
  void setupJob(string);
  void initializePileupInfo(string, string);  
  void read_info(string);
  void setupGeneral(TTree*, string);
  void setCutNeeds();

  void smearLepton(Lepton&, CUTS, const PartStats&);
  void smearJet(const PartStats&);

  bool JetMatchesLepton(const Lepton&, const TLorentzVector&, double, CUTS);
  TLorentzVector matchLeptonToGen(const TLorentzVector&, const PartStats&, CUTS);
  TLorentzVector matchTauToGen(const TLorentzVector&, double);

  void getGoodTauNu();  
  void getGoodGen(int, int, CUTS, const PartStats&);
  void getGoodRecoLeptons(const Lepton&, const CUTS, const CUTS, const PartStats&);
  void getGoodRecoJets(CUTS, const PartStats&);

  void getGoodMetTopologyLepton(const Lepton&, CUTS,CUTS, const PartStats&);
  void getGoodLeptonCombos(Lepton&, Lepton&, CUTS,CUTS,CUTS, const PartStats&);
  void getGoodDiJets(const PartStats&);

  void VBFTopologyCut();
  void TriggerCuts(vector<int>&, const vector<string>&, CUTS);

  void SVFit(const Lepton&, const Lepton&, CUTS, svFitStandalone::kDecayType, svFitStandalone::kDecayType, string, int, double);
  pair<svFitStandalone::kDecayType, svFitStandalone::kDecayType> getTypePair(CUTS ePos);

  double calculateLeptonMetMt(const TLorentzVector&);
  double diParticleMass(const TLorentzVector&, const TLorentzVector&, string);
  bool passDiParticleApprox(const TLorentzVector&, const TLorentzVector&, string);
  bool isZdecay(const TLorentzVector&, const Lepton&);

  bool isOverlaping(const TLorentzVector&, Lepton&, CUTS, double);
  bool passProng(string, int);
  bool isInTheCracks(float);
  bool passedLooseJetID(int);

  double getPZeta(const TLorentzVector&, const TLorentzVector&);
  double getPZetaVis(const TLorentzVector&, const TLorentzVector&);
  double normPhi(double);
  double absnormPhi(double);

  void updateMet();
  double getPileupWeight(float);

  TFile* f;
  TTree* BOOM;
  TH1F *hPUmc;
  TH1F *hPUdata;

  Generated* _Gen;
  Electron* _Electron;
  Muon* _Muon;
  Taus* _Tau;
  Jet* _Jet;
  Histogramer histo;

  unordered_map<string, PartStats> distats;

  vector<int>* trigPlace[nTrigReq];
  vector<string>* trigName[nTrigReq];

  PartStats genStat;
  unordered_map<string, double> genMap;
  std::array<std::vector<int>, static_cast<int>(CUTS::enumSize)> goodParts;
  
  vector<int> cuts_per, cuts_cumul;

  TLorentzVector theMETVector;
  double deltaMEx, deltaMEy, sumpxForMht, sumpyForMht, sumptForHt, phiForMht;
  double againstElectron, againstMuon, maxIso, minIso;
  int leadIndex, maxCut;
  bool isData, CalculatePUSystematics;

  vector<double>* Trigger_decision = 0;
  vector<string>* Trigger_names = 0;
  float nTruePU = 0;
  int bestVertices = 0;
  double gen_weight = 0;
  double Met_px = 0;
  double Met_py = 0;
  double Met_pz = 0;
  TMatrixD MetCov;

  double pu_weight, wgt;
  unordered_map<CUTS, bool, EnumHash> need_cut;

  unordered_map<string, CUTS> fill_num = { {"FillVertices", CUTS::eRVertex}, {"FillTauJet1", CUTS::eRTau1}, {"FillTauJet2", CUTS::eRTau2}, {"FillElectron1", CUTS::eRElec1}, {"FillElectron2", CUTS::eRElec2}, {"FillMuon1", CUTS::eRMuon1}, {"FillMuon2", CUTS::eRMuon2}, {"FillJet1", CUTS::eRJet1}, {"FillJet2", CUTS::eRJet2}, {"FillBJet", CUTS::eRBJet}, {"FillCentralJet", CUTS::eRCenJet}, {"FillSusyCuts", CUTS::eSusyCom}, {"FillDiMuon", CUTS::eDiMuon}, {"FillDiTau", CUTS::eDiTau}, {"FillMuon1Tau1", CUTS::eMuon1Tau1}, {"FillMuon1Tau2", CUTS::eMuon1Tau2}, {"FillMuon2Tau1", CUTS::eMuon2Tau1}, {"FillMuon2Tau2", CUTS::eMuon2Tau2}, {"FillElectron1Tau1", CUTS::eElec1Tau1}, {"FillElectron1Tau2", CUTS::eElec1Tau2}, {"FillElectron2Tau1", CUTS::eElec2Tau1}, {"FillElectron2Tau2", CUTS::eElec2Tau2}, {"FillMuon1Electron1", CUTS::eMuon1Elec1}, {"FillMuon1Electron2", CUTS::eMuon1Elec2}, {"FillMuon2Electron1", CUTS::eMuon2Elec1}, {"FillMuon2Electron2", CUTS::eMuon2Elec2} };
  
  std::unordered_map<string, CUTS> cut_num = { {"NGenTau", CUTS::eGTau}, {"NGenTop", CUTS::eGTop}, {"NGenElectron", CUTS::eGElec}, \
    {"NGenMuon", CUTS::eGMuon}, {"NGenZ", CUTS::eGZ}, {"NGenW", CUTS::eGW}, {"NGenHiggs", CUTS::eGHiggs}, \
    {"NRecoVertex", CUTS::eRVertex}, {"NRecoMuon1", CUTS::eRMuon1}, {"NRecoMuon2", CUTS::eRMuon2}, \
    {"NRecoElectron1", CUTS::eRElec1}, {"NRecoElectron2",CUTS::eRElec2}, {"NRecoTau1", CUTS::eRTau1},  \
    {"NRecoTau2", CUTS::eRTau2}, {"NRecoJet1", CUTS::eRJet1}, {"NRecoJet2", CUTS::eRJet2}, \
    {"NRecoCentralJet", CUTS::eRCenJet}, {"NRecoBJet", CUTS::eRBJet}, {"NRecoTriggers1", CUTS::eRTrig1}, 
    {"NRecoTriggers2", CUTS::eRTrig2}, {"NRecoFirstLeadingJet", CUTS::eR1stJet}, {"NRecoSecondLeadingJet", CUTS::eR2ndJet},
    {"NRecoMuon1MetTopology", CUTS::eTMuon1}, {"NRecoMuon2MetTopology", CUTS::eTMuon2}, 
    {"NRecoElectron1MetTopology", CUTS::eTElec1}, {"NRecoElectron2MetTopology", CUTS::eTElec2}, 
    {"NRecoTau1MetTopology", CUTS::eTTau1}, {"NRecoTau2MetTopology", CUTS::eTTau2}, {"NDiMuonCombinations", CUTS::eDiMuon},
    {"NDiElectronCombinations", CUTS::eDiElec}, {"NDiTauCombinations", CUTS::eDiTau}, {"NDiJetCombinations", CUTS::eDiJet},
    {"NMuon1Tau1Combinations", CUTS::eMuon1Tau1}, {"NMuon1Tau2Combinations", CUTS::eMuon1Tau2}, 
    {"NMuon2Tau1Combinations", CUTS::eMuon2Tau1}, {"NMuon2Tau2Combinations", CUTS::eMuon2Tau2},
    {"NElectron1Tau1Combinations", CUTS::eElec1Tau1}, {"NElectron1Tau2Combinations", CUTS::eElec1Tau2},
    {"NElectron2Tau1Combinations", CUTS::eElec2Tau1}, {"NElectron2Tau2Combinations", CUTS::eElec2Tau2},
    {"NMuon1Electron1Combinations", CUTS::eMuon1Elec1}, {"NMuon1Electron2Combinations", CUTS::eMuon1Elec2},
    {"NMuon2Electron1Combinations", CUTS::eMuon2Elec1}, {"NMuon2Electron2Combinations", CUTS::eMuon2Elec2},
    {"NSusyCombinations", CUTS::eSusyCom}, {"METCut", CUTS::eMET} };


};


#endif
