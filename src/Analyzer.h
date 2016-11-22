#ifndef Analyzer_h
#define Analyzer_h

struct CRTester;

// system include files
#include <memory>

// user include files
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>


#include <TDirectory.h>
#include <TLorentzVector.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>

#include "Particle.h"
#include "Histo.h"
#include "Cut_enum.h"
#include "FillInfo.h"

double normPhi(double phi);
double absnormPhi(double phi);

//#define const
using namespace std;

static const int nTrigReq = 2;

class Analyzer {

 public:
  Analyzer(string, string, bool setCR = false);
  ~Analyzer();
  void clear_values();
  void preprocess(int);
  void fillCuts();
  void printCuts();
  void writeout();
  int nentries;
  void fill_histogram();
  void setControlRegions() { histo.setControlRegions();}

  vector<int>* getList(CUTS ePos) {return goodParts[ePos];}
  double getMet() {return theMETVector.Pt();}
  double getHT() {return sumptForHt;}
  double getMHT() {return sqrt((sumpxForMht * sumpxForMht) + (sumpyForMht * sumpyForMht));}


 private:
  
  void CRfillCuts();
  ///// Functions /////
  void fill_Folder(string, const int);

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


  double calculateLeptonMetMt(const TLorentzVector&);
  double diParticleMass(const TLorentzVector&, const TLorentzVector&, string);
  bool passDiParticleApprox(const TLorentzVector&, const TLorentzVector&, string);
  bool isZdecay(const TLorentzVector&, const Lepton&);

  bool isOverlaping(const TLorentzVector&, Lepton&, CUTS, double);
  bool passProng(string, int);
  bool isInTheCracks(float);
  bool passedLooseJetID(int);

  pair<double, double> getPZeta(const TLorentzVector&, const TLorentzVector&);
  void create_fillInfo();

  inline bool passCutRange(string, double, const PartStats&);  

  void updateMet();
  double getPileupWeight(float);
  unordered_map<CUTS, vector<int>*, EnumHash> getArray();

  double getCRVal(string);
  void setupCR(string, double);

  ///// values /////

  TFile* f;
  TTree* BOOM;
  double hPU[100];

  Generated* _Gen;
  Electron* _Electron;
  Muon* _Muon;
  Taus* _Tau;
  Jet* _Jet;
  Histogramer histo;
  PartStats genStat;

  unordered_map<string, PartStats> distats;
  unordered_map<string, FillVals*> fillInfo;
  unordered_map<string, double> genMap;
  unordered_map<CUTS, vector<int>*, EnumHash> goodParts;
  unordered_map<CUTS, bool, EnumHash> need_cut;

  static const unordered_map<string, CUTS> cut_num;
  static const unordered_map<CUTS, vector<CUTS>, EnumHash> adjList;




  vector<int>* trigPlace[nTrigReq];
  bool setTrigger = false;
  vector<string>* trigName[nTrigReq];

  vector<int> cuts_per, cuts_cumul;

  TLorentzVector theMETVector;
  double deltaMEx, deltaMEy, sumpxForMht, sumpyForMht, sumptForHt, phiForMht;

  double maxIso, minIso;
  int leadIndex, maxCut;
  bool isData, CalculatePUSystematics;

  vector<double>* Trigger_decision = 0;
  vector<string>* Trigger_names = 0;
  float nTruePU = 0;
  int bestVertices = 0;
  double gen_weight = 0;
  double Met[3] = {0, 0, 0};


  static vector<CUTS> genCuts;
  double pu_weight, wgt;


  vector<CRTester*> testVec;
};

vector<CUTS> Analyzer::genCuts = {CUTS::eGTau, CUTS::eNuTau, CUTS::eGTop, CUTS::eGElec, CUTS::eGMuon, CUTS::eGZ, CUTS::eGW, CUTS::eGHiggs};   

struct CRTester {
    
  const FillVals* info;
  const string variable;
  const double cutVal;

  CRTester(FillVals* _info, string var, double val) : info(_info), variable(var), cutVal(val) {}
  bool test(Analyzer* analyzer);
};




#endif
