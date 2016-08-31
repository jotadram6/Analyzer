#include "Analyzer.h"
#define ival(x) static_cast<int>(x)
#define BIG_NUM 46340
/*For speed, uncomment these commands and the define statement in the h file.  It will lead to 
  unchecked bounds on values, but will improve speed */
//#define const
//#define at(x) operator[](x)
typedef vector<int>::iterator vec_iter;

//Filespace that has all of the .in files
const string FILESPACE = "PartDet/";
const string PUSPACE = "Pileup/";
//////////PUBLIC FUNCTIONS////////////////////

///Constructor
Analyzer::Analyzer(string infile, string outfile) : hPUmc(new TH1F("hPUmc", "hPUmc", 100, 0, 100)), hPUdata(new TH1F("hPUdata", "hPUdata", 100, 0, 100)), MetCov(2,2) {
  cout << "setup start" << endl;
  f = TFile::Open(infile.c_str());
  f->cd("TNT");
  BOOM = (TTree*)f->Get("TNT/BOOM");
  nentries = (int) BOOM->GetEntries();
  BOOM->SetBranchStatus("*", 0);
  std::cout << "TOTAL EVENTS: " << nentries << std::endl;

  MetCov[0][0] = 0;
  MetCov[0][1] = 0;
  MetCov[1][0] = 0;
  MetCov[1][1] = 0;

  for(int i=0; i < nTrigReq; i++) {
    vector<int>* tmpi = new vector<int>();
    vector<string>* tmps = new vector<string>();
    trigPlace[i] = tmpi;
    trigName[i] = tmps;
  }

  setupGeneral(BOOM,infile);

  isData = distats["Run"].bmap.at("isData");
  CalculatePUSystematics = distats["Run"].bmap.at("CalculatePUSystematics");
  histo = Histogramer(1, FILESPACE+"Hist_entries.in", FILESPACE+"Cuts.in", outfile, isData);
  setCutNeeds();
  for(unordered_map<string,PartStats>::iterator it = distats.begin(); it != distats.end(); it++) {
    if(it->second.bmap["MassBySVFit"]) {
      histo.setupSVFit(it->first, it->second.smap["SVHistname"], it->second.dmap["SVbins"], it->second.dmap["SVmin"], it->second.dmap["SVmax"]);
    }
  }

  initializePileupInfo(distats["Run"].smap.at("MCHistos"), distats["Run"].smap.at("DataHistos"));

  //////need to initialize histo and get values for cut arrays

  cuts_per.resize(histo.get_cuts()->size());
  cuts_cumul.resize(histo.get_cuts()->size());

  if(!isData) {
    _Gen = new Generated(BOOM, FILESPACE + "Gen_info.in");
    genStat = _Gen->pstats["Gen"];
    genMap = genStat.dmap;
  }
  _Electron = new Electron(BOOM, FILESPACE + "Electron_info.in");
  _Muon = new Muon(BOOM, FILESPACE + "Muon_info.in");
  _Tau = new Taus(BOOM, FILESPACE + "Tau_info.in");
  _Jet = new Jet(BOOM, FILESPACE + "Jet_info.in");

  std::cout << "setup complete" << std::endl << endl;

}

////destructor
Analyzer::~Analyzer() {
  delete f;
  delete _Electron;
  delete _Muon;
  delete _Tau;
  delete _Jet;
  if(!isData) delete _Gen;
  
  for(int i=0; i < nTrigReq; i++) {
    delete trigPlace[i];
    delete trigName[i];
  }
}


///resets values so analysis can start
void Analyzer::clear_values() {
  for(int i=0; i < (int)goodParts.size(); i++) {
    goodParts[i].clear();
  }
  deltaMEx=0;
  deltaMEy=0;
  sumpxForMht = 0.0;
  sumpyForMht = 0.0;
  sumptForHt  = 0.0;
  leadIndex=-1;
  maxCut = 0;
}

///Function that does most of the work.  Calculates the number of each particle
void Analyzer::preprocess(int event) {
  BOOM->GetEntry(event);

  //TODO: add in pdf vector(set to 1 for now);
  
  theMETVector.SetPxPyPzE(Met_px, Met_py, Met_pz, sqrt(pow(Met_px,2) + pow(Met_py,2)));
  pu_weight = (!isData && CalculatePUSystematics) ? getPileupWeight(nTruePU) : 1.0;

  // SET NUMBER OF GEN PARTICLES
  // TODOGeneralize to remove magic numbers
  if(!isData){
    getGoodGen(genMap.at("TauID"), genMap.at("TauStatus"), CUTS::eGTau, genStat);
    getGoodGen(genMap.at("TopID"), genMap.at("TopStatus"), CUTS::eGTop, genStat);
    getGoodGen(genMap.at("ElectronID"), genMap.at("ElectronStatus"), CUTS::eGElec, genStat);
    getGoodGen(genMap.at("MuonID"), genMap.at("MuonStatus"), CUTS::eGMuon, genStat);
    getGoodGen(genMap.at("ZID"), genMap.at("ZStatus"), CUTS::eGZ, genStat);
    getGoodGen(genMap.at("WID"), genMap.at("WStatus"), CUTS::eGW, genStat);
    getGoodGen(genMap.at("HiggsID"), genMap.at("HiggsStatus"), CUTS::eGHiggs, genStat);
    getGoodTauNu();
  }

  //////Smearing  
  smearLepton(*_Electron, CUTS::eGElec, _Electron->pstats["Smear"]);
  smearLepton(*_Muon, CUTS::eGMuon, _Muon->pstats["Smear"]);
  smearLepton(*_Tau, CUTS::eGTau, _Tau->pstats["Smear"]);
  smearJet(_Jet->pstats["Smear"]);

  //////Triggers and Vertices
  goodParts[ival(CUTS::eRVertex)].resize(bestVertices);
  TriggerCuts(*(trigPlace[0]), *(trigName[0]), CUTS::eRTrig1);
  TriggerCuts(*(trigPlace[1]), *(trigName[1]), CUTS::eRTrig2);

  // // SET NUMBER OF RECO PARTICLES
  // // MUST BE IN ORDER: Muon/Electron, Tau, Jet
  getGoodRecoLeptons(*_Electron, CUTS::eRElec1, CUTS::eGElec, _Electron->pstats["Elec1"]);
  getGoodRecoLeptons(*_Electron, CUTS::eRElec2, CUTS::eGElec, _Electron->pstats["Elec2"]);
  getGoodRecoLeptons(*_Muon, CUTS::eRMuon1, CUTS::eGMuon, _Muon->pstats["Muon1"]);
  getGoodRecoLeptons(*_Muon, CUTS::eRMuon2, CUTS::eGMuon, _Muon->pstats["Muon2"]);
  getGoodRecoLeptons(*_Tau, CUTS::eRTau1, CUTS::eGTau, _Tau->pstats["Tau1"]);
  getGoodRecoLeptons(*_Tau, CUTS::eRTau2, CUTS::eGTau, _Tau->pstats["Tau2"]);

  getGoodRecoJets(CUTS::eRJet1, _Jet->pstats["Jet1"]);
  getGoodRecoJets(CUTS::eRJet2, _Jet->pstats["Jet2"]);
  getGoodRecoJets(CUTS::eRCenJet, _Jet->pstats["CentralJet"]);
  getGoodRecoJets(CUTS::eRBJet, _Jet->pstats["BJet"]);

  getGoodRecoJets(CUTS::eR1stJet, _Jet->pstats["FirstLeadingJet"]);
  leadIndex = goodParts[ival(CUTS::eR1stJet)].at(0); 
  getGoodRecoJets(CUTS::eR2ndJet, _Jet->pstats["SecondLeadingJet"]);

  ////Updates Met and does MET cut
  updateMet();

  /////  SET NUMBER OF RECO MET TOPOLOGY PARTICLES
  getGoodMetTopologyLepton(*_Electron, CUTS::eRElec1, CUTS::eTElec1, _Electron->pstats["Elec1"]);
  getGoodMetTopologyLepton(*_Electron, CUTS::eRElec2, CUTS::eTElec2, _Electron->pstats["Elec2"]);
  getGoodMetTopologyLepton(*_Muon, CUTS::eRMuon1, CUTS::eTMuon1, _Muon->pstats["Muon1"]);
  getGoodMetTopologyLepton(*_Muon, CUTS::eRMuon2, CUTS::eTMuon2, _Muon->pstats["Muon2"]);
  getGoodMetTopologyLepton(*_Tau, CUTS::eRTau1, CUTS::eTTau1, _Tau->pstats["Tau1"]);
  getGoodMetTopologyLepton(*_Tau, CUTS::eRTau2, CUTS::eTTau2, _Tau->pstats["Tau2"]);

  ///VBF Susy cut on leadin jets
  if(goodParts[ival(CUTS::eR1stJet)].at(0) != -1 && goodParts[ival(CUTS::eR2ndJet)].at(0) != -1) VBFTopologyCut();

  /////lepton lepton topology cuts
  getGoodLeptonCombos(*_Electron, *_Tau, CUTS::eRElec1,CUTS::eRTau1, CUTS::eElec1Tau1, distats["Electron1Tau1"]);
  getGoodLeptonCombos(*_Electron, *_Tau, CUTS::eRElec2, CUTS::eRTau1, CUTS::eElec2Tau1, distats["Electron2Tau1"]);
  getGoodLeptonCombos(*_Electron, *_Tau, CUTS::eRElec1, CUTS::eRTau2, CUTS::eElec1Tau2, distats["Electron1Tau2"]);
  getGoodLeptonCombos(*_Electron, *_Tau, CUTS::eRElec2, CUTS::eRTau2, CUTS::eElec2Tau2, distats["Electron2Tau2"]);

  getGoodLeptonCombos(*_Muon, *_Tau, CUTS::eRMuon1, CUTS::eRTau1, CUTS::eMuon1Tau1, distats["Muon1Tau1"]);
  getGoodLeptonCombos(*_Muon, *_Tau, CUTS::eRMuon1, CUTS::eRTau2, CUTS::eMuon1Tau2, distats["Muon1Tau2"]);
  getGoodLeptonCombos(*_Muon, *_Tau, CUTS::eRMuon2, CUTS::eRTau1, CUTS::eMuon2Tau1, distats["Muon2Tau1"]);
  getGoodLeptonCombos(*_Muon, *_Tau, CUTS::eRMuon2, CUTS::eRTau2, CUTS::eMuon2Tau2, distats["Muon2Tau2"]);

  getGoodLeptonCombos(*_Muon, *_Electron, CUTS::eRMuon1, CUTS::eRElec1, CUTS::eMuon1Elec1, distats["Muon1Electron1"]);
  getGoodLeptonCombos(*_Muon, *_Electron, CUTS::eRMuon1, CUTS::eRElec2, CUTS::eMuon1Elec2, distats["Muon1Electron2"]);
  getGoodLeptonCombos(*_Muon, *_Electron, CUTS::eRMuon2, CUTS::eRElec1, CUTS::eMuon2Elec1, distats["Muon2Electron1"]);
  getGoodLeptonCombos(*_Muon, *_Electron, CUTS::eRMuon2, CUTS::eRElec2, CUTS::eMuon2Elec2, distats["Muon2Electron2"]);

  ////DIlepton topology cuts
  getGoodLeptonCombos(*_Tau, *_Tau, CUTS::eRTau1, CUTS::eRTau2, CUTS::eDiTau, distats["DiTau"]);
  getGoodLeptonCombos(*_Electron, *_Electron, CUTS::eRElec1, CUTS::eRElec2, CUTS::eDiElec, distats["DiElectron"]);
  getGoodLeptonCombos(*_Muon, *_Muon, CUTS::eRMuon1, CUTS::eRMuon2, CUTS::eDiMuon, distats["DiMuon"]);

  ////Dijet cuts
  getGoodDiJets(distats["DiJet"]);

  if(event % 50000 == 0) {
    cout << "Event #" << event << endl;
  }
}


////Reads cuts from Cuts.in file and see if the event has enough particles
void Analyzer::fillCuts() {
  unordered_map<string,pair<int,int> >* cut_info = histo.get_cuts();
  vector<string>* cut_order = histo.get_order();

  string cut;
  int min, max;
  bool prevTrue = true;
  int nparticles, i=0;
  maxCut=0;

  

  for(vector<string>::iterator it=cut_order->begin(); it != cut_order->end(); it++, i++) {
    if(isData && it->find("Gen") != string::npos) continue;
    cut = *it;
    min= cut_info->at(cut).first;
    max= cut_info->at(cut).second;
    nparticles = goodParts[ival(cut_num[cut])].size();
    if( (nparticles >= min) && (nparticles <= max || max == -1)) {
      if((cut_num[cut] == CUTS::eR1stJet || cut_num[cut] == CUTS::eR2ndJet) && goodParts[ival(cut_num[cut])].at(0) == -1 ) {
	prevTrue = false;
	continue;  ////dirty dirty hack
      }
      cuts_per[i]++;
      cuts_cumul[i] += (prevTrue) ? 1 : 0;
      maxCut += (prevTrue) ? 1 : 0;
    } else prevTrue = false;
  }

}


///Prints the number of events that passed each cut per event and cumulatively
//done at the end of the analysis
void Analyzer::printCuts() {
  vector<string>* cut_order = histo.get_order();
  int i =0;

  cout.setf(ios::floatfield,ios::fixed);
  cout<<setprecision(3);
  cout << "\n";
  cout << "Selection Efficiency " << "\n";
  cout << "Total events: " << nentries << "\n";
  cout << "               Name                 Indiv.         Cumulative\n";
  cout << "---------------------------------------------------------------------------\n";
  for(vector<string>::iterator it=cut_order->begin(); it != cut_order->end(); it++, i++) {
    cout << setw(28) << *it << "    ";
    if(isData && it->find("Gen") != string::npos) cout << "Skipped" << endl;
    else cout << setw(5) << cuts_per.at(i) << "  ( " << setw(5) << ((float)cuts_per.at(i)) / nentries << ") "
	      << setw(5) << cuts_cumul.at(i) << "  ( " << setw(5) << ((float)cuts_cumul.at(i)) / nentries << ") " << endl;
  }
  cout << "---------------------------------------------------------------------------\n";  
  histo.fill_histogram();
}

/////////////PRIVATE FUNCTIONS////////////////



///Calculates met from values from each file plus smearing and treating muons as neutrinos
void Analyzer::updateMet() {
  ////// Neutrino update before calculation
  if(distats["Run"].bmap.at("TreatMuonsAsNeutrinos")) {
    for(vec_iter it=goodParts[ival(CUTS::eRMuon1)].begin(); it!=goodParts[ival(CUTS::eRMuon1)].end(); it++) {
      if(find(goodParts[ival(CUTS::eRMuon2)].begin(), goodParts[ival(CUTS::eRMuon2)].end(), (*it)) != goodParts[ival(CUTS::eRMuon2)].end() ) continue;
      deltaMEx += _Muon->smearP.at(*it).Px();
      deltaMEy += _Muon->smearP.at(*it).Py();
    }    
    for(vec_iter it=goodParts[ival(CUTS::eRMuon2)].begin(); it!=goodParts[ival(CUTS::eRMuon2)].end(); it++) {
      deltaMEx += _Muon->smearP.at(*it).Px();
      deltaMEy += _Muon->smearP.at(*it).Py();
    }
  }
  ///---MHT and HT calculations----////
  int i=0;
  for(vector<TLorentzVector>::iterator it=_Jet->smearP.begin(); it!=_Jet->smearP.end(); it++, i++) {
    if( (it->Pt() > distats["Run"].dmap.at("JetPtForMhtAndHt")) && (fabs(it->Eta()) < distats["Run"].dmap.at("JetEtaForMhtAndHt")) ) {
      if(distats["Run"].bmap.at("ApplyJetLooseIDforMhtAndHt") && !passedLooseJetID(i) ) continue;
      
      sumpxForMht -= it->Px();
      sumpyForMht -= it->Py();
      sumptForHt  += it->Pt();
    }
  }
  phiForMht = atan2(sumpyForMht,sumpxForMht);
  //  if(sumpxForMht < 0) phiForMht += (sumpyForMht >= 0) ? TMath::Pi() : -TMath::Pi();

  theMETVector.SetPxPyPzE(theMETVector.Px()+deltaMEx, theMETVector.Py()+deltaMEy, theMETVector.Pz(), 
  			  TMath::Sqrt(pow(theMETVector.Px()+deltaMEx,2) + pow(theMETVector.Py()+deltaMEy,2)));

  /////MET CUTS

  if(distats["Run"].bmap.at("DiscrByMet")) {
    if(theMETVector.Pt() < distats["Run"].pmap.at("RecoMetCut").first) return;
    if(theMETVector.Pt() > distats["Run"].pmap.at("RecoMetCut").second) return;
  }
  if(distats["Run"].bmap.at("DiscrByMHT") && sqrt(pow(sumpxForMht,2) + pow(sumpyForMht,2)) < distats["Run"].dmap.at("MhtCut")) return;

  if(distats["Run"].bmap.at("DiscrByHT") && sumptForHt < distats["Run"].dmap.at("HtCut")) return; 
  
  goodParts[ival(CUTS::eMET)].push_back(1);
}


/////sets up other values needed for analysis that aren't particle specific
void Analyzer::setupGeneral(TTree* BOOM, string infile) {
  BOOM->SetBranchStatus("Trigger_decision", 1);
  BOOM->SetBranchStatus("Trigger_names", 1);
  BOOM->SetBranchStatus("nTruePUInteractions", 1);
  BOOM->SetBranchStatus("bestVertices", 1);
  BOOM->SetBranchStatus("weightevt", 1);
  BOOM->SetBranchStatus("Met_type1PF_px", 1);
  BOOM->SetBranchStatus("Met_type1PF_py", 1);
  BOOM->SetBranchStatus("Met_type1PF_pz", 1);
  BOOM->SetBranchStatus("Met_type1PF_cov00", 1);
  BOOM->SetBranchStatus("Met_type1PF_cov01", 1);
  BOOM->SetBranchStatus("Met_type1PF_cov10", 1);
  BOOM->SetBranchStatus("Met_type1PF_cov11", 1);

  BOOM->SetBranchAddress("Trigger_decision", &Trigger_decision);
  BOOM->SetBranchAddress("Trigger_names", &Trigger_names);
  BOOM->SetBranchAddress("nTruePUInteractions", &nTruePU);
  BOOM->SetBranchAddress("bestVertices", &bestVertices);
  BOOM->SetBranchAddress("weightevt", &gen_weight);
  BOOM->SetBranchAddress("Met_type1PF_px", &Met_px);
  BOOM->SetBranchAddress("Met_type1PF_py", &Met_py);
  BOOM->SetBranchAddress("Met_type1PF_pz", &Met_pz);
  BOOM->SetBranchAddress("Met_type1PF_cov00", &MetCov[0][0]);
  BOOM->SetBranchAddress("Met_type1PF_cov01", &MetCov[0][1]);
  BOOM->SetBranchAddress("Met_type1PF_cov10", &MetCov[1][0]);
  BOOM->SetBranchAddress("Met_type1PF_cov11", &MetCov[1][1]);

  read_info(FILESPACE + "ElectronTau_info.in");
  read_info(FILESPACE + "MuonTau_info.in");
  read_info(FILESPACE + "MuonElectron_info.in");
  read_info(FILESPACE + "DiParticle_info.in");
  read_info(FILESPACE + "VBFCuts_info.in");
  read_info(FILESPACE + "Run_info.in");
}


///parsing method that gets info on diparts and basic run info
//put in map called "distats"
void Analyzer::read_info(string filename) {
  typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
  ifstream info_file(filename);
  boost::char_separator<char> sep(", \t");

  if(!info_file) {
    std::cout << "could not open file " << filename <<std::endl;
    exit(1);
  }

  vector<string> stemp;
  string group, line;
  while(getline(info_file, line)) {
    tokenizer tokens(line, sep);
    stemp.clear();
    for(tokenizer::iterator iter = tokens.begin();iter != tokens.end(); iter++) {
      if( ((*iter)[0] == '/' && (*iter)[0] == '/') || ((*iter)[0] == '#') ) break;
      stemp.push_back(*iter);
    }
    if(stemp.size() == 0) continue;
    else if(stemp.size() == 1) {
      group = stemp[0];
      continue;
    } else if(group == "") {
      cout << "error in " << filename << "; no groups specified for data" << endl;
      exit(1);
    } else if(stemp.size() == 2) {
      if(stemp.at(0).find("Trigger") != string::npos) {
	int ntrig = (stemp.at(0).find("1") != string::npos) ? 0 : 1;
	trigName[ntrig]->push_back(stemp.at(1));
	trigPlace[ntrig]->push_back(0);
	continue;
      }
	
      char* p;
      strtod(stemp[1].c_str(), &p);
      if(stemp[1] == "1" || stemp[1] == "true") distats[group].bmap[stemp[0]]=true;
      else if(stemp[1] == "0" || stemp[1] == "false") distats[group].bmap[stemp[0]]=false; 
      else if(*p) distats[group].smap[stemp[0]] = stemp[1];
      else  distats[group].dmap[stemp[0]]=stod(stemp[1]);

    } else if(stemp.size() == 3) distats[group].pmap[stemp[0]] = make_pair(stod(stemp[1]), stod(stemp[2]));
    else if(stemp.size() == 4) {
      distats[group].smap["SVHistname"] = stemp[0];
      distats[group].dmap["SVbins"] = stod(stemp[1]);
      distats[group].dmap["SVmin"] = stod(stemp[2]);
      distats[group].dmap["SVmax"] = stod(stemp[3]);
    }
  }
  info_file.close();
}

void Analyzer::setCutNeeds() {
  vector<string>* cuts = histo.get_order();
  vector<string>* fillgroups = histo.get_groups();
  for(vector<string>::iterator it = cuts->begin(); it != cuts->end(); ++it) {
    if(cut_num.find(*it) != cut_num.end()) need_cut[cut_num.at(*it)] = true;
  }
  cout << endl;
  for(vector<string>::iterator it = fillgroups->begin(); it != fillgroups->end(); ++it) {
    if(fill_num.find(*it) != fill_num.end()) need_cut[fill_num.at(*it)] = true;
  }
}


///Smears lepton only if specified and not a data file.  Otherwise, just filles up lorentz vectors
//of the data into the vector container smearP with is in each lepton object.
void Analyzer::smearLepton(Lepton& lepton, CUTS eGenPos, const PartStats& stats) {
  lepton.smearP.clear();

  double smearedPt;
  double smearedEta;
  double smearedPhi;
  double smearedEnergy;


  for(int i = 0; i < (int)lepton.pt->size(); i++) {
    TLorentzVector tmpSmear;
    tmpSmear.SetPtEtaPhiE(lepton.pt->at(i), lepton.eta->at(i), lepton.phi->at(i), lepton.energy->at(i));

    if(isData || !stats.bmap.at("SmearTheParticle")) {
      lepton.smearP.push_back(tmpSmear);
      continue;
    }

    TLorentzVector genVec =  matchLeptonToGen(tmpSmear, lepton.pstats["Smear"],eGenPos);
    if(genVec == TLorentzVector(0,0,0,0)) {      
      lepton.smearP.push_back(tmpSmear);
      continue;
    }

    smearedPt = (genVec.Pt()*stats.dmap.at("PtScaleOffset")) + (tmpSmear.Pt() - genVec.Pt())*stats.dmap.at("PtSigmaOffset");
    smearedEta =(genVec.Eta()*stats.dmap.at("EtaScaleOffset")) + (tmpSmear.Eta() - genVec.Eta())*stats.dmap.at("EtaSigmaOffset");
    smearedPhi = (genVec.Phi() * stats.dmap.at("PhiScaleOffset")) + (tmpSmear.Phi() - genVec.Phi())*stats.dmap.at("PhiSigmaOffset");
    smearedEnergy = (genVec.Energy()*stats.dmap.at("EnergyScaleOffset")) + (tmpSmear.Energy() - genVec.Energy())*stats.dmap.at("EnergySigmaOffset");
    
    TLorentzVector final;
    final.SetPtEtaPhiE(smearedPt, smearedEta, smearedPhi, smearedEnergy);
    lepton.smearP.push_back(final);
    deltaMEx += tmpSmear.Px() - final.Px();
    deltaMEy += tmpSmear.Py() - final.Py();
  }
}

///Same as smearlepton, just jet specific
void Analyzer::smearJet(const PartStats& stats) {
  _Jet->smearP.clear();
  TLorentzVector jetV;

  for(int i=0; i< (int)_Jet->pt->size(); i++) {
    jetV.SetPtEtaPhiE(_Jet->pt->at(i), _Jet->eta->at(i), _Jet->phi->at(i), _Jet->energy->at(i));

    if(isData || !stats.bmap.at("SmearTheJet")) {
      _Jet->smearP.push_back(jetV);
      continue;
    }
    
    if(JetMatchesLepton(*_Muon, jetV, stats.dmap.at("MuonMatchingDeltaR"), CUTS::eGMuon) ||
       JetMatchesLepton(*_Tau, jetV, stats.dmap.at("TauMatchingDeltaR"), CUTS::eGTau) ||       
       JetMatchesLepton(*_Electron, jetV,stats.dmap.at("ElectronMatchingDeltaR"), CUTS::eGElec)){

      _Jet->smearP.push_back(jetV);

      continue;
    }

    _Jet->smearP.push_back(stats.dmap.at("JetEnergyScaleOffset") * jetV);
    deltaMEx += (1 - stats.dmap.at("JetEnergyScaleOffset"))*jetV.Px();
    deltaMEy += (1 -stats.dmap.at("JetEnergyScaleOffset"))*jetV.Py();
  }
}

/////checks if jet is close to a lepton and the lepton is a gen particle, then the jet is a lepton object, so
//this jet isn't smeared
bool Analyzer::JetMatchesLepton(const Lepton& lepton, const TLorentzVector& jetV, double partDeltaR, CUTS eGenPos) {
  TLorentzVector tempV;
  for(int j = 0; j < (int)lepton.pt->size(); j++) {
    tempV.SetPtEtaPhiE(lepton.pt->at(j), lepton.eta->at(j), lepton.phi->at(j), lepton.energy->at(j));
    if(jetV.DeltaR(tempV) < partDeltaR && matchLeptonToGen(tempV, lepton.pstats.at("Smear"), eGenPos) != TLorentzVector(0,0,0,0)) return true;
  }
  return false;
}


////checks if reco object matchs a gen object.  If so, then reco object is for sure a correctly identified particle
TLorentzVector Analyzer::matchLeptonToGen(const TLorentzVector& lvec, const PartStats& stats, CUTS ePos) {
  if(ePos == CUTS::eGTau) {
    return matchTauToGen(lvec, stats.dmap.at("GenMatchingDeltaR"));
  }
  TLorentzVector genVec = TLorentzVector(0,0,0,0);
  
  for(vec_iter it=goodParts[ival(ePos)].begin(); it !=goodParts[ival(ePos)].end();it++) {
    genVec.SetPtEtaPhiE(_Gen->pt->at(*it), _Gen->eta->at(*it), _Gen->phi->at(*it), _Gen->energy->at(*it));
    if(lvec.DeltaR(genVec) <= stats.dmap.at("GenMatchingDeltaR")) {
      unordered_map<string,bool>::const_iterator mother = stats.bmap.find("UseMotherID");
      if(mother != stats.bmap.end() && mother->second && abs(_Gen->motherpdg_id->at(*it)) != stats.dmap.at("MotherID")) continue; 
      return genVec;
    }
  }
  
  return TLorentzVector(0,0,0,0);
}


///Tau specific matching fucntion.  Works by seeing if a tau doesn't decay into a muon/electron and has
//a matching tau neutrino showing that the tau decayed and decayed hadronically
TLorentzVector Analyzer::matchTauToGen(const TLorentzVector& lvec, double lDeltaR) {
  TLorentzVector genVec(0,0,0,0);
  int i = 0;
  for(vec_iter it=goodParts[ival(CUTS::eGTau)].begin(); it !=goodParts[ival(CUTS::eGTau)].end();it++, i++) {
    int nu = goodParts[ival(CUTS::eNuTau)].at(i);
    if(nu == -1) continue;

    TLorentzVector tmp1, tmp2;
    tmp1.SetPtEtaPhiE(_Gen->pt->at(*it), _Gen->eta->at(*it), _Gen->phi->at(*it), _Gen->energy->at(*it));
    tmp2.SetPtEtaPhiE(_Gen->pt->at(nu), _Gen->eta->at(nu), _Gen->phi->at(nu), _Gen->energy->at(nu));
    genVec = tmp1 - tmp2;
    if(lvec.DeltaR(genVec) <= lDeltaR) {
      return genVec;
    }
  }
  return TLorentzVector(0,0,0,0);

}


////Calculates the number of gen particles.  Based on id number and status of each particle
void Analyzer::getGoodGen(int particle_id, int particle_status, CUTS ePos, const PartStats& stats) {
  for(int j = 0; j < (int)_Gen->pt->size(); j++) {
    if(particle_id == 15 && (_Gen->pt->at(j) < stats.pmap.at("TauPtCut").first || _Gen->pt->at(j) > stats.pmap.at("TauPtCut").second || abs(_Gen->eta->at(j)) > stats.dmap.at("TauEtaCut"))) continue;
    
    if((abs(_Gen->pdg_id->at(j)) == particle_id) && (_Gen->status->at(j) == particle_status)) {
      goodParts[ival(ePos)].push_back(j);
    }
  }
}

////Tau neutrino specific function used for calculating the number of hadronic taus
void Analyzer::getGoodTauNu() {
  for(vec_iter it=goodParts[ival(CUTS::eGTau)].begin(); it !=goodParts[ival(CUTS::eGTau)].end();it++) {
    bool leptonDecay = false;
    int nu = -1;
    for(int j = 0; j < (int)_Gen->pt->size(); j++) {
      if(abs(_Gen->BmotherIndex->at(j)) == (*it)) {
	if( (abs(_Gen->pdg_id->at(j)) == 16) && (abs(_Gen->motherpdg_id->at(j)) == 15) && (_Gen->status->at(_Gen->BmotherIndex->at(j)) == 2) ) nu = j;
	else if( (abs(_Gen->pdg_id->at(j)) == 12) || (abs(_Gen->pdg_id->at(j)) == 14) ) leptonDecay = true;
      }
    }
    nu = (leptonDecay) ? -1 : nu;
    goodParts[ival(CUTS::eNuTau)].push_back(nu);
  }
}

///Function used to find the number of reco leptons that pass the various cuts.
///Divided into if blocks for the different lepton requirements.
void Analyzer::getGoodRecoLeptons(const Lepton& lep, const CUTS ePos, const CUTS eGenPos, const PartStats& stats) {
  int i = 0;

  for(vector<TLorentzVector>::const_iterator it=lep.smearP.begin(); it != lep.smearP.end(); it++, i++) {
    TLorentzVector lvec = (*it);

    if (fabs(lvec.Eta()) > stats.dmap.at("EtaCut")) continue;
    if (lvec.Pt() < stats.pmap.at("PtCut").first || lvec.Pt() > stats.pmap.at("PtCut").second) continue;

    if((lep.pstats.at("Smear").bmap.at("MatchToGen")) && (!isData)) {   /////check
      if(matchLeptonToGen(lvec, lep.pstats.at("Smear") ,eGenPos) == TLorentzVector(0,0,0,0)) continue;
    }

    if(ePos == CUTS::eRMuon1 || ePos == CUTS::eRMuon2) {      ////////////////MUON CUTS/////////////
      const Muon& partM = static_cast<const Muon&>(lep);
    
      if(stats.bmap.at("DoDiscrByTightID") && (partM.tight->at(i) == 0)) continue;
      if(stats.bmap.at("DoDiscrBySoftID") && (partM.soft->at(i) == 0)) continue;
      
      if (stats.bmap.at("DoDiscrByIsolation")) {
	double maxIsoval = max(0.0, partM.isoNeutralHadron->at(i) + partM.isoPhoton->at(i) - 0.5 * partM.isoPU->at(i));
	double isoSum = (partM.isoCharged->at(i) + maxIsoval) / lvec.Pt();
	if(isoSum < stats.pmap.at("IsoSumPtCutValue").first || isoSum >= stats.pmap.at("IsoSumPtCutValue").second) continue;
      }
    } else if(ePos == CUTS::eRElec1 || ePos == CUTS::eRElec2) {    ///////////////ELECTRON CUT///////////
      const Electron& partE = static_cast<const Electron&>(lep);

      //----Require electron to pass ID discriminators
      if(stats.bmap.at("DoDiscrByVetoID") && (partE.isPassVeto->at(i) == 0)) continue;
      if(stats.bmap.at("DoDiscrByLooseID") && (partE.isPassLoose->at(i) == 0)) continue;
      if(stats.bmap.at("DoDiscrByMediumID") && (partE.isPassMedium->at(i) == 0)) continue;
      if(stats.bmap.at("DoDiscrByTightID") && (partE.isPassTight->at(i) == 0)) continue;
      if(stats.bmap.at("DoDiscrByHEEPID") && (partE.isPassHEEPId->at(i) == 0)) continue;

      if (stats.bmap.at("DoDiscrByIsolation")) {
	double maxIsoval = std::max(0.0 , partE.isoNeutralHadrons->at(i) + partE.isoPhotons->at(i) - 0.5 * partE.isoPU->at(i) );
	double isoSum = (partE.isoChargedHadrons->at(i) + maxIsoval) / lvec.Pt();
	if(isoSum < stats.pmap.at("IsoSumPtCutValue").first || isoSum >= stats.pmap.at("IsoSumPtCutValue").second) continue;
      }

    } else if(ePos == CUTS::eRTau1 || ePos == CUTS::eRTau2) {   /////////////TAU CUT/////////////////
      const Taus& partT = static_cast<const Taus&>(lep);
      if (stats.bmap.at("DoDiscrByLeadTrack")) {
	if(partT.leadChargedCandPt->at(i) < stats.dmap.at("LeadTrackThreshold")) continue;
      }

      // ----Isolation requirement
      if (stats.bmap.at("DoDiscrByIsolation")) {
	//--- max isolation requirement
	maxIso = (ePos == CUTS::eRTau1) ? partT.maxIso.first->at(i) : partT.maxIso.second->at(i);
	if(maxIso < 0.5) continue;
	if(stats.smap.at("DiscrByMinIsolation") != "ZERO") {
	  minIso = (ePos == CUTS::eRTau1) ? partT.minIso.first->at(i) : partT.minIso.second->at(i);
	  if (minIso > 0.5) continue;
	}
      }

      // ----Require 1 or 3 prongs
      if(stats.smap.at("DiscrByProngType").find("hps") != string::npos && partT.decayModeFindingNewDMs->at(i) < 0.5) continue;
      if(!passProng(stats.smap.at("DiscrByProngType"), partT.nProngs->at(i))) continue;

      // ----Electron and Muon vetos
      if (stats.bmap.at("DoDiscrAgainstElectron") && !stats.bmap.at("SelectTausThatAreElectrons")) {
	againstElectron = (ePos == CUTS::eRTau1) ? partT.againstElectron.first->at(i) : partT.againstElectron.second->at(i);
	if (againstElectron < 0.5) continue;
      } else if (stats.bmap.at("SelectTausThatAreElectrons")) {
	againstElectron = (ePos == CUTS::eRTau1) ? partT.againstElectron.first->at(i) : partT.againstElectron.second->at(i);
	if (againstElectron > 0.5) continue;
      }

      if (stats.bmap.at("DoDiscrAgainstMuon") && !stats.bmap.at("SelectTausThatAreMuons")) {
	againstMuon = (ePos == CUTS::eRTau1) ? partT.againstMuon.first->at(i) : partT.againstMuon.second->at(i);
	if (againstMuon < 0.5) continue;
      } else if (stats.bmap.at("SelectTausThatAreMuons")) {
	againstMuon = (ePos == CUTS::eRTau1) ? partT.againstMuon.first->at(i) : partT.againstMuon.second->at(i);
	if (againstMuon > 0.5) continue;
      }

      if (stats.bmap.at("DoDiscrByCrackCut") && isInTheCracks(lvec.Eta())) continue;

      // ----anti-overlap requirements
      if (stats.bmap.at("RemoveOverlapWithMuon1s") && isOverlaping(lvec, *_Muon, CUTS::eRMuon1, stats.dmap.at("Muon1MatchingDeltaR"))) continue;
      if (stats.bmap.at("RemoveOverlapWithMuon2s") && isOverlaping(lvec, *_Muon, CUTS::eRMuon2, stats.dmap.at("Muon2MatchingDeltaR"))) continue;
      if (stats.bmap.at("RemoveOverlapWithElectron1s") && isOverlaping(lvec, *_Electron, CUTS::eRElec1, stats.dmap.at("Electron1MatchingDeltaR"))) continue;
      if (stats.bmap.at("RemoveOverlapWithElectron2s") && isOverlaping(lvec, *_Electron, CUTS::eRElec2, stats.dmap.at("Electron2MatchingDeltaR"))) continue;
    }
    goodParts[ival(ePos)].push_back(i);    
  }

}

////Jet specific function for finding the number of jets that pass the cuts.
//used to find the nubmer of good jet1, jet2, central jet, 1st and 2nd leading jets and bjet.
void Analyzer::getGoodRecoJets(CUTS ePos, const PartStats& stats) {
  int i=0;

  for(vector<TLorentzVector>::iterator it=_Jet->smearP.begin(); it != _Jet->smearP.end(); it++, i++) {
    TLorentzVector lvec = (*it);
    ///if else loop for central jet requirements

    if( ePos == CUTS::eRCenJet) {
      if(fabs(lvec.Eta()) > 2.5) continue;
    } else if (fabs(lvec.Eta()) < stats.pmap.at("EtaCut").first || fabs(lvec.Eta()) > stats.pmap.at("EtaCut").second) continue;
    if (lvec.Pt() < stats.dmap.at("PtCut")) continue;

    /// BJet specific
    if(ePos == CUTS::eRBJet) {
      if(stats.bmap.at("ApplyJetBTagging") && _Jet->bDiscriminator->at(i) <= stats.dmap.at("JetBTaggingCut")) continue;
      if((stats.bmap.at("MatchBToGen")) && !isData && abs(_Jet->partonFlavour->at(i)) != 5) continue;
    } else if (stats.bmap.at("ApplyLooseID") && !passedLooseJetID(i)) continue;

  // ----anti-overlap requirements
    if(stats.bmap.at("RemoveOverlapWithMuon1s") && isOverlaping(lvec, *_Muon, CUTS::eRMuon1, stats.dmap.at("Muon1MatchingDeltaR"))) continue;
    if(stats.bmap.at("RemoveOverlapWithMuon2s") && isOverlaping(lvec, *_Muon, CUTS::eRMuon2, stats.dmap.at("Muon2MatchingDeltaR"))) continue;
    if(stats.bmap.at("RemoveOverlapWithElectron1s") && isOverlaping(lvec, *_Electron, CUTS::eRElec1, stats.dmap.at("Electron1MatchingDeltaR"))) continue;
    if(stats.bmap.at("RemoveOverlapWithElectron2s") && isOverlaping(lvec, *_Electron, CUTS::eRElec2, stats.dmap.at("Electron2MatchingDeltaR"))) continue;
    if(stats.bmap.at("RemoveOverlapWithTau1s") && isOverlaping(lvec, *_Tau, CUTS::eRTau1, stats.dmap.at("Tau1MatchingDeltaR"))) continue;
    if(stats.bmap.at("RemoveOverlapWithTau2s") && isOverlaping(lvec, *_Tau, CUTS::eRTau2, stats.dmap.at("Tau2MatchingDeltaR"))) continue;

    /////fill up array
    goodParts[ival(ePos)].push_back(i);    
  }
  
  if(ePos == CUTS::eR1stJet || ePos == CUTS::eR2ndJet) {
    int potential = -1;
    double prevPt = -1;
    for(vec_iter leadit = goodParts[ival(ePos)].begin(); leadit != goodParts[ival(ePos)].end(); ++leadit) {
      if(((ePos == CUTS::eR2ndJet && (*leadit) != leadIndex) || ePos == CUTS::eR1stJet) && _Jet->smearP.at(*leadit).Pt() > prevPt) {
	potential = (*leadit);
	prevPt = _Jet->smearP.at(*leadit).Pt();
      }
    }
    goodParts[ival(ePos)].clear();
    goodParts[ival(ePos)].push_back(potential);
  }

} 

///function to see if a lepton is overlapping with another particle.  Used to tell if jet or tau
//came ro decayed into those leptons
bool Analyzer::isOverlaping(const TLorentzVector& lvec, Lepton& overlapper, CUTS ePos, double MatchingDeltaR) {
  for(vec_iter it=goodParts[ival(ePos)].begin(); it < goodParts[ival(ePos)].end(); it++) {
    if(lvec.DeltaR(overlapper.smearP.at(*it)) < MatchingDeltaR) return true;
  }
  return false;
}

///Tests if tau decays into the specified number of jet prongs.
bool Analyzer::passProng(string prong, int value) {
  return ( (prong.find("1") != string::npos && value == 1) ||
	   (prong.find("2") != string::npos && value == 2) ||
	   (prong.find("3") != string::npos && value == 3) );
}


////Tests if tau is within the cracks of the detector (the specified eta ranges)
bool Analyzer::isInTheCracks(float etaValue){
  return (fabs(etaValue) < 0.018 ||
          (fabs(etaValue)>0.423 && fabs(etaValue)<0.461) ||
          (fabs(etaValue)>0.770 && fabs(etaValue)<0.806) ||
          (fabs(etaValue)>1.127 && fabs(etaValue)<1.163) ||
          (fabs(etaValue)>1.460 && fabs(etaValue)<1.558));
}
 

//Tests if a jet meets a litany of different tests
bool Analyzer::passedLooseJetID(int nobj) {
  if (_Jet->neutralHadEnergyFraction->at(nobj) >= 0.99) return false;
  if (_Jet->neutralEmEmEnergyFraction->at(nobj) >= 0.99) return false;
  if (_Jet->numberOfConstituents->at(nobj) <= 1) return false;
  if (_Jet->muonEnergyFraction->at(nobj) >= 0.80) return false;
  if ( (fabs(_Jet->smearP.at(nobj).Eta()) < 2.4) && 
           ((_Jet->chargedHadronEnergyFraction->at(nobj) <= 0.0) || 
	    (_Jet->chargedMultiplicity->at(nobj) <= 0.0) || 
	    (_Jet->chargedEmEnergyFraction->at(nobj) >= 0.99) )) return false;
  return true;
}


///sees if the event passed one of the two cuts provided
void Analyzer::TriggerCuts(vector<int>& prevTrig, const vector<string>& trigvec, CUTS ePos) {
  if(! need_cut[ePos]) return;
  for(int i = 0; i < (int)trigvec.size(); i++) {
    if(prevTrig[i] >= (int)Trigger_names->size() || trigvec.at(i) != Trigger_names->at(prevTrig.at(i)) ) {
      for(int j = 0; j < (int)Trigger_names->size(); j++) {
	if(Trigger_names->at(j).find(trigvec.at(i)) != string::npos) {
	  prevTrig.at(i) = j;
	  break;
	}
      }
    }
    if(prevTrig.at(i) < (int)Trigger_names->size() && Trigger_decision->at(prevTrig.at(i)) == 1) {
      goodParts[ival(ePos)].push_back(0);
      return;
    }
  }
}


////VBF specific cuts dealing with the leading jets.
void Analyzer::VBFTopologyCut() {
  if(! need_cut[CUTS::eSusyCom]) return;
  PartStats stats = distats["VBFSUSY"];
  TLorentzVector ljet1 = _Jet->smearP.at(goodParts[ival(CUTS::eR1stJet)].at(0));
  TLorentzVector ljet2 = _Jet->smearP.at(goodParts[ival(CUTS::eR2ndJet)].at(0));
  
  if(stats.bmap.at("DiscrByMass")) {
    if((ljet1 + ljet2).M() < stats.pmap.at("MassCut").first) return;
    if((ljet1 + ljet2).M() > stats.pmap.at("MassCut").second) return;
  }

  if(stats.bmap.at("DiscrByPt")) {
    if((ljet1 + ljet2).Pt() < stats.pmap.at("PtCut").first) return;
    if((ljet1 + ljet2).Pt() > stats.pmap.at("PtCut").second) return;
  }
  
  if(stats.bmap.at("DiscrByDeltaEta")) {
    if(fabs(ljet1.Eta() - ljet2.Eta()) < stats.pmap.at("DeltaEtaCut").first) return;
    if(fabs(ljet1.Eta() - ljet2.Eta()) > stats.pmap.at("DeltaEtaCut").second) return;
  }

  if(stats.bmap.at("DiscrByDeltaPhi")) {
    if(absnormPhi(ljet1.Phi() - ljet2.Phi()) < stats.pmap.at("DeltaPhiCut").first) return;
    if(absnormPhi(ljet1.Phi() - ljet2.Phi()) > stats.pmap.at("DeltaPhiCut").second) return;
  }
  
  if(stats.bmap.at("DiscrByOSEta")) {
    if((ljet1.Eta() * ljet2.Eta()) >= 0) return;
  }

  double dphi1 = normPhi(ljet1.Phi() - theMETVector.Phi());
  double dphi2 = normPhi(ljet2.Phi() - theMETVector.Phi());
  double r1, r2, alpha;
  
  if(stats.bmap.at("DiscrByR1")) {
    r1 = sqrt( pow(dphi1,2.0) + pow((TMath::Pi() - dphi2),2.0) );
    if(r1 < stats.pmap.at("R1Cut").first || r1 > stats.pmap.at("R1Cut").second) return;

  }
  if(stats.bmap.at("DiscrByR2")) {
    r2 = sqrt( pow(dphi2,2.0) + pow((TMath::Pi() - dphi1),2.0) );
    if(r2 < stats.pmap.at("R2Cut").first || r2 > stats.pmap.at("R2Cut").second) return;
  }
  if(stats.bmap.at("DiscrByAlpha")) {
    TLorentzVector addVec = ljet1 + ljet2;
    alpha = (addVec.M() > 0) ? ljet2.Pt() / addVec.M() : -1;
    if(alpha < stats.pmap.at("AlphaCut").first || alpha > stats.pmap.at("AlphaCut").second) return;
  }
  if(stats.bmap.at("DiscrByDphi1")) {
    if(abs(dphi1) < stats.pmap.at("Dphi1Cut").first || abs(dphi1) > stats.pmap.at("Dphi1Cut").second) return;
  }
  if(stats.bmap.at("DiscrByDphi2")) {
    if(abs(dphi2) < stats.pmap.at("Dphi2Cut").first || abs(dphi2) > stats.pmap.at("Dphi2Cut").second) return;
  }

  goodParts[ival(CUTS::eSusyCom)].push_back(0);
}


void Analyzer::getGoodMetTopologyLepton(const Lepton& lep, CUTS eReco, CUTS ePos, const PartStats& stats) {
  if(! need_cut[ePos]) return;
  for(vec_iter it=goodParts[ival(eReco)].begin(); it != goodParts[ival(eReco)].end(); it++) {
    if ((ePos != CUTS::eTTau1 && ePos != CUTS::eTTau2) && stats.bmap.at("DiscrIfIsZdecay")) { 
      if(isZdecay(lep.smearP.at(*it), lep)) continue;
    }

    if (stats.bmap.at("DiscrByMetDphi")) {
      double Dphi = absnormPhi(lep.smearP.at(*it).Phi() - theMETVector.Phi());
      if(Dphi < stats.pmap.at("MetDphiCut").first || Dphi > stats.pmap.at("MetDphiCut").second) continue;
    }

    if (stats.bmap.at("DiscrByMetMt")) {
      double LeptonMetMt = calculateLeptonMetMt(lep.smearP.at(*it));
      if( LeptonMetMt < stats.pmap.at("MetMtCut").first || LeptonMetMt > stats.pmap.at("MetMtCut").second) continue;
    }
    goodParts[ival(ePos)].push_back(*it);
  }
}


//-----Calculate lepton+met transverse mass
double Analyzer::calculateLeptonMetMt(const TLorentzVector& Tobj) {
  double px = Tobj.Px() + theMETVector.Px();
  double py = Tobj.Py() + theMETVector.Py();
  double et = Tobj.Et() + theMETVector.Energy(); //TMath::Sqrt((theMETVector.Px() * theMETVector.Px()) + (theMETVector.Py() * theMETVector.Py()));
  double mt2 = et*et - (px*px + py*py);
  return (mt2 >= 0) ? sqrt(mt2) : -1;
}

/////all it does is add lorentz vectors. 
/////keep in case needed later
// TLorentzVector Analyzer::CalculateTheDiJet4Momentum(TLorentzVector* Tobj1, TLorentzVector* Tobj2) {
//   return (*Tobj1) + (*Tobj2);
// }


/////Calculate the diparticle mass based on how to calculate it
///can use Collinear Approximation, which can fail (number failed available in a histogram)
///can use VectorSumOfVisProductAndMet which is sum of particles and met
///Other which is adding without met
double Analyzer::diParticleMass(const TLorentzVector& Tobj1, const TLorentzVector& Tobj2, string howCalc) {
  bool ratioNotInRange = false;
  TLorentzVector The_LorentzVect;


  //////check this equation/////
  if(howCalc == "CollinearApprox") {
    double x1_numerator = (Tobj1.Px() * Tobj2.Py()) - (Tobj2.Px() * Tobj1.Py());
    double x1_denominator = (Tobj2.Py() * (Tobj1.Px() + theMETVector.Px())) - (Tobj2.Px() * (Tobj1.Py() + theMETVector.Py()));
    double x1 = ( x1_denominator != 0. ) ? x1_numerator/x1_denominator : -1.;
    double x2_numerator = x1_numerator;
    double x2_denominator = (Tobj1.Px() * (Tobj2.Py() + theMETVector.Py())) - (Tobj1.Py() * (Tobj2.Px() + theMETVector.Px()));
    double x2 = ( x2_denominator != 0. ) ? x2_numerator/x2_denominator : -1.;
    ratioNotInRange=!((x1 > 0. && x1 < 1.) && (x2 > 0. && x2 < 1.));
    if (!ratioNotInRange) {
      The_LorentzVect.SetPxPyPzE( (Tobj1.Px() / x1) + (Tobj2.Px() / x2), (Tobj1.Py() / x1) + (Tobj2.Py() / x2), (Tobj1.Pz() / x1) + (Tobj2.Pz() / x2), (Tobj1.Energy() / x1) + (Tobj2.Energy() / x2) );
      return The_LorentzVect.M();
    }
  } 

  if(howCalc == "VectorSumOfVisProductsAndMet" || ratioNotInRange) {
    double px = Tobj1.Px() + Tobj2.Px() + theMETVector.Px();
    double py = Tobj1.Py() + Tobj2.Py() + theMETVector.Py();
    double pz = Tobj1.Pz() + Tobj2.Pz();
    double e  = Tobj1.Energy() + Tobj2.Energy() + theMETVector.Energy();///TMath::Sqrt(pow(theMETVector.Px(),2) + pow(theMETVector.Py(),2));
    The_LorentzVect.SetPxPyPzE(px, py, pz, e);
    return The_LorentzVect.M();
  }

  return (Tobj1 + Tobj2).M();
}

////Tests if the CollinearApproximation works for finding the mass of teh particles
bool Analyzer::passDiParticleApprox(const TLorentzVector& Tobj1, const TLorentzVector& Tobj2, string howCalc) {
  if(howCalc == "CollinearApprox") {
    double x1_numerator = (Tobj1.Px() * Tobj2.Py()) - (Tobj2.Px() * Tobj1.Py());
    double x1_denominator = (Tobj2.Py() * (Tobj1.Px() + theMETVector.Px())) - (Tobj2.Px() * (Tobj1.Py() + theMETVector.Py()));
    double x1 = ( x1_denominator != 0. ) ? x1_numerator/x1_denominator : -1.;
    double x2_numerator = x1_numerator;
    double x2_denominator = (Tobj1.Px() * (Tobj2.Py() + theMETVector.Py())) - (Tobj1.Py() * (Tobj2.Px() + theMETVector.Px()));
    double x2 = ( x2_denominator != 0. ) ? x2_numerator/x2_denominator : -1.;
    return (x1 > 0. && x1 < 1.) && (x2 > 0. && x2 < 1.);
  } else {
    return true;
  }
}




/////abs for values
///Find the number of lepton combos that pass the dilepton cuts
void Analyzer::getGoodLeptonCombos(Lepton& lep1, Lepton& lep2, CUTS ePos1, CUTS ePos2, CUTS ePosFin, const PartStats& stats) {
  if(! need_cut[ePosFin]) return;
  bool sameParticle = (&lep1 == &lep2);
  TLorentzVector part1, part2;

  for(vec_iter i1=goodParts[ival(ePos1)].begin(); i1 != goodParts[ival(ePos1)].end(); i1++) {
    for(vec_iter i2=goodParts[ival(ePos2)].begin(); i2 != goodParts[ival(ePos2)].end(); i2++) {
      if(sameParticle && (*i2) <= (*i1)) continue;
      part1 = lep1.smearP.at(*i1);
      part2 = lep2.smearP.at(*i2);

      if(stats.bmap.at("DiscrByDeltaR") && (part1.DeltaR(part2)) < stats.dmap.at("DeltaRCut")) continue;
   
      if(stats.smap.at("DiscrByOSLSType") == "LS" && (lep1.charge->at(*i1) * lep2.charge->at(*i2) <= 0)) continue;
      else if(stats.smap.at("DiscrByOSLSType") == "OS" && (lep1.charge->at(*i1) * lep2.charge->at(*i2) >= 0)) continue;

      if(stats.bmap.at("DiscrByCosDphi")) {
	double Dphi = absnormPhi( part1.Phi() - part2.Phi());
	if(cos(Dphi) < stats.pmap.at("CosDphiCut").first || cos(Dphi) > stats.pmap.at("CosDphiCut").second) continue;
      }
  // ----Mass window requirement
      
      if (stats.bmap.at("DiscrByMassReco")) {
      	double diMass = diParticleMass(part1,part2, stats.smap.at("HowCalculateMassReco"));
      	if( diMass < stats.pmap.at("MassCut").first || diMass > stats.pmap.at("MassCut").second) continue;
      }

      if (stats.bmap.at("DiscrByCDFzeta2D")) {
      	double CDFzeta = stats.dmap.at("PZetaCutCoefficient") * getPZeta(part1, part2) 
	  + stats.dmap.at("PZetaVisCutCoefficient") * getPZetaVis(part1, part2);
      	if( CDFzeta < stats.pmap.at("CDFzeta2DCutValue").first || CDFzeta > stats.pmap.at("CDFzeta2DCutValue").second ) continue;
      }

      //////////abs on the difference????
      ///////////////////
      if(stats.bmap.find("DeltaPtAndMet") != stats.bmap.end() && stats.bmap.at("DiscrByCosDphi_DeltaPtAndMet")) {
      	double DPhi = absnormPhi(atan2(part1.Py() + part2.Py(), part1.Px() + part2.Px()) - theMETVector.Phi());
      	if( cos(DPhi) < stats.pmap.at("CosDphi_DeltaPtMetCut").first || cos(DPhi) > stats.pmap.at("CosDphi_DeltaPtMetCut").second) continue;
      }
      if (stats.bmap.at("DiscrByDeltaPtDivSumPt")) {
	double ptDiv = (lep1.smearP.at(*i1).Pt() - lep2.smearP.at(*i2).Pt()) / (lep1.smearP.at(*i1).Pt() + lep2.smearP.at(*i2).Pt());
	if( ptDiv < stats.pmap.at("DeltaPtDivSumPtCutValue").first || ptDiv > stats.pmap.at("DeltaPtDivSumPtCutValue").second) continue;
      }

      if (stats.bmap.at("DiscrByDeltaPt")) {
	double deltaPt = lep1.smearP.at(*i1).Pt() - lep2.smearP.at(*i2).Pt(); 
	if(deltaPt < stats.pmap.at("DeltaPtCutValue").first || deltaPt > stats.pmap.at("DeltaPtCutValue").second) continue;
      }
      ///Particlesp that lead to good combo are nGen * part1 + part2
      /// final / nGen = part1 (make sure is integer)
      /// final % nGen = part2 
      goodParts[ival(ePosFin)].push_back((*i1)*BIG_NUM + (*i2));
    }
  }
}


//////////////LOOK INTO DIJET PICKING
///////HOW TO GET RID OF REDUNCENCIES??

/////Same as gooddilepton, just jet specific
void Analyzer::getGoodDiJets(const PartStats& stats) {
  if(! need_cut[CUTS::eDiJet]) return;
  // ----Separation cut between jets (remove overlaps)
  for(vec_iter ij2=goodParts[ival(CUTS::eRJet2)].begin(); ij2 != goodParts[ival(CUTS::eRJet2)].end(); ij2++) {
    for(vec_iter ij1=goodParts[ival(CUTS::eRJet1)].begin(); ij1 != goodParts[ival(CUTS::eRJet1)].end() && (*ij1) < (*ij2); ij1++) {
      if (stats.bmap.at("DiscrByDeltaR")) {
	if(_Jet->smearP.at(*ij1).DeltaR(_Jet->smearP.at(*ij2)) < stats.dmap.at("DeltaRCut")) continue;
      }

      if (stats.bmap.at("DiscrByDeltaEta")) {
	if(fabs(_Jet->smearP.at(*ij1).Eta() - _Jet->smearP.at(*ij2).Eta()) < stats.pmap.at("DeltaEtaCut").first) continue;
	if(fabs(_Jet->smearP.at(*ij1).Eta() - _Jet->smearP.at(*ij2).Eta()) > stats.pmap.at("DeltaEtaCut").second) continue;
      }

      if (stats.bmap.at("DiscrByDeltaPhi")) {
	if(absnormPhi(_Jet->smearP.at(*ij1).Phi() - _Jet->smearP.at(*ij2).Phi()) < stats.pmap.at("DeltaPhiCut").first) continue;
	if(absnormPhi(_Jet->smearP.at(*ij1).Phi() - _Jet->smearP.at(*ij2).Phi()) > stats.pmap.at("DeltaPhiCut").second) continue;
      } 

      if (stats.bmap.at("DiscrByOSEta")) {
	if((_Jet->smearP.at(*ij1).Eta() * _Jet->smearP.at(*ij2).Eta()) >= 0) continue;
      }
      // ----Require both legs to be almost back-to-back in phi
      if (stats.bmap.at("DiscrByCosDphi")) {
	if(cos(absnormPhi(_Jet->smearP.at(*ij1).Phi() - _Jet->smearP.at(*ij2).Phi())) < stats.pmap.at("CosDphiCut").first) continue;
	if(cos(absnormPhi(_Jet->smearP.at(*ij1).Phi() - _Jet->smearP.at(*ij2).Phi())) > stats.pmap.at("CosDphiCut").second) continue;
      }
      // ----Mass window requirement
      if (stats.bmap.at("DiscrByMassReco")) {
	if( ((_Jet->smearP.at(*ij1) + _Jet->smearP.at(*ij2)).M() < stats.pmap.at("MassCut").first) || ((_Jet->smearP.at(*ij1) + _Jet->smearP.at(*ij2)).M() > stats.pmap.at("MassCut").second) ) continue;
      }
      ///Particlesp that lead to good combo are totjet * part1 + part2
      /// final / totjet = part1 (make sure is integer)
      /// final % totjet = part2 
      goodParts[ival(CUTS::eDiJet)].push_back((*ij1)*_Jet->smearP.size() + (*ij2));
    }
  }
}


///////Only tested for if is Zdecay, can include massptasymmpair later?
/////Tests to see if a light lepton came form a zdecay
bool Analyzer::isZdecay(const TLorentzVector& theObject, const Lepton& lep) {
  bool eventIsZdecay = false;
  const float zMass = 90.1876;
  const float zWidth = 2.4952;
  float zmmPtAsymmetry = -10.;

  // if mass is within 3 sigmas of z or pt asymmetry is small set to true.
  for(vector<TLorentzVector>::const_iterator lepit= lep.smearP.begin(); lepit != lep.smearP.end(); lepit++) {
    if(theObject.DeltaR(*lepit) < 0.3) continue;
    if(theObject == (*lepit)) continue;

    TLorentzVector The_LorentzVect = theObject + (*lepit);
    zmmPtAsymmetry = (theObject.Pt() - lepit->Pt()) / (theObject.Pt() + lepit->Pt());

    if( (abs(The_LorentzVect.M() - zMass) < 3.*zWidth) || (fabs(zmmPtAsymmetry) < 0.20) ) {
      eventIsZdecay = true;
      break;
    }
  }

  return eventIsZdecay;
}


///Calculates the Pzeta value
double Analyzer::getPZeta(const TLorentzVector& Tobj1, const TLorentzVector& Tobj2) {
  double zetaX = cos(Tobj1.Phi()) + cos(Tobj2.Phi());
  double zetaY = sin(Tobj1.Phi()) + sin(Tobj2.Phi());
  double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
  if ( zetaR > 0. ) { zetaX /= zetaR; zetaY /= zetaR; }
  double visPx = Tobj1.Px() + Tobj2.Px();
  double visPy = Tobj1.Py() + Tobj2.Py();
  double px = visPx + theMETVector.Px();
  double py = visPy + theMETVector.Py();
  double pZeta = px*zetaX + py*zetaY;
  return pZeta;
}


///calculates the visible pzeta value
double Analyzer::getPZetaVis(const TLorentzVector& Tobj1, const TLorentzVector& Tobj2) {
  double zetaX = cos(Tobj1.Phi()) + cos(Tobj2.Phi());
  double zetaY = sin(Tobj1.Phi()) + sin(Tobj2.Phi());
  double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
  if ( zetaR > 0. ) { zetaX /= zetaR; zetaY /= zetaR; }
  double visPx = Tobj1.Px() + Tobj2.Px();
  double visPy = Tobj1.Py() + Tobj2.Py();
  double pZetaVis = visPx*zetaX + visPy*zetaY;
  return pZetaVis;
}

///Normalizes phi to be between -PI and PI
double Analyzer::normPhi(double phi) {
  static double const TWO_PI = TMath::Pi() * 2;
  while ( phi <= -TMath::Pi() ) phi += TWO_PI;
  while ( phi >  TMath::Pi() ) phi -= TWO_PI;
  return phi;
}

///Takes the absolute value of of normPhi (made because constant use)
double Analyzer::absnormPhi(double phi) {
  return abs(normPhi(phi));
}

void Analyzer::SVFit(const Lepton& lep1, const Lepton& lep2, CUTS ePos, svFitStandalone::kDecayType l1Type, svFitStandalone::kDecayType l2Type, string group, int max, double wgt) {


  TLorentzVector part1, part2;
  string keyname = distats[group.substr(4)].smap["SVHistname"];
  
  for(vec_iter it=goodParts[ival(ePos)].begin(); it != goodParts[ival(ePos)].end(); it++) {
    int p1= (*it) / BIG_NUM;
    int p2= (*it) % BIG_NUM;

    part1 = lep1.smearP.at(p1);
    part2 = lep2.smearP.at(p2);

    vector<svFitStandalone::MeasuredTauLepton> leptons;
    leptons.push_back(svFitStandalone::MeasuredTauLepton(l1Type,part1.Pt(), part1.Eta(), part1.Phi(), part1.M()));
    leptons.push_back(svFitStandalone::MeasuredTauLepton(l2Type,part2.Pt(), part2.Eta(), part2.Phi(), part2.M()));
    SVfitStandaloneAlgorithm algo(leptons, theMETVector.Px(), theMETVector.Py(), MetCov , 0);
    //    algo.maxObjFunctionCalls(5);
    algo.integrateMarkovChain();
    
    if ( algo.isValidSolution() ) {
      histo.addVal(algo.mass(), group,max, keyname, wgt);
      std::cout << "... m svfit : " << algo.mass() << " +/- " << algo.massUncert()  << std::endl; // return value is in units of GeV
    } else {
      std::cout << "... m svfit : ---" << std::endl;
    }

  }
}


////Grabs a list of the groups of histograms to be filled and asked Fill_folder to fill up the histograms
void Analyzer::fill_histogram() {
  if(distats["Run"].bmap["ApplyGenWeight"] && gen_weight == 0.0) return;
  fillCuts();
  vector<string> groups = *histo.get_groups();
  wgt = pu_weight;
  if(distats["Run"].bmap["ApplyGenWeight"]) wgt *= (gen_weight > 0) ? 1.0 : -1.0;

  for(vector<string>::iterator it = groups.begin(); it!=groups.end(); it++) {
    fill_Folder(*it, maxCut);
  }
}

///Function that fills up the histograms
void Analyzer::fill_Folder(string group, int max) {
  if(group == "FillRun") {
    histo.addVal(false, group,histo.get_order()->size(), "Events", wgt);
    histo.addVal(true, group,max, "Events", wgt);
    histo.addVal(bestVertices, group,max, "NVertices", wgt);

  } else if(!isData && group == "FillGen") {

    int nhadtau = 0;
    TLorentzVector genVec;
    int i = 0;
    for(vec_iter it=goodParts[ival(CUTS::eGTau)].begin(); it!=goodParts[ival(CUTS::eGTau)].end(); it++, i++) {

      int nu = goodParts[ival(CUTS::eNuTau)].at(i);
      if(nu != -1) {
	TLorentzVector tmp1, tmp2;
	tmp1.SetPtEtaPhiE(_Gen->pt->at(*it), _Gen->eta->at(*it), _Gen->phi->at(*it), _Gen->energy->at(*it));
	tmp2.SetPtEtaPhiE(_Gen->pt->at(nu), _Gen->eta->at(nu), _Gen->phi->at(nu), _Gen->energy->at(nu));
	genVec = tmp1 - tmp2;
	histo.addVal(genVec.Pt(), group,max, "HadTauPt", wgt);
	histo.addVal(genVec.Eta(), group,max, "HadTauEta", wgt);
	nhadtau++;

      }
      histo.addVal(_Gen->energy->at(*it), group,max, "TauEnergy", wgt);
      histo.addVal(_Gen->pt->at(*it), group,max, "TauPt", wgt);
      histo.addVal(_Gen->eta->at(*it), group,max, "TauEta", wgt);
      histo.addVal(_Gen->phi->at(*it), group,max, "TauPhi", wgt);
      for(vec_iter it2=it+1; it2!=goodParts[ival(CUTS::eGTau)].end(); it2++) {
	TLorentzVector genObjt1;
	genObjt1.SetPtEtaPhiE(_Gen->pt->at(*it), _Gen->eta->at(*it), _Gen->phi->at(*it), _Gen->energy->at(*it));
	TLorentzVector genObjt2;
	genObjt2.SetPtEtaPhiE(_Gen->pt->at(*it2), _Gen->eta->at(*it2), _Gen->phi->at(*it2), _Gen->energy->at(*it2));
	histo.addVal(diParticleMass(genObjt1,genObjt2, "none"), group,max, "DiTauMass", wgt);
      }
    }
    histo.addVal(goodParts[ival(CUTS::eGTau)].size(), group,max, "NTau", wgt);
    histo.addVal(nhadtau, group,max, "NHadTau", wgt);

    for(vec_iter it=goodParts[ival(CUTS::eGMuon)].begin(); it!=goodParts[ival(CUTS::eGMuon)].end(); it++) {
      histo.addVal(_Gen->energy->at(*it), group,max, "MuonEnergy", wgt);
      histo.addVal(_Gen->pt->at(*it), group,max, "MuonPt", wgt);
      histo.addVal(_Gen->eta->at(*it), group,max, "MuonEta", wgt);
      histo.addVal(_Gen->phi->at(*it), group,max, "MuonPhi", wgt);
    }
    histo.addVal(goodParts[ival(CUTS::eGMuon)].size(), group,max, "NMuon", wgt);

          

  } else if(group == "FillTauJet1" || group == "FillTauJet2" || group == "FillMuon1" || group == "FillMuon2" || group == "FillJet1" || group == "FillJet2" || group == "FillBJet" || group == "FillCentralJet" || group == "FillElectron1" || group == "FillElectron2") {
    Particle* part;
    if(group == "FillTauJet1" || group == "FillTauJet2") part=_Tau;
    else if(group == "FillMuon1" || group == "FillMuon2") part=_Muon;
    else if(group == "FillElectron1" || group == "FillElectron2") part=_Electron; 
    else part = _Jet;
    CUTS ePos = fill_num[group];

    for(vec_iter it=goodParts[ival(ePos)].begin(); it!=goodParts[ival(ePos)].end(); it++) {
      histo.addVal(part->smearP.at(*it).Energy(), group,max, "Energy", wgt);
      histo.addVal(part->smearP.at(*it).Pt(), group,max, "Pt", wgt);
      histo.addVal(part->smearP.at(*it).Eta(), group,max, "Eta", wgt);
      histo.addVal(part->smearP.at(*it).Phi(), group,max, "Phi", wgt);
      if(dynamic_cast<Taus*>(part) != NULL) {
	histo.addVal(_Tau->nProngs->at(*it), group,max, "NumSignalTracks", wgt);
  	histo.addVal(_Tau->charge->at(*it), group,max, "Charge", wgt);
	histo.addVal(_Tau->leadChargedCandPt->at(*it), group,max, "SeedTrackPt", wgt);
      } else if(dynamic_cast<Muon*>(part) != NULL) {
	histo.addVal(calculateLeptonMetMt(_Muon->smearP.at(*it)), group,max, "MetMt", wgt);  
      }
        else if(dynamic_cast<Electron*>(part) != NULL) {
          histo.addVal(calculateLeptonMetMt(_Electron->smearP.at(*it)), group,max, "MetMt", wgt); 
      }
    }
    
    if((ePos == CUTS::eRMuon1 || ePos == CUTS::eRMuon2 || ePos == CUTS::eRTau1 || ePos == CUTS::eRTau2 || ePos == CUTS::eRElec1 || ePos == CUTS::eRElec2 ) && goodParts[ival(ePos)].size() > 0) {
      double leadpt = 0;
      double leadeta = 0;
      for(vec_iter it=goodParts[ival(ePos)].begin(); it!=goodParts[ival(ePos)].end(); it++) {
	if(part->smearP.at(*it).Pt() >= leadpt) {
	  leadpt = part->smearP.at(*it).Pt();
	  leadeta = part->smearP.at(*it).Eta();
	}
      }

      histo.addVal(leadpt, group, max, "FirstLeadingPt", wgt);
      histo.addVal(leadeta, group, max, "FirstLeadingEta", wgt);
    }

    histo.addVal(goodParts[ival(ePos)].size(), group,max, "N", wgt);


  } else if(group == "FillSusyCuts") {

    histo.addVal(sqrt((sumpxForMht * sumpxForMht) + (sumpyForMht * sumpyForMht)),group,max, "MHT", wgt);
    histo.addVal(sumptForHt,group,max, "HT", wgt);  
    histo.addVal(sumptForHt + sqrt((sumpxForMht * sumpxForMht) + (sumpyForMht * sumpyForMht)),group,max, "Meff", wgt);
    histo.addVal(theMETVector.Pt(), group,max, "Met", wgt);
    if(goodParts[ival(CUTS::eR1stJet)].at(0) !=-1 && goodParts[ival(CUTS::eR2ndJet)].at(0) != -1) {
      TLorentzVector DiJet = _Jet->smearP.at(goodParts[ival(CUTS::eR1stJet)].at(0)) + _Jet->smearP.at(goodParts[ival(CUTS::eR2ndJet)].at(0));
      histo.addVal(absnormPhi(theMETVector.Phi() - DiJet.Phi()), group,max, "MetDiJetDeltaPhi", wgt);
    }
    
  } else if(group == "FillLeadingJet" && goodParts[ival(CUTS::eSusyCom)].size() == 0) {
    double eta1 = -100, eta2 = -100;
    double pt1 = 0, pt2 = 0;
    if(goodParts[ival(CUTS::eR1stJet)].at(0) != -1) {
      pt1 = _Jet->smearP.at(goodParts[ival(CUTS::eR1stJet)].at(0)).Pt();
      eta1 = _Jet->smearP.at(goodParts[ival(CUTS::eR1stJet)].at(0)).Eta();
    }
    if(goodParts[ival(CUTS::eR2ndJet)].at(0) != -1) {
      pt2 = _Jet->smearP.at(goodParts[ival(CUTS::eR2ndJet)].at(0)).Pt();
      eta2 = _Jet->smearP.at(goodParts[ival(CUTS::eR2ndJet)].at(0)).Eta();
    }
    histo.addVal(pt1,group,max, "FirstPt", wgt);
    histo.addVal(eta1,group,max, "FirstEta", wgt);

    histo.addVal(pt2,group,max, "SecondPt", wgt);
    histo.addVal(eta2,group,max, "SecondEta", wgt);

  } else if(group == "FillLeadingJet" && goodParts[ival(CUTS::eSusyCom)].size() != 0) {
    TLorentzVector first = _Jet->smearP.at(goodParts[ival(CUTS::eR1stJet)].at(0));
    TLorentzVector second = _Jet->smearP.at(goodParts[ival(CUTS::eR2ndJet)].at(0));
    
    histo.addVal(first.Pt(),group,max, "FirstPt", wgt);
    histo.addVal(second.Pt(),group,max, "SecondPt", wgt);

    histo.addVal(first.Eta(),group,max, "FirstEta", wgt);
    histo.addVal(second.Eta(),group,max, "SecondEta", wgt);
    
    TLorentzVector LeadDiJet = first + second;
    
    histo.addVal(LeadDiJet.M(), group,max, "Mass", wgt); 
    histo.addVal(LeadDiJet.Pt(), group,max, "Pt", wgt);  
    histo.addVal(fabs(first.Eta() - second.Eta()), group,max, "DeltaEta", wgt); 
    histo.addVal(first.DeltaR(second),group,max, "DeltaR", wgt);  

    double dphiDijets = absnormPhi(first.Phi() - second.Phi());
    double dphi1 = normPhi(first.Phi() - theMETVector.Phi());
    double dphi2 = normPhi(second.Phi() - theMETVector.Phi());
    double alpha = (LeadDiJet.M() > 0) ? second.Pt() / LeadDiJet.M() : 999999999.0;

    histo.addVal(dphiDijets,group,max, "LeadSublDijetDphi", wgt); 
    histo.addVal(theMETVector.Pt(),dphiDijets, group,max, "MetVsDiJetDeltaPhiLeadSubl", wgt);
    histo.addVal(fabs(first.Eta()-second.Eta()), dphiDijets, group,max, "DeltaEtaVsDeltaPhiLeadSubl", wgt);

    histo.addVal(sqrt( pow(dphi1,2.0) + pow((TMath::Pi() - dphi2),2.0) ), group,max, "R1", wgt);
    histo.addVal(sqrt( pow(dphi2,2.0) + pow((TMath::Pi() - dphi1),2.0)), group,max, "R2", wgt);
    histo.addVal(normPhi(first.Phi() - phiForMht), group,max, "Dphi1MHT", wgt); 
    histo.addVal(normPhi(second.Phi() - phiForMht), group,max, "Dphi2MHT", wgt);
    histo.addVal(dphi1,group,max, "Dphi1", wgt);
    histo.addVal(dphi2,group,max, "Dphi2", wgt);
    histo.addVal(dphi1,dphi2,group,max, "Dphi1VsDphi2", wgt);
    histo.addVal(alpha,group,max, "Alpha", wgt);


    //dijet info
  } else if(group == "FillDiJet") {
    double leaddijetmass = 0;
    double leaddijetpt = 0;
    double leaddijetdeltaR = 0;
    double leaddijetdeltaEta = 0;
    double etaproduct = 0;
    for(vec_iter it=goodParts[ival(CUTS::eDiJet)].begin(); it!=goodParts[ival(CUTS::eDiJet)].end(); it++) {
      int p1 = (*it) / _Jet->smearP.size();
      int p2 = (*it) % _Jet->smearP.size();
      TLorentzVector jet1 = _Jet->smearP.at(p1);
      TLorentzVector jet2 = _Jet->smearP.at(p2);
      TLorentzVector DiJet = jet1 + jet2;
	  
      if(DiJet.M() > leaddijetmass) {
	leaddijetmass = DiJet.M();
	etaproduct = (jet1.Eta() * jet2.Eta() > 0) ? 1 : -1;
      }
      if(DiJet.Pt() > leaddijetpt) leaddijetpt = DiJet.Pt();
      if(fabs(jet1.Eta() - jet2.Eta()) > leaddijetdeltaEta) leaddijetdeltaEta = fabs(jet1.Eta() - jet2.Eta());
      if(jet1.DeltaR(jet2) > leaddijetdeltaR) leaddijetdeltaR = jet1.DeltaR(jet2);

      histo.addVal(DiJet.M(), group,max, "Mass", wgt);
      histo.addVal(DiJet.Pt(), group,max, "Pt", wgt);
      histo.addVal(fabs(jet1.Eta() - jet2.Eta()), group,max, "DeltaEta", wgt);
      histo.addVal(absnormPhi(jet1.Phi() - jet2.Phi()), group,max, "DeltaPhi", wgt);
      histo.addVal(jet1.DeltaR(jet2), group,max, "DeltaR", wgt);
    }

    histo.addVal(leaddijetmass, group,max, "LeadMass", wgt);
    histo.addVal(leaddijetpt, group,max, "LeadPt", wgt);  
    histo.addVal(leaddijetdeltaEta, group,max, "LeadDeltaEta", wgt);
    histo.addVal(leaddijetdeltaR, group,max, "LeadDeltaR", wgt);
    histo.addVal(etaproduct, group,max, "LeadEtaProduct", wgt);


    ////diparticle stuff
  } else if(group == "FillDiMuon" || group == "FillDiTau" || group == "FillMuon1Tau1" || group == "FillMuon1Tau2" || group == "FillMuon2Tau1" || group == "FillMuon2Tau2"  || group == "FillElectron1Tau1" || group == "FillElectron1Tau2" || group == "FillElectron2Tau1" || group == "FillElectron2Tau2" || group == "FillMuon1Electron1" || group == "FillMuon1Electron2" || group == "FillMuon2Electron1" || group == "FillMuon2Electron2") { 
    Lepton* lep1 = NULL;
    Lepton* lep2 = NULL;
    CUTS ePos = fill_num[group];
    string digroup = group;
    digroup.erase(0,4);
    if(ePos == CUTS::eMuon1Tau1 || ePos == CUTS::eMuon1Tau2 || ePos == CUTS::eMuon2Tau1 || ePos == CUTS::eMuon2Tau2) {
      lep1 = _Muon; lep2 = _Tau;
    } else if(ePos == CUTS::eElec1Tau1 || ePos == CUTS::eElec1Tau2 || ePos == CUTS::eElec2Tau1 || ePos == CUTS::eElec2Tau2) {
      lep1 = _Electron; lep2 = _Tau;
    } else if(ePos == CUTS::eMuon1Elec1 || ePos == CUTS::eMuon1Elec2 || ePos == CUTS::eMuon2Elec1 || ePos == CUTS::eMuon2Elec2) {
      lep1 = _Muon; lep2 = _Electron;

    } else if(ePos == CUTS::eDiMuon) {
      lep1 = _Muon; lep2 = _Muon;
    } else if(ePos == CUTS::eDiTau) { lep1 = _Tau; lep2 = _Tau; 
    } else if (ePos == CUTS::eDiElec) { lep1 = _Electron; lep2 = _Electron; }

    if(distats[digroup].bmap["MassBySVFit"] && maxCut == (int)histo.get_cuts()->size()) {
      pair<svFitStandalone::kDecayType, svFitStandalone::kDecayType> typePair = getTypePair(ePos);
      SVFit(*lep1, *lep2, ePos, typePair.first, typePair.second, group, max, wgt);
    }

    TLorentzVector part1;
    TLorentzVector part2;
    
    for(vec_iter it=goodParts[ival(ePos)].begin(); it!=goodParts[ival(ePos)].end(); it++) {
      int p1= (*it) / BIG_NUM;
      int p2= (*it) % BIG_NUM;

      part1 = lep1->smearP.at(p1);
      part2 = lep2->smearP.at(p2);


      histo.addVal(part1.Pt(),part2.Pt(), group,max, "Part1PtVsPart2Pt", wgt);
      histo.addVal(part1.DeltaR(part2), group,max, "DeltaR", wgt); 
      if(group.find("Di") != string::npos) {
	histo.addVal((part1.Pt() - part2.Pt()) / (part1.Pt() + part2.Pt()), group,max, "DeltaPtDivSumPt", wgt);  
	histo.addVal(part1.Pt() - part2.Pt(), group,max, "DeltaPt", wgt);
      } else {
	histo.addVal((part2.Pt() - part1.Pt()) / (part1.Pt() + part2.Pt()), group,max, "DeltaPtDivSumPt", wgt);  
	histo.addVal(part2.Pt() - part1.Pt(), group,max, "DeltaPt", wgt);
      }
      histo.addVal(cos(absnormPhi(part2.Phi() - part1.Phi())), group,max, "CosDphi", wgt);
      histo.addVal(absnormPhi(part1.Phi() - theMETVector.Phi()), group,max, "Part1MetDeltaPhi", wgt);
      histo.addVal(absnormPhi(part1.Phi() - theMETVector.Phi()), cos(absnormPhi(part2.Phi() - part1.Phi())), group,max, "Part1MetDeltaPhiVsCosDphi", wgt);
      histo.addVal(absnormPhi(part2.Phi() - theMETVector.Phi()), group,max, "Part2MetDeltaPhi", wgt);
      histo.addVal(cos(absnormPhi(atan2(part1.Py() - part2.Py(), part1.Px() - part2.Px()) - theMETVector.Phi())), group,max, "CosDphi_DeltaPtAndMet", wgt);

      double diMass = diParticleMass(part1,part2, distats[digroup].smap.at("HowCalculateMassReco"));
      if(passDiParticleApprox(part1,part2, distats[digroup].smap.at("HowCalculateMassReco"))) {
	histo.addVal(diMass, group,max, "ReconstructableMass", wgt);
      } else {
      	histo.addVal(diMass, group,max, "NotReconstructableMass", wgt);
      }

      double PZeta = getPZeta(part1,part2);
      double PZetaVis = getPZetaVis(part1,part2);
      histo.addVal(calculateLeptonMetMt(part1), group,max, "Part1MetMt", wgt);
      histo.addVal(calculateLeptonMetMt(part2), group,max, "Part2MetMt", wgt); 
      histo.addVal(lep2->charge->at(p2) * lep1->charge->at(p1), group,max, "OSLS", wgt);  
      histo.addVal(PZeta, group,max, "PZeta", wgt); 
      histo.addVal(PZetaVis, group,max, "PZetaVis", wgt);  
      histo.addVal(PZetaVis,PZeta, group,max, "Zeta2D", wgt);  
      histo.addVal((distats.at(digroup).dmap.at("PZetaCutCoefficient") * PZeta) + (distats.at(digroup).dmap.at("PZetaVisCutCoefficient") * PZetaVis), group, max, "Zeta1D", wgt);

      if ((goodParts[ival(CUTS::eR1stJet)].at(0) != -1) && (goodParts[ival(CUTS::eR2ndJet)].at(0) != -1)) {
	TLorentzVector TheLeadDiJetVect = _Jet->smearP.at(goodParts[ival(CUTS::eR1stJet)].at(0)) + _Jet->smearP.at(goodParts[ival(CUTS::eR2ndJet)].at(0));

	histo.addVal(absnormPhi(part1.Phi() - TheLeadDiJetVect.Phi()), group,max, "Part1DiJetDeltaPhi", wgt);
	histo.addVal(absnormPhi(part2.Phi() - TheLeadDiJetVect.Phi()), group,max, "Part2DiJetDeltaPhi", wgt);
	histo.addVal(diParticleMass(TheLeadDiJetVect, part1+part2, "VectorSumOfVisProductsAndMet"), group, max, "DiJetReconstructableMass", wgt); 
      }
      
      if(dynamic_cast<Taus*>(lep1) == NULL) {
	histo.addVal(isZdecay(part1, *lep1), group,max, "Part1IsZdecay", wgt); 
      }
      if(dynamic_cast<Taus*>(lep2) == NULL){ 
	histo.addVal(isZdecay(part2, *lep2), group,max, "Part2IsZdecay", wgt); 
      }
    }
  }
}

pair<svFitStandalone::kDecayType, svFitStandalone::kDecayType> Analyzer::getTypePair(CUTS ePos) {
  using namespace svFitStandalone;
  if(ePos == CUTS::eElec1Tau1 || ePos == CUTS::eElec1Tau2 || ePos == CUTS::eElec2Tau1 || ePos == CUTS::eElec2Tau2)
    return make_pair(kTauToElecDecay, kTauToHadDecay);
  else if(ePos == CUTS::eMuon1Tau1 || ePos == CUTS::eMuon1Tau2 || ePos == CUTS::eMuon2Tau1 || ePos == CUTS::eMuon2Tau2)
    return make_pair(kTauToMuDecay, kTauToHadDecay);
  else if(ePos == CUTS::eMuon1Elec1 || ePos == CUTS::eMuon1Elec2 || ePos == CUTS::eMuon2Elec1 || ePos == CUTS::eMuon2Elec2)
    return make_pair(kTauToMuDecay, kTauToElecDecay);
  else if(ePos == CUTS::eDiTau)
    return make_pair(kTauToHadDecay, kTauToHadDecay);
  else if(ePos == CUTS::eDiMuon)
    return make_pair(kTauToMuDecay, kTauToMuDecay);
  else 
    return make_pair(kTauToElecDecay, kTauToElecDecay);
}

void Analyzer::initializePileupInfo(string MCHisto, string DataHisto) {
  // Filenames must be c_strings below. Here is the conversion from strings to c_strings
  // As you can see above cstr1 corresponds to MC and cstr2 corresponds to data.

  TFile *file1 = new TFile((PUSPACE+MCHisto).c_str());
  TH1* histmc = static_cast<TH1*>(file1->Get("analyzeHiMassTau/NVertices_0"));
  if(!histmc) {throw std::runtime_error("failed to extract histogram");}
  for(int bin=0; bin<=(histmc->GetXaxis()->GetNbins() + 1); bin++) {
    hPUmc->SetBinContent(bin,histmc->GetBinContent(bin));
  }
  file1->Close();

  TFile* file2 = new TFile((PUSPACE+DataHisto).c_str());
  TH1* histdata = static_cast<TH1*>(file2->Get("analyzeHiMassTau/NVertices_0"));
  if(!histdata) {throw std::runtime_error("failed to extract histogram");}
  for(int bin=0; bin<=(histdata->GetXaxis()->GetNbins() + 1); bin++) {
    hPUdata->SetBinContent(bin,histdata->GetBinContent(bin));
  }
  file2->Close();
}

double Analyzer::getPileupWeight(float ntruePUInt) {
  int bin;
  double MCintegral;
  double MCvalue;
  double Dataintegral;
  double Datavalue;

  // The probability that data (or MC) has N pileup interactions is value / integral
  // The ratio of the data and MC probability density functions gives us our pileup weights

  //std::cout << "Grabbing pileup info. " << std::endl;
  bin = hPUmc->GetBin(ntruePUInt+1);
  MCvalue = hPUmc->GetBinContent(bin);
  MCintegral = hPUmc->Integral();
  Datavalue = hPUdata->GetBinContent(bin);
  Dataintegral = hPUdata->Integral();

  return ((MCvalue * Dataintegral) != 0) ? (Datavalue * MCintegral) / (MCvalue * Dataintegral) : 1.0;
}





