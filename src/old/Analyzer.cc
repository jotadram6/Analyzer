#include "Analyzer.h"
#define ival(x) static_cast<int>(x)
#define BIG_NUM 46340

#define histAddVal2(val1, val2, name) histo.addVal(val1, val2, group, max, name, wgt)
#define histAddVal(val, name) histo.addVal(val, group, max, name, wgt)
#define SetBranch(name, variable) BOOM->SetBranchStatus(name, 1);  BOOM->SetBranchAddress(name, &variable);

typedef vector<int>::iterator vec_iter;

//Filespace that has all of the .in files
const string FILESPACE = "PartDet/";
const string PUSPACE = "Pileup/";
//////////PUBLIC FUNCTIONS////////////////////

array<vector<int>*, static_cast<int>(CUTS::enumSize)> getArray();

///Constructor
Analyzer::Analyzer(string infile, string outfile) : hPU(new TH1D("hPU", "hPU", 100, 0, 100)), goodParts(getArray()) {
  cout << "setup start" << endl;
  f = TFile::Open(infile.c_str());
  f->cd("TNT");
  BOOM = (TTree*)f->Get("TNT/BOOM");
  nentries = (int) BOOM->GetEntries();
  BOOM->SetBranchStatus("*", 0);
  std::cout << "TOTAL EVENTS: " << nentries << std::endl;

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

unordered_map<CUTS, vector<int>*, EnumHash> getArray() {
  unordered_map<CUTS, vector<int>*, EnumHash> rmap;
  for(auto e: Enum<CUTS>()) {
    rmap[e] = new vector<int>();
  }
  return rmap;
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
    goodParts[i]->clear();
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
  
  theMETVector.SetPxPyPzE(Met[0], Met[1], 0, sqrt(pow(Met[0],2) + pow(Met[1],2)));
  // theMETVector.SetPxPyPzE(Met_px, Met_py, Met_pz, sqrt(pow(Met_px,2) + pow(Met_py,2)));
  pu_weight = (!isData && CalculatePUSystematics) ? hPU->GetBinContent(nTruePU+1) : 1.0;

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
  goodParts[ival(CUTS::eRVertex)]->resize(bestVertices);
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
  leadIndex = goodParts[ival(CUTS::eR1stJet)]->at(0); 
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
  if(goodParts[ival(CUTS::eR1stJet)]->at(0) != -1 && goodParts[ival(CUTS::eR2ndJet)]->at(0) != -1) VBFTopologyCut();

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
    nparticles = goodParts[ival(cut_num[cut])]->size();
    if( (nparticles >= min) && (nparticles <= max || max == -1)) {
      if((cut_num[cut] == CUTS::eR1stJet || cut_num[cut] == CUTS::eR2ndJet) && goodParts[ival(cut_num[cut])]->at(0) == -1 ) {
	prevTrue = false;
	continue;  ////dirty dirty hack
      }
      cuts_per[i]++;
      cuts_cumul[i] += (prevTrue) ? 1 : 0;
      maxCut += (prevTrue) ? 1 : 0;
    } else prevTrue = false;
  }

}

/////// maxcut made -1 if doesn't pass all of the cuts
// void Analyzer::CRfillCuts() {
//   unordered_map<string,pair<int,int> >* cut_info = histo.get_cuts();
//   vector<string>* cut_order = histo.get_order();

//   string cut;
//   int min, max;
//   int nparticles, i=0;
//   maxCut=0;

  

//   for(vector<string>::iterator it=cut_order->begin(); it != cut_order->end(); it++, i++) {
//     if(isData && it->find("Gen") != string::npos) continue;
//     cut = *it;
//     min= cut_info->at(cut).first;
//     max= cut_info->at(cut).second;
//     nparticles = goodParts[ival(cut_num[cut])].size();
//     if( (nparticles > min) || (nparticles > max && max != -1)) {
//       maxCut = -1;
//     } else if((cut_num[cut] == CUTS::eR1stJet || cut_num[cut] == CUTS::eR2ndJet) && goodParts[ival(cut_num[cut])].at(0) == -1 ) maxCut = -1;

//   }
// }

// void Analyzer::CRfolderMax() {
//   string variable1, variable2;

  
// }


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
    for(vec_iter it=goodParts[ival(CUTS::eRMuon1)]->begin(); it!=goodParts[ival(CUTS::eRMuon1)]->end(); it++) {
      if(find(goodParts[ival(CUTS::eRMuon2)]->begin(), goodParts[ival(CUTS::eRMuon2)]->end(), (*it)) != goodParts[ival(CUTS::eRMuon2)]->end() ) continue;
      deltaMEx += _Muon->smearP.at(*it).Px();
      deltaMEy += _Muon->smearP.at(*it).Py();
    }    
    for(vec_iter it=goodParts[ival(CUTS::eRMuon2)]->begin(); it!=goodParts[ival(CUTS::eRMuon2)]->end(); it++) {
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

  theMETVector.SetPxPyPzE(theMETVector.Px()+deltaMEx, theMETVector.Py()+deltaMEy, theMETVector.Pz(), 
  			  TMath::Sqrt(pow(theMETVector.Px()+deltaMEx,2) + pow(theMETVector.Py()+deltaMEy,2)));

  /////MET CUTS

  if(!passCutRange("Met", theMETVector.Pt(), distats["Run"])) return;
  if(distats["Run"].bmap.at("DiscrByHT") && sumptForHt < distats["Run"].dmap.at("HtCut")) return; 
  
  goodParts[ival(CUTS::eMET)]->push_back(1);
}


/////sets up other values needed for analysis that aren't particle specific
void Analyzer::setupGeneral(TTree* BOOM, string infile) {
  SetBranch("Trigger_decision", Trigger_decision);
  SetBranch("Trigger_names", Trigger_names);
  SetBranch("nTruePUInteractions", nTruePU);
  SetBranch("bestVertices", bestVertices);
  SetBranch("weightevt", gen_weight);
  SetBranch("Met_type1PF_px", Met[0]);
  SetBranch("Met_type1PF_py", Met[1]);
  SetBranch("Met_type1PF_pz", Met[2]);

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

    if(cut_num.find(*it) != cut_num.end()) {// cout << *it << endl;
      need_cut[cut_num.at(*it)] = true;}
  }
  cout << endl;
  for(vector<string>::iterator it = fillgroups->begin(); it != fillgroups->end(); ++it) {
    if(fill_num.find(*it) != fill_num.end()) {// cout << *it << endl;
      need_cut[fill_num.at(*it)] = true;}
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
  
  for(vec_iter it=goodParts[ival(ePos)]->begin(); it !=goodParts[ival(ePos)]->end();it++) {
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
  for(vec_iter it=goodParts[ival(CUTS::eGTau)]->begin(); it !=goodParts[ival(CUTS::eGTau)]->end();it++, i++) {
    int nu = goodParts[ival(CUTS::eNuTau)]->at(i);
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
      goodParts[ival(ePos)]->push_back(j);
    }
  }
}

////Tau neutrino specific function used for calculating the number of hadronic taus
void Analyzer::getGoodTauNu() {
  for(vec_iter it=goodParts[ival(CUTS::eGTau)]->begin(); it !=goodParts[ival(CUTS::eGTau)]->end();it++) {
    bool leptonDecay = false;
    int nu = -1;
    for(int j = 0; j < (int)_Gen->pt->size(); j++) {
      if(abs(_Gen->BmotherIndex->at(j)) == (*it)) {
	if( (abs(_Gen->pdg_id->at(j)) == 16) && (abs(_Gen->motherpdg_id->at(j)) == 15) && (_Gen->status->at(_Gen->BmotherIndex->at(j)) == 2) ) nu = j;
	else if( (abs(_Gen->pdg_id->at(j)) == 12) || (abs(_Gen->pdg_id->at(j)) == 14) ) leptonDecay = true;
      }
    }
    nu = (leptonDecay) ? -1 : nu;
    goodParts[ival(CUTS::eNuTau)]->push_back(nu);
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

    if (stats.bmap.at("DoDiscrByIsolation")) {
      double firstIso = (stats.pmap.find("IsoSumPtCutValue") != stats.pmap.end()) ? stats.pmap.at("IsoSumPtCutValue").first : (ival(ePos) - ival(CUTS::eRTau1) + 1);
      double secondIso = (stats.pmap.find("IsoSumPtCutValue") != stats.pmap.end()) ? stats.pmap.at("IsoSumPtCutValue").second : 0;
      if(!lep.get_Iso(i, firstIso, secondIso)) continue;
    }

    if(lep.type == PType::Muon) {      ////////////////MUON CUTS/////////////
      if(stats.bmap.at("DoDiscrByTightID") && (_Muon->tight->at(i) == 0)) continue;
      if(stats.bmap.at("DoDiscrBySoftID") && (_Muon->soft->at(i) == 0)) continue;
      
    } else if(lep.type == PType::Electron) {    ///////////////ELECTRON CUT///////////

      //----Require electron to pass ID discriminators
      if(stats.bmap.at("DoDiscrByVetoID") && (_Electron->isPassVeto->at(i) == 0)) continue;
      if(stats.bmap.at("DoDiscrByLooseID") && (_Electron->isPassLoose->at(i) == 0)) continue;
      if(stats.bmap.at("DoDiscrByMediumID") && (_Electron->isPassMedium->at(i) == 0)) continue;
      if(stats.bmap.at("DoDiscrByTightID") && (_Electron->isPassTight->at(i) == 0)) continue;
      if(stats.bmap.at("DoDiscrByHEEPID") && (_Electron->isPassHEEPId->at(i) == 0)) continue;

     
    } else if(lep.type == PType::Tau) {   /////////////TAU CUT/////////////////
      if (stats.bmap.at("DoDiscrByLeadTrack")) {
	if(_Tau->leadChargedCandPt->at(i) < stats.dmap.at("LeadTrackThreshold")) continue;
      }
     
      // ----Require 1 or 3 prongs
      if(stats.smap.at("DiscrByProngType").find("hps") != string::npos && _Tau->decayModeFindingNewDMs->at(i) < 0.5) continue;
      if(!passProng(stats.smap.at("DiscrByProngType"), _Tau->nProngs->at(i))) continue;

      // ----Electron and Muon vetos
      double against = (ePos == CUTS::eRTau1) ? _Tau->againstElectron.first->at(i) : _Tau->againstElectron.second->at(i);
      if (stats.bmap.at("DoDiscrAgainstElectron") && against < 0.5) continue;
      else if (stats.bmap.at("SelectTausThatAreElectrons") && against > 0.5) continue;
    
      against = (ePos == CUTS::eRTau1) ? _Tau->againstMuon.first->at(i) : _Tau->againstMuon.second->at(i);
      if (stats.bmap.at("DoDiscrAgainstMuon") && against < 0.5) continue;
      else if (stats.bmap.at("SelectTausThatAreMuons") && against > 0.5) continue;

      if (stats.bmap.at("DoDiscrByCrackCut") && isInTheCracks(lvec.Eta())) continue;

      // ----anti-overlap requirements
      if (stats.bmap.at("RemoveOverlapWithMuon1s") && isOverlaping(lvec, *_Muon, CUTS::eRMuon1, stats.dmap.at("Muon1MatchingDeltaR"))) continue;
      if (stats.bmap.at("RemoveOverlapWithMuon2s") && isOverlaping(lvec, *_Muon, CUTS::eRMuon2, stats.dmap.at("Muon2MatchingDeltaR"))) continue;
      if (stats.bmap.at("RemoveOverlapWithElectron1s") && isOverlaping(lvec, *_Electron, CUTS::eRElec1, stats.dmap.at("Electron1MatchingDeltaR"))) continue;
      if (stats.bmap.at("RemoveOverlapWithElectron2s") && isOverlaping(lvec, *_Electron, CUTS::eRElec2, stats.dmap.at("Electron2MatchingDeltaR"))) continue;
    }
    goodParts[ival(ePos)]->push_back(i);    
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
    goodParts[ival(ePos)]->push_back(i);    
  }
  
  if(ePos == CUTS::eR1stJet || ePos == CUTS::eR2ndJet) {
    int potential = -1;
    double prevPt = -1;
    for(vec_iter leadit = goodParts[ival(ePos)]->begin(); leadit != goodParts[ival(ePos)]->end(); ++leadit) {
      if(((ePos == CUTS::eR2ndJet && (*leadit) != leadIndex) || ePos == CUTS::eR1stJet) && _Jet->smearP.at(*leadit).Pt() > prevPt) {
	potential = (*leadit);
	prevPt = _Jet->smearP.at(*leadit).Pt();
      }
    }
    goodParts[ival(ePos)]->clear();
    goodParts[ival(ePos)]->push_back(potential);
  }

} 

///function to see if a lepton is overlapping with another particle.  Used to tell if jet or tau
//came ro decayed into those leptons
bool Analyzer::isOverlaping(const TLorentzVector& lvec, Lepton& overlapper, CUTS ePos, double MatchingDeltaR) {
  for(vec_iter it=goodParts[ival(ePos)]->begin(); it < goodParts[ival(ePos)]->end(); it++) {
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
      goodParts[ival(ePos)]->push_back(0);
      return;
    }
  }
}


////VBF specific cuts dealing with the leading jets.
void Analyzer::VBFTopologyCut() {
  if(! need_cut[CUTS::eSusyCom]) return;
  PartStats stats = distats["VBFSUSY"];
  TLorentzVector ljet1 = _Jet->smearP.at(goodParts[ival(CUTS::eR1stJet)]->at(0));
  TLorentzVector ljet2 = _Jet->smearP.at(goodParts[ival(CUTS::eR2ndJet)]->at(0));
  
  if(!passCutRange("Mass", (ljet1 + ljet2).M(), stats)) return;
  if(!passCutRange("Pt", (ljet1 + ljet2).Pt(), stats)) return;
  if(!passCutRange("DeltaEta", abs(ljet1.Eta() - ljet2.Eta()), stats)) return;
  if(!passCutRange("DeltaEta", absnormPhi(ljet1.Phi() - ljet2.Phi()), stats)) return;
  
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
  if( !passCutRange("Dphi1", abs(dphi1), stats)) return;
  if( !passCutRange("Dphi2", abs(dphi2), stats)) return;

  goodParts[ival(CUTS::eSusyCom)]->push_back(0);
}

inline bool Analyzer::passCutRange(string CutName, double value, const PartStats& stats) {
  return ( !(stats.bmap.at("DiscrBy" + CutName)) || (value > stats.pmap.at(CutName + "Cut").first && value < stats.pmap.at(CutName + "Cut").second) );

}

void Analyzer::getGoodMetTopologyLepton(const Lepton& lep, CUTS eReco, CUTS ePos, const PartStats& stats) {
  if(! need_cut[ePos]) return;
  for(vec_iter it=goodParts[ival(eReco)]->begin(); it != goodParts[ival(eReco)]->end(); it++) {
    if ((ePos != CUTS::eTTau1 && ePos != CUTS::eTTau2) && stats.bmap.at("DiscrIfIsZdecay")) { 
      if(isZdecay(lep.smearP.at(*it), lep)) continue;
    }
    if(!passCutRange("MetDphi", absnormPhi(lep.smearP.at(*it).Phi() - theMETVector.Phi()), stats)) continue;
    if(!passCutRange("MetMt", calculateLeptonMetMt(lep.smearP.at(*it)), stats)) continue;

    goodParts[ival(ePos)]->push_back(*it);
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
    double denominator = (Tobj1.Px() * Tobj2.Py()) - (Tobj2.Px() * Tobj1.Py());
    double x1 = (Tobj2.Py()*theMETVector.Px() - Tobj2.Px()*theMETVector.Py())/denominator;
    double x2 = (Tobj1.Px()*theMETVector.Py() - Tobj1.Py()*theMETVector.Px())/denominator;
    ratioNotInRange=!((x1 < 0.) && (x2 < 0.));
    if (!ratioNotInRange) {
      The_LorentzVect.SetPxPyPzE( (Tobj1.Px()*(1 + x1) + Tobj2.Px()*(1+x2)), (Tobj1.Py()*(1+x1) + Tobj2.Py()*(1+x2)), (Tobj1.Pz()*(1+x1) + Tobj2.Pz()*(1+x2)), (Tobj1.Energy()*(1+x1) + Tobj2.Energy()*(1+x2)) );
      return The_LorentzVect.M();
    }
  } 

  if(howCalc == "VectorSumOfVisProductsAndMet" || ratioNotInRange) {
    return (Tobj1 + Tobj2 + theMETVector).M();
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

  for(vec_iter i1=goodParts[ival(ePos1)]->begin(); i1 != goodParts[ival(ePos1)]->end(); i1++) {
    for(vec_iter i2=goodParts[ival(ePos2)]->begin(); i2 != goodParts[ival(ePos2)]->end(); i2++) {
      if(sameParticle && (*i2) <= (*i1)) continue;
      part1 = lep1.smearP.at(*i1);
      part2 = lep2.smearP.at(*i2);

      if(stats.bmap.at("DiscrByDeltaR") && (part1.DeltaR(part2)) < stats.dmap.at("DeltaRCut")) continue;
   
      if(stats.smap.at("DiscrByOSLSType") == "LS" && (lep1.charge->at(*i1) * lep2.charge->at(*i2) <= 0)) continue;
      else if(stats.smap.at("DiscrByOSLSType") == "OS" && (lep1.charge->at(*i1) * lep2.charge->at(*i2) >= 0)) continue;

      if( !passCutRange("CosDphi", absnormPhi( part1.Phi() - part2.Phi()), stats)) continue;

  // ----Mass window requirement
      
      if (stats.bmap.at("DiscrByMassReco")) {
      	double diMass = diParticleMass(part1,part2, stats.smap.at("HowCalculateMassReco"));
      	if( diMass < stats.pmap.at("MassCut").first || diMass > stats.pmap.at("MassCut").second) continue;
      }

      if (stats.bmap.at("DiscrByCDFzeta2D")) {
      	double CDFzeta = stats.dmap.at("PZetaCutCoefficient") * getPZeta(part1, part2).first 
	  + stats.dmap.at("PZetaVisCutCoefficient") * getPZeta(part1, part2).second;
      	if( CDFzeta < stats.pmap.at("CDFzeta2DCutValue").first || CDFzeta > stats.pmap.at("CDFzeta2DCutValue").second ) continue;
      }

      ///      optional cut
      // if(stats.bmap.find("DeltaPtAndMet") != stats.bmap.end() && stats.bmap.at("DiscrByCosDphi_DeltaPtAndMet")) {
      // 	double DPhi = absnormPhi(atan2(part1.Py() + part2.Py(), part1.Px() + part2.Px()) - theMETVector.Phi());
      // 	if( cos(DPhi) < stats.pmap.at("CosDphi_DeltaPtMetCut").first || cos(DPhi) > stats.pmap.at("CosDphi_DeltaPtMetCut").second) continue;
      // }

      //////////abs on the difference????
      ///////////////////

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
      goodParts[ival(ePosFin)]->push_back((*i1)*BIG_NUM + (*i2));
    }
  }
}


//////////////LOOK INTO DIJET PICKING
///////HOW TO GET RID OF REDUNCENCIES??

/////Same as gooddilepton, just jet specific
void Analyzer::getGoodDiJets(const PartStats& stats) {
  if(! need_cut[CUTS::eDiJet]) return;
   TLorentzVector jet1, jet2;
  // ----Separation cut between jets (remove overlaps)
  for(vec_iter ij2=goodParts[ival(CUTS::eRJet2)]->begin(); ij2 != goodParts[ival(CUTS::eRJet2)]->end(); ij2++) {
    jet2 = _Jet->smearP.at(*ij2);
    for(vec_iter ij1=goodParts[ival(CUTS::eRJet1)]->begin(); ij1 != goodParts[ival(CUTS::eRJet1)]->end() && (*ij1) < (*ij2); ij1++) {
      jet1 = _Jet->smearP.at(*ij1);

      if (stats.bmap.at("DiscrByDeltaR")) {
	if(jet1.DeltaR(jet2) < stats.dmap.at("DeltaRCut")) continue;
      }

      if( !passCutRange("DeltaEta", abs(jet1.Eta() - jet2.Eta()), stats) ) continue;
      if( !passCutRange("DeltaPhi", abs(jet1.Phi() - jet2.Phi()), stats) ) continue;

      if (stats.bmap.at("DiscrByOSEta")) {
	if((jet1.Eta() * jet2.Eta()) >= 0) continue;
      }
      // ----Require both legs to be almost back-to-back in phi
      if( !passCutRange("CosDphi", cos(absnormPhi(jet1.Phi() - jet2.Phi())), stats) ) continue;

      // ----Mass window requirement
      if (stats.bmap.at("DiscrByMassReco")) {
	if( ((jet1 + jet2).M() < stats.pmap.at("MassCut").first) || ((jet1 + jet2).M() > stats.pmap.at("MassCut").second) ) continue;
      }
      ///Particlesp that lead to good combo are totjet * part1 + part2
      /// final / totjet = part1 (make sure is integer)
      /// final % totjet = part2 
      goodParts[ival(CUTS::eDiJet)]->push_back((*ij1)*_Jet->smearP.size() + (*ij2));
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
pair<double, double> Analyzer::getPZeta(const TLorentzVector& Tobj1, const TLorentzVector& Tobj2) {
  double zetaX = cos(Tobj1.Phi()) + cos(Tobj2.Phi());
  double zetaY = sin(Tobj1.Phi()) + sin(Tobj2.Phi());
  double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
  if ( zetaR > 0. ) { zetaX /= zetaR; zetaY /= zetaR; }
  double visPx = Tobj1.Px() + Tobj2.Px();
  double visPy = Tobj1.Py() + Tobj2.Py();
  double px = visPx + theMETVector.Px();
  double py = visPy + theMETVector.Py();
  return make_pair(px*zetaX + py*zetaY, visPx*zetaX + visPy*zetaY);
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

    histo.addVal(false, group,histo.get_groups()->size(), "Events", 1);
    histAddVal(true, "Events");
    histAddVal(bestVertices, "NVertices");
  } else if(!isData && group == "FillGen") {

    int nhadtau = 0;
    TLorentzVector genVec;
    int i = 0;
    for(vec_iter it=goodParts[ival(CUTS::eGTau)]->begin(); it!=goodParts[ival(CUTS::eGTau)]->end(); it++, i++) {

      int nu = goodParts[ival(CUTS::eNuTau)]->at(i);
      if(nu != -1) {
	TLorentzVector tmp1, tmp2;
	tmp1.SetPtEtaPhiE(_Gen->pt->at(*it), _Gen->eta->at(*it), _Gen->phi->at(*it), _Gen->energy->at(*it));
	tmp2.SetPtEtaPhiE(_Gen->pt->at(nu), _Gen->eta->at(nu), _Gen->phi->at(nu), _Gen->energy->at(nu));
	genVec = tmp1 - tmp2;
	histAddVal(genVec.Pt(), "HadTauPt");
	histAddVal(genVec.Eta(), "HadTauEta");
	nhadtau++;

      }
      histAddVal(_Gen->energy->at(*it), "TauEnergy");
      histAddVal(_Gen->pt->at(*it), "TauPt");
      histAddVal(_Gen->eta->at(*it), "TauEta");
      histAddVal(_Gen->phi->at(*it), "TauPhi");
      for(vec_iter it2=it+1; it2!=goodParts[ival(CUTS::eGTau)]->end(); it2++) {
	TLorentzVector genObjt1;
	genObjt1.SetPtEtaPhiE(_Gen->pt->at(*it), _Gen->eta->at(*it), _Gen->phi->at(*it), _Gen->energy->at(*it));
	TLorentzVector genObjt2;
	genObjt2.SetPtEtaPhiE(_Gen->pt->at(*it2), _Gen->eta->at(*it2), _Gen->phi->at(*it2), _Gen->energy->at(*it2));
	histAddVal(diParticleMass(genObjt1,genObjt2, "none"), "DiTauMass");
      }
    }
    histAddVal(goodParts[ival(CUTS::eGTau)]->size(), "NTau");
    histAddVal(nhadtau, "NHadTau");

    for(vec_iter it=goodParts[ival(CUTS::eGMuon)]->begin(); it!=goodParts[ival(CUTS::eGMuon)]->end(); it++) {
      histAddVal(_Gen->energy->at(*it), "MuonEnergy");
      histAddVal(_Gen->pt->at(*it), "MuonPt");
      histAddVal(_Gen->eta->at(*it), "MuonEta");
      histAddVal(_Gen->phi->at(*it), "MuonPhi");
    }
    histAddVal(goodParts[ival(CUTS::eGMuon)]->size(), "NMuon");

          

  } else if(group == "FillTauJet1" || group == "FillTauJet2" || group == "FillMuon1" || group == "FillMuon2" || group == "FillJet1" || group == "FillJet2" || group == "FillBJet" || group == "FillCentralJet" || group == "FillElectron1" || group == "FillElectron2") {
    Particle* part;
    if(group == "FillTauJet1" || group == "FillTauJet2") part=_Tau;
    else if(group == "FillMuon1" || group == "FillMuon2") part=_Muon;
    else if(group == "FillElectron1" || group == "FillElectron2") part=_Electron; 
    else part = _Jet;
    CUTS ePos = fill_num[group];

    for(vec_iter it=goodParts[ival(ePos)]->begin(); it!=goodParts[ival(ePos)]->end(); it++) {
      histAddVal(part->smearP.at(*it).Energy(), "Energy");
      histAddVal(part->smearP.at(*it).Pt(), "Pt");
      histAddVal(part->smearP.at(*it).Eta(), "Eta");
      histAddVal(part->smearP.at(*it).Phi(), "Phi");
      if(part->type == PType::Tau) {
	histAddVal(_Tau->nProngs->at(*it), "NumSignalTracks");
  	histAddVal(_Tau->charge->at(*it), "Charge");
	histAddVal(_Tau->leadChargedCandPt->at(*it), "SeedTrackPt");
      } else if(part->type == PType::Muon) {
	histAddVal(calculateLeptonMetMt(_Muon->smearP.at(*it)), "MetMt");  
      } else if(part->type == PType::Electron) {
          histAddVal(calculateLeptonMetMt(_Electron->smearP.at(*it)), "MetMt"); 
      }
    }
    
    if((ePos == CUTS::eRMuon1 || ePos == CUTS::eRMuon2 || ePos == CUTS::eRTau1 || ePos == CUTS::eRTau2 || ePos == CUTS::eRElec1 || ePos == CUTS::eRElec2 ) && goodParts[ival(ePos)]->size() > 0) {
      double leadpt = 0;
      double leadeta = 0;
      for(vec_iter it=goodParts[ival(ePos)]->begin(); it!=goodParts[ival(ePos)]->end(); it++) {
	if(part->smearP.at(*it).Pt() >= leadpt) {
	  leadpt = part->smearP.at(*it).Pt();
	  leadeta = part->smearP.at(*it).Eta();
	}
      }

      histAddVal(leadpt, "FirstLeadingPt");
      histAddVal(leadeta, "FirstLeadingEta");
    }

    histAddVal(goodParts[ival(ePos)]->size(), "N");


  } else if(group == "FillSusyCuts") {

    histAddVal(sqrt((sumpxForMht * sumpxForMht) + (sumpyForMht * sumpyForMht)), "MHT");
    histAddVal(sumptForHt, "HT");  
    histAddVal(sumptForHt + sqrt((sumpxForMht * sumpxForMht) + (sumpyForMht * sumpyForMht)), "Meff");
    histAddVal(theMETVector.Pt(), "Met");
    if(goodParts[ival(CUTS::eR1stJet)]->at(0) !=-1 && goodParts[ival(CUTS::eR2ndJet)]->at(0) != -1) {
      TLorentzVector DiJet = _Jet->smearP.at(goodParts[ival(CUTS::eR1stJet)]->at(0)) + _Jet->smearP.at(goodParts[ival(CUTS::eR2ndJet)]->at(0));
      histAddVal(absnormPhi(theMETVector.Phi() - DiJet.Phi()), "MetDiJetDeltaPhi");
    }
    
  } else if(group == "FillLeadingJet" && goodParts[ival(CUTS::eSusyCom)]->size() == 0) {
    double eta1 = -100, eta2 = -100;
    double pt1 = 0, pt2 = 0;
    if(goodParts[ival(CUTS::eR1stJet)]->at(0) != -1) {
      pt1 = _Jet->smearP.at(goodParts[ival(CUTS::eR1stJet)]->at(0)).Pt();
      eta1 = _Jet->smearP.at(goodParts[ival(CUTS::eR1stJet)]->at(0)).Eta();
    }
    if(goodParts[ival(CUTS::eR2ndJet)]->at(0) != -1) {
      pt2 = _Jet->smearP.at(goodParts[ival(CUTS::eR2ndJet)]->at(0)).Pt();
      eta2 = _Jet->smearP.at(goodParts[ival(CUTS::eR2ndJet)]->at(0)).Eta();
    }
    histAddVal(pt1, "FirstPt");
    histAddVal(eta1, "FirstEta");

    histAddVal(pt2, "SecondPt");
    histAddVal(eta2, "SecondEta");

  } else if(group == "FillLeadingJet" && goodParts[ival(CUTS::eSusyCom)]->size() != 0) {
    TLorentzVector first = _Jet->smearP.at(goodParts[ival(CUTS::eR1stJet)]->at(0));
    TLorentzVector second = _Jet->smearP.at(goodParts[ival(CUTS::eR2ndJet)]->at(0));
    
    histAddVal(first.Pt(), "FirstPt");
    histAddVal(second.Pt(), "SecondPt");

    histAddVal(first.Eta(), "FirstEta");
    histAddVal(second.Eta(), "SecondEta");
    
    TLorentzVector LeadDiJet = first + second;
    
    histAddVal(LeadDiJet.M(), "Mass"); 
    histAddVal(LeadDiJet.Pt(), "Pt");  
    histAddVal(fabs(first.Eta() - second.Eta()), "DeltaEta"); 
    histAddVal(first.DeltaR(second), "DeltaR");  

    double dphiDijets = absnormPhi(first.Phi() - second.Phi());
    double dphi1 = normPhi(first.Phi() - theMETVector.Phi());
    double dphi2 = normPhi(second.Phi() - theMETVector.Phi());
    double alpha = (LeadDiJet.M() > 0) ? second.Pt() / LeadDiJet.M() : 999999999.0;

    histAddVal(dphiDijets, "LeadSublDijetDphi"); 
    histAddVal2(theMETVector.Pt(),dphiDijets, "MetVsDiJetDeltaPhiLeadSubl");
    histAddVal2(fabs(first.Eta()-second.Eta()), dphiDijets, "DeltaEtaVsDeltaPhiLeadSubl");

    histAddVal(sqrt( pow(dphi1,2.0) + pow((TMath::Pi() - dphi2),2.0) ), "R1");
    histAddVal(sqrt( pow(dphi2,2.0) + pow((TMath::Pi() - dphi1),2.0)), "R2");
    histAddVal(normPhi(first.Phi() - phiForMht), "Dphi1MHT"); 
    histAddVal(normPhi(second.Phi() - phiForMht), "Dphi2MHT");
    histAddVal(dphi1, "Dphi1");
    histAddVal(dphi2, "Dphi2");
    histAddVal2(dphi1,dphi2, "Dphi1VsDphi2");
    histAddVal(alpha, "Alpha");


    //dijet info
  } else if(group == "FillDiJet") {
    double leaddijetmass = 0;
    double leaddijetpt = 0;
    double leaddijetdeltaR = 0;
    double leaddijetdeltaEta = 0;
    double etaproduct = 0;
    for(vec_iter it=goodParts[ival(CUTS::eDiJet)]->begin(); it!=goodParts[ival(CUTS::eDiJet)]->end(); it++) {
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

      histAddVal(DiJet.M(), "Mass");
      histAddVal(DiJet.Pt(), "Pt");
      histAddVal(fabs(jet1.Eta() - jet2.Eta()), "DeltaEta");
      histAddVal(absnormPhi(jet1.Phi() - jet2.Phi()), "DeltaPhi");
      histAddVal(jet1.DeltaR(jet2), "DeltaR");
    }

    histAddVal(leaddijetmass, "LeadMass");
    histAddVal(leaddijetpt, "LeadPt");  
    histAddVal(leaddijetdeltaEta, "LeadDeltaEta");
    histAddVal(leaddijetdeltaR, "LeadDeltaR");
    histAddVal(etaproduct, "LeadEtaProduct");


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


    TLorentzVector part1;
    TLorentzVector part2;
    
    for(vec_iter it=goodParts[ival(ePos)]->begin(); it!=goodParts[ival(ePos)]->end(); it++) {
      int p1= (*it) / BIG_NUM;
      int p2= (*it) % BIG_NUM;

      part1 = lep1->smearP.at(p1);
      part2 = lep2->smearP.at(p2);


      histAddVal2(part1.Pt(),part2.Pt(), "Part1PtVsPart2Pt");
      histAddVal(part1.DeltaR(part2), "DeltaR"); 
      if(group.find("Di") != string::npos) {
	histAddVal((part1.Pt() - part2.Pt()) / (part1.Pt() + part2.Pt()), "DeltaPtDivSumPt");  
	histAddVal(part1.Pt() - part2.Pt(), "DeltaPt");
      } else {
	histAddVal((part2.Pt() - part1.Pt()) / (part1.Pt() + part2.Pt()), "DeltaPtDivSumPt");  
	histAddVal(part2.Pt() - part1.Pt(), "DeltaPt");
      }
      histAddVal(cos(absnormPhi(part2.Phi() - part1.Phi())), "CosDphi");
      histAddVal(absnormPhi(part1.Phi() - theMETVector.Phi()), "Part1MetDeltaPhi");
      histAddVal2(absnormPhi(part1.Phi() - theMETVector.Phi()), cos(absnormPhi(part2.Phi() - part1.Phi())), "Part1MetDeltaPhiVsCosDphi");
      histAddVal(absnormPhi(part2.Phi() - theMETVector.Phi()), "Part2MetDeltaPhi");
      histAddVal(cos(absnormPhi(atan2(part1.Py() - part2.Py(), part1.Px() - part2.Px()) - theMETVector.Phi())), "CosDphi_DeltaPtAndMet");

      double diMass = diParticleMass(part1,part2, distats[digroup].smap.at("HowCalculateMassReco"));
      if(passDiParticleApprox(part1,part2, distats[digroup].smap.at("HowCalculateMassReco"))) {
	histAddVal(diMass, "ReconstructableMass");
      } else {
      	histAddVal(diMass, "NotReconstructableMass");
      }

      double PZeta = getPZeta(part1,part2).first;
      double PZetaVis = getPZeta(part1,part2).second;
      histAddVal(calculateLeptonMetMt(part1), "Part1MetMt");
      histAddVal(calculateLeptonMetMt(part2), "Part2MetMt"); 
      histAddVal(lep2->charge->at(p2) * lep1->charge->at(p1), "OSLS");  
      histAddVal(PZeta, "PZeta"); 
      histAddVal(PZetaVis, "PZetaVis");  
      histAddVal2(PZetaVis,PZeta, "Zeta2D");  
      histAddVal((distats.at(digroup).dmap.at("PZetaCutCoefficient") * PZeta) + (distats.at(digroup).dmap.at("PZetaVisCutCoefficient") * PZetaVis), "Zeta1D");

      if ((goodParts[ival(CUTS::eR1stJet)]->at(0) != -1) && (goodParts[ival(CUTS::eR2ndJet)]->at(0) != -1)) {
	TLorentzVector TheLeadDiJetVect = _Jet->smearP.at(goodParts[ival(CUTS::eR1stJet)]->at(0)) + _Jet->smearP.at(goodParts[ival(CUTS::eR2ndJet)]->at(0));

	histAddVal(absnormPhi(part1.Phi() - TheLeadDiJetVect.Phi()), "Part1DiJetDeltaPhi");
	histAddVal(absnormPhi(part2.Phi() - TheLeadDiJetVect.Phi()), "Part2DiJetDeltaPhi");
	histAddVal(diParticleMass(TheLeadDiJetVect, part1+part2, "VectorSumOfVisProductsAndMet"), "DiJetReconstructableMass"); 
      }
      
      if(lep1->type != PType::Tau) {
	histAddVal(isZdecay(part1, *lep1), "Part1IsZdecay"); 
      }
      if(lep2->type != PType::Tau){ 
	histAddVal(isZdecay(part2, *lep2), "Part2IsZdecay"); 
      }
    }
  }
}


void Analyzer::initializePileupInfo(string MCHisto, string DataHisto) {
  // Filenames must be c_strings below. Here is the conversion from strings to c_strings
  // As you can see above cstr1 corresponds to MC and cstr2 corresponds to data.

  TFile *file1 = new TFile((PUSPACE+MCHisto).c_str());
  TH1D* histmc = (TH1D*)file1->FindObjectAny("NVertices_0");
  if(!histmc) throw std::runtime_error("failed to extract histogram");

  TFile* file2 = new TFile((PUSPACE+DataHisto).c_str());
  TH1D* histdata = (TH1D*)file2->FindObjectAny("NVertices_0");
  if(!histdata) throw std::runtime_error("failed to extract histogram");

  double factor = histmc->Integral() / histdata->Integral();
  double value = 1;
  for(int bin=0; bin<=(histdata->GetXaxis()->GetNbins() + 1); bin++) {
    if(histmc->GetBinContent(bin) == 0) value = 1;
    else value = factor*histdata->GetBinContent(bin) / histmc->GetBinContent(bin);
    hPU->SetBinContent(bin,value);
  }

  file1->Close();
  file2->Close();

}

