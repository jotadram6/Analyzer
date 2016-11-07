#include "Particle.h"

#define SetBranch(name, variable) BOOM->SetBranchStatus(name, 1);  BOOM->SetBranchAddress(name, &variable);

Particle::Particle(TTree* BOOM, string GenName, string filename) {
  type = PType::None;
  getPartStats(filename);

  SetBranch((GenName+"_pt").c_str(), pt);
  SetBranch((GenName+"_eta").c_str(), eta);
  SetBranch((GenName+"_phi").c_str(), phi);
  SetBranch((GenName+"_energy").c_str(), energy);

}

Generated::Generated(TTree* BOOM, string filename) : Particle(BOOM, "Gen", filename) {
  SetBranch("Gen_pdg_id", pdg_id);
  SetBranch("Gen_motherpdg_id", motherpdg_id);
  SetBranch("Gen_status", status);
  SetBranch("Gen_BmotherIndex", BmotherIndex);
}



Jet::Jet(TTree* BOOM, string filename) : Particle(BOOM, "Jet", filename) {
  type = PType::Jet;
  SetBranch("Jet_neutralHadEnergyFraction", neutralHadEnergyFraction);
  SetBranch("Jet_neutralEmEmEnergyFraction", neutralEmEmEnergyFraction);
  SetBranch("Jet_numberOfConstituents", numberOfConstituents);
  SetBranch("Jet_muonEnergyFraction", muonEnergyFraction);
  SetBranch("Jet_chargedHadronEnergyFraction", chargedHadronEnergyFraction);
  SetBranch("Jet_chargedMultiplicity", chargedMultiplicity);
  SetBranch("Jet_chargedEmEnergyFraction", chargedEmEnergyFraction);
  SetBranch("Jet_partonFlavour", partonFlavour);
  SetBranch("Jet_bDiscriminator_pfCISVV2", bDiscriminator);
}
    
//template <typename T>
Lepton::Lepton(TTree* BOOM, string GenName, string EndName) : Particle(BOOM, GenName, EndName) {
  SetBranch((GenName+"_charge").c_str(), charge);
}


Electron::Electron(TTree* BOOM, string filename) : Lepton(BOOM, "patElectron", filename) {
  type = PType::Electron;
  if(pstats["Elec1"].bmap["DoDiscrByIsolation"] || pstats["Elec2"].bmap["DoDiscrByIsolation"]) {  
    SetBranch("patElectron_isoChargedHadrons", isoChargedHadrons);
    SetBranch("patElectron_isoNeutralHadrons", isoNeutralHadrons);
    SetBranch("patElectron_isoPhotons", isoPhotons);
    SetBranch("patElectron_isoPU", isoPU);
  }
  if(pstats["Elec1"].bmap["DoDiscrByVetoID"] || pstats["Elec2"].bmap["DoDiscrByVetoID"]) {
    SetBranch("patElectron_isPassVeto", isPassVeto);
  }
  if(pstats["Elec1"].bmap["DoDiscrByLooseID"] || pstats["Elec2"].bmap["DoDiscrByLooseID"]) {
    SetBranch("patElectron_isPassLoose", isPassLoose);
  }
  if(pstats["Elec1"].bmap["DoDiscrByMediumID"] || pstats["Elec2"].bmap["DoDiscrByMediumID"]) {
    SetBranch("patElectron_isPassMedium", isPassMedium);
  }
  if(pstats["Elec1"].bmap["DoDiscrByTightID"] || pstats["Elec2"].bmap["DoDiscrByTightID"]) {
    SetBranch("patElectron_isPassTight", isPassTight);
  }
  if(pstats["Elec1"].bmap["DoDiscrByHEEPID"] || pstats["Elec2"].bmap["DoDiscrByHEEPID"]) {
    SetBranch("patElectron_isPassHEEPId", isPassHEEPId);
  }
}


Muon::Muon(TTree* BOOM, string filename) : Lepton(BOOM, "Muon", filename) {
  type = PType::Muon;

  if(pstats["Muon1"].bmap["DoDiscrByTightID"] || pstats["Muon2"].bmap["DoDiscrByTightID"]) {
    SetBranch("Muon_tight", tight);
  }
  if(pstats["Muon1"].bmap["DoDiscrBySoftID"] || pstats["Muon2"].bmap["DoDiscrBySoftID"]) {
    SetBranch("Muon_soft", soft);
  }
  if(pstats["Muon1"].bmap["DoDiscrByIsolation"] || pstats["Muon2"].bmap["DoDiscrByIsolation"]) {
    SetBranch("Muon_isoCharged", isoCharged);
    SetBranch("Muon_isoNeutralHadron", isoNeutralHadron);
    SetBranch("Muon_isoPhoton", isoPhoton);
    SetBranch("Muon_isoPU", isoPU);
  }
}


///////fix against stuff
Taus::Taus(TTree* BOOM, string filename) : Lepton(BOOM, "Tau", filename) {
  type = PType::Tau;

  ////Electron discrimination
  if((pstats["Tau1"].bmap["DoDiscrAgainstElectron"] || pstats["Tau1"].bmap["SelectTausThatAreElectrons"]) &&
     (pstats["Tau2"].bmap["DoDiscrAgainstElectron"] || pstats["Tau2"].bmap["SelectTausThatAreElectrons"]) &&
     (pstats["Tau1"].smap["DiscrAgainstElectron"] == pstats["Tau2"].smap["DiscrAgainstElectron"]) ) {
    SetBranch(("Tau_"+pstats["Tau1"].smap["DiscrAgainstElectron"]).c_str(), againstElectron.first);
    againstElectron.second = againstElectron.first;
  } else {
    if(pstats["Tau1"].bmap["DoDiscrAgainstElectron"] || pstats["Tau1"].bmap["SelectTausThatAreElectrons"]) {
      SetBranch(("Tau_"+pstats["Tau1"].smap["DiscrAgainstElectron"]).c_str(), againstElectron.first);
    }
    if(pstats["Tau2"].bmap["DoDiscrAgainstElectron"] || pstats["Tau2"].bmap["SelectTausThatAreElectrons"]) {
      SetBranch(("Tau_"+pstats["Tau2"].smap["DiscrAgainstElectron"]).c_str(), againstElectron.second);
    }
  }
  ////Muon discrimination
  if((pstats["Tau1"].bmap["DoDiscrAgainstMuon"] || pstats["Tau1"].bmap["SelectTausThatAreMuons"]) &&
     (pstats["Tau2"].bmap["DoDiscrAgainstMuon"] || pstats["Tau2"].bmap["SelectTausThatAreMuons"]) &&
     (pstats["Tau1"].smap["DiscrAgainstMuon"] == pstats["Tau2"].smap["DiscrAgainstMuon"]) ) {
    SetBranch(("Tau_"+pstats["Tau1"].smap["DiscrAgainstMuon"]).c_str(), againstMuon.first);
    againstMuon.second = againstMuon.first;
  } else {
    if(pstats["Tau1"].bmap["DoDiscrAgainstMuon"] || pstats["Tau1"].bmap["SelectTausThatAreMuons"]) {
      SetBranch(("Tau_"+pstats["Tau1"].smap["DiscrAgainstMuon"]).c_str(), againstMuon.first);
    }
    if(pstats["Tau2"].bmap["DoDiscrAgainstMuon"] || pstats["Tau2"].bmap["SelectTausThatAreMuons"]) {
      SetBranch(("Tau_"+pstats["Tau2"].smap["DiscrAgainstMuon"]).c_str(), againstMuon.second);
    }
  }

  /////Isolation discrimination
  if(pstats["Tau1"].bmap["DoDiscrByIsolation"] && pstats["Tau2"].bmap["DoDiscrByIsolation"] &&
     pstats["Tau1"].smap["DiscrByMaxIsolation"] == pstats["Tau2"].smap["DiscrByMaxIsolation"]) {
    SetBranch(("Tau_"+pstats["Tau1"].smap["DiscrByMaxIsolation"]).c_str(), (maxIso.first));      
    maxIso.second = maxIso.first;
  } else {
    if(pstats["Tau1"].bmap["DoDiscrByIsolation"]) {
      SetBranch(("Tau_"+pstats["Tau1"].smap["DiscrByMaxIsolation"]).c_str(), (maxIso.first));      
    }
    if(pstats["Tau2"].bmap["DoDiscrByIsolation"]) {
      SetBranch(("Tau_"+pstats["Tau2"].smap["DiscrByMaxIsolation"]).c_str(), (maxIso.second));
    }
  }      ////min stuff
  if(pstats["Tau1"].bmap["DoDiscrByIsolation"] && pstats["Tau2"].bmap["DoDiscrByIsolation"] &&
     pstats["Tau1"].smap["DiscrByMinIsolation"] == pstats["Tau2"].smap["DiscrByMinIsolation"] &&
     pstats["Tau1"].smap["DiscrByMinIsolation"] != "ZERO") {
    SetBranch(("Tau_"+pstats["Tau1"].smap["DiscrByMinIsolation"]).c_str(), (minIso.first));      
    minIso.second = minIso.first;
  } else {
    if(pstats["Tau1"].bmap["DoDiscrByIsolation"] && pstats["Tau1"].smap["DiscrByMinIsolation"] != "ZERO") {

      SetBranch(("Tau_"+pstats["Tau1"].smap["DiscrByMinIsolation"]).c_str(), (minIso.first));      
    }
    if(pstats["Tau2"].bmap["DoDiscrByIsolation"] && pstats["Tau2"].smap["DiscrByMinIsolation"] != "ZERO") {
      SetBranch(("Tau_"+pstats["Tau2"].smap["DiscrByMaxIsolation"]).c_str(), (minIso.second));
    }
  }      


  SetBranch("Tau_decayModeFindingNewDMs", decayModeFindingNewDMs);
  SetBranch("Tau_nProngs", nProngs);
  SetBranch("Tau_leadChargedCandPt", leadChargedCandPt);

  
  //////NOT USED BRANCHES/////
  //  SetBranchStatus("Tau_decayModeFinding", 1);
  //  SetBranch("Tau_leadChargedCandCharge", &leadChargedCandCharge);
  //  SetBranch("Tau_leadChargedCandEta", &leadChargedCandEta);
  //  SetBranch("Tau_leadChargedCandPhi", &leadChargedCandPhi);
  //  SetBranch("Tau_chargedIsoPtSum", &chargedIsoPtSum);
  //  SetBranch("Tau_neutralIsoPtSum", &neutralIsoPtSum);  
  //  SetBranch("Tau_puCorrPtSum", &puCorrPtSum);
}

void Particle::getPartStats(string filename) {
  typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
  ifstream info_file(filename);
  boost::char_separator<char> sep(", \t");

  if(!info_file) {
    std::cout << "could not open file " << filename <<std::endl;
    return;
  }

  vector<string> stemp;
  string group,line;
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

      if(stemp[1] == "1" || stemp[1] == "true" ) pstats[group].bmap[stemp[0]] = true;
      else if(stemp[1] == "0"  || stemp[1] == "false" ) pstats[group].bmap[stemp[0]]=false; 

      else if(stemp[1].find_first_not_of("0123456789+-.") == string::npos) pstats[group].dmap[stemp[0]]=stod(stemp[1]);
      else pstats[group].smap[stemp[0]] = stemp[1];
    } else  pstats[group].pmap[stemp[0]] = make_pair(stod(stemp[1]), stod(stemp[2]));
  }
  info_file.close();

}




bool Electron::get_Iso(int index, double min, double max) const {
  double maxIsoval = std::max(0.0, isoNeutralHadrons->at(index) + isoPhotons->at(index) - 0.5 * isoPU->at(index));
  double isoSum = (isoChargedHadrons->at(index) + maxIsoval) / smearP.at(index).Pt();
  return (isoSum > min && isoSum < max);
}

bool Muon::get_Iso(int index, double min, double max) const {
  double maxIsoval = std::max(0.0, isoNeutralHadron->at(index) + isoPhoton->at(index) - 0.5 * isoPU->at(index));
  double isoSum = (isoCharged->at(index) + maxIsoval) / smearP.at(index).Pt();
  return (isoSum > min && isoSum < max);
}

bool Taus::get_Iso(int index, double onetwo, double max) const {
  // double maxIsoval = (onetwo == 1) ? maxIso.first->at(i) : maxIso.second->at(i);

  // if
  //   minIso = (ePos == CUTS::eRTau1) ? minIso.first->at(i) : minIso.second->at(i);
  // if (minIso > 0.5) continue;

  // return (maxIso > 0.5 && minIso > 0.5);
  return false;
}
