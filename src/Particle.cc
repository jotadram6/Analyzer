#include "Particle.h"

Particle::Particle(TTree* BOOM, string GenName, string filename) {
  getPartStats(filename);

  BOOM->SetBranchStatus((GenName+"_pt").c_str(), 1);
  BOOM->SetBranchStatus((GenName+"_eta").c_str(), 1);
  BOOM->SetBranchStatus((GenName+"_phi").c_str(), 1);
  BOOM->SetBranchStatus((GenName+"_energy").c_str(), 1);
  
  BOOM->SetBranchAddress((GenName+"_pt").c_str(), &pt);
  BOOM->SetBranchAddress((GenName+"_eta").c_str(), &eta);
  BOOM->SetBranchAddress((GenName+"_phi").c_str(), &phi);
  BOOM->SetBranchAddress((GenName+"_energy").c_str(), &energy);

}

Generated::Generated(TTree* BOOM, string filename) : Particle(BOOM, "Gen", filename) {
  BOOM->SetBranchStatus("Gen_pdg_id", 1);
  BOOM->SetBranchStatus("Gen_motherpdg_id", 1);
  BOOM->SetBranchStatus("Gen_status", 1);
  BOOM->SetBranchStatus("Gen_BmotherIndex", 1);

  BOOM->SetBranchAddress("Gen_pdg_id", &pdg_id);
  BOOM->SetBranchAddress("Gen_motherpdg_id", &motherpdg_id);
  BOOM->SetBranchAddress("Gen_status", &status);
  BOOM->SetBranchAddress("Gen_BmotherIndex", &BmotherIndex);
}



Jet::Jet(TTree* BOOM, string filename) : Particle(BOOM, "Jet", filename) {
  BOOM->SetBranchStatus("Jet_neutralHadEnergyFraction", 1);
  BOOM->SetBranchStatus("Jet_neutralEmEmEnergyFraction", 1);
  BOOM->SetBranchStatus("Jet_numberOfConstituents", 1);
  BOOM->SetBranchStatus("Jet_muonEnergyFraction", 1);
  BOOM->SetBranchStatus("Jet_chargedHadronEnergyFraction", 1);
  BOOM->SetBranchStatus("Jet_chargedMultiplicity", 1);
  BOOM->SetBranchStatus("Jet_chargedEmEnergyFraction", 1);
  BOOM->SetBranchStatus("Jet_partonFlavour", 1);
  BOOM->SetBranchStatus("Jet_bDiscriminator_pfCISVV2", 1);

  BOOM->SetBranchAddress("Jet_neutralHadEnergyFraction", &neutralHadEnergyFraction);
  BOOM->SetBranchAddress("Jet_neutralEmEmEnergyFraction", &neutralEmEmEnergyFraction);
  BOOM->SetBranchAddress("Jet_numberOfConstituents", &numberOfConstituents);
  BOOM->SetBranchAddress("Jet_muonEnergyFraction", &muonEnergyFraction);
  BOOM->SetBranchAddress("Jet_chargedHadronEnergyFraction", &chargedHadronEnergyFraction);
  BOOM->SetBranchAddress("Jet_chargedMultiplicity", &chargedMultiplicity);
  BOOM->SetBranchAddress("Jet_chargedEmEnergyFraction", &chargedEmEnergyFraction);
  BOOM->SetBranchAddress("Jet_partonFlavour", &partonFlavour);
  BOOM->SetBranchAddress("Jet_bDiscriminator_pfCISVV2", &bDiscriminator);
}
    
//template <typename T>
Lepton::Lepton(TTree* BOOM, string GenName, string EndName) : Particle(BOOM, GenName, EndName) {
  BOOM->SetBranchStatus((GenName+"_charge").c_str(), 1);

  BOOM->SetBranchAddress((GenName+"_charge").c_str(), &charge);
}




Electron::Electron(TTree* BOOM, string filename) : Lepton(BOOM, "patElectron", filename) {
  if(pstats["Elec1"].bmap["DoDiscrByIsolation"] || pstats["Elec2"].bmap["DoDiscrByIsolation"]) {  
    BOOM->SetBranchStatus("patElectron_isoChargedHadrons", 1);
    BOOM->SetBranchStatus("patElectron_isoNeutralHadrons", 1);
    BOOM->SetBranchStatus("patElectron_isoPhotons", 1);
    BOOM->SetBranchStatus("patElectron_isoPU", 1);

    BOOM->SetBranchAddress("patElectron_isoChargedHadrons", &isoChargedHadrons);
    BOOM->SetBranchAddress("patElectron_isoNeutralHadrons", &isoNeutralHadrons);
    BOOM->SetBranchAddress("patElectron_isoPhotons", &isoPhotons);
    BOOM->SetBranchAddress("patElectron_isoPU", &isoPU);
  }
  if(pstats["Elec1"].bmap["DoDiscrByVetoID"] || pstats["Elec2"].bmap["DoDiscrByVetoID"]) {
    BOOM->SetBranchStatus("patElectron_isPassVeto", 1);
    BOOM->SetBranchAddress("patElectron_isPassVeto", &isPassVeto);
  }
  if(pstats["Elec1"].bmap["DoDiscrByLooseID"] || pstats["Elec2"].bmap["DoDiscrByLooseID"]) {
    BOOM->SetBranchStatus("patElectron_isPassLoose", 1);
    BOOM->SetBranchAddress("patElectron_isPassLoose", &isPassLoose);
  }
  if(pstats["Elec1"].bmap["DoDiscrByMediumID"] || pstats["Elec2"].bmap["DoDiscrByMediumID"]) {
    BOOM->SetBranchStatus("patElectron_isPassMedium", 1);
    BOOM->SetBranchAddress("patElectron_isPassMedium", &isPassMedium);
  }
  if(pstats["Elec1"].bmap["DoDiscrByTightID"] || pstats["Elec2"].bmap["DoDiscrByTightID"]) {
    BOOM->SetBranchStatus("patElectron_isPassTight", 1);
    BOOM->SetBranchAddress("patElectron_isPassTight", &isPassTight);
  }
  if(pstats["Elec1"].bmap["DoDiscrByHEEPID"] || pstats["Elec2"].bmap["DoDiscrByHEEPID"]) {
    BOOM->SetBranchStatus("patElectron_isPassHEEPId", 1);
    BOOM->SetBranchAddress("patElectron_isPassHEEPId", &isPassHEEPId);
  }
}


Muon::Muon(TTree* BOOM, string filename) : Lepton(BOOM, "Muon", filename) {
  if(pstats["Muon1"].bmap["DoDiscrByTightID"] || pstats["Muon2"].bmap["DoDiscrByTightID"]) {
    BOOM->SetBranchStatus("Muon_tight", 1);
    BOOM->SetBranchAddress("Muon_tight", &tight);
  }
  if(pstats["Muon1"].bmap["DoDiscrBySoftID"] || pstats["Muon2"].bmap["DoDiscrBySoftID"]) {
    BOOM->SetBranchStatus("Muon_soft", 1);
    BOOM->SetBranchAddress("Muon_soft", &soft);
  }
  if(pstats["Muon1"].bmap["DoDiscrByIsolation"] || pstats["Muon2"].bmap["DoDiscrByIsolation"]) {
    BOOM->SetBranchStatus("Muon_isoCharged", 1);
    BOOM->SetBranchStatus("Muon_isoNeutralHadron", 1);
    BOOM->SetBranchStatus("Muon_isoPhoton", 1);
    BOOM->SetBranchStatus("Muon_isoPU", 1);

    BOOM->SetBranchAddress("Muon_isoCharged", &isoCharged);
    BOOM->SetBranchAddress("Muon_isoNeutralHadron", &isoNeutralHadron);
    BOOM->SetBranchAddress("Muon_isoPhoton", &isoPhoton);
    BOOM->SetBranchAddress("Muon_isoPU", &isoPU);
  }
}


///////fix against stuff
Taus::Taus(TTree* BOOM, string filename) : Lepton(BOOM, "Tau", filename) {
  ////Electron discrimination
  if((pstats["Tau1"].bmap["DoDiscrAgainstElectron"] || pstats["Tau1"].bmap["SelectTausThatAreElectrons"]) &&
     (pstats["Tau2"].bmap["DoDiscrAgainstElectron"] || pstats["Tau2"].bmap["SelectTausThatAreElectrons"]) &&
     (pstats["Tau1"].smap["DiscrAgainstElectron"] == pstats["Tau2"].smap["DiscrAgainstElectron"]) ) {
    BOOM->SetBranchStatus(("Tau_"+pstats["Tau1"].smap["DiscrAgainstElectron"]).c_str(), 1);    
    BOOM->SetBranchAddress(("Tau_"+pstats["Tau1"].smap["DiscrAgainstElectron"]).c_str(), &againstElectron.first);
    againstElectron.second = againstElectron.first;
  } else {
    if(pstats["Tau1"].bmap["DoDiscrAgainstElectron"] || pstats["Tau1"].bmap["SelectTausThatAreElectrons"]) {
      BOOM->SetBranchStatus(("Tau_"+pstats["Tau1"].smap["DiscrAgainstElectron"]).c_str(), 1);    
      BOOM->SetBranchAddress(("Tau_"+pstats["Tau1"].smap["DiscrAgainstElectron"]).c_str(), &againstElectron.first);
    }
    if(pstats["Tau2"].bmap["DoDiscrAgainstElectron"] || pstats["Tau2"].bmap["SelectTausThatAreElectrons"]) {
      BOOM->SetBranchStatus(("Tau_"+pstats["Tau2"].smap["DiscrAgainstElectron"]).c_str(), 1);    
      BOOM->SetBranchAddress(("Tau_"+pstats["Tau2"].smap["DiscrAgainstElectron"]).c_str(), &againstElectron.second);
    }
  }
  ////Muon discrimination
  if((pstats["Tau1"].bmap["DoDiscrAgainstMuon"] || pstats["Tau1"].bmap["SelectTausThatAreMuons"]) &&
     (pstats["Tau2"].bmap["DoDiscrAgainstMuon"] || pstats["Tau2"].bmap["SelectTausThatAreMuons"]) &&
     (pstats["Tau1"].smap["DiscrAgainstMuon"] == pstats["Tau2"].smap["DiscrAgainstMuon"]) ) {
    BOOM->SetBranchStatus(("Tau_"+pstats["Tau1"].smap["DiscrAgainstMuon"]).c_str(), 1);    
    BOOM->SetBranchAddress(("Tau_"+pstats["Tau1"].smap["DiscrAgainstMuon"]).c_str(), &againstMuon.first);
    againstMuon.second = againstMuon.first;
  } else {
    if(pstats["Tau1"].bmap["DoDiscrAgainstMuon"] || pstats["Tau1"].bmap["SelectTausThatAreMuons"]) {
      BOOM->SetBranchStatus(("Tau_"+pstats["Tau1"].smap["DiscrAgainstMuon"]).c_str(), 1);    
      BOOM->SetBranchAddress(("Tau_"+pstats["Tau1"].smap["DiscrAgainstMuon"]).c_str(), &againstMuon.first);
    }
    if(pstats["Tau2"].bmap["DoDiscrAgainstMuon"] || pstats["Tau2"].bmap["SelectTausThatAreMuons"]) {
      BOOM->SetBranchStatus(("Tau_"+pstats["Tau2"].smap["DiscrAgainstMuon"]).c_str(), 1);    
      BOOM->SetBranchAddress(("Tau_"+pstats["Tau2"].smap["DiscrAgainstMuon"]).c_str(), &againstMuon.second);
    }
  }

  /////Isolation discrimination
  if(pstats["Tau1"].bmap["DoDiscrByIsolation"] && pstats["Tau2"].bmap["DoDiscrByIsolation"] &&
     pstats["Tau1"].smap["DiscrByMaxIsolation"] == pstats["Tau2"].smap["DiscrByMaxIsolation"]) {
    BOOM->SetBranchStatus(("Tau_"+pstats["Tau1"].smap["DiscrByMaxIsolation"]).c_str(), 1);
    BOOM->SetBranchAddress(("Tau_"+pstats["Tau1"].smap["DiscrByMaxIsolation"]).c_str(), &(maxIso.first));      
    maxIso.second = maxIso.first;
  } else {
    if(pstats["Tau1"].bmap["DoDiscrByIsolation"]) {
      BOOM->SetBranchStatus(("Tau_"+pstats["Tau1"].smap["DiscrByMaxIsolation"]).c_str(), 1);
      BOOM->SetBranchAddress(("Tau_"+pstats["Tau1"].smap["DiscrByMaxIsolation"]).c_str(), &(maxIso.first));      
    }
    if(pstats["Tau2"].bmap["DoDiscrByIsolation"]) {
      BOOM->SetBranchStatus(("Tau_"+pstats["Tau2"].smap["DiscrByMaxIsolation"]).c_str(), 1);
      BOOM->SetBranchAddress(("Tau_"+pstats["Tau2"].smap["DiscrByMaxIsolation"]).c_str(), &(maxIso.second));
    }
  }      ////min stuff
  if(pstats["Tau1"].bmap["DoDiscrByIsolation"] && pstats["Tau2"].bmap["DoDiscrByIsolation"] &&
     pstats["Tau1"].smap["DiscrByMinIsolation"] == pstats["Tau2"].smap["DiscrByMinIsolation"] &&
     pstats["Tau1"].smap["DiscrByMinIsolation"] != "ZERO") {
    BOOM->SetBranchStatus(("Tau_"+pstats["Tau1"].smap["DiscrByMinIsolation"]).c_str(), 1);
    BOOM->SetBranchAddress(("Tau_"+pstats["Tau1"].smap["DiscrByMinIsolation"]).c_str(), &(minIso.first));      
    minIso.second = minIso.first;
  } else {
    if(pstats["Tau1"].bmap["DoDiscrByIsolation"] && pstats["Tau1"].smap["DiscrByMinIsolation"] != "ZERO") {
      BOOM->SetBranchStatus(("Tau_"+pstats["Tau1"].smap["DiscrByMinIsolation"]).c_str(), 1);
      BOOM->SetBranchAddress(("Tau_"+pstats["Tau1"].smap["DiscrByMinIsolation"]).c_str(), &(minIso.first));      
    }
    if(pstats["Tau2"].bmap["DoDiscrByIsolation"] && pstats["Tau2"].smap["DiscrByMinIsolation"] != "ZERO") {
      BOOM->SetBranchStatus(("Tau_"+pstats["Tau2"].smap["DiscrByMaxIsolation"]).c_str(), 1);
      BOOM->SetBranchAddress(("Tau_"+pstats["Tau2"].smap["DiscrByMaxIsolation"]).c_str(), &(minIso.second));
    }
  }      



  BOOM->SetBranchStatus("Tau_decayModeFindingNewDMs", 1);
  BOOM->SetBranchStatus("Tau_nProngs", 1);
  BOOM->SetBranchStatus("Tau_leadChargedCandPt", 1);

  BOOM->SetBranchAddress("Tau_decayModeFindingNewDMs", &decayModeFindingNewDMs);
  BOOM->SetBranchAddress("Tau_nProngs", &nProngs);
  BOOM->SetBranchAddress("Tau_leadChargedCandPt", &leadChargedCandPt);

  
  //////NOT USED BRANCHES/////
  //  BOOM->SetBranchStatus("Tau_decayModeFinding", 1);
  //  BOOM->SetBranchAddress("Tau_leadChargedCandCharge", &leadChargedCandCharge);
  //  BOOM->SetBranchAddress("Tau_leadChargedCandEta", &leadChargedCandEta);
  //  BOOM->SetBranchAddress("Tau_leadChargedCandPhi", &leadChargedCandPhi);
  //  BOOM->SetBranchAddress("Tau_chargedIsoPtSum", &chargedIsoPtSum);
  //  BOOM->SetBranchAddress("Tau_neutralIsoPtSum", &neutralIsoPtSum);  
  //  BOOM->SetBranchAddress("Tau_puCorrPtSum", &puCorrPtSum);
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




