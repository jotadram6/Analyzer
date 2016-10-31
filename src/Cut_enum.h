#ifndef CUTS_ENUM_H_
#define CUTS_ENUM_H_  

struct EnumHash
{
    template <typename T>
    std::size_t operator()(T t) const
  {
    return static_cast<std::size_t>(t);
  }
};

enum class CUTS { eGTau, eGTop, eGElec, eGMuon, eGZ, eGW, eGHiggs, eRVertex, eRMuon1, eRMuon2, eRElec1, eRElec2, eRTau1, eRTau2, eRJet1, eRJet2, eRCenJet, eR1stJet, eR2ndJet, eRBJet, eRTrig1, eRTrig2, eTMuon1, eTMuon2, eTElec1, eTElec2, eTTau1, eTTau2, eDiElec, eDiMuon, eDiTau, eDiJet, eMuon1Tau1, eMuon1Tau2, eMuon2Tau1, eMuon2Tau2, eElec1Tau1, eElec1Tau2, eElec2Tau1,eElec2Tau2, eMuon1Elec1, eMuon1Elec2, eMuon2Elec1, eMuon2Elec2, eSusyCom, eMET, eNuTau, enumSize };

#endif
