#include <functional>
#include <set>
#include <algorithm>

class CUTS {
private:
  // Constructors
  explicit CUTS(int Value);
  // Predicate for finding the corresponding instance
  struct Enum_Predicate_Corresponds:
    public std::unary_function<const CUTS*, bool> {
    Enum_Predicate_Corresponds(int Value): m_value(Value) { }
    bool operator()(const CUTS* E)
    { return E->Get_Value() == m_value; }
  private:
    const int m_value;
  };
  // Comparison functor for the set of instances
  struct Enum_Ptr_Less:
    public std::binary_function<const CUTS*, const CUTS*, bool> {
    bool operator()(const CUTS* E_1, const CUTS* E_2)
    { return E_1->Get_Value() < E_2->Get_Value(); }
  };
public:
  // Compiler-generated copy constructor and operator= are OK.
  typedef std::set<const CUTS*, Enum_Ptr_Less> instances_list;
  typedef typename instances_list::const_iterator const_iterator;
  // Access to int value
  int Get_Value(void) const { return m_value; }
  static int Min(void) { return (*s_instances.begin())->m_value; }
  static int Max(void) { return (*s_instances.rbegin())->m_value; }
  static const CUTS* Corresponding_Enum(int Value)
  { const_iterator it = std::find_if(s_instances.begin(), s_instances.end(), 
				Enum_Predicate_Corresponds(Value));
    return (it != s_instances.end()) ? *it : NULL; }
  static bool Is_Valid_Value(int Value) { return Corresponding_Enum(Value) != NULL; }
  // Number of elements
  static typename instances_list::size_type size(void) { return s_instances.size();
  }
  // Iteration
  static const_iterator begin(void) { return s_instances.begin(); }
  static const_iterator end(void) { return s_instances.end(); }
private:
  int m_value;
  static instances_list s_instances;

 public:
  static const CUTS eGTau;
  static const CUTS eGTop;
  static const CUTS eGElec;
  static const CUTS eGMuon;
  static const CUTS eGZ;
  static const CUTS eGW;
  static const CUTS eGHiggs;

  static const CUTS eRVertex;
  static const CUTS eRMuon1;
  static const CUTS eRMuon2;
  static const CUTS eRElec1;
  static const CUTS eRElec2;
  static const CUTS eRTau1;
  static const CUTS eRTau2;
  static const CUTS eRJet1;
  static const CUTS eRJet2;
  static const CUTS eRCenJet;
  static const CUTS eR1stJet;
  static const CUTS eR2ndJet;
  static const CUTS eRBJet;
  static const CUTS eRTrig1;
  static const CUTS eRTrig2;

  static const CUTS eTMuon1;
  static const CUTS eTMuon2;
  static const CUTS eTElec1;
  static const CUTS eTElec2;
  static const CUTS eTTau1;
  static const CUTS eTTau2;
  
  static const CUTS eDiElec;
  static const CUTS eDiMuon;
  static const CUTS eDiTau;
  static const CUTS eDiJet;
  
  static const CUTS eMuon1Tau1;
  static const CUTS eMuon1Tau2;
  static const CUTS eMuon2Tau1;
  static const CUTS eMuon2Tau2;
  static const CUTS eElec1Tau1;
  static const CUTS eElec1Tau2;
  static const CUTS eElec2Tau1;
  static const CUTS eElec2Tau2;
  static const CUTS eMuon1Elec1;
  static const CUTS eMuon1Elec2;
  static const CUTS eMuon2Elec1;
  static const CUTS eMuon2Elec2;

  static const CUTS eSusyCom;
  static const CUTS eMET;
  static const CUTS eNuTau;

};

inline CUTS::CUTS(int Value):
  m_value(Value)
{
  s_instances.insert(this);
}
