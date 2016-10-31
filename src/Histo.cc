#include "Histo.h"

Histogramer::Histogramer() {
}

Histogramer::Histogramer(int _Npdf, string histname, string cutname, string outfilename, bool _isData): outname(outfilename), Npdf(_Npdf), isData(_isData) {
  outfile = NULL;
  read_cuts(cutname);
  NFolders = folders.size();

  read_hist(histname);

}

Histogramer& Histogramer::operator=(const Histogramer& rhs) {
  if(this == &rhs) return *this;

  outname = rhs.outname;
  NFolders = rhs.NFolders;
  isData = rhs.isData;
  
  cuts = rhs.cuts;
  cut_order = rhs.cut_order;
  folders = rhs.folders;
  folder_num = rhs.folder_num;
  data_order.reserve(rhs.data_order.size());
  data_order = rhs.data_order;

  for(unordered_map<string, DataBinner*>::const_iterator mit = rhs.data.begin(); mit != rhs.data.end(); mit++) {
    data[mit->first] = new DataBinner(*(mit->second));
  }
  return *this;
}

Histogramer::Histogramer(const Histogramer& rhs) {

  NFolders = rhs.NFolders;
  outname = rhs.outname;
  isData = rhs.isData;
  cuts = rhs.cuts;
  cut_order = rhs.cut_order;
  folders = rhs.folders;
  folder_num = rhs.folder_num;
  data_order = rhs.data_order;  
  for(unordered_map<string, DataBinner*>::const_iterator mit = rhs.data.begin(); mit != rhs.data.end(); mit++) {
    data[mit->first] = new DataBinner(*(mit->second));
  }
}



Histogramer::~Histogramer() {
  if(outfile != NULL)
    outfile->Close();
  
  for(vector<string>::iterator it = data_order.begin(); it != data_order.end(); it++) {
    delete data[*it];
    data[*it] = NULL;
  }
}


void Histogramer::read_hist(string filename) {
  typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
  ifstream info_file(filename);
  boost::char_separator<char> sep(", \t");  

  if(!info_file) {
    cout << "ERROR: Didn't Read Histo File!" << endl;
    cout << filename << endl;
    exit(1);
  }

  vector<string> stemp;
  string group,line;
  bool accept = false;
  while(getline(info_file, line)) {
    tokenizer tokens(line, sep);
    stemp.clear();
    for(tokenizer::iterator iter = tokens.begin();iter != tokens.end(); iter++) {
      if( ((*iter)[0] == '/' && (*iter)[0] == '/') || ((*iter)[0] == '#') ) break;
      stemp.push_back(*iter);
    }

    if(stemp.size() == 0) continue;
    else if(stemp.size() == 2) {
      group = stemp[0];
      accept = stoi(stemp[1]) && !(isData && group.find("Gen") != string::npos);
      if(accept) {
	data[group] = new DataBinner();
	data_order.push_back(group); 
      }
    } else if(!accept) continue;
    else if(stemp.size() == 4) {
      string name = extractHistname(group, stemp[0]);
      data[group]->Add_Hist(name, stemp[0], stod(stemp[1]), stod(stemp[2]), stod(stemp[3]), NFolders);
    } else if(stemp.size() == 7) {
      string name = extractHistname(group, stemp[0]);
      data[group]->Add_Hist(name, stemp[0], stod(stemp[1]), stod(stemp[2]), stod(stemp[3]),stod(stemp[4]), stod(stemp[5]), stod(stemp[6]), NFolders);
    }
  }

  info_file.close(); 
}

string Histogramer::extractHistname(string group, string histo) {
  regex reg ("((Tau|Muon|Electron)+(1|2)+)");
  smatch m;

  string stringkey = group.erase(0,4);
  regex first (stringkey+"(_)?");
  histo = regex_replace(histo,first, "");
  if(stringkey.find("Di") != string::npos) {
    stringkey=stringkey+"1"+stringkey+"2";
  }

  int i=1;
  while(regex_search(stringkey,m,reg)) {
    regex key(m[0].str());
    histo = regex_replace(histo,key,"Part"+to_string(i));
    stringkey = m.suffix().str();
    i++;
  }
  return histo;
}

void Histogramer::setupSVFit(string group, string name,  double bins, double min, double max) {
  data["Fill" + group]->Add_Hist(name, group + "SVFit_Mass", bins, min, max, NFolders);
}


void Histogramer::read_cuts(string filename) {
  typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
  ifstream info_file(filename);
  boost::char_separator<char> sep(", \t");  

  if(!info_file) {
    cout << "ERROR: Didn't Read Histo File!" << endl;
    cout << filename << endl;
    exit(1);
  }

  vector<string> stemp;
  string name,line;
  int i = 0;

  while(getline(info_file, line)) {
    tokenizer tokens(line, sep);
    stemp.clear();
    for(tokenizer::iterator iter = tokens.begin();iter != tokens.end(); iter++) {
      if( ((*iter)[0] == '/' && (*iter)[0] == '/') || ((*iter)[0] == '#') ) break;
      stemp.push_back(*iter);
    }

    if(stemp.size() == 3) {
      name = stemp[0];
      if(name[0]=='*' && name[1]=='*' && name[2]=='*') {
	name.erase(0,3);
	folders.push_back(name);
	folder_num.push_back(i);
      } else if(stemp[1] == "0" && stemp[2] == "-1") continue;   ////remove unnecessary cuts
      cuts[name] = std::make_pair(stoi(stemp[1]),stoi(stemp[2]));
      cut_order.push_back(name);
      i++;
    }
  }

  if(folders.size() == 0 && cut_order.size() == 0){
    folders.push_back(name);
    folder_num.push_back(i-1);
  } else if(cut_order.size() != 0 && ( folders.size() == 0 || folders.back() != cut_order.back())) {
    folders.push_back(cut_order.back());
    folder_num.push_back(cut_order.size() -1);
  }

  info_file.close(); 
 }

void Histogramer::fill_histogram() {

  outfile = new TFile(outname.c_str(), "RECREATE");
  for(vector<string>::iterator it = folders.begin(); it != folders.end(); it++) {
    outfile->mkdir( it->c_str() );	
  }
  
  for(vector<string>::iterator it = data_order.begin(); it != data_order.end(); ++it) {
    data[*it]->write_histogram(outfile, folders);
  }
  outfile->Close();
}  


unordered_map<string,pair<int,int>>* Histogramer::get_cuts() {
  return &cuts;
}

vector<string>* Histogramer::get_order() {
  return &cut_order;
}

vector<string>* Histogramer::get_groups() {
  return &data_order;
}


void Histogramer::addVal(double value, string group, int maxcut, string histn, double weight) {
  int maxFolder=0;

  for(int i = 0; i < NFolders; i++) {
    if(maxcut > folder_num[i]) maxFolder++;
    else break;
  }

  data[group]->AddPoint(histn, maxFolder, value, weight);
}

void Histogramer::addVal(double valuex, double valuey, string group, int maxcut, string histn, double weight) {
  int maxFolder=0;

  for(int i = 0; i < NFolders; i++) {
    if(maxcut > folder_num[i]) maxFolder++;
    else break;
  }
  data[group]->AddPoint(histn, maxFolder, valuex, valuey, weight);
}


    
