#ifndef Histo_h
#define Histo_h

// system include files
#include <memory>

// user include files
#include <Math/VectorUtil.h>
#include <fstream>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <string>
#include <vector>
#include <unordered_map>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <regex>
#include "DataBinner.h"
#include "tokenizer.hpp"


using namespace std;

class Histogramer {

 public:
  Histogramer();
  Histogramer(int, string, string, string, bool);
  Histogramer(const Histogramer&);
  Histogramer& operator=(const Histogramer&);
  ~Histogramer();
  
  vector<int> get_folders();
  
  unordered_map<string,pair<int,int>>* get_cuts();
  vector<string>* get_order();
  vector<string>* get_groups();
  void setupSVFit(string, string, double, double, double);
  void addVal(double, string, int, string, double);
  void addVal(double, double, string, int, string, double);
  void fill_histogram();
  void setControlRegions();

 private:
  TFile * outfile;
  string outname;
  int NFolders;
  int Npdf;
  bool isData;
  unordered_map<string, pair<int,int>> cuts;
  vector<string> cut_order;
  vector<string> folders;
  vector<int> folder_num;
  unordered_map<string, DataBinner*> data;
  vector<string> data_order;

  void read_hist(string);
  void read_cuts(string);

  
  string extractHistname(string, string);
};

#endif

