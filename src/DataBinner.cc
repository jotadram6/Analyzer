#include "DataBinner.h"

using namespace std;

Piece1D::Piece1D(string _name, int _bins, double _begin, double _end, int _Nfold) :  
  data(_Nfold*(_bins+2), 0),  begin(_begin), end(_end), bins(_bins), Nfold(_Nfold),fold_width(_bins+2), name(_name) {
  
  width = (end - begin)/bins;
  
}

Piece1D::Piece1D(const Piece1D& rhs) {

  data = rhs.data;
  begin = rhs.begin; 
  end = rhs.end;
  width = rhs.width;
  bins = rhs.bins;
  Nfold = rhs.Nfold;
  fold_width = rhs.fold_width;
  name = rhs.name;
}

Piece1D::~Piece1D() {}
  

int Piece1D::get_bin(double y) {
  return (int)((y-begin)/width);
}

void Piece1D::bin(int folder, double y, double weight) {
  int by = (y < begin) ? -1 : get_bin(y);
  by = (y > end) ? bins+1 : by+1;
  data[fold_width*folder + by] += weight;
}

void Piece1D::write_histogram(vector<string>& folders, TFile* outfile) {
  double entries;
  TH1D histogram(name.c_str(), name.c_str(), bins, begin, end);
  for(int i = 0; i<Nfold; i++) {
    outfile->cd(folders.at(i).c_str());
    entries = 0;
    for(int j = 0; j < (bins+2); j++) {
      histogram.SetBinContent(j, data[i*fold_width + j]);
      entries += data[i*fold_width + j];
    }
    histogram.SetEntries((int)(entries+0.5));
    histogram.Write();

  }
}

/*------------------------------------------------------------------------------------------*/

Piece2D::Piece2D(string _name, int _binx, double _beginx, double _endx, int _biny, double _beginy, double _endy, int _Nfold) :  
  data(_Nfold*(_binx+2)*(_biny+2), 0), beginx(_beginx), endx(_endx), beginy(_beginy), endy(_endy), binx(_binx), biny(_biny), Nfold(_Nfold),
  fold_width((_binx+2)*(_biny+2)), name(_name) {
  
  widthx = (endx - beginx)/binx;
  widthy = (endy - beginy)/biny;
  
}

Piece2D::Piece2D(const Piece2D& rhs) {
  data = rhs.data;
  beginx = rhs.beginx; 
  endx = rhs.endx;
  beginy = rhs.beginy; 
  endy = rhs.endy; 
  widthx = rhs.widthx;
  widthy = rhs.widthy;
  binx = rhs.binx;
  biny = rhs.biny;
  Nfold = rhs.Nfold;
  fold_width = rhs.fold_width;
  name = rhs.name;
 }

Piece2D::~Piece2D() {}
  

int Piece2D::get_bin(double val, bool isX) {
  if(isX) return (int)((val-beginx)/widthx);
  else return (int)((val-beginy)/widthy);
  return -1;
}

void Piece2D::bin(int folder, double x, double y, double weight) {
  int bx = (x < beginx) ? -1 : get_bin(x, true);
  bx = (x > endx) ? binx+1 : bx+1;
  int by = (y < beginy) ? -1 : get_bin(y, false);
  by = (y > endy) ? biny+1 : by+1;

  //  cout << x << ", " << y << " : " << bx << ", " << by << endl;

  data[fold_width*folder + bx + by*(binx+2)] += weight;
  
}

void Piece2D::write_histogram(vector<string>& folders, TFile* outfile) {
  double entries;
  TH2D histogram(name.c_str(), name.c_str(), binx, beginx, endx, biny, beginy, endy);

  for(int i = 0; i<Nfold; i++) {
    
    outfile->cd(folders.at(i).c_str());
    entries = 0;
    for(int j = 0; j < (binx+2)*(biny+2); j++) {
      histogram.SetBinContent(j, data[i*fold_width + j]);
      entries += data[i*fold_width + j];
    }

    histogram.SetEntries((int)(entries+0.5));
    histogram.Write();
  }
}

/*---------------------------------------------------------------------------------------*/			      
 
DataBinner::DataBinner(){}

DataBinner::DataBinner(const DataBinner& rhs) {
  order = rhs.order;

  for(unordered_map<string, DataPiece*>::const_iterator it = rhs.datamap.begin(); it!=rhs.datamap.end(); it++) {
    if(dynamic_cast<Piece1D*>(it->second) != NULL) {
      datamap[it->first] = new Piece1D(*dynamic_cast<Piece1D*>(it->second));
    } else if(dynamic_cast<Piece2D*>(it->second) != NULL) {
      datamap[it->first] = new Piece2D(*dynamic_cast<Piece2D*>(it->second));
    }
  }

}

DataBinner::~DataBinner() {
  for(unordered_map<string, DataPiece*>::iterator it = datamap.begin(); it!=datamap.end(); it++) {
    delete it->second;
    it->second = NULL;

  }
}

void DataBinner::Add_Hist(string shortname, string fullname, int bin, double left, double right, int Nfolder) {
  datamap[shortname] = new Piece1D(fullname, bin, left, right, Nfolder);
  order.push_back(shortname);
}

void DataBinner::Add_Hist(string shortname, string fullname, int binx, double leftx, double rightx, int biny, double lefty, double righty, int Nfolder) {
  datamap[shortname] = new Piece2D(fullname, binx, leftx, rightx, biny, lefty, righty, Nfolder);
  order.push_back(shortname);
}


void DataBinner::AddPoint(string name, int maxfolder, double value, double weight) {
  if(datamap.find(name) == datamap.end())  return;

  for(int i=0; i < maxfolder; i++) {
    datamap[name]->bin(i,value, weight);
  }
}

void DataBinner::AddPoint(string name, int maxfolder, double valuex, double valuey, double weight) {
  if(datamap.find(name) == datamap.end()) return;

  for(int i=0; i < maxfolder; i++) {
    datamap[name]->bin(i,valuex, valuey, weight);
  }
}

void DataBinner::write_histogram(TFile* outfile, vector<string>& folders) {
  for(vector<string>::iterator it = order.begin(); it != order.end(); it++) {
    datamap[*it]->write_histogram(folders, outfile);
  }
}

