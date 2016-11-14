#include "Analyzer.h"

int main (int argc, char* argv[]) {
  
  bool setCR = false;

  if(argc < 3) {
    std::cout << "You have entered too little arguments, please type:" << std::endl;
    std::cout << "./Analyzer infile.root outfile.root" << std::endl;
    exit(EXIT_FAILURE);
    ////TODO, add ability to give multiple input files
  }

  string inputname, outputname;
  for(int i = 1; i < argc; ++i) {
    if(argv[i][0] == '-') {
      if(strcmp(argv[i], "-CR") == 0) setCR = true;
      else {
	cout << "wrong option, exiting" << endl;
	exit(0);
      }
    } else if(inputname == "") inputname = argv[i];
    else outputname = argv[i];
  }

  ifstream ifile(inputname);
  if ( !ifile && inputname.find("root://") == string::npos && inputname.find("root\\://") == string::npos) {
    std::cout << "The file '" << inputname << "' doesn't exist" << std::endl;
    exit(EXIT_FAILURE);
  }

  Analyzer testing(inputname, outputname, setCR);
  

  for(int i=0; i < testing.nentries; i++) {
    testing.clear_values();
    testing.preprocess(i);
    testing.fill_histogram();
  }
  testing.printCuts();
  return 0;
}
