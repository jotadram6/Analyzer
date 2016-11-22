
void compare() {
  TFile* oldFile = new TFile("oldfile.root");
  TFile* newFile = new TFile("newfile.root");

  TDirectory* dir1 = (TDirectory*)oldFile->Get("METCut");
  TDirectory* dir2 = (TDirectory*)newFile->Get("METCut");
  if(dir1 == NULL) {
    cout << "bad directory" << endl;
    exit(1);
  }

  int count =0, all = 0;
  TIter iter(dir1->GetListOfKeys());
  TKey* key;
  while((key = (TKey*)iter())) {
    all++;
    string histname = key->GetName();

    TH1F* histo1 = (TH1F*)( dir1->FindObjectAny(histname.c_str()) );
    TH1F* histo2 = (TH1F*)( dir2->FindObjectAny(histname.c_str()) );
    
    if(histo1 != 0 && histo2 != 0) {
      bool passed = true;
      if(histo1->GetNbinsX() != histo2->GetNbinsX()) {
	cout << histname << " has " << histo1->GetNbinsX() << " and " << histo2->GetNbinsX() << " bins" << endl;
		count++;
	continue;
      }
      
      for(int i =0; i < histo1->GetNbinsX()+2; i++) {
	if(abs(histo1->GetBinContent(i) - histo2->GetBinContent(i)) > 1) {
	  

	  cout << histname << " has different entries " << endl;
	  count++;
	  passed = false;
	  break;
	}
      }
      
    } else {
      count++;
      cout << histname << " not found" << endl;
    }
  }

  cout << "Didn't match " << count << " out of " << all << " graphs" << endl;

}
