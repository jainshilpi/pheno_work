#define Analyse_cxx
#include "Analyse.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include "TSystemFile.h"
#include "TSystemDirectory.h"

void Analyse::Loop(string fname, string fnameout)
{
  
  cout<<"Now inside loop, fname : fout : "<<fname<<" "<<fnameout<<endl;
  

  //bool debug = true;
  bool debug = false;
  
  TTree *tin = list_files(fname.c_str());
  Init(tin);

  TFile *fout = new TFile(fnameout.c_str(), "RECREATE");
  map<string,TH1F*> hmap;
  hmap["phoPt"] = new TH1F("phoPt","",200,0,500);
  hmap["jet1Pt"] = new TH1F("jet1Pt","",200,0,500);
  hmap["jet2Pt"] = new TH1F("jet2Pt","",200,0,500);

  hmap["phoEta"] = new TH1F("phoEta","",100,-2.5,2.5);
  hmap["jet1Eta"] = new TH1F("jet1Eta","",200,-5,5);
  hmap["jet2Eta"] = new TH1F("jet2Eta","",200,-5,5);

  hmap["eta1TimesEta2"] = new TH1F("eta1TimesEta2", "", 200, -25, 25);

  hmap["dEta"] = new TH1F("dEta","",200, -10, 10);

  hmap["jetMass"] = new TH1F("jetMass","",500, 0, 3000);

  hmap["eta_pho_jetSys"] = new TH1F("eta_pho_jetSys", "", 200, 10,10);
  hmap["eta_pho_jet1"] = new TH1F("eta_pho_jet1", "", 200, 10,10);
  hmap["eta_pho_jet2"] = new TH1F("eta_pho_jet2", "", 200, 10,10);
  
  hmap["sumTrkPt"] = new TH1F("sumTrkPt","",500,0,200);
  hmap["ntracks"] = new TH1F("ntracks","",200,0,200);
  
  hmap["trkPt"] = new TH1F("trkPt","",500,0,50);
  hmap["trkEta"] = new TH1F("trkEta","",200,-3,3);

  for(map<string,TH1F*>::iterator it = hmap.begin(); it != hmap.end(); ++it) {
    
    hmap[it->first]->Sumw2();
  }
  
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      int phoInd = selPho();
      int jet1Ind = -99, jet2Ind = -99;
      
      selVBFJets(jet1Ind, jet2Ind);


      ////this is the event selection so do all here now
      if(phoInd>=0 && jet1Ind>=0 && jet2Ind>=0){

	if(debug){
	  cout<<"Found pho and jets, ind "<<phoInd<<" "<<jet1Ind<<" "<<jet2Ind<<endl;
	}

	hmap["phoPt"]->Fill(Photon_PT[phoInd]);
	hmap["phoEta"]->Fill(Photon_Eta[phoInd]);

	hmap["jet1Pt"]->Fill(Jet_PT[jet1Ind]);
	hmap["jet2Pt"]->Fill(Jet_PT[jet2Ind]);

	double jet1eta = Jet_Eta[jet1Ind];
	double jet2eta = Jet_Eta[jet2Ind];

	hmap["jet1Eta"]->Fill(jet1eta);
	hmap["jet2Eta"]->Fill(jet2eta);
	
	hmap["eta1TimesEta2"]->Fill(jet1eta*jet2eta);
	hmap["dEta"]->Fill(jet1eta-jet2eta);

	///jet mass
	TLorentzVector jet1, jet2;
	jet1.SetPtEtaPhiM(Jet_PT[jet1Ind], Jet_Eta[jet1Ind], Jet_Phi[jet1Ind], Jet_Mass[jet1Ind]);
	jet2.SetPtEtaPhiM(Jet_PT[jet2Ind], Jet_Eta[jet2Ind], Jet_Phi[jet2Ind], Jet_Mass[jet2Ind]);

	double jetmass = (jet1 + jet2).M();
	hmap["jetMass"]->Fill(jetmass); 
	double eta_pho_jetSys = Photon_Eta[phoInd] - (jet1eta+jet2eta)/2.;
	hmap["eta_pho_jetSys"]->Fill(eta_pho_jetSys); 

	double eta_pho_jet1 = Photon_Eta[phoInd] - jet1eta;
	double eta_pho_jet2 = Photon_Eta[phoInd] - jet2eta;
	
	hmap["eta_pho_jet1"]->Fill(eta_pho_jet1); 
	hmap["eta_pho_jet2"]->Fill(eta_pho_jet2); 

	///nTracks between the two jets
	double sumTrkPt = 0;
	int ntracks = selTracks(jet1Ind, jet2Ind, sumTrkPt, hmap);
	hmap["ntracks"]->Fill(ntracks);
	hmap["sumTrkPt"]->Fill(sumTrkPt);
      }
      
   }

  for(map<string,TH1F*>::iterator it = hmap.begin(); it != hmap.end(); ++it) {
    hmap[it->first]->Write();
  }

}


int Analyse::selPho(){

  int phoind = -99;
  double maxPt = -99;
  
  for(int ipho=0; ipho<Photon_; ipho++){

    if( fabs(Photon_Eta[ipho]) > 2.5 ) continue;
    //if( Photon_EhadOverEem[ipho] > 0.05 ) continue;
    if( Photon_PT[ipho] < 10. ) continue;
    
    if(maxPt < Photon_PT[ipho])
      {
	maxPt = Photon_PT[ipho];
	phoind = ipho;
      }
  }

  return phoind;
}


void Analyse::selVBFJets(int &jet1ind, int &jet2ind){

  jet1ind = -99;
  jet2ind = -99;
  
  double max1Pt = -99.;
  double max2Pt = -99.;

  ///1st jet selection
  for(int ijet=0; ijet<Jet_; ijet++){

    if( fabs(Jet_Eta[ijet]) > 4.7 ) continue;
    if( Jet_PT[ijet] < 30 ) continue;
      
    if(max1Pt < Jet_PT[ijet])
      {
	max1Pt = Jet_PT[ijet];
	jet1ind = ijet;
      }
  }///for(int ijet=0; ijet<Jet_; ijet++)

  ///2nd jet selection
  for(int ijet=0; ijet<Jet_; ijet++){

    if(ijet==jet1ind) continue;
    
    if( fabs(Jet_Eta[ijet]) > 4.7 ) continue;
    if( Jet_PT[ijet] < 30 ) continue;
    
    if(max2Pt < Jet_PT[ijet])
      {
	max2Pt = Jet_PT[ijet];
	jet2ind = ijet;
      }
  }///for(int ijet=0; ijet<Jet_; ijet++)

  
}

int Analyse::selTracks(int jet1Ind, int jet2Ind, double &sumTrkPt, map<string, TH1F*> hmap){

  int nTrks = 0;

  sumTrkPt = 0;
  
  for(int itrk=0; itrk<Track_; itrk++){
    
    bool sel_midTracks = (Track_Eta[itrk]>Jet_Eta[jet1Ind] && Track_Eta[itrk]<Jet_Eta[jet2Ind]) || (Track_Eta[itrk]>Jet_Eta[jet2Ind] && Track_Eta[itrk]<Jet_Eta[jet1Ind]);
    
    if( !sel_midTracks ) continue;
    //if( Photon_EhadOverEem[itrk] > 0.05 ) continue;
    if( Track_PT[itrk] < 0.5 ) continue;

    nTrks++;
    sumTrkPt += Track_PT[itrk];

    hmap["trkEta"]->Fill(Track_Eta[itrk]);
    hmap["trkPt"]->Fill(Track_PT[itrk]);

  }

  return nTrks;
}


TChain* Analyse::list_files(const char *dirname, const char *ext){
  
  TChain *t = new TChain("Delphes");
  TSystemDirectory dir(dirname, dirname); 
  TList *files = dir.GetListOfFiles(); 
  if (files) 
    { 
      TSystemFile *file; 
      TString fname; 
      TIter next(files); 
      while ((file=(TSystemFile*)next())) { 
	fname = file->GetName(); 
	if (!file->IsDirectory() && fname.EndsWith(ext)) 
	  { 
	    cout << fname.Data() << endl; 
	    t->Add(Form("%s/%s",dirname,fname.Data()));
	  } 
      } 
    } 
  return t;
}
