#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TH1F.h"
#include "TChain.h"
#include "TROOT.h"

#include <iostream>
#include <vector>

#define EPSILON 1E-2

using namespace TMath;
using namespace std;

const Double_t GeV=1e3;
const Double_t MeV=1e0;
const Int_t nLayers = 6, nPixEta = 64, nPixPhi = 64;
const Double_t calorSizeXY = 144 * 16; // in mm

TString file_ref="_vector";

void Matrix2Vector(TString infile = "PFlowNtupleFile_DiffPos.root")
{

  // old variables
  float cell_Energy[nLayers][nPixPhi][nPixEta];
  float cellCh_Energy[nLayers][nPixPhi][nPixEta];
  float cellNu_Energy[nLayers][nPixPhi][nPixEta];
  float Total_Ch_Energy;
  float Total_Nu_Energy;
  float True_Ch_Energy;
  float True_Nu_Energy;
  
  // new variables
  vector<float>  cell_e;
  vector<float>  cell_x;
  vector<float>  cell_y;
  vector<float>  cell_dx;
  vector<float>  cell_dy;
  vector<int>  cell_l;
  vector<float>  cell_chfrac;
  float ECAL1_e;
  float ECAL2_e;
  float ECAL3_e;
  float HCAL1_e;
  float HCAL2_e;
  float HCAL3_e;
  float Cal_e;
  float Caltotal_e;
  
  TString directory = getenv("PWD");
  TChain * oldtree = new TChain("EventTree");   
  TString outfile = infile; outfile.ReplaceAll(".root",file_ref+".root");
  cout << "Will create "<<directory <<"/"<<outfile<< endl;
  oldtree->Add(directory+"/"+infile);

  oldtree->SetBranchAddress("cell_Energy", cell_Energy);
  oldtree->SetBranchAddress("cellCh_Energy", cellCh_Energy);
  oldtree->SetBranchAddress("cellNu_Energy", cellNu_Energy);
  oldtree->SetBranchAddress("Total_Ch_Energy", &Total_Ch_Energy);
  oldtree->SetBranchAddress("Total_Nu_Energy", &Total_Nu_Energy);
  oldtree->SetBranchAddress("True_Ch_Energy", &True_Ch_Energy);
  oldtree->SetBranchAddress("True_Nu_Energy", &True_Nu_Energy);


  TFile * outputfile = new TFile(directory+"/"+outfile,"recreate");
  TTree *newtree = new TTree("physics",Form("Edep"));
  newtree->Branch("cell_e",&cell_e);
  newtree->Branch("cell_x",&cell_x);
  newtree->Branch("cell_y",&cell_y);
  //newtree->Branch("cell_dx",&cell_dx);
  //newtree->Branch("cell_dy",&cell_dy);
  newtree->Branch("cell_l",&cell_l);
  newtree->Branch("cell_chfrac",&cell_chfrac);
  newtree->Branch("ECAL1_e",&ECAL1_e,"ECAL1_e/F");
  newtree->Branch("ECAL2_e",&ECAL2_e,"ECAL2_e/F");
  newtree->Branch("ECAL3_e",&ECAL3_e,"ECAL3_e/F");
  newtree->Branch("HCAL1_e",&HCAL1_e,"HCAL1_e/F");
  newtree->Branch("HCAL2_e",&HCAL2_e,"HCAL2_e/F");
  newtree->Branch("HCAL3_e",&HCAL3_e,"HCAL3_e/F");
  newtree->Branch("Cal_e",&Cal_e,"Cal_e/F");
  newtree->Branch("Caltotal_e",&Caltotal_e,"Caltotal_e/F");
 
  Int_t nentries = (Int_t)oldtree->GetEntries();
  printf("%d entries to be processed\n",nentries);
  for(Int_t entry=0;entry<nentries;entry++) {
	oldtree->GetEntry(entry);
	if ((entry%(nentries/10))==0) printf("%f complete\n",100*(double)entry/nentries);
	
	// Reset variables in the new tree:
	cell_e.clear();
	cell_x.clear();
	cell_y.clear();
	cell_dx.clear();
	cell_dy.clear();
	cell_l.clear();
	cell_chfrac.clear();
	ECAL1_e=0; ECAL2_e=0; ECAL3_e=0;
	HCAL1_e=0; HCAL2_e=0; HCAL3_e=0;
	Cal_e=Total_Ch_Energy+Total_Nu_Energy;
	Caltotal_e=True_Ch_Energy+True_Nu_Energy;
	
	for (int l = 0; l < nLayers; l++){
	for(int bin_eta=0; bin_eta<nPixEta;bin_eta++){
	for(int bin_phi=0; bin_phi<nPixPhi;bin_phi++){
		double e = (double)cell_Energy[l][bin_phi][bin_eta];
		if(e<EPSILON) continue;
		cell_e.push_back(e);
		cell_x.push_back(bin_eta*double(calorSizeXY/nPixEta)-calorSizeXY/2.);
		cell_y.push_back(bin_phi*double(calorSizeXY/nPixPhi)-calorSizeXY/2.);
		cell_dx.push_back(double(calorSizeXY/nPixEta));
		cell_dy.push_back(double(calorSizeXY/nPixPhi));
		cell_l.push_back(l+1);
		cell_chfrac.push_back(cellCh_Energy[l][bin_phi][bin_eta]/e);
		if(l==0) ECAL1_e+=e;
		if(l==1) ECAL2_e+=e;
		if(l==2) ECAL3_e+=e;
		if(l==3) HCAL1_e+=e;
		if(l==4) HCAL2_e+=e;
		if(l==5) HCAL3_e+=e;
	}}}
	
	// Fill output tree:
	newtree->Fill();
  }

  outputfile->cd();
  printf("Writes %s\n",outputfile->GetName());
  newtree->Write();
  outputfile->Write();
  outputfile->Close();
  delete outputfile;
  
}

