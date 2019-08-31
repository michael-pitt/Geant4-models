#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TH1F.h"
#include "TChain.h"
#include "TROOT.h"

#include <iostream>
#include <vector>

#define EPSILON 1E-4

using namespace TMath;
using namespace std;

// number of layers
const Int_t nLayers = 6;

// Old layer dimentions
const nPixEta = 64, nPixPhi = 64;

// New layer dimentions
TString file_ref="ATLAS_resolution"; 
const nPixEta_new[nLayers] = {64, 64, 32, 16, 16, 8};
const nPixPhi_new[nLayers] = {32, 64, 32, 16, 16, 8};

// common detector parameters
const Double_t GeV=1e3;
const Double_t MeV=1e0;
const Double_t calorSizeXY = 144 * 16; // in mm
const Double_t cell_dEta = calorSizeXY / nPixEta;
const Double_t cell_dPhi = calorSizeXY / nPixPhi;
const int kNaxPrimaries = 2;

void Matrix2Matrix(TString infile = "events_6D64x64.root")
{

// old branches
  float cell_Energy[nLayers][nPixPhi][nPixEta];
  float cellCh_Energy[nLayers][nPixPhi][nPixEta];
  float cellNu_Energy[nLayers][nPixPhi][nPixEta];
  
// nea branches
  float cell_Energy_L1[nPixPhi_new[0]][nPixEta_new[0]];
  float cell_Energy_L2[nPixPhi_new[1]][nPixEta_new[1]];
  float cell_Energy_L3[nPixPhi_new[2]][nPixEta_new[2]];
  float cell_Energy_L4[nPixPhi_new[3]][nPixEta_new[3]];
  float cell_Energy_L5[nPixPhi_new[4]][nPixEta_new[4]];
  float cell_Energy_L6[nPixPhi_new[5]][nPixEta_new[5]];
  
  float cellCh_Energy_L1[nPixPhi_new[0]][nPixEta_new[0]];
  float cellCh_Energy_L2[nPixPhi_new[1]][nPixEta_new[1]];
  float cellCh_Energy_L3[nPixPhi_new[2]][nPixEta_new[2]];
  float cellCh_Energy_L4[nPixPhi_new[3]][nPixEta_new[3]];
  float cellCh_Energy_L5[nPixPhi_new[4]][nPixEta_new[4]];
  float cellCh_Energy_L6[nPixPhi_new[5]][nPixEta_new[5]];

  float cellNu_Energy_L1[nPixPhi_new[0]][nPixEta_new[0]];
  float cellNu_Energy_L2[nPixPhi_new[1]][nPixEta_new[1]];
  float cellNu_Energy_L3[nPixPhi_new[2]][nPixEta_new[2]];
  float cellNu_Energy_L4[nPixPhi_new[3]][nPixEta_new[3]];
  float cellNu_Energy_L5[nPixPhi_new[4]][nPixEta_new[4]];
  float cellNu_Energy_L6[nPixPhi_new[5]][nPixEta_new[5]];
 
  TString directory = getenv("PWD");
  TChain * oldtree = new TChain("EventTree");   
  TString outfile = infile; outfile.ReplaceAll(".root",file_ref+".root");
  cout << "Will create "<<directory <<"/"<<outfile<< endl;
  oldtree->Add(directory+"/"+infile);

  oldtree->SetBranchAddress("cell_Energy", cell_Energy);
  oldtree->SetBranchAddress("cellCh_Energy", cellCh_Energy);
  oldtree->SetBranchAddress("cellNu_Energy", cellNu_Energy);

  TFile * outputfile = new TFile(directory+"/"+outfile,"recreate");
  oldtree->LoadTree(0);
  TTree *newtree = oldtree->GetTree()->CloneTree(0);
  newtree->SetBranchStatus("*",0);
  newtree->SetBranchStatus("n_primaries",1);
  newtree->SetBranchStatus("par_*",1);
  
  newtree->Branch("cell_Energy_L1",cell_Energy_L1,Form("cell_Energy_L1[%d][%d]/F",nPixPhi_new[0],nPixEta_new[0]));
  newtree->Branch("cell_Energy_L2",cell_Energy_L2,Form("cell_Energy_L2[%d][%d]/F",nPixPhi_new[1],nPixEta_new[1]));
  newtree->Branch("cell_Energy_L3",cell_Energy_L3,Form("cell_Energy_L3[%d][%d]/F",nPixPhi_new[2],nPixEta_new[2]));
  newtree->Branch("cell_Energy_L4",cell_Energy_L4,Form("cell_Energy_L4[%d][%d]/F",nPixPhi_new[3],nPixEta_new[3]));
  newtree->Branch("cell_Energy_L5",cell_Energy_L5,Form("cell_Energy_L5[%d][%d]/F",nPixPhi_new[4],nPixEta_new[4]));
  newtree->Branch("cell_Energy_L6",cell_Energy_L6,Form("cell_Energy_L6[%d][%d]/F",nPixPhi_new[5],nPixEta_new[5]));
  
  newtree->Branch("cellCh_Energy_L1",cellCh_Energy_L1,Form("cellCh_Energy_L1[%d][%d]/F",nPixPhi_new[0],nPixEta_new[0]));
  newtree->Branch("cellCh_Energy_L2",cellCh_Energy_L2,Form("cellCh_Energy_L2[%d][%d]/F",nPixPhi_new[1],nPixEta_new[1]));
  newtree->Branch("cellCh_Energy_L3",cellCh_Energy_L3,Form("cellCh_Energy_L3[%d][%d]/F",nPixPhi_new[2],nPixEta_new[2]));
  newtree->Branch("cellCh_Energy_L4",cellCh_Energy_L4,Form("cellCh_Energy_L4[%d][%d]/F",nPixPhi_new[3],nPixEta_new[3]));
  newtree->Branch("cellCh_Energy_L5",cellCh_Energy_L5,Form("cellCh_Energy_L5[%d][%d]/F",nPixPhi_new[4],nPixEta_new[4]));
  newtree->Branch("cellCh_Energy_L6",cellCh_Energy_L6,Form("cellCh_Energy_L6[%d][%d]/F",nPixPhi_new[5],nPixEta_new[5]));
  
  newtree->Branch("cellNu_Energy_L1",cellNu_Energy_L1,Form("cellNu_Energy_L1[%d][%d]/F",nPixPhi_new[0],nPixEta_new[0]));
  newtree->Branch("cellNu_Energy_L2",cellNu_Energy_L2,Form("cellNu_Energy_L2[%d][%d]/F",nPixPhi_new[1],nPixEta_new[1]));
  newtree->Branch("cellNu_Energy_L3",cellNu_Energy_L3,Form("cellNu_Energy_L3[%d][%d]/F",nPixPhi_new[2],nPixEta_new[2]));
  newtree->Branch("cellNu_Energy_L4",cellNu_Energy_L4,Form("cellNu_Energy_L4[%d][%d]/F",nPixPhi_new[3],nPixEta_new[3]));
  newtree->Branch("cellNu_Energy_L5",cellNu_Energy_L5,Form("cellNu_Energy_L5[%d][%d]/F",nPixPhi_new[4],nPixEta_new[4]));
  newtree->Branch("cellNu_Energy_L6",cellNu_Energy_L6,Form("cellNu_Energy_L6[%d][%d]/F",nPixPhi_new[5],nPixEta_new[5]));

  Int_t nentries = (Int_t)oldtree->GetEntries();
  printf("%d entries to be processed\n",nentries);
  for(Int_t entry=0;entry<nentries;entry++) {
	oldtree->GetEntry(entry);
	if ((entry%(nentries/10))==0) printf("%f complete\n",100*(double)entry/nentries);
	
	// Reset variables in the new tree:
	for(int i=0;i<nPixPhi_new[0];i++){for(int j=0;j<nPixEta_new[0];j++){
		cellCh_Energy_L1[i][j] = 0; cellNu_Energy_L1[i][j] = 0;
	}}
	for(int i=0;i<nPixPhi_new[1];i++){for(int j=0;j<nPixEta_new[1];j++){
		cellCh_Energy_L2[i][j] = 0; cellNu_Energy_L2[i][j] = 0;
	}}
	for(int i=0;i<nPixPhi_new[2];i++){for(int j=0;j<nPixEta_new[2];j++){
		cellCh_Energy_L3[i][j] = 0; cellNu_Energy_L3[i][j] = 0;
	}}
	for(int i=0;i<nPixPhi_new[3];i++){for(int j=0;j<nPixEta_new[3];j++){
		cellCh_Energy_L4[i][j] = 0; cellNu_Energy_L4[i][j] = 0;
	}}
	for(int i=0;i<nPixPhi_new[4];i++){for(int j=0;j<nPixEta_new[4];j++){
		cellCh_Energy_L5[i][j] = 0; cellNu_Energy_L5[i][j] = 0;
	}}
	for(int i=0;i<nPixPhi_new[5];i++){for(int j=0;j<nPixEta_new[5];j++){
		cellCh_Energy_L6[i][j] = 0; cellNu_Energy_L6[i][j] = 0;
	}}


/*



	
	int n_Cells = cell_e->size();
	for (int i = 0; i < n_Cells; i++){
		float e = (float)cell_e->at(i), f = (float)cell_chfrac->at(i);
		int layer = cell_l->at(i);
		double x = cell_x->at(i), y = cell_y->at(i);	
		double dx = cell_dx->at(i), dy = cell_dy->at(i);
		int BinEtaMin = (x-dx/2 + calorSizeXY/2) / cell_dEta;
		int BinEtaMax = (x+dx/2 + calorSizeXY/2) / cell_dEta;
		int BinPhiMin = (y-dy/2 + calorSizeXY/2) / cell_dPhi;
		int BinPhiMax = (y+dy/2 + calorSizeXY/2) / cell_dPhi;
		if(BinEtaMax==BinEtaMin) BinEtaMax++;
		if(BinPhiMax==BinPhiMin) BinPhiMax++;
		float denom = (BinEtaMax-BinEtaMin)*(BinPhiMax-BinPhiMin); //democratic
		for(int bin_eta=BinEtaMin;bin_eta<BinEtaMax;bin_eta++){
		for(int bin_phi=BinPhiMin;bin_phi<BinPhiMax;bin_phi++){
		  cellCh_Energy[layer-1][bin_phi][bin_eta] += e*f/denom;
          cellNu_Energy[layer-1][bin_phi][bin_eta] += e*(1.0-f)/denom;
		}}
	}
	
	// Fill output tree:
	for(int i=0;i<nPixPhi;i++){for(int j=0;j<nPixEta;j++){for(int k=0;k<nLayers;k++){
		cell_Energy[k][i][j] = cellCh_Energy[k][i][j]+cellNu_Energy[k][i][j];
	}}}
*/
	newtree->Fill();
  }

  outputfile->cd();
  printf("Writes %s\n",outputfile->GetName());
  newtree->Write();
  outputfile->Write();
  outputfile->Close();
  delete outputfile;
  
}

