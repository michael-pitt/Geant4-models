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
const Int_t nPixEta = 64, nPixPhi = 64;

// New layer dimentions
TString file_ref="ATLAS_resolution"; 
const Int_t nPixEta_new[nLayers] = {64, 64, 32, 16, 16, 8};
const Int_t nPixPhi_new[nLayers] = {32, 64, 32, 16, 16, 8};

// common detector parameters
const Double_t GeV=1e3;
const Double_t MeV=1e0;
const Double_t calorSizeXY = 144 * 16; // in mm
const Double_t cell_dEta = calorSizeXY / nPixEta;
const Double_t cell_dPhi = calorSizeXY / nPixPhi;
const int kNaxPrimaries = 2;

void Matrix2Matrix6L(TString infile = "events_6D64x64.root")
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
  TString outfile = infile; outfile.ReplaceAll(".root","_"+file_ref+".root");
  cout << "Will create "<<directory <<"/"<<outfile<< endl;
  oldtree->Add(directory+"/"+infile);

  oldtree->SetBranchAddress("cell_Energy", cell_Energy);
  oldtree->SetBranchAddress("cellCh_Energy", cellCh_Energy);
  oldtree->SetBranchAddress("cellNu_Energy", cellNu_Energy);

  TFile * outputfile = new TFile(directory+"/"+outfile,"recreate");
  //oldtree->LoadTree(0);
  //TTree *newtree = oldtree->GetTree()->CloneTree(0);
  //newtree->SetBranchStatus("*",0);
  //newtree->SetBranchStatus("n_primaries",1);
  //newtree->SetBranchStatus("par_*",1);
  TTree *newtree = new TTree("EventTree","Detector images with varying resolution per layer");

  
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
	
	// Fill all variables, depending on the output resolution
	// NOTE: output resolution is always lower than the input!!!
	for(int i=0;i<nPixPhi_new[0];i++){for(int j=0;j<nPixEta_new[0];j++){
		cellCh_Energy_L1[i][j] = 0; cellNu_Energy_L1[i][j] = 0;
		int phi_bins = nPixPhi / nPixPhi_new[0];
		int eta_bins = nPixEta / nPixEta_new[0];
		for(int i_phi = i*phi_bins; i_phi < (i+1)*phi_bins; i_phi++){
		for(int j_eta = j*eta_bins; j_eta < (j+1)*eta_bins; j_eta++){
			cellCh_Energy_L1[i][j] += cellCh_Energy[0][i_phi][j_eta];
			cellNu_Energy_L1[i][j] += cellNu_Energy[0][i_phi][j_eta];
		}}
		cell_Energy_L1[i][j] = cellCh_Energy_L1[i][j] + cellNu_Energy_L1[i][j];
	}}

	for(int i=0;i<nPixPhi_new[1];i++){for(int j=0;j<nPixEta_new[1];j++){
		cellCh_Energy_L2[i][j] = 0; cellNu_Energy_L2[i][j] = 0;
		int phi_bins = nPixPhi / nPixPhi_new[1];
		int eta_bins = nPixEta / nPixEta_new[1];
		for(int i_phi = i*phi_bins; i_phi < (i+1)*phi_bins; i_phi++){
		for(int j_eta = j*eta_bins; j_eta < (j+1)*eta_bins; j_eta++){
			cellCh_Energy_L2[i][j] += cellCh_Energy[0][i_phi][j_eta];
			cellNu_Energy_L2[i][j] += cellNu_Energy[0][i_phi][j_eta];
		}}
		cell_Energy_L2[i][j] = cellCh_Energy_L2[i][j] + cellNu_Energy_L2[i][j];
	}}
	for(int i=0;i<nPixPhi_new[2];i++){for(int j=0;j<nPixEta_new[2];j++){
		cellCh_Energy_L3[i][j] = 0; cellNu_Energy_L3[i][j] = 0;
		int phi_bins = nPixPhi / nPixPhi_new[2];
		int eta_bins = nPixEta / nPixEta_new[2];
		for(int i_phi = i*phi_bins; i_phi < (i+1)*phi_bins; i_phi++){
		for(int j_eta = j*eta_bins; j_eta < (j+1)*eta_bins; j_eta++){
			cellCh_Energy_L3[i][j] += cellCh_Energy[0][i_phi][j_eta];
			cellNu_Energy_L3[i][j] += cellNu_Energy[0][i_phi][j_eta];
		}}
		cell_Energy_L3[i][j] = cellCh_Energy_L3[i][j] + cellNu_Energy_L3[i][j];
	}}
	for(int i=0;i<nPixPhi_new[3];i++){for(int j=0;j<nPixEta_new[3];j++){
		cellCh_Energy_L4[i][j] = 0; cellNu_Energy_L4[i][j] = 0;
		int phi_bins = nPixPhi / nPixPhi_new[3];
		int eta_bins = nPixEta / nPixEta_new[3];
		for(int i_phi = i*phi_bins; i_phi < (i+1)*phi_bins; i_phi++){
		for(int j_eta = j*eta_bins; j_eta < (j+1)*eta_bins; j_eta++){
			cellCh_Energy_L4[i][j] += cellCh_Energy[0][i_phi][j_eta];
			cellNu_Energy_L4[i][j] += cellNu_Energy[0][i_phi][j_eta];
		}}
		cell_Energy_L4[i][j] = cellCh_Energy_L4[i][j] + cellNu_Energy_L4[i][j];
	}}
	for(int i=0;i<nPixPhi_new[4];i++){for(int j=0;j<nPixEta_new[4];j++){
		cellCh_Energy_L5[i][j] = 0; cellNu_Energy_L5[i][j] = 0;
		int phi_bins = nPixPhi / nPixPhi_new[4];
		int eta_bins = nPixEta / nPixEta_new[4];
		for(int i_phi = i*phi_bins; i_phi < (i+1)*phi_bins; i_phi++){
		for(int j_eta = j*eta_bins; j_eta < (j+1)*eta_bins; j_eta++){
			cellCh_Energy_L5[i][j] += cellCh_Energy[0][i_phi][j_eta];
			cellNu_Energy_L5[i][j] += cellNu_Energy[0][i_phi][j_eta];
		}}
		cell_Energy_L5[i][j] = cellCh_Energy_L5[i][j] + cellNu_Energy_L5[i][j];
	}}
	for(int i=0;i<nPixPhi_new[5];i++){for(int j=0;j<nPixEta_new[5];j++){
		cellCh_Energy_L6[i][j] = 0; cellNu_Energy_L6[i][j] = 0;
		int phi_bins = nPixPhi / nPixPhi_new[5];
		int eta_bins = nPixEta / nPixEta_new[5];
		for(int i_phi = i*phi_bins; i_phi < (i+1)*phi_bins; i_phi++){
		for(int j_eta = j*eta_bins; j_eta < (j+1)*eta_bins; j_eta++){
			cellCh_Energy_L6[i][j] += cellCh_Energy[0][i_phi][j_eta];
			cellNu_Energy_L6[i][j] += cellNu_Energy[0][i_phi][j_eta];
		}}
		cell_Energy_L6[i][j] = cellCh_Energy_L6[i][j] + cellNu_Energy_L6[i][j];
	}}
	
	newtree->Fill();
  }

  outputfile->cd();
  printf("Writes %s\n",outputfile->GetName());
  newtree->Write();
  outputfile->Write();
  outputfile->Close();
  delete outputfile;
  
}

