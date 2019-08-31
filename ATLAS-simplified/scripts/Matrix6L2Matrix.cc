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

// Old dimentions of layers
TString file_ref="6D64x64"; 
const Int_t nPixEta_new[nLayers] = {64, 64, 32, 16, 16, 8};
const Int_t nPixPhi_new[nLayers] = {32, 64, 32, 16, 16, 8};

// New output dimentions
const Int_t nPixEta = 64, nPixPhi = 64;

// common detector parameters
const Double_t GeV=1e3;
const Double_t MeV=1e0;
const Double_t calorSizeXY = 144 * 16; // in mm
const Double_t cell_dEta = calorSizeXY / nPixEta;
const Double_t cell_dPhi = calorSizeXY / nPixPhi;
const int kNaxPrimaries = 2;

void Matrix6L2Matrix(TString infile = "events_6D64x64_ATLAS_resolution.root")
{

  
// old branches
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

// new branches
  float cell_Energy[nLayers][nPixPhi][nPixEta];
  float cellCh_Energy[nLayers][nPixPhi][nPixEta];
  float cellNu_Energy[nLayers][nPixPhi][nPixEta];

 
  TString directory = getenv("PWD");
  TChain * oldtree = new TChain("EventTree");   
  TString outfile = infile; outfile.ReplaceAll(".root","_"+file_ref+".root");
  cout << "Will create "<<directory <<"/"<<outfile<< endl;
  oldtree->Add(directory+"/"+infile);

  oldtree->SetBranchAddress("cell_Energy_L1",cell_Energy_L1);
  oldtree->SetBranchAddress("cell_Energy_L2",cell_Energy_L2);
  oldtree->SetBranchAddress("cell_Energy_L3",cell_Energy_L3);
  oldtree->SetBranchAddress("cell_Energy_L4",cell_Energy_L4);
  oldtree->SetBranchAddress("cell_Energy_L5",cell_Energy_L5);
  oldtree->SetBranchAddress("cell_Energy_L6",cell_Energy_L6);
  
  oldtree->SetBranchAddress("cellCh_Energy_L1",cellCh_Energy_L1);
  oldtree->SetBranchAddress("cellCh_Energy_L2",cellCh_Energy_L2);
  oldtree->SetBranchAddress("cellCh_Energy_L3",cellCh_Energy_L3);
  oldtree->SetBranchAddress("cellCh_Energy_L4",cellCh_Energy_L4);
  oldtree->SetBranchAddress("cellCh_Energy_L5",cellCh_Energy_L5);
  oldtree->SetBranchAddress("cellCh_Energy_L6",cellCh_Energy_L6);
  
  oldtree->SetBranchAddress("cellNu_Energy_L1",cellNu_Energy_L1);
  oldtree->SetBranchAddress("cellNu_Energy_L2",cellNu_Energy_L2);
  oldtree->SetBranchAddress("cellNu_Energy_L3",cellNu_Energy_L3);
  oldtree->SetBranchAddress("cellNu_Energy_L4",cellNu_Energy_L4);
  oldtree->SetBranchAddress("cellNu_Energy_L5",cellNu_Energy_L5);
  oldtree->SetBranchAddress("cellNu_Energy_L6",cellNu_Energy_L6);

  TFile * outputfile = new TFile(directory+"/"+outfile,"recreate");
  //oldtree->LoadTree(0);
  //TTree *newtree = oldtree->GetTree()->CloneTree(0);
  //newtree->SetBranchStatus("*",0);
  //newtree->SetBranchStatus("n_primaries",1);
  //newtree->SetBranchStatus("par_*",1);
  TTree *newtree = new TTree("EventTree","Detector upscaling to 6x64x64 with democratic approach");

  newtree->Branch("cell_Energy",cell_Energy,Form("cell_Energy[%d][%d][%d]/F",nLayers,nPixPhi,nPixEta));
  newtree->Branch("cellCh_Energy",cellCh_Energy,Form("cellCh_Energy[%d][%d][%d]/F",nLayers,nPixPhi,nPixEta));
  newtree->Branch("cellNu_Energy",cellNu_Energy,Form("cellNu_Energy[%d][%d][%d]/F",nLayers,nPixPhi,nPixEta));


  Int_t nentries = (Int_t)oldtree->GetEntries();
  printf("%d entries to be processed\n",nentries);
  for(Int_t entry=0;entry<nentries;entry++) {
	oldtree->GetEntry(entry);
	if ((entry%(nentries/10))==0) printf("%f complete\n",100*(double)entry/nentries);
	
	// Fill the variables using democratic upscaling
	
	// NOTE: output resolution is always lower than the input!!!
	for(int i=0;i<nPixPhi;i++){for(int j=0;j<nPixEta;j++){
		int phi_bins = nPixPhi / nPixPhi_new[0];
		int eta_bins = nPixEta / nPixEta_new[0];
		
		cellCh_Energy[0][i][j] = cellCh_Energy_L1[i/phi_bins][j/eta_bins] / float(eta_bins*phi_bins);
		cellNu_Energy[0][i][j] = cellNu_Energy_L1[i/phi_bins][j/eta_bins] / float(eta_bins*phi_bins);
		cell_Energy[0][i][j] = cellCh_Energy[0][i][j] + cellNu_Energy[0][i][j];
	}}

	for(int i=0;i<nPixPhi;i++){for(int j=0;j<nPixEta;j++){
		int phi_bins = nPixPhi / nPixPhi_new[1];
		int eta_bins = nPixEta / nPixEta_new[1];
		
		cellCh_Energy[1][i][j] = cellCh_Energy_L2[i/phi_bins][j/eta_bins] / float(eta_bins*phi_bins);
		cellNu_Energy[1][i][j] = cellNu_Energy_L2[i/phi_bins][j/eta_bins] / float(eta_bins*phi_bins);
		cell_Energy[1][i][j] = cellCh_Energy[1][i][j] + cellNu_Energy[1][i][j];
	}}
	
	for(int i=0;i<nPixPhi;i++){for(int j=0;j<nPixEta;j++){
		int phi_bins = nPixPhi / nPixPhi_new[2];
		int eta_bins = nPixEta / nPixEta_new[2];
		
		cellCh_Energy[2][i][j] = cellCh_Energy_L3[i/phi_bins][j/eta_bins] / float(eta_bins*phi_bins);
		cellNu_Energy[2][i][j] = cellNu_Energy_L3[i/phi_bins][j/eta_bins] / float(eta_bins*phi_bins);
		cell_Energy[2][i][j] = cellCh_Energy[2][i][j] + cellNu_Energy[2][i][j];
	}}

	for(int i=0;i<nPixPhi;i++){for(int j=0;j<nPixEta;j++){
		int phi_bins = nPixPhi / nPixPhi_new[3];
		int eta_bins = nPixEta / nPixEta_new[3];
		
		cellCh_Energy[3][i][j] = cellCh_Energy_L4[i/phi_bins][j/eta_bins] / float(eta_bins*phi_bins);
		cellNu_Energy[3][i][j] = cellNu_Energy_L4[i/phi_bins][j/eta_bins] / float(eta_bins*phi_bins);
		cell_Energy[3][i][j] = cellCh_Energy[3][i][j] + cellNu_Energy[3][i][j];
	}}
	
	for(int i=0;i<nPixPhi;i++){for(int j=0;j<nPixEta;j++){
		int phi_bins = nPixPhi / nPixPhi_new[4];
		int eta_bins = nPixEta / nPixEta_new[4];
		
		cellCh_Energy[4][i][j] = cellCh_Energy_L5[i/phi_bins][j/eta_bins] / float(eta_bins*phi_bins);
		cellNu_Energy[4][i][j] = cellNu_Energy_L5[i/phi_bins][j/eta_bins] / float(eta_bins*phi_bins);
		cell_Energy[4][i][j] = cellCh_Energy[4][i][j] + cellNu_Energy[4][i][j];
	}}
	
	for(int i=0;i<nPixPhi;i++){for(int j=0;j<nPixEta;j++){
		int phi_bins = nPixPhi / nPixPhi_new[5];
		int eta_bins = nPixEta / nPixEta_new[5];
		
		cellCh_Energy[5][i][j] = cellCh_Energy_L6[i/phi_bins][j/eta_bins] / float(eta_bins*phi_bins);
		cellNu_Energy[5][i][j] = cellNu_Energy_L6[i/phi_bins][j/eta_bins] / float(eta_bins*phi_bins);
		cell_Energy[5][i][j] = cellCh_Energy[5][i][j] + cellNu_Energy[5][i][j];
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

