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

const Double_t GeV=1e3;
const Double_t MeV=1e0;
const Int_t nLayers = 6, nPixEta = 64, nPixPhi = 64; // eta = x and phi = yes
const Double_t calorSizeXY = 144 * 16; // in mm
TString file_ref=Form("_%dD%dx%d",nLayers,nPixPhi,nPixEta); 
const Double_t cell_dEta = calorSizeXY / nPixEta;
const Double_t cell_dPhi = calorSizeXY / nPixPhi;
const int kNaxPrimaries = 2;

void CreateInputMatrix(TString infile = "PFlowNtuple.root")
{

  float cell_Energy[nLayers][nPixPhi][nPixEta];
  float cellCh_Energy[nLayers][nPixPhi][nPixEta];
  float cellNu_Energy[nLayers][nPixPhi][nPixEta];
  vector<double>  *cell_e;
  vector<double>  *cell_x;
  vector<double>  *cell_y;
  vector<double>  *cell_dx;
  vector<double>  *cell_dy;
  vector<int>  *cell_l;
  vector<double>  *cell_chfrac;
  int n_primaries;
  int  par_ID[kNaxPrimaries];
  float  par_dx0[kNaxPrimaries];
  float  par_dy0[kNaxPrimaries];
  float  par_dx1[kNaxPrimaries];
  float  par_dy1[kNaxPrimaries];
  vector<double>  *particle_x;
  vector<double>  *particle_y;
  vector<double>  *particle_px;
  vector<double>  *particle_py;
  vector<double>  *particle_e;
  vector<int>  *particle_pdgId;
  
  TString directory = "/mnt/Lustre/agrp/pitt/ML/Simulation/Analysis";
  TChain * oldtree = new TChain("physics");   
  TString outfile = infile; outfile.ReplaceAll(".root",file_ref+".root");
  cout << "Will create "<<directory <<"/"<<outfile<< endl;
  oldtree->Add(directory+"/"+infile);

  oldtree->SetBranchAddress("cell_e", &cell_e);
  oldtree->SetBranchAddress("cell_x", &cell_x);
  oldtree->SetBranchAddress("cell_dx", &cell_dx);
  oldtree->SetBranchAddress("cell_y", &cell_y);
  oldtree->SetBranchAddress("cell_dy", &cell_dy);
  oldtree->SetBranchAddress("cell_l", &cell_l);
  oldtree->SetBranchAddress("cell_chfrac", &cell_chfrac);
  oldtree->SetBranchAddress("particle_x", &particle_x);
  oldtree->SetBranchAddress("particle_y", &particle_y);
  oldtree->SetBranchAddress("particle_px", &particle_px);
  oldtree->SetBranchAddress("particle_py", &particle_py);
  oldtree->SetBranchAddress("particle_e", &particle_e);
  oldtree->SetBranchAddress("particle_pdgId", &particle_pdgId);


  TFile * outputfile = new TFile(directory+"/"+outfile,"recreate");
  TTree *newtree = new TTree("EventTree",Form("Detector images of dim. of %d x %d x %d",nLayers,nPixPhi,nPixEta));
  newtree->Branch("cell_Energy",cell_Energy,Form("cell_Energy[%d][%d][%d]/F",nLayers,nPixPhi,nPixEta));
  newtree->Branch("cellCh_Energy",cellCh_Energy,Form("cellCh_Energy[%d][%d][%d]/F",nLayers,nPixPhi,nPixEta));
  newtree->Branch("cellNu_Energy",cellNu_Energy,Form("cellNu_Energy[%d][%d][%d]/F",nLayers,nPixPhi,nPixEta));
  newtree->Branch("n_primaries", &n_primaries,"n_primaries/I");
  newtree->Branch("par_ID", par_ID,"par_ID[n_primaries]/I");
  newtree->Branch("par_dx0", par_dx0,"par_dx0[n_primaries]/F");
  newtree->Branch("par_dy0", par_dy0,"par_dy0[n_primaries]/F");
  newtree->Branch("par_dx1", par_dx1,"par_dx1[n_primaries]/F");
  newtree->Branch("par_dy1", par_dy1,"par_dy1[n_primaries]/F");


  Int_t nentries = (Int_t)oldtree->GetEntries();
  printf("%d entries to be processed\n",nentries);
  for(Int_t entry=0;entry<nentries;entry++) {
	oldtree->GetEntry(entry);
	if ((entry%(nentries/10))==0) printf("%f complete\n",100*(double)entry/nentries);
	
	// Reset variables in the new tree:
	for(int i=0;i<nPixPhi;i++){for(int j=0;j<nPixEta;j++){for(int k=0;k<nLayers;k++){
		cellCh_Energy[k][i][j] = 0; cellNu_Energy[k][i][j] = 0;
	}}}
	
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
	n_primaries = int(particle_pdgId->size());
	for (int i=0;i<n_primaries;i++){
		par_ID[i] = (int(particle_pdgId->at(i)));
		par_dx0[i] = (float)particle_x->at(i);
		par_dy0[i] = (float)particle_y->at(i);
		par_dx1[i] = (float)(particle_x->at(i) + particle_px->at(i)*1500./particle_e->at(i));
		par_dy1[i] = (float)(particle_y->at(i) + particle_py->at(i)*1500./particle_e->at(i));
	}
	newtree->Fill();
  }

  outputfile->cd();
  printf("Writes %s\n",outputfile->GetName());
  newtree->Write();
  outputfile->Write();
  outputfile->Close();
  delete outputfile;
  
}

