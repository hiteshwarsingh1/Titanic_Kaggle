// g++ `root-config --cflags` chineu1.C -o chineu1 `root-config --glibs`
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <vector>

#include "TCanvas.h"
#include "TPad.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TTree.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TGraph.h"
#include "TMinuit.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TRandom3.h"
#include "TStyle.h"



using namespace std;
int main(){
  //void chineu(){
  Double_t mass_41[15]={1.0e-8,5.0e-8,1.0e-7,5.0e-7,1.0e-6,5.0e-6,1.0e-5,5.0e-5,1.0e-4,5.0e-4,1.0e-3,5.0e-3,1.0e-2,5.0e-2,0.1};
  Double_t sin2_th14[7]={3.0e-4,13.0e-4,25.0e-4,127.0e-4,257.0e-4,1464.0e-4,0.5};
  
  double NH_pdf, IH_pdf;
  double NH_data, IH_data, NN_pdf, IN_pdf, NI_pdf, II_pdf, cv1;
  double chi_nn_nh, chi_nn_ih, chi_in_nh, chi_in_ih;
  double chi_ni_nh, chi_ni_ih, chi_ii_nh, chi_ii_ih;
  double mnn_nh_s1, mnn_ih_s1, min_nh_s1, min_ih_s1;
  double mni_nh_s1, mni_ih_s1, mii_nh_s1, mii_ih_s1;
  double chi_nn, chi_ni, chi_in, chi_ii;
  
  
  double mnn, mni, min, mii;
  char spectrafile[200], outputfile[200];
  char hist1[200], hist2[200], hist3[200], hist4[200];
   
  gStyle->SetOptStat(0);
  ofstream out,out1;
  out.open("Chisquareval_stat.txt");
  out1.open("chi_interim.txt");
  ifstream in1;
  in1.open("spectra_file.txt");//Opening input file containing name of input root files
  
  TFile *f1 = new TFile("spectrafile.root","RECREATE");//File to cross-check spectra
  TH1D *wo_bkg_nh = new TH1D("wo_bkg_hist_nh","NH_wo/bkg_spectra",560,0.8,12.0);
  //TH1D *w_bkg_nh = new TH1D("w_bkg_hist_nh","NH_w/bkg_spectra",560,0.8,12.0);
  TH1D *wo_bkg_ih = new TH1D("wo_bkg_hist_ih","IH_wo/bkg_spectra",560,0.8,12.0);
  //TH1D *w_bkg_ih = new TH1D("w_bkg_hist_ih","IH_w/bkg_spectra",560,0.8,12.0);
  
  
  while(!in1.eof()){
    memset(spectrafile,0,sizeof(spectrafile));
    in1>>spectrafile;
    out<<spectrafile<<endl;
    if(!in1.eof()){
      TFile *f = new TFile(spectrafile,"READ");
      TH1D *nh = (TH1D*)f->Get("three_fav_nh");
      TH1D *ih = (TH1D*)f->Get("three_fav_ih");

      TH1D *tnh = (TH1D*)f->Get("test_nh");
      TH1D *tih = (TH1D*)f->Get("test_ih");
      int nbin= nh->GetXaxis()->GetNbins();
      cout<<nbin<<endl;

      for(int ij=0; ij<15; ij++){//loop for mass
	for(int jk=0; jk<7; jk++){//loop for sin_theta
	  memset(hist1,0,sizeof(hist1));
	  sprintf(hist1, "Hist_%d_%d_NH_NH", ij, jk);
	  memset(hist2,0,sizeof(hist2));
	  sprintf(hist2, "Hist_%d_%d_IH_NH", ij, jk);
	  memset(hist3,0,sizeof(hist3));
	  sprintf(hist3, "Hist_%d_%d_NH_IH", ij, jk);
	  memset(hist4,0,sizeof(hist4));
	  sprintf(hist4, "Hist_%d_%d_IH_IH", ij, jk);

	  TH1D *nh_nh = (TH1D*)f->Get(hist1);
	  TH1D *ih_nh = (TH1D*)f->Get(hist2);
	  TH1D *nh_ih = (TH1D*)f->Get(hist3);
	  TH1D *ih_ih = (TH1D*)f->Get(hist4);
	  	    
	  ////// Initiation of minimization ////////
	  mnn_nh_s1=mnn_ih_s1=min_nh_s1=min_ih_s1=1000000.0;
	  mni_nh_s1=mni_ih_s1=mii_nh_s1=mii_ih_s1=1000000.0;
	  
	  for(int pq=0; pq<1000; pq++){//loop for statistical fluctutation

	    TRandom3 *gaus = new TRandom3(time(0));
	    
	    chi_nn_nh=chi_nn_ih=0.0;
	    chi_ni_nh=chi_ni_ih=0.0;
	    chi_in_nh=chi_in_ih=0.0;
	    chi_ii_nh=chi_ii_ih=0.0;
	    chi_nn=chi_ni=0.0;
	    chi_in=chi_ii=0.0;

	    for(int kl=1;kl<=nbin; kl++){//loop over histogram bins

	      ///////////// 3+1 flavour histograms///////////////
	      double cnn = nh_nh->GetBinContent(kl);
	      double cin = ih_nh->GetBinContent(kl);
	      double cni = nh_ih->GetBinContent(kl);
	      double cii = ih_ih->GetBinContent(kl);

	      ////////////// 3 flavour histograms ///////////////
	      double dn = nh->GetBinContent(kl);
	      double di = ih->GetBinContent(kl);

	      ////////////// 3 flavour test histograms ///////////
	      double tn = tnh->GetBinContent(kl);
	      double ti = tih->GetBinContent(kl);
	      
	      
	      /////// Initiation ////
	      NH_pdf=IH_pdf=0.0;
	      NN_pdf=IN_pdf=NI_pdf=II_pdf=0.0;
	      NH_data=IH_data=0.0;
	      	      
	      ////////// Definition //////
	      NH_data=dn;//// Experiment data for NH
	      IH_data=di;//// Experiment data for IH

	      NN_pdf=gaus->Gaus(cnn, sqrt(cnn));//// Theory 3+1 model NH-NH
	      IN_pdf=gaus->Gaus(cin, sqrt(cin));//// Theory 3+1 model IH-NH
	      NI_pdf=gaus->Gaus(cni, sqrt(cni));//// Theory 3+1 model NH-IH
	      II_pdf=gaus->Gaus(cii, sqrt(cii));//// Theory 3+1 model IH-IH

	      NH_pdf=gaus->Gaus(tn, sqrt(tn));////Theory 3 neutrino model NH
	      IH_pdf=gaus->Gaus(ti, sqrt(ti));////Theory 3 neutrino model IH

	      //cout<<NN_pdf<<" "<<IN_pdf<<" "<<NI_pdf<<" "<<II_pdf<<endl;

	      //////////////// Filling up histograms for background//////////
	      if(pq==0 && ij==0 && jk==0){
		wo_bkg_nh->SetBinContent(kl,dn);
		wo_bkg_ih->SetBinContent(kl,di);}
		
	      
	      ////////////CHI-square Calculation for 3+1 model //////////////
	      cv1=0.0;
	      ////////////////////// NH-NH ////////////////
	      if(NH_data>0 && NN_pdf>0){
		cv1 = 2.0*((NN_pdf - NH_data) - NH_data*log((NN_pdf/NH_data)));
		chi_nn_nh = chi_nn_nh + cv1;}
	      if(NH_data==0 || NN_pdf==0){chi_nn_nh = chi_nn_nh + 2.0*NN_pdf;}
	      if(NN_pdf<0){chi_nn_nh = chi_nn_nh+0.0;}
	      
	      cv1=0.0;
	      if(IH_data>0 && NN_pdf>0){
		cv1 = 2.0*((NN_pdf - IH_data) - IH_data*log((NN_pdf/IH_data)));
		chi_nn_ih = chi_nn_ih + cv1;}
	      if(IH_data==0 || NN_pdf==0){chi_nn_ih = chi_nn_ih + 2.0*NN_pdf;}
	      if(NN_pdf<0){chi_nn_ih = chi_nn_ih+0.0;}


	      cv1=0.0;
	      ////////////////////// IH-NH /////////////////////
	      if(NH_data>0 && IN_pdf>0){
		cv1 = 2.0*((IN_pdf - NH_data) - NH_data*log((IN_pdf/NH_data)));
		chi_in_nh = chi_in_nh + cv1;}
	      if(NH_data==0 || IN_pdf==0){chi_in_nh = chi_in_nh + 2.0*IN_pdf;}
	      if(IN_pdf<0){chi_in_nh = chi_in_nh + 0.0;}
	      
	      cv1=0.0;
	      if(IH_data>0 && IN_pdf>0){
		cv1 = 2.0*((IN_pdf - IH_data) - IH_data*log((IN_pdf/IH_data)));
		chi_in_ih = chi_in_ih + cv1;}
	      if(IH_data==0 || IN_pdf==0){chi_in_ih = chi_in_ih + 2.0*IN_pdf;}
	      if(IN_pdf<0){chi_in_ih = chi_in_ih + 0.0;}


	      cv1=0.0;
	      ////////////////////////NH-IH//////////////////////////
	      if(NH_data>0 && NI_pdf>0){
		cv1 = 2.0*((NI_pdf - NH_data) - NH_data*log((NI_pdf/NH_data)));
		chi_ni_nh = chi_ni_nh + cv1;}
	      if(NH_data==0 || NI_pdf==0){chi_ni_nh = chi_ni_nh + 2.0*NI_pdf;}
	      if(NI_pdf<0){chi_ni_nh = chi_ni_nh + 0.0;}
	      
	      cv1=0.0;
	      if(IH_data>0 && NI_pdf>0){
		cv1 = 2.0*((NI_pdf - IH_data) - IH_data*log((NI_pdf/IH_data)));
		chi_ni_ih = chi_ni_ih + cv1;}
	      if(IH_data==0 || NI_pdf==0){chi_ni_ih = chi_ni_ih + 2.0*NI_pdf;}
	      if(NI_pdf<0){chi_ni_ih = chi_ni_ih + 0.0;}

	      
	      cv1=0.0;
	      ////////////////////////IH-IH//////////////////////////
	      if(NH_data>0 && II_pdf>0){
		cv1 = 2.0*((II_pdf - NH_data) - NH_data*log((II_pdf/NH_data)));
		chi_ii_nh = chi_ii_nh + cv1;}
	      if(NH_data==0 || II_pdf==0){chi_ii_nh = chi_ii_nh + 2.0*II_pdf;}
	      if(II_pdf<0){chi_ii_nh = chi_ii_nh + 0.0;}
	      cv1=0.0;
	      if(IH_data>0 && II_pdf>0){
		cv1 = 2.0*((II_pdf - IH_data) - IH_data*log((II_pdf/IH_data)));
		chi_ii_ih = chi_ii_ih + cv1;}
	      if(IH_data==0 || II_pdf==0){chi_ii_ih = chi_ii_ih + 2.0*II_pdf;}
	      if(II_pdf<0){chi_ii_ih = chi_ii_ih + 0.0;}


	      
	      /////////// CHI-square calculation for 3 neutrino model /////
	      
	      cv1=0.0;
	      //////////////////////// Experiment_NH /////////////////////////
	      if(NH_pdf>0 && NH_data>0){
		cv1 = 2.0*((NH_pdf-NH_data)-NH_data*log((NH_pdf/NH_data)));
		chi_nn=chi_nn+cv1;}
	      if(NH_data==0 || NH_pdf==0){chi_nn = chi_nn + 2.0*NH_pdf;}
	      if(NH_pdf<0){chi_nn = chi_nn + 0.0;}
		
	      cv1=0.0;
	      if(IH_pdf>0 && NH_data>0){
		cv1 = 2.0*((IH_pdf-NH_data)-NH_data*log((IH_pdf/NH_data)));
		chi_ni=chi_ni+cv1;}
	      if(NH_data==0 || IH_pdf==0){chi_ni = chi_ni + 2.0*IH_pdf;}
	      if(IH_pdf<0){chi_ni = chi_ni + 0.0;}		

	      cv1=0.0;
	      ////////////////////////Experiment_IH /////////////////////////
	      if(NH_pdf>0 && IH_data>0){
		cv1 = 2.0*((NH_pdf-IH_data)-IH_data*log((NH_pdf/IH_data)));
		chi_in=chi_in+cv1;}
	      if(IH_data==0 || NH_pdf==0){chi_in = chi_in + 2.0*NH_pdf;}
	      if(NH_pdf<0){chi_in = chi_in + 0.0;}

	      cv1=0.0;
	      if(IH_pdf>0 && IH_data>0){
		cv1 = 2.0*((IH_pdf-IH_data)-IH_data*log((IH_pdf/IH_data)));
		chi_ii=chi_ii+cv1;}
	      if(IH_data==0 || IH_pdf==0){chi_ii = chi_ii + 2.0*IH_pdf;}
	      if(IH_pdf<0){chi_ii = chi_ii + 0.0;}
	      
	      
		  
	    }//for(int kl=1; kl<=nbin; kl++)
	    delete gaus;
	    gaus = NULL;
	    out1<<chi_nn<<" "<<chi_ni<<" "<<chi_in<<" "<<chi_ii<<endl;
	    /////////// minimization over statistical fluctutation/////
	    if(mnn_nh_s1>=abs(chi_nn_nh-chi_nn)){mnn_nh_s1=chi_nn_nh-chi_nn;}
	    if(mnn_ih_s1>=abs(chi_nn_ih-chi_in)){mnn_ih_s1=chi_nn_ih-chi_in;}

	    if(min_nh_s1>=abs(chi_in_nh-chi_ni)){min_nh_s1=chi_in_nh-chi_ni;}
	    if(min_ih_s1>=abs(chi_in_ih-chi_ii)){min_ih_s1=chi_in_ih-chi_ii;}

	    if(mni_nh_s1>=abs(chi_ni_nh-chi_nn)){mni_nh_s1=chi_ni_nh-chi_nn;}
	    if(mni_ih_s1>=abs(chi_ni_ih-chi_ni)){mni_ih_s1=chi_ni_ih-chi_ni;}

	    if(mii_nh_s1>=abs(chi_ii_nh-chi_in)){mii_nh_s1=chi_ii_nh-chi_in;}
	    if(mii_ih_s1>=abs(chi_ii_ih-chi_ii)){mii_ih_s1=chi_ii_ih-chi_ii;}
	  }//for(int pq=0; pq<1000;pq++)

	  /////////////// minimization over flavour mass hierarchy//////////
	  if(abs(mnn_nh_s1)<=abs(mnn_ih_s1)){mnn=mnn_nh_s1;}
	  else{mnn=mnn_ih_s1;}
	  if(abs(min_nh_s1)<=abs(min_ih_s1)){min=min_nh_s1;}
	  else{min=min_ih_s1;}
	  if(abs(mni_nh_s1)<=abs(mni_ih_s1)){mni=mni_nh_s1;}
	  else{mni=mni_ih_s1;}
	  if(abs(mii_nh_s1)<=abs(mii_ih_s1)){mii=mii_nh_s1;}
	  else{mii=mii_ih_s1;}

	  out<<mass_41[ij]<<" "<<sin2_th14[jk]<<" "<<mnn<<endl;
	  // out<<mass_41[ij]<<" "<<sin2_th14[jk]<<" "<<chi_in_nh<<" "<<min<<endl;
	  // out<<mass_41[ij]<<" "<<sin2_th14[jk]<<" "<<chi_ni_nh<<" "<<mni<<endl;
	  // out<<mass_41[ij]<<" "<<sin2_th14[jk]<<" "<<chi_ii_nh<<" "<<mii<<endl;
	  
	}//for(int jk=0; jk<7; jk++)  
      }//for(int ij=0; ij<15; ij++)
      f1->Write();
      f1->Close();
      f->Close();

    }//if(!in1.eof())
  }//while(!in1.eof())
  out.close();
  out1.close(); 
  in1.close();
  return 0;
}
