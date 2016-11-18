#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include <boost/algorithm/string.hpp>

#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TStyle.h"

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TNtuple.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TF1.h"
#include "TString.h"
#include "TLatex.h"
#include "TMath.h"
#include "TROOT.h"

#include "HGCSSSimHit.hh"
#include "HGCSSSamplingSection.hh"

//===============================================================================
//Declaration here of any helpful function definition after main
std::string IntToString (int number1);

int main(){//main

  //We work with config: v_HGCALEE_v5=12
  const int numberofconfigs = 1; // 11
  std::vector<int> configs;
  configs.clear();
  configs.push_back(12);//v_HGCALEE_v5
  //======================================================================================
  //Step1 only or not
  bool step1 = true; 
  bool step2 = true; 
  //======================================================================================
  //The root files with the individual plots for all test beam configs 
  const int numberoffiles = 1;
  //======================================================================================
  //3 layers times 4 Si Pads
  // const int numberofpads = 12; 
  //======================================================================================
  //Particle type and energy
  TString particle = "e-"; //mu+
  const int numberofenergies = 12; //15 30 50 80 100 120 150 180 200 250 300 500
  //const int numberofenergies = 1; // 15 30 50 80 100 150
 // int particleenergy = 15;
  std::vector<int> particleenergies;
  particleenergies.clear();
  particleenergies.push_back(15);
  particleenergies.push_back(30);
  particleenergies.push_back(50);
  particleenergies.push_back(80);
  particleenergies.push_back(100);
  particleenergies.push_back(120);
  particleenergies.push_back(150);
  particleenergies.push_back(180);
  particleenergies.push_back(200);
  particleenergies.push_back(250);
  particleenergies.push_back(300);
  particleenergies.push_back(500);

  //======================================================================================
  //Read all files and tree
  TString filename[numberofconfigs][numberofenergies][numberoffiles];
  for (int k=0; k<numberofconfigs; k++){
    for (int l=0; l<numberofenergies; l++){

      for (int j=0; j<numberoffiles; j++){
      TString fp = "PFCal_" + IntToString(configs[k]);	
      TString sp = "_" + particle;
      TString tp = "_" + IntToString(particleenergies[l]);
      // TString fop = "GeV_" + IntToString(j);
      TString fop = "GeV";
      filename[k][l][j] = fp + sp + tp + fop;
      //std::cout << filename[k] << std::endl;
      }
    }
  }

  //The tree in the files that we want to get
  // TTree *lTree[numberofconfigs][numberofenergies]; 
  // ...and chains in the cases of multiple files
  TChain *lChain[numberofconfigs][numberofenergies];

  // TFile* files[numberofconfigs][numberofenergies][numberoffiles];
  for (int k=0; k<numberofconfigs; k++){
    for (int l=0; l<numberofenergies; l++){
      lChain[k][l] = new TChain("HGCSSTree"); 
      for (int j=0; j<numberoffiles; j++){
	// lChain[k][l]->Add("root://eoscms//eos/cms/store/group/phys_b2g/apsallid/PFCal/DetectorConfigurations/Config_12/"+filename[k][l][j]+".root");
	lChain[k][l]->Add("/afs/cern.ch/work/a/apsallid/CMS/Geant4/PFCal/PFCalEE/analysis/data/"+filename[k][l][j]+".root");
	// lChain[k][l]->Add("/tmp/apsallid/"+filename[k][l][j]+".root");

      }
    }
  }
  //Set the branches here
  std::vector<HGCSSSimHit> * simhitvec = 0;
  std::vector<HGCSSSamplingSection> * ssvec = 0;
  for (int k=0; k<numberofconfigs; k++){
    for (int l=0; l<numberofenergies; l++){
      lChain[k][l]->SetBranchAddress("HGCSSSimHitVec",&simhitvec);
      lChain[k][l]->SetBranchAddress("HGCSSSamplingSectionVec",&ssvec);  
    }
  }
  //======================================================================================
  // Histos
  //dE vs MIPs
  std::map< int, std::vector<TH2F*> > h_dEvsMIPs;
  h_dEvsMIPs.clear();
  //Esilicon vs (Epassive+Esilicon) for the sampling fraction
  std::map< int, std::vector<TH2F*> > h_MIPSvsdEplusMIPs;
  h_MIPSvsdEplusMIPs.clear();
  //Shower depth
  std::map< int, std::vector<TH1F*> > h_Showerprofile;
  h_Showerprofile.clear();
  //Shower depth true
  std::vector<TH1F*> h_ShowerDepth_reco;
  h_ShowerDepth_reco.clear();
  //Total sampling fraction
  std::vector<TH1F*> h_TotalSamplingFraction;
  h_TotalSamplingFraction.clear();
  //Energy loss
  std::map< int, std::vector<TH1F*> > h_loss;
  h_loss.clear();
  //SF for layers defined as: SF = (dE + Esilicon)/dE = (Epassive+Eactive)/Eactive
  std::map< int, std::vector<TH1F*> > h_SF;
  h_SF.clear();
  //======================================================================================
  // For the fit
  std::map< int, std::vector<TGraphErrors*> > dEvsMIPs_calib;
  dEvsMIPs_calib.clear();
  std::map< int, std::vector<TGraphErrors*> > MIPSvsdEplusMIPs_calib;
  MIPSvsdEplusMIPs_calib.clear();
  //Average energy lost as a function of the transversed thickness
  std::vector<TGraphErrors*> dEvsXos;
  dEvsXos.clear();
  //======================================================================================
  //Weights
  std::map< int, std::vector<double> > wXos;
  wXos.clear();
  std::map< int, std::vector<double> > wdEdx;
  wdEdx.clear();
  std::map< int, std::vector<double> > wlamdba;
  wlamdba.clear();
  std::map< int, std::vector<double> > wdevsmips;
  wdevsmips.clear();
  std::map< int, std::vector<double> > wdevsmipserror;
  wdevsmipserror.clear();
  std::map< int, std::vector<double> > wdevsmips2;
  wdevsmips2.clear();
  std::map< int, std::vector<double> > wdevsmipserror2;
  wdevsmipserror2.clear();

  lChain[0][0]->GetEntry(0);
  TGraphErrors *calib = new TGraphErrors();
  for(unsigned iL(0); iL<(*ssvec).size(); iL++){
    for (Int_t k=0; k<numberofenergies; k++){
      std::string auxnam1 = "h_dEvsMIPs_" + IntToString((int) iL);
      std::string auxnam1_sf = "h_SF_" + IntToString((int) iL);
      std::string auxnam1_sf2 = "h_MIPSvsdEplusMIPs_" + IntToString((int) iL);
      std::string auxnam1_tgr = "dEvsMIPs_" + IntToString((int) iL);
      std::string auxnam1_tgr1 = "dEvsXos_" + IntToString((int) iL);
      std::string auxnam1_tgr2 = "MIPSvsdEplusMIPs_" + IntToString((int) iL);
      std::string auxnam1_sh = "Showerprofile_" + IntToString((int) iL);
      std::string auxnam2 = "_" + IntToString(particleenergies[k]);
      std::string auxnam = auxnam1 + auxnam2;
      std::string auxnam_sf = auxnam1_sf + auxnam2;
      std::string auxnam_sf2 = auxnam1_sf2 + auxnam2;
      std::string auxnam_tgr = auxnam1_tgr + auxnam2;
      std::string auxnam_tgr1 = auxnam1_tgr1 + auxnam2;
      std::string auxnam_tgr2 = auxnam1_tgr2 + auxnam2;
      std::string auxnam_sh = auxnam1_sh + auxnam2;
      h_dEvsMIPs[(int) iL].push_back(new TH2F( auxnam.c_str(),";E (MIPs);SimHits;E (GeV);SimHits",1000,0.,1000.,520,0.,520.));
      h_MIPSvsdEplusMIPs[(int) iL].push_back(new TH2F( auxnam_sf2.c_str(),";E_{active} + E_{passive};E_{active};E (GeV);SimHits",520,0.,520.,1000,0.,1000.));
      h_Showerprofile[(int) iL].push_back(new TH1F( auxnam_sh.c_str(),";E (MIPs);Events",1000,0.,1000.));
      h_loss[(int) iL].push_back(new TH1F( auxnam_tgr1.c_str(),";dE (MIPs);Events",10000,0.,10000.));
      h_SF[(int) iL].push_back(new TH1F( auxnam_sf.c_str(),";Scale factors;Events",1000,0.,0.1));

      dEvsMIPs_calib[(int) iL].push_back( (TGraphErrors *) calib->Clone( auxnam_tgr.c_str() )  );
      MIPSvsdEplusMIPs_calib[(int) iL].push_back( (TGraphErrors *) calib->Clone( auxnam_tgr2.c_str() )  );
      //don't forget 0 is before calo, 31 is for leakage
      if (iL==0 || iL==1){
	wXos[(int) iL].push_back( 1. );wdEdx[(int) iL].push_back( 1. );wlamdba[(int) iL].push_back( 1. );
      } else {
	wXos[(int) iL].push_back( (*ssvec)[iL].volX0trans()/( (*ssvec)[0].volX0trans() + (*ssvec)[1].volX0trans() )  );
	wdEdx[(int) iL].push_back( (*ssvec)[iL].voldEdx()/( (*ssvec)[0].voldEdx() + (*ssvec)[1].voldEdx() ) );
	wlamdba[(int) iL].push_back( (*ssvec)[iL].volLambdatrans()/( (*ssvec)[0].volLambdatrans() + (*ssvec)[1].volLambdatrans() ) );
      }

    }
  }
  
  for (Int_t k=0; k<numberofenergies; k++){
    std::string auxnam_tgr1 = "dEvsXos_" + IntToString(particleenergies[k]);
    h_ShowerDepth_reco.push_back(new TH1F( ("ShowerDepth_reco_"+ IntToString(particleenergies[k])).c_str(),";Shower depth (1/X_{0});Events",300,0.,30.) );
    h_TotalSamplingFraction.push_back(new TH1F( ("TotalSF_"+ IntToString(particleenergies[k])).c_str(),";Total sampling fraction;Events",100,0.,0.1) );
    dEvsXos.push_back( (TGraphErrors *) calib->Clone( auxnam_tgr1.c_str() )  );
  }



 
  //======================================================================================
  //The files that we will store the results of this analysis
  TString res[numberofconfigs];
  for (int k=0; k<numberofconfigs; k++){
    TString fp = "SiWEcal_" + IntToString(configs[k]);	
    res[k] = fp + "_final.root";
    //std::cout << res[k] << std::endl;
  }
  
  TFile* results[numberofconfigs];

  if (step1){

    for (int k=0; k<numberofconfigs; k++){
      results[k]= new TFile(res[k],"recreate");
      std::cout << "Results file " << res[k] << std::endl;

      //======================================================================================
      //Loop on energies
      for (int l=0; l<numberofenergies; l++){
	//======================================================================================
	// Loop on entries
	for (Int_t ievt=0; ievt<lChain[k][l]->GetEntries(); ievt++) {
	  // 	// if (ievt != 0) {continue;}
	  lChain[k][l]->GetEntry(ievt);
	  if (ievt%(1000)==0) std::cout << "Entry " << ievt << std::endl;
	  if (ievt==10){break;}
	  // std::cout << "Entry " << ievt << std::endl;
	  // double lossbefcalo = 0.;
	  double totalscalefactor = 0.;
	  double lossinsensors = 0.;
	  double totalloss = 0.;
	  //For the xotransfered add the before calo xo's
	  std::vector<double> xotrans;
	  xotrans.clear();
	  double currentxotrans = 0.;
	  int layermax = -1.;
	  double currmaxenergy = 0.;

	  //Loop on sampling sections
	  for(unsigned iL(0); iL<(*ssvec).size(); iL++){
	    // HGCSSSamplingSection lSamSec = (*ssvec)[iL];
	    // double energyabsorber = lSamSec.absorberE();

	    // std::cout << "Sampling Section " << iL << " energyabsorber " << energyabsorber  << std::endl; 
	  
	    //dE vs MIPs 
	    //For the energy lost start with material just after the previous Silicon layer (e.g. layer1) and 
	    //stop just before the specific layer (layer 2). Be careful with layer 0 before calorimeter. 
	    //For the conversion we use the simulation value that we calculated: 0.081 MeV/MIP
	    h_dEvsMIPs[(int) (*ssvec)[iL].volNb()][l]->Fill( ( (*ssvec)[iL].measuredE()/(0.081) ) , (*ssvec)[iL].absorberE()/1000. ); //x is MIPs, y is dE in GeV. 
	    Int_t np=dEvsMIPs_calib[(int) iL][l]->GetN();
	    // std::cout << "np " << np << std::endl;
	    dEvsMIPs_calib[(int) iL][l]->SetPoint(np, ( (*ssvec)[iL].measuredE()/(0.081) )  , (*ssvec)[iL].absorberE()/1000.  );
	    // dEvsMIPs_calib[(int) iL][l]->SetPointError();

	    //Esilicon vs (Epassive+Esilicon) for the sampling fraction
	    h_MIPSvsdEplusMIPs[(int) iL][l]->Fill( ( (*ssvec)[iL].measuredE() + (*ssvec)[iL].absorberE() )/(1000.) , (*ssvec)[iL].measuredE()/1000. ); //x is GeV, y is GeV. 


	    Int_t np2=MIPSvsdEplusMIPs_calib[(int) iL][l]->GetN();
	    // std::cout << "np " << np << std::endl;
	    MIPSvsdEplusMIPs_calib[(int) iL][l]->SetPoint(np2, ( (*ssvec)[iL].measuredE() + (*ssvec)[iL].absorberE() )/(1000.)  , (*ssvec)[iL].measuredE()/1000.  );
	    // MIPSvsdEplusMIPs_calib[(int) iL][l]->SetPointError();
	    


	    h_Showerprofile[(int) iL][l]->Fill( ( (*ssvec)[iL].measuredE()/(0.081) ) );//in MIPs
	    h_loss[(int) iL][l]->Fill( (*ssvec)[iL].absorberE()/(0.081) );//in MIPs
	    
	    if(iL!=0){
	      //Count from layer 1 as start of calo and not include before calorimeter.
	      currentxotrans = currentxotrans + (*ssvec)[iL].volX0trans();
	      xotrans.push_back(currentxotrans);
	      if ( currmaxenergy < (*ssvec)[iL].measuredE() ) { 
		layermax = (int) iL;
		currmaxenergy = (*ssvec)[iL].measuredE();	      
	      }
	    }
	    //double w = (*ssvec)[iL].volX0trans()/(*ssvec)[0].volX0trans();

	    if ( (iL != 0) && (iL != 31)){
	      totalloss = totalloss + (*ssvec)[iL].measuredE() + (*ssvec)[iL].absorberE(); //in MeV
	      lossinsensors = lossinsensors + (*ssvec)[iL].measuredE(); //in MeV
	    }
	    
	    //SF = Esilicon/(Epassive+Esilicon)
	    // totalscalefactor = totalscalefactor + (*ssvec)[iL].measuredE() / ( (*ssvec)[iL].measuredE() + (*ssvec)[iL].absorberE() );//in MeV
	    // scalefactor = ( ( (lossbefcalo * 1000.) + (*ssvec)[iL].measuredE() ) / (lossbefcalo * 1000.) ); //both in MeV here
	    // h_SF[(int) iL][l]->Fill( scalefactor ); 
	    // std::cout << scalefactor << std::endl;
	    // }	  


	  } //loop on sampling sections

	  //For the total sampling fraction per event
	  //Sum_1^30(Esilicon_i) / Sum_1^30(Epassive+Esilicon)
	  totalscalefactor = lossinsensors/totalloss;
	  h_TotalSamplingFraction[l]->Fill( totalscalefactor );

	  //layer max 9 is the 8 element (0,...8) of vector
	  h_ShowerDepth_reco[l]->Fill( xotrans[layermax-1]  );
	  // std::cout << "layermax " << layermax << " xo " << xotrans[layermax-1] << std::endl;

	} //  Loop on entries

	//Loop again on sampling sections for the average plots
	//0 is before calo, -1 for the leakage
	for(unsigned iL(1); iL<(*ssvec).size()-1; iL++){
	  //For the average energy lost as a function of the transversed thickness
	  Int_t nnp=dEvsXos[l]->GetN();
	  dEvsXos[l]->SetPoint( nnp, (*ssvec)[iL].volX0trans() , h_loss[(int) iL][l]->GetMean() ); // y is in MIPs, x is in radiation length
	  dEvsXos[l]->SetPointError( nnp, 0 , h_loss[(int) iL][l]->GetRMS() ); // y is in MIPs, x is in radiation length
	}
	


	//======================================================================================
	//Write histos 
	dEvsXos[l]->Write();
	h_TotalSamplingFraction[l]->Write();
	h_ShowerDepth_reco[l]->Write();
	lChain[0][0]->GetEntry(0);
	for(unsigned iL(0); iL<(*ssvec).size(); iL++){
	  // h_dEvsMIPs[(int) iL][l]->Write();
	  h_MIPSvsdEplusMIPs[(int) iL][l]->Write();
	  // h_SF[(int) iL][l]->Write();
	  dEvsMIPs_calib[(int) iL][l]->Write();
	  MIPSvsdEplusMIPs_calib[(int) iL][l]->Write();
	}


	//======================================================================================
	//We should here clear the histograms because we want them empty for the next file. 
	//Reset histos
	for(unsigned iL(0); iL<(*ssvec).size(); iL++){
	  h_dEvsMIPs[(int) iL][l]->Reset();
	  h_MIPSvsdEplusMIPs[(int) iL][l]->Reset();
	  h_SF[(int) iL][l]->Reset();
	  // dEvsMIPs_calib[(int) iL][l]->Reset();
	}

	h_TotalSamplingFraction[l]->Reset();
	h_ShowerDepth_reco[l]->Reset();
       
      }//Loop on energies
   
      results[k]->Close();

    } //Loop on files

  }// if to avoid running again on input files

  //======================================================================================
  //Make plots
  TCanvas *c1[(*ssvec).size()-1];
  TCanvas *c2[(*ssvec).size()-1];
  TCanvas *c4[(*ssvec).size()-1][numberofenergies];
  TCanvas *c9[(*ssvec).size()-1];
  TCanvas *c10[(*ssvec).size()-1][numberofenergies];
  for(unsigned iL(0); iL<(*ssvec).size()-1; iL++){
    c1[(int) iL] = new TCanvas(("c1_"+ IntToString((int) iL)).c_str(), "  ");
    c9[(int) iL] = new TCanvas(("c9_"+ IntToString((int) iL)).c_str(), "  ");
    c2[(int) iL] = new TCanvas(("c2_"+ IntToString((int) iL)).c_str(), "  ");
    for (int l=0; l<numberofenergies; l++){
      std::string can_title = "c4_"+ IntToString((int) iL);
      std::string can_title2 = "c10_"+ IntToString((int) iL);
      c4[(int) iL][l] = new TCanvas( (can_title+"_"+IntToString(l)) .c_str(), "  ");
      c10[(int) iL][l] = new TCanvas( (can_title2+"_"+IntToString(l)) .c_str(), "  ");
    }
  }
  TCanvas *c3 = new TCanvas("c3", "  ");
  TCanvas *c5 = new TCanvas("c5", "  ");
  TCanvas *c6 = new TCanvas("c6", "  ");
  TCanvas *c7 = new TCanvas("c7", "  ");
  TCanvas *c8 = new TCanvas("c8", "  ");
  TCanvas *c11 = new TCanvas("c11", "  ");
  TCanvas *c12 = new TCanvas("c12", "  ");
  //For the legend
  // TLegend* leg[numberoffiles];
  TLegend* leg[(*ssvec).size()-1];
  for(unsigned iL(0); iL<(*ssvec).size()-1; iL++){
    leg[(int) iL] = new TLegend(0.5, 0.7, 0.8, 0.9);
    leg[(int) iL]->SetHeader("Energy");
    leg[(int) iL]->SetFillColor(17);
  }
  TLegend* leg_sf2[(*ssvec).size()-1];
  for(unsigned iL(0); iL<(*ssvec).size()-1; iL++){
    leg_sf2[(int) iL] = new TLegend(0.5, 0.7, 0.8, 0.9);
    leg_sf2[(int) iL]->SetHeader("Energy");
    leg_sf2[(int) iL]->SetFillColor(17);
  }
  TLegend* leg_sf[(*ssvec).size()-1];
  for(unsigned iL(0); iL<(*ssvec).size()-1; iL++){
    leg_sf[(int) iL] = new TLegend(0.5, 0.7, 0.8, 0.9);
    leg_sf[(int) iL]->SetHeader("Energy");
    leg_sf[(int) iL]->SetFillColor(17);
  }
  TLegend* leg_devsmipswithfit[(*ssvec).size()-1];
  for(unsigned iL(0); iL<(*ssvec).size()-1; iL++){
    leg_devsmipswithfit[(int) iL] = new TLegend(0.5, 0.7, 0.8, 0.9);
    leg_devsmipswithfit[(int) iL]->SetHeader("Energy");
    leg_devsmipswithfit[(int) iL]->SetFillColor(17);
  }
  TLegend* leg_fit = new TLegend(0.5, 0.7, 0.8, 0.9);
  leg_fit->SetHeader("Layers");
  leg_fit->SetFillColor(17);
  TLegend* leg_fit_c = new TLegend(0.5, 0.7, 0.8, 0.9);
  leg_fit_c->SetHeader("Layers");
  leg_fit_c->SetFillColor(17);
  TLegend* leg_fit2 = new TLegend(0.5, 0.7, 0.8, 0.9);
  leg_fit2->SetHeader("Layers");
  leg_fit2->SetFillColor(17);
  TLegend* leg_fit2_c = new TLegend(0.5, 0.7, 0.8, 0.9);
  leg_fit2_c->SetHeader("Layers");
  leg_fit2_c->SetFillColor(17);
  TLegend* leg_devsxos = new TLegend(0.5, 0.7, 0.8, 0.9);
  leg_devsxos->SetHeader("Layers");
  leg_devsxos->SetFillColor(17);
  TLegend* leg_w = new TLegend(0.5, 0.7, 0.8, 0.9);
  leg_w->SetHeader("Layers");
  leg_w->SetFillColor(17);
  TLegend* leg_allw = new TLegend(0.5, 0.7, 0.8, 0.9);
  leg_allw->SetHeader("Layers");
  leg_allw->SetFillColor(17);

  TString procNa[numberofenergies];
  // procNa[0] = "15 GeV";
  // procNa[1] = "30 GeV";
  // // procNa[2] = "50 GeV";
  // procNa[2] = "80 GeV";
  // procNa[3] = "100 GeV";
  // // procNa[5] = "120 GeV";
  // procNa[4] = "150 GeV";
  // // procNa[7] = "180 GeV";
  // procNa[5] = "200 GeV";
  // // procNa[9] = "250 GeV";
  // procNa[6] = "300 GeV";
  // procNa[7] = "500 GeV";
  procNa[0] = "15 GeV";
  procNa[1] = "30 GeV";
  procNa[2] = "50 GeV";
  procNa[3] = "80 GeV";
  procNa[4] = "100 GeV";
  procNa[5] = "120 GeV";
  procNa[6] = "150 GeV";
  procNa[7] = "180 GeV";
  procNa[8] = "200 GeV";
  procNa[9] = "250 GeV";
  procNa[10] = "300 GeV";
  procNa[11] = "500 GeV";

  TString procNa_weights[4];
  procNa_weights[0] = "devsMIPs";
  procNa_weights[1] = "X_{0}";
  procNa_weights[2] = "dEdx";
  procNa_weights[3] = "#lamdba";


  TString procNa_layers[(*ssvec).size()-1];
  procNa_layers[0] = "Before calo";
  procNa_layers[1] = "Layer 1";
  procNa_layers[2] = "Layer 2";
  procNa_layers[3] = "Layer 3";
  procNa_layers[4] = "Layer 4";
  procNa_layers[5] = "Layer 5";
  procNa_layers[6] = "Layer 6";
  procNa_layers[7] = "Layer 7";
  procNa_layers[8] = "Layer 8";
  procNa_layers[9] = "Layer 9";
  procNa_layers[10] = "Layer 10";
  procNa_layers[11] = "Layer 11";
  procNa_layers[12] = "Layer 12";
  procNa_layers[13] = "Layer 13";
  procNa_layers[14] = "Layer 14";
  procNa_layers[15] = "Layer 15";
  procNa_layers[16] = "Layer 16";
  procNa_layers[17] = "Layer 17";
  procNa_layers[18] = "Layer 18";
  procNa_layers[19] = "Layer 19";
  procNa_layers[20] = "Layer 20";
  procNa_layers[21] = "Layer 21";
  procNa_layers[22] = "Layer 22";
  procNa_layers[23] = "Layer 23";
  procNa_layers[24] = "Layer 24";
  procNa_layers[25] = "Layer 25";
  procNa_layers[26] = "Layer 26";
  procNa_layers[27] = "Layer 27";
  procNa_layers[28] = "Layer 28";
  procNa_layers[29] = "Layer 29";
  procNa_layers[30] = "Layer 30";

  //======================================================================================
  //The files that we will store the results of this analysis for the combined plot
  TString res_com[numberofconfigs];
  TFile* results_com[numberofconfigs];
  for (int k=0; k<numberofconfigs; k++){
    TString fp = "SiWEcal_" + IntToString(configs[k]);	
    res_com[k] = fp + "_combinedplots.root";
    //std::cout << res[k] << std::endl;
  }
  TH2F *hist[(*ssvec).size()-1];
  TH2F *hist2[(*ssvec).size()-1];
  TH1F *hist_sf[(*ssvec).size()-1];
  TGraphErrors *gr[(*ssvec).size()-1];
  TGraphErrors *gr2[(*ssvec).size()-1];
  // //use a TProfile to convert the 2-d to 1-d fit problem
  // TProfile *prof[(*ssvec).size()-1];
  //For accessing the fit results
  TF1 *fit[(*ssvec).size()-1][numberofenergies];
  TF1 *fit2[(*ssvec).size()-1][numberofenergies];
  // TF1 *line0 = new TF1("line0","[0]*x+[1]",0.,20);
  // line0->SetParameter(0, 1);
  //line0->SetParameter(1, 0);
  TGraphErrors *gr_dEvsMIPs_fit[(*ssvec).size()-1];
  TGraphErrors *gr_dEvsMIPs_fit_c[(*ssvec).size()-1];
  TGraphErrors *gr_MIPSvsdEplusMIPs_fit[(*ssvec).size()-1];
  TGraphErrors *gr_MIPSvsdEplusMIPs_fit_c[(*ssvec).size()-1];
  // TMultiGraph* mg = new TMultiGraph();
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  TString titleofplot1,titleofplot2,titleofplot3,titleofplot4,titleofplot; 

  if (step2){

    for (int k=0; k<numberofconfigs; k++){
      results[k]= TFile::Open(res[k],"read");
      std::cout << "Results file " << res[k] << std::endl;
    
      //======================================================================================
      //Loop on energies
      for (int l=0; l<numberofenergies; l++){
	// for (int l=0; l<1; l++){

	//======================================================================================
	//Loop on layers
	//-1 in layers for the leakage layer
	//Starting from 1 since 0 is before calorimeter
	for(unsigned iL(1); iL<(*ssvec).size()-1; iL++){
	  // for(unsigned iL(0); iL<2; iL++){

	  std::string auxnam1 = "h_dEvsMIPs_" + IntToString((int) iL);
	  std::string auxnam1_tgr = "dEvsMIPs_" + IntToString((int) iL);
	  std::string auxnam1_tgr2 = "MIPSvsdEplusMIPs_" + IntToString((int) iL);
	  std::string auxnam1_sf2 = "h_MIPSvsdEplusMIPs_" + IntToString((int) iL);
	  std::string auxnam2 = "_" + IntToString(particleenergies[l]);
	  std::string auxnam = auxnam1 + auxnam2;
	  std::string auxnam_tgr = auxnam1_tgr + auxnam2;
	  std::string auxnam_tgr2 = auxnam1_tgr2 + auxnam2;
	  std::string auxnam_sf2 = auxnam1_sf2 + auxnam2;


	  //For the dE vs MIPs
	  //------------------------------------------------------------------------------------------------
	  // std::cout << auxnam << std::endl;
	  // hist[(int) iL] = (TH2F*)results[k]->Get(auxnam.c_str());
	  // titleofplot1 = "Energy lost in materials in front of sensor vs measured energy for layer "; 
	  // titleofplot2 = IntToString( ((int) iL) ) + " and ";
	  // titleofplot3 =  particle +  " beam particle gun"; 
	  // titleofplot = titleofplot1 + titleofplot2 + titleofplot3;
	  // hist[(int) iL]->SetTitle(titleofplot); 
	  // hist[(int) iL]->GetXaxis()->SetTitle("Measured energy (MIPs)"); 
	  // hist[(int) iL]->GetYaxis()->SetTitle("dE (GeV)");
	  // // if (l==5) hist[1]->GetXaxis()->SetRangeUser(0.,(hist[(int) iL]->FindLastBinAbove(0.000001,1)) + 100);
	  // // hist[(int) iL]->SetMarkerSize(1.5);
	  // // c1[(int) iL]->SetLogy();
	  // if ( l == 0){hist[(int) iL]->SetMarkerColor(4);hist[(int) iL]->SetMarkerStyle(20);}
	  // if ( l == 1){hist[(int) iL]->SetMarkerColor(2);hist[(int) iL]->SetMarkerStyle(21);}
	  // if ( l == 2){hist[(int) iL]->SetMarkerColor(1);hist[(int) iL]->SetMarkerStyle(22);}
	  // if ( l == 3){hist[(int) iL]->SetMarkerColor(3);hist[(int) iL]->SetMarkerStyle(23);}
	  // if ( l == 4){hist[(int) iL]->SetMarkerColor(5);hist[(int) iL]->SetMarkerStyle(28);}
	  // if ( l == 5){hist[(int) iL]->SetMarkerColor(6);hist[(int) iL]->SetMarkerStyle(29);}
	  // if ( l == 6){hist[(int) iL]->SetMarkerColor(7);hist[(int) iL]->SetMarkerStyle(30);}
	  // if ( l == 7){hist[(int) iL]->SetMarkerColor(8);hist[(int) iL]->SetMarkerStyle(31);}
	  // if ( l == 8){hist[(int) iL]->SetMarkerColor(9);hist[(int) iL]->SetMarkerStyle(33);}
	  // if ( l == 9){hist[(int) iL]->SetMarkerColor(12);hist[(int) iL]->SetMarkerStyle(34);}
	  // if ( l == 10){hist[(int) iL]->SetMarkerColor(46);hist[(int) iL]->SetMarkerStyle(25);}
	  // if ( l == 11){hist[(int) iL]->SetMarkerColor(49);hist[(int) iL]->SetMarkerStyle(26);}
	  // //hist[(int) iL]->SetLineColor(l+1);
	  // leg[(int) iL]->AddEntry(hist[(int) iL], procNa[l], "P");
	  // // leg[(int) iL]->AddEntry(gr[(int) iL], procNa[l], "P");
	  // hist[(int) iL]->GetYaxis()->SetRangeUser(0.,40.);
	  // hist[(int) iL]->GetXaxis()->SetRangeUser(0.,3200.); 
	  // c1[(int) iL]->cd();
	  // c1[(int) iL]->Update(); 
	  // // if (l==5){mg->Draw("ap same");}
	  // l == 0 ? hist[(int) iL]->Draw("HIST") : hist[(int) iL]->Draw("HISTsame");
	  // l == 0 ? leg[(int) iL]->Draw() : leg[(int) iL]->Draw("same");
	  // // l == 0 ? fit[(int) iL]->Draw() : fit[(int) iL]->Draw("same");
	  // c1[(int) iL]->Update(); 
	  // // c1[(int) iL]->SaveAs("dEvsMIPs_"+particle+"_layer"+ IntToString( ((int) iL) ) + ".png");
	


	  //For the Esilicon vs (Epassive+Esilicon), that is the sampling fraction
	  //------------------------------------------------------------------------------------------------
	  std::cout << auxnam_sf2 << std::endl;
	  hist2[(int) iL] = (TH2F*)results[k]->Get(auxnam_sf2.c_str());
	  titleofplot1 = "Measured energy vs absorbed plus measured energy for layer "; 
	  titleofplot2 = IntToString( ((int) iL) ) + " and ";
	  titleofplot3 =  particle +  " beam particle gun"; 
	  titleofplot = titleofplot1 + titleofplot2 + titleofplot3;
	  hist2[(int) iL]->SetTitle(titleofplot); 
	  hist2[(int) iL]->GetXaxis()->SetTitle("E_{active} + E_{passive} (GeV)"); 
	  hist2[(int) iL]->GetYaxis()->SetTitle("E_{active} (GeV)");
	  // if (l==5) hist2[1]->GetXaxis()->SetRangeUser(0.,(hist2[(int) iL]->FindLastBinAbove(0.000001,1)) + 100);
	  // hist2[(int) iL]->SetMarkerSize(1.5);
	  // c1[(int) iL]->SetLogy();
	  if ( l == 0){hist2[(int) iL]->SetMarkerColor(4);hist2[(int) iL]->SetMarkerStyle(20);}
	  if ( l == 1){hist2[(int) iL]->SetMarkerColor(2);hist2[(int) iL]->SetMarkerStyle(21);}
	  if ( l == 2){hist2[(int) iL]->SetMarkerColor(1);hist2[(int) iL]->SetMarkerStyle(22);}
	  if ( l == 3){hist2[(int) iL]->SetMarkerColor(3);hist2[(int) iL]->SetMarkerStyle(23);}
	  if ( l == 4){hist2[(int) iL]->SetMarkerColor(5);hist2[(int) iL]->SetMarkerStyle(28);}
	  if ( l == 5){hist2[(int) iL]->SetMarkerColor(6);hist2[(int) iL]->SetMarkerStyle(29);}
	  if ( l == 6){hist2[(int) iL]->SetMarkerColor(7);hist2[(int) iL]->SetMarkerStyle(30);}
	  if ( l == 7){hist2[(int) iL]->SetMarkerColor(8);hist2[(int) iL]->SetMarkerStyle(31);}
	  if ( l == 8){hist2[(int) iL]->SetMarkerColor(9);hist2[(int) iL]->SetMarkerStyle(33);}
	  if ( l == 9){hist2[(int) iL]->SetMarkerColor(12);hist2[(int) iL]->SetMarkerStyle(34);}
	  if ( l == 10){hist2[(int) iL]->SetMarkerColor(46);hist2[(int) iL]->SetMarkerStyle(25);}
	  if ( l == 11){hist2[(int) iL]->SetMarkerColor(49);hist2[(int) iL]->SetMarkerStyle(26);}
	  //hist2[(int) iL]->SetLineColor(l+1);
	  leg_sf2[(int) iL]->AddEntry(hist2[(int) iL], procNa[l], "P");
	  // leg_sf2[(int) iL]->AddEntry(gr[(int) iL], procNa[l], "P");
	  hist2[(int) iL]->GetYaxis()->SetRangeUser(0.,40.);
	  hist2[(int) iL]->GetXaxis()->SetRangeUser(0.,520.); 
	  c9[(int) iL]->cd();
	  c9[(int) iL]->Update(); 
	  // if (l==5){mg->Draw("ap same");}
	  l == 0 ? hist2[(int) iL]->Draw("HIST") : hist2[(int) iL]->Draw("HISTsame");
	  l == 0 ? leg_sf2[(int) iL]->Draw() : leg_sf2[(int) iL]->Draw("same");
	  // l == 0 ? fit[(int) iL]->Draw() : fit[(int) iL]->Draw("same");
	  c9[(int) iL]->Update(); 
	  c9[(int) iL]->SaveAs("MIPSvsdEplusMIPs_"+particle+"_layer"+ IntToString( ((int) iL) ) + ".png");
	
	  //For the scale factors
	  //------------------------------------------------------------------------------------------------
	  // std::string auxnam1_sf = "h_SF_" + IntToString((int) iL);
	  // std::string auxnam_sf = auxnam1_sf + auxnam2;
	  // hist_sf[(int) iL] = (TH1F*)results[k]->Get(auxnam_sf.c_str());
	  // titleofplot1 = "Scale factors for layer "; 
	  // titleofplot2 = IntToString( ((int) iL) ) + " and ";
	  // titleofplot3 =  particle +  " beam particle gun"; 
	  // titleofplot = titleofplot1 + titleofplot2 + titleofplot3;
	  // hist_sf[(int) iL]->SetTitle(titleofplot); 
	  // hist_sf[(int) iL]->GetXaxis()->SetTitle("Scale factor = E_{active}/(E_{passive} + E_{active})"); 
	  // hist_sf[(int) iL]->GetYaxis()->SetTitle("Events/0.0001");
	  // // c1->SetLogy();
	  // if ( l == 0){hist_sf[(int) iL]->SetLineColor(4);}
	  // if ( l == 1){hist_sf[(int) iL]->SetLineColor(2);}
	  // if ( l == 2){hist_sf[(int) iL]->SetLineColor(1);}
	  // if ( l == 3){hist_sf[(int) iL]->SetLineColor(3);}
	  // if ( l == 4){hist_sf[(int) iL]->SetLineColor(5);}
	  // if ( l == 5){hist_sf[(int) iL]->SetLineColor(6);}
	  // if ( l == 6){hist_sf[(int) iL]->SetLineColor(7);}
	  // if ( l == 7){hist_sf[(int) iL]->SetLineColor(8);}
	  // if ( l == 8){hist_sf[(int) iL]->SetLineColor(9);}
	  // if ( l == 9){hist_sf[(int) iL]->SetLineColor(12);}
	  // if ( l == 10){hist_sf[(int) iL]->SetLineColor(46);}
	  // if ( l == 11){hist_sf[(int) iL]->SetLineColor(49);}
	  // //hist->SetLineColor(l+1);
	  // leg_sf[(int) iL]->AddEntry(hist_sf[(int) iL], procNa[l], "L");
	  // //hist->GetYaxis()->SetRangeUser(0.,900.);//  21: 100. 22: 900.
	  // c2[(int) iL]->cd();
	  // c2[(int) iL]->Update(); 
	  // l == 0 ? hist_sf[(int) iL]->Draw("HIST") : hist_sf[(int) iL]->Draw("HISTsame");
	  // l == 0 ? leg_sf[(int) iL]->Draw() : leg_sf[(int) iL]->Draw("same");
	  // c2[(int) iL]->Update(); 

	  //For the fit 
	  //------------------------------------------------------------------------------------------------
	  // gr[(int) iL] = (TGraphErrors*)results[k]->Get(auxnam_tgr.c_str());
	  // gr[(int) iL]->GetXaxis()->SetRangeUser(0.,150000.);
	  // gr[(int) iL]->Fit("pol1");
	
	  // // prof[(int) iL] = hist[(int) iL]->ProfileX();
	  // // prof[(int) iL]->Fit("pol1");
	  // fit[(int) iL][l] = gr[(int) iL]->GetFunction("pol1");
	  // //Results of the fit
	  // std::cout << "First parameter: " << fit[(int) iL][l]->GetParameter(0) << " Second parameter: " << fit[(int) iL][l]->GetParameter(1) << std::endl;
	
	  // titleofplot1 = "Energy lost in materials in front of sensor vs measured energy for layer "; 
	  // titleofplot2 = IntToString( ((int) iL) ) + " and ";
	  // titleofplot3 =  particle +  " "; 
	  // titleofplot4 =  IntToString(particleenergies[l]) +  " GeV beam particle gun"; 
	  // titleofplot = titleofplot1 + titleofplot2 + titleofplot3 + titleofplot4;
	  // gr[(int) iL]->SetTitle(titleofplot); 
	  // gr[(int) iL]->GetHistogram()->SetXTitle("Measured energy (MIPs)"); 
	  // gr[(int) iL]->GetHistogram()->SetYTitle("dE in passive material (GeV)");
	  // if ( l == 0){fit[(int) iL][l]->SetLineColor(4);gr[(int) iL]->SetMarkerColor(4);gr[(int) iL]->SetMarkerStyle(20);}
	  // if ( l == 1){fit[(int) iL][l]->SetLineColor(2);gr[(int) iL]->SetMarkerColor(2);gr[(int) iL]->SetMarkerStyle(21);}
	  // if ( l == 2){fit[(int) iL][l]->SetLineColor(1);gr[(int) iL]->SetMarkerColor(1);gr[(int) iL]->SetMarkerStyle(22);}
	  // if ( l == 3){fit[(int) iL][l]->SetLineColor(3);gr[(int) iL]->SetMarkerColor(3);gr[(int) iL]->SetMarkerStyle(23);}
	  // if ( l == 4){fit[(int) iL][l]->SetLineColor(5);gr[(int) iL]->SetMarkerColor(5);gr[(int) iL]->SetMarkerStyle(28);}
	  // if ( l == 5){fit[(int) iL][l]->SetLineColor(6);gr[(int) iL]->SetMarkerColor(6);gr[(int) iL]->SetMarkerStyle(29);}
	  // if ( l == 6){fit[(int) iL][l]->SetLineColor(7);gr[(int) iL]->SetMarkerColor(7);gr[(int) iL]->SetMarkerStyle(30);}
	  // if ( l == 7){fit[(int) iL][l]->SetLineColor(8);gr[(int) iL]->SetMarkerColor(8);gr[(int) iL]->SetMarkerStyle(31);}
	  // if ( l == 8){fit[(int) iL][l]->SetLineColor(9);gr[(int) iL]->SetMarkerColor(9);gr[(int) iL]->SetMarkerStyle(33);}
	  // if ( l == 9){fit[(int) iL][l]->SetLineColor(12);gr[(int) iL]->SetMarkerColor(12);gr[(int) iL]->SetMarkerStyle(34);}
	  // if ( l == 10){fit[(int) iL][l]->SetLineColor(46);gr[(int) iL]->SetMarkerColor(46);gr[(int) iL]->SetMarkerStyle(25);}
	  // if ( l == 11){fit[(int) iL][l]->SetLineColor(49);gr[(int) iL]->SetMarkerColor(49);gr[(int) iL]->SetMarkerStyle(26);}
	  // //hist[(int) iL]->SetLineColor(l+1);
	  // // leg_devsmipswithfit[(int) iL]->AddEntry(hist[(int) iL], procNa[l], "P");
	  // // leg_devsmipswithfit[(int) iL]->AddEntry(gr[(int) iL], procNa[l], "P");
	  // //hist[(int) iL]->GetYaxis()->SetRangeUser(0.,900.);//  21: 100. 22: 900.
	  // c4[(int) iL][l]->cd();
	  // c4[(int) iL][l]->Update(); 
	  // gr[(int) iL]->Draw("AP");
	  // // l == 0 ? gr[(int) iL]->Draw("AP") : gr[(int) iL]->Draw("APS");
	  // // leg_devsmipswithfit[(int) iL]->Draw("PS");
	  // // l == 0 ? leg_devsmipswithfit[(int) iL]->Draw() : leg_devsmipswithfit[(int) iL]->Draw("same");
	  // // l == 0 ? fit[(int) iL]->Draw() : fit[(int) iL]->Draw("same");
	  // gr[(int) iL]->GetHistogram()->GetXaxis()->SetRangeUser(0.,150000.);

	  char buf[500];
	  TLatex lat;
	  double latx = 0.;
	  double laty = 0.;
	  // // if (l==0 || l==2 || l==4 || l==7){
	  // latx = gr[(int) iL]->GetHistogram()->GetXaxis()->GetXmin()+(gr[(int) iL]->GetHistogram()->GetXaxis()->GetXmax()-gr[(int) iL]->GetHistogram()->GetXaxis()->GetXmin())/20.;
	  // laty = gr[(int) iL]->GetHistogram()->GetMaximum();
	  // sprintf(buf,"dE = p0 + p1 * MIPs");
	  // lat.DrawLatex(latx,laty*0.8,buf);
	  // sprintf(buf,"p0 = %3.3f +/- %3.3f",fit[(int) iL][l]->GetParameter(0),fit[(int) iL][l]->GetParError(0));
	  // lat.DrawLatex(latx,laty*0.7,buf);
	  // sprintf(buf,"p1 = %3.3f +/- %3.3f",fit[(int) iL][l]->GetParameter(1),fit[(int) iL][l]->GetParError(1));
	  // lat.DrawLatex(latx,laty*0.6,buf);
	  // sprintf(buf,"#chi^{2}/N = %3.3f/%d = %3.3f",fit[(int) iL][l]->GetChisquare(),fit[(int) iL][l]->GetNDF(),fit[(int) iL][l]->GetChisquare()/fit[(int) iL][l]->GetNDF());
	  // // }
	  // lat.DrawLatex(latx,laty*0.5,buf);
 

	  // c4[(int) iL][l]->Update(); 
	  TString cantopng1 = "dEvsMIPs_plusfit_"+ particle;
	  TString cantopng2 = "_" + IntToString(particleenergies[l]);
	  TString cantopng3 = "GeV_layer" + IntToString((int) iL);
	  TString cantopng = cantopng1 + cantopng2 + cantopng3 + ".png";
	
	  // if (l==0 || l==2 || l==4 || l==7){
  	  // c4[(int) iL][l]->SaveAs(cantopng);
	  // }
	  
	  
	  //For the second fit 
	  //------------------------------------------------------------------------------------------------
	  gr2[(int) iL] = (TGraphErrors*)results[k]->Get(auxnam_tgr2.c_str());
	  // gr2[(int) iL]->GetXaxis()->SetRangeUser(0.,150000.);
	  gr2[(int) iL]->Fit("pol1");
	
	  // prof[(int) iL] = hist[(int) iL]->ProfileX();
	  // prof[(int) iL]->Fit("pol1");
	  fit2[(int) iL][l] = gr2[(int) iL]->GetFunction("pol1");
	  //Results of the fit
	  std::cout << "First parameter: " << fit2[(int) iL][l]->GetParameter(0) << " Second parameter: " << fit2[(int) iL][l]->GetParameter(1) << std::endl;
	
	  titleofplot1 = "Measured energy vs absorbed plus measured energy for layer "; 
	  titleofplot2 = IntToString( ((int) iL) ) + " and ";
	  titleofplot3 =  particle +  " "; 
	  titleofplot4 =  IntToString(particleenergies[l]) +  " GeV beam particle gun"; 
	  titleofplot = titleofplot1 + titleofplot2 + titleofplot3 + titleofplot4;
	  gr2[(int) iL]->SetTitle(titleofplot); 
	  gr2[(int) iL]->GetHistogram()->SetXTitle("E_{active} + E_{passive} (GeV)"); 
	  gr2[(int) iL]->GetHistogram()->SetYTitle("E_{active} (GeV)");
	  if ( l == 0){fit2[(int) iL][l]->SetLineColor(4);gr2[(int) iL]->SetMarkerColor(4);gr2[(int) iL]->SetMarkerStyle(20);}
	  if ( l == 1){fit2[(int) iL][l]->SetLineColor(2);gr2[(int) iL]->SetMarkerColor(2);gr2[(int) iL]->SetMarkerStyle(21);}
	  if ( l == 2){fit2[(int) iL][l]->SetLineColor(1);gr2[(int) iL]->SetMarkerColor(1);gr2[(int) iL]->SetMarkerStyle(22);}
	  if ( l == 3){fit2[(int) iL][l]->SetLineColor(3);gr2[(int) iL]->SetMarkerColor(3);gr2[(int) iL]->SetMarkerStyle(23);}
	  if ( l == 4){fit2[(int) iL][l]->SetLineColor(5);gr2[(int) iL]->SetMarkerColor(5);gr2[(int) iL]->SetMarkerStyle(28);}
	  if ( l == 5){fit2[(int) iL][l]->SetLineColor(6);gr2[(int) iL]->SetMarkerColor(6);gr2[(int) iL]->SetMarkerStyle(29);}
	  if ( l == 6){fit2[(int) iL][l]->SetLineColor(7);gr2[(int) iL]->SetMarkerColor(7);gr2[(int) iL]->SetMarkerStyle(30);}
	  if ( l == 7){fit2[(int) iL][l]->SetLineColor(8);gr2[(int) iL]->SetMarkerColor(8);gr2[(int) iL]->SetMarkerStyle(31);}
	  if ( l == 8){fit2[(int) iL][l]->SetLineColor(9);gr2[(int) iL]->SetMarkerColor(9);gr2[(int) iL]->SetMarkerStyle(33);}
	  if ( l == 9){fit2[(int) iL][l]->SetLineColor(12);gr2[(int) iL]->SetMarkerColor(12);gr2[(int) iL]->SetMarkerStyle(34);}
	  if ( l == 10){fit2[(int) iL][l]->SetLineColor(46);gr2[(int) iL]->SetMarkerColor(46);gr2[(int) iL]->SetMarkerStyle(25);}
	  if ( l == 11){fit2[(int) iL][l]->SetLineColor(49);gr2[(int) iL]->SetMarkerColor(49);gr2[(int) iL]->SetMarkerStyle(26);}
	  //hist[(int) iL]->SetLineColor(l+1);
	  // leg_devsmipswithfit[(int) iL]->AddEntry(hist[(int) iL], procNa[l], "P");
	  // leg_devsmipswithfit[(int) iL]->AddEntry(gr2[(int) iL], procNa[l], "P");
	  //hist[(int) iL]->GetYaxis()->SetRangeUser(0.,900.);//  21: 100. 22: 900.
	  c10[(int) iL][l]->cd();
	  c10[(int) iL][l]->Update(); 
	  gr2[(int) iL]->Draw("AP");
	  // l == 0 ? gr2[(int) iL]->Draw("AP") : gr2[(int) iL]->Draw("APS");
	  // leg_devsmipswithfit[(int) iL]->Draw("PS");
	  // l == 0 ? leg_devsmipswithfit[(int) iL]->Draw() : leg_devsmipswithfit[(int) iL]->Draw("same");
	  // l == 0 ? fit[(int) iL]->Draw() : fit[(int) iL]->Draw("same");
	  gr2[(int) iL]->GetHistogram()->GetXaxis()->SetRangeUser(0.,520.);


	  latx = 0.;
	  laty = 0.;
	  latx = gr2[(int) iL]->GetHistogram()->GetXaxis()->GetXmin()+(gr2[(int) iL]->GetHistogram()->GetXaxis()->GetXmax()-gr2[(int) iL]->GetHistogram()->GetXaxis()->GetXmin())/20.;
	  laty = gr2[(int) iL]->GetHistogram()->GetMaximum();
	  sprintf(buf,"E_{active} = p0 + p1 * (E_{active} + E_{passive}) ");
	  lat.DrawLatex(latx,laty*0.8,buf);
	  sprintf(buf,"p0 = %3.3f +/- %3.3f",fit2[(int) iL][l]->GetParameter(0),fit2[(int) iL][l]->GetParError(0));
	  lat.DrawLatex(latx,laty*0.7,buf);
	  sprintf(buf,"p1 = %3.3f +/- %3.3f",fit2[(int) iL][l]->GetParameter(1),fit2[(int) iL][l]->GetParError(1));
	  lat.DrawLatex(latx,laty*0.6,buf);
	  sprintf(buf,"#chi^{2}/N = %3.3f/%d = %3.3f",fit2[(int) iL][l]->GetChisquare(),fit2[(int) iL][l]->GetNDF(),fit2[(int) iL][l]->GetChisquare()/fit2[(int) iL][l]->GetNDF());
	  lat.DrawLatex(latx,laty*0.5,buf);
 

	  c10[(int) iL][l]->Update(); 
	  cantopng1 = "MIPSvsdEplusMIPs_plusfit_"+ particle;
	  cantopng2 = "_" + IntToString(particleenergies[l]);
	  cantopng3 = "GeV_layer" + IntToString((int) iL);
	  cantopng = cantopng1 + cantopng2 + cantopng3 + ".png";
	
  	  c10[(int) iL][l]->SaveAs(cantopng);







	}//Loop on layers
      
      }//Loop on energies

      results_com[k]= new TFile(res_com[k],"recreate");
    
      //======================================================================================
      //Write canvas with combined plot 
      // c1->Print("Leakage.pdf",".pdf");
      for(unsigned iL(1); iL<(*ssvec).size()-1; iL++){
	c1[(int) iL]->Write();
	c9[(int) iL]->Write();
	c2[(int) iL]->Write();
	for (int l=0; l<numberofenergies; l++){
	  c4[(int) iL][l]->Write();
	  c10[(int) iL][l]->Write();
	}
	// gr[(int) iL]->Write();
	gr2[(int) iL]->Write();
      
      }
      //Reset histos
      // hist->Reset();

      //======================================================================================
      //Save the fit results for faster running
      // std::ofstream myfile;
      // TString titleoffile = "dEvsMIPs_" + particle + ".txt";
      // myfile.open(titleoffile);

      // //Here print the p0 and p1 that are going to be used for Ereco
      // for (int l=0; l<numberofenergies; l++){
      // 	std::cout << "=========================================================" << std::endl;
      // 	std::cout << "Particle " <<  particle << std::endl;
      // 	std::cout << "Energy of the beam " <<  particleenergies[l] << std::endl;
      // 	for(unsigned iL(1); iL<(*ssvec).size()-1; iL++){
      // 	  std::cout << "---------------------------------------------------------" << std::endl;
      // 	  std::cout << "Layer " << iL << std::endl;
      // 	  std::cout << "dEvsMIPs slope: " <<  fit[(int) iL][l]->GetParameter(1) << " error slope " << fit[(int) iL][l]->GetParError(1) << std::endl;
      // 	  std::cout << "dEvsMIPs constant: " <<  fit[(int) iL][l]->GetParameter(0) << " error constant " << fit[(int) iL][l]->GetParError(0) << std::endl;

      // 	  //particle particleenergy layer constant slope errorconstant errorslope
      // 	  myfile << particle << " " << particleenergies[l] << " " << iL << " " << fit[(int) iL][l]->GetParameter(0) << " " << fit[(int) iL][l]->GetParameter(1) << " " << fit[(int) iL][l]->GetParError(0) <<  " " << fit[(int) iL][l]->GetParError(1) << "\n";

      // 	  //the weights in our scheme
      // 	  wdevsmips[(int) iL].push_back( fit[(int) iL][l]->GetParameter(1) );
      // 	  wdevsmipserror[(int) iL].push_back( fit[(int) iL][l]->GetParError(1) );

      // 	}
      // }

      //======================================================================================
      //Save the second fit results for faster running
      std::ofstream myfile;
      TString titleoffile = "MIPSvsdEplusMIPs_" + particle + ".txt";
      myfile.open(titleoffile);

      //Here print the p0 and p1 that are going to be used for Ereco
      for (int l=0; l<numberofenergies; l++){
	std::cout << "=========================================================" << std::endl;
	std::cout << "Particle " <<  particle << std::endl;
	std::cout << "Energy of the beam " <<  particleenergies[l] << std::endl;
	for(unsigned iL(1); iL<(*ssvec).size()-1; iL++){
	  std::cout << "---------------------------------------------------------" << std::endl;
	  std::cout << "Layer " << iL << std::endl;
	  std::cout << "MIPSvsdEplusMIPs slope: " <<  fit2[(int) iL][l]->GetParameter(1) << " error slope " << fit2[(int) iL][l]->GetParError(1) << std::endl;
	  std::cout << "MIPSvsdEplusMIPs constant: " <<  fit2[(int) iL][l]->GetParameter(0) << " error constant " << fit2[(int) iL][l]->GetParError(0) << std::endl;

	  //particle particleenergy layer constant slope errorconstant errorslope
	  myfile << particle << " " << particleenergies[l] << " " << iL << " " << fit2[(int) iL][l]->GetParameter(0) << " " << fit2[(int) iL][l]->GetParameter(1) << " " << fit2[(int) iL][l]->GetParError(0) <<  " " << fit2[(int) iL][l]->GetParError(1) << "\n";

	  //the weights in our scheme
	  wdevsmips2[(int) iL].push_back( fit2[(int) iL][l]->GetParameter(1) );
	  wdevsmipserror2[(int) iL].push_back( fit2[(int) iL][l]->GetParError(1) );

	}
      }

      //=========================================================================================
      //Here we will make the plot for the fit results slope
      // c3->cd();

      // for(unsigned iL(1); iL<(*ssvec).size()-1; iL++){
      // 	// for(unsigned iL(0); iL<2; iL++){
      // 	std::string auxnam1_tgr = "dEoverMIPs_" + IntToString((int) iL);
      // 	// gr_dEvsMIPs_fit[(int) iL]= new TGraphErrors(auxnam1_tgr.c_str());
      // 	gr_dEvsMIPs_fit[(int) iL]= (TGraphErrors *) calib->Clone( auxnam1_tgr.c_str() );

      
      // 	// gr_dEvsMIPs_fit[(int) iL]->SetMarkerStyle( 19+((int) iL)); 
      
      // 	// gr_dEvsMIPs_fit[(int) iL]->SetMarkerSize(0.2); 
      // 	gr_dEvsMIPs_fit[(int) iL]->SetMarkerColor(((int) iL) + 1);  
      // 	gr_dEvsMIPs_fit[(int) iL]->SetLineColor(((int) iL) + 1);

      // 	if ((((int) iL) == 0)){
      // 	  gr_dEvsMIPs_fit[(int) iL]->SetMarkerColor(28);  
      // 	  gr_dEvsMIPs_fit[(int) iL]->SetLineColor(28);
      // 	}
      // 	if (((int) iL) == 10){
      // 	  gr_dEvsMIPs_fit[(int) iL]->SetMarkerColor(46);  
      // 	  gr_dEvsMIPs_fit[(int) iL]->SetLineColor(46);
      // 	}
 
      // 	c3->Update(); 

      // 	// Int_t np=gr_dEvsMIPs_fit[(int) iL]->GetN();
      // 	for (int l=0; l<numberofenergies; l++){
      // 	  // for (int l=0; l<1; l++){
      // 	  gr_dEvsMIPs_fit[(int) iL]->SetPoint( l , particleenergies[l] , fit[(int) iL][l]->GetParameter(1) );
      // 	  gr_dEvsMIPs_fit[(int) iL]->SetPointError(l , 0. , fit[(int) iL][l]->GetParError(1) );
      // 	}
      // 	leg_fit->AddEntry(gr_dEvsMIPs_fit[(int) iL], procNa_layers[(int) iL], "LP");
      // 	//Be careful here if you include before calo layer !!!!!!!!!!!!!!!!!!!!!!!!!!
      // 	((int) iL) == 1 ? gr_dEvsMIPs_fit[(int) iL]->Draw("APL") : gr_dEvsMIPs_fit[(int) iL]->Draw("PSL");
      // }
      // leg_fit->Draw("PS");
      // c3->Update(); 

      // gr_dEvsMIPs_fit[1]->SetTitle( "Energy lost vs measured MIPs fit results for " + particle + " for all beam energies" );
      // gr_dEvsMIPs_fit[1]->GetHistogram()->SetXTitle("Beam Energy (GeV)"); 
      // gr_dEvsMIPs_fit[1]->GetHistogram()->SetYTitle("E_{passive} = f(E_{active}): Slope (GeV/MIP)");
      // gr_dEvsMIPs_fit[1]->GetXaxis()->SetRangeUser(0.,520.);  
      // gr_dEvsMIPs_fit[1]->GetYaxis()->SetRangeUser(-1.,1.);  

      // c3->Update(); 
      // c3->SaveAs("dEvsMIPs_fitresults_slope" + particle + ".png"); 
      // c3->Write();

      //=========================================================================================
      //Here we will make the plot for the fit results constant term
      // c8->cd();

      // for(unsigned iL(1); iL<(*ssvec).size()-1; iL++){
      // 	// for(unsigned iL(0); iL<2; iL++){
      // 	std::string auxnam1_tgr = "dEoverMIPs_const" + IntToString((int) iL);
      // 	// gr_dEvsMIPs_fit[(int) iL]= new TGraphErrors(auxnam1_tgr.c_str());
      // 	gr_dEvsMIPs_fit_c[(int) iL]= (TGraphErrors *) calib->Clone( auxnam1_tgr.c_str() );

      
      // 	// gr_dEvsMIPs_fit_c[(int) iL]->SetMarkerStyle( 19+((int) iL)); 
      
      // 	// gr_dEvsMIPs_fit_c[(int) iL]->SetMarkerSize(0.2); 
      // 	gr_dEvsMIPs_fit_c[(int) iL]->SetMarkerColor(((int) iL) + 1);  
      // 	gr_dEvsMIPs_fit_c[(int) iL]->SetLineColor(((int) iL) + 1);

      // 	if ((((int) iL) == 0)){
      // 	  gr_dEvsMIPs_fit_c[(int) iL]->SetMarkerColor(28);  
      // 	  gr_dEvsMIPs_fit_c[(int) iL]->SetLineColor(28);
      // 	}
      // 	if (((int) iL) == 10){
      // 	  gr_dEvsMIPs_fit_c[(int) iL]->SetMarkerColor(46);  
      // 	  gr_dEvsMIPs_fit_c[(int) iL]->SetLineColor(46);
      // 	}
 
      // 	c8->Update(); 

      // 	// Int_t np=gr_dEvsMIPs_fit_c[(int) iL]->GetN();
      // 	for (int l=0; l<numberofenergies; l++){
      // 	  // for (int l=0; l<1; l++){
      // 	  gr_dEvsMIPs_fit_c[(int) iL]->SetPoint( l , particleenergies[l] , fit[(int) iL][l]->GetParameter(0) );
      // 	  gr_dEvsMIPs_fit_c[(int) iL]->SetPointError(l , 0. , fit[(int) iL][l]->GetParError(0) );
      // 	}
      // 	leg_fit_c->AddEntry(gr_dEvsMIPs_fit_c[(int) iL], procNa_layers[(int) iL], "LP");
      // 	//Be careful here if you include before calo layer !!!!!!!!!!!!!!!!!!!!!!!!!!
      // 	((int) iL) == 1 ? gr_dEvsMIPs_fit_c[(int) iL]->Draw("APL") : gr_dEvsMIPs_fit_c[(int) iL]->Draw("PSL");
      // }
      // leg_fit_c->Draw("PS");
      // c8->Update(); 

      // gr_dEvsMIPs_fit_c[1]->SetTitle( "Energy lost vs measured MIPs fit results for " + particle + " for all beam energies" );
      // gr_dEvsMIPs_fit_c[1]->GetHistogram()->SetXTitle("Beam Energy (GeV)"); 
      // gr_dEvsMIPs_fit_c[1]->GetHistogram()->SetYTitle("E_{passive} = f(E_{active}): Constant term (GeV)");
      // gr_dEvsMIPs_fit_c[1]->GetXaxis()->SetRangeUser(0.,520.);  
      // gr_dEvsMIPs_fit_c[1]->GetYaxis()->SetRangeUser(-10.,10.);  

      // c8->Update(); 
      // c8->SaveAs("dEvsMIPs_fitresults_constantterm" + particle + ".png"); 
      // c8->Write();

      
      //=========================================================================================
      //Here we will make the plot for the second fit results slope
      c11->cd();

      for(unsigned iL(1); iL<(*ssvec).size()-1; iL++){
	// for(unsigned iL(0); iL<2; iL++){
	std::string auxnam1_tgr2 = "MIPSoverdEplusMIPs_" + IntToString((int) iL);
	// gr_MIPSvsdEplusMIPs_fit[(int) iL]= new TGraphErrors(auxnam1_tgr.c_str());
	gr_MIPSvsdEplusMIPs_fit[(int) iL]= (TGraphErrors *) calib->Clone( auxnam1_tgr2.c_str() );

      
	// gr_MIPSvsdEplusMIPs_fit[(int) iL]->SetMarkerStyle( 19+((int) iL)); 
      
	// gr_MIPSvsdEplusMIPs_fit[(int) iL]->SetMarkerSize(0.2); 
	gr_MIPSvsdEplusMIPs_fit[(int) iL]->SetMarkerColor(((int) iL) + 1);  
	gr_MIPSvsdEplusMIPs_fit[(int) iL]->SetLineColor(((int) iL) + 1);

	if ((((int) iL) == 0)){
	  gr_MIPSvsdEplusMIPs_fit[(int) iL]->SetMarkerColor(28);  
	  gr_MIPSvsdEplusMIPs_fit[(int) iL]->SetLineColor(28);
	}
	if (((int) iL) == 10){
	  gr_MIPSvsdEplusMIPs_fit[(int) iL]->SetMarkerColor(46);  
	  gr_MIPSvsdEplusMIPs_fit[(int) iL]->SetLineColor(46);
	}
 
	c11->Update(); 

	// Int_t np=gr_MIPSvsdEplusMIPs_fit[(int) iL]->GetN();
	for (int l=0; l<numberofenergies; l++){
	  // for (int l=0; l<1; l++){
	  gr_MIPSvsdEplusMIPs_fit[(int) iL]->SetPoint( l , particleenergies[l] , fit2[(int) iL][l]->GetParameter(1) );
	  gr_MIPSvsdEplusMIPs_fit[(int) iL]->SetPointError(l , 0. , fit2[(int) iL][l]->GetParError(1) );
	}
	leg_fit2->AddEntry(gr_MIPSvsdEplusMIPs_fit[(int) iL], procNa_layers[(int) iL], "LP");
	//Be careful here if you include before calo layer !!!!!!!!!!!!!!!!!!!!!!!!!!
	((int) iL) == 1 ? gr_MIPSvsdEplusMIPs_fit[(int) iL]->Draw("APL") : gr_MIPSvsdEplusMIPs_fit[(int) iL]->Draw("PSL");
      }
      leg_fit2->Draw("PS");
      c11->Update(); 

      gr_MIPSvsdEplusMIPs_fit[1]->SetTitle( "Measured energy vs absorbed plus measured energy fit results for " + particle + " for all beam energies" );
      gr_MIPSvsdEplusMIPs_fit[1]->GetHistogram()->SetXTitle("Beam Energy (GeV)"); 
      gr_MIPSvsdEplusMIPs_fit[1]->GetHistogram()->SetYTitle("E_{active} = f(E_{active} + E_{passive}): Slope ");
      gr_MIPSvsdEplusMIPs_fit[1]->GetXaxis()->SetRangeUser(0.,520.);  
      gr_MIPSvsdEplusMIPs_fit[1]->GetYaxis()->SetRangeUser(-1.,1.);  

      c11->Update(); 
      c11->SaveAs("MIPSvsdEplusMIPs_fitresults_slope" + particle + ".png"); 
      c11->Write();


      //=========================================================================================
      //Here we will make the plot for the second fit results constant term
      c12->cd();

      for(unsigned iL(1); iL<(*ssvec).size()-1; iL++){
	// for(unsigned iL(0); iL<2; iL++){
	std::string auxnam1_tgr = "MIPSoverdEplusMIPs_const" + IntToString((int) iL);
	// gr_MIPSvsdEplusMIPs_fit[(int) iL]= new TGraphErrors(auxnam1_tgr.c_str());
	gr_MIPSvsdEplusMIPs_fit_c[(int) iL]= (TGraphErrors *) calib->Clone( auxnam1_tgr.c_str() );

      
	// gr_MIPSvsdEplusMIPs_fit_c[(int) iL]->SetMarkerStyle( 19+((int) iL)); 
      
	// gr_MIPSvsdEplusMIPs_fit_c[(int) iL]->SetMarkerSize(0.2); 
	gr_MIPSvsdEplusMIPs_fit_c[(int) iL]->SetMarkerColor(((int) iL) + 1);  
	gr_MIPSvsdEplusMIPs_fit_c[(int) iL]->SetLineColor(((int) iL) + 1);

	if ((((int) iL) == 0)){
	  gr_MIPSvsdEplusMIPs_fit_c[(int) iL]->SetMarkerColor(28);  
	  gr_MIPSvsdEplusMIPs_fit_c[(int) iL]->SetLineColor(28);
	}
	if (((int) iL) == 10){
	  gr_MIPSvsdEplusMIPs_fit_c[(int) iL]->SetMarkerColor(46);  
	  gr_MIPSvsdEplusMIPs_fit_c[(int) iL]->SetLineColor(46);
	}
 
	c12->Update(); 

	// Int_t np=gr_MIPSvsdEplusMIPs_fit_c[(int) iL]->GetN();
	for (int l=0; l<numberofenergies; l++){
	  // for (int l=0; l<1; l++){
	  gr_MIPSvsdEplusMIPs_fit_c[(int) iL]->SetPoint( l , particleenergies[l] , fit2[(int) iL][l]->GetParameter(0) );
	  gr_MIPSvsdEplusMIPs_fit_c[(int) iL]->SetPointError(l , 0. , fit2[(int) iL][l]->GetParError(0) );
	}
	leg_fit2_c->AddEntry(gr_MIPSvsdEplusMIPs_fit_c[(int) iL], procNa_layers[(int) iL], "LP");
	//Be careful here if you include before calo layer !!!!!!!!!!!!!!!!!!!!!!!!!!
	((int) iL) == 1 ? gr_MIPSvsdEplusMIPs_fit_c[(int) iL]->Draw("APL") : gr_MIPSvsdEplusMIPs_fit_c[(int) iL]->Draw("PSL");
      }
      leg_fit2_c->Draw("PS");
      c12->Update(); 

      gr_MIPSvsdEplusMIPs_fit_c[1]->SetTitle( "Measured energy vs absorbed plus measured energy fit results for " + particle + " for all beam energies" );
      gr_MIPSvsdEplusMIPs_fit_c[1]->GetHistogram()->SetXTitle("Beam Energy (GeV)"); 
      gr_MIPSvsdEplusMIPs_fit_c[1]->GetHistogram()->SetYTitle("E_{active} = f(E_{active} + E_{passive}): Constant term (GeV)");
      gr_MIPSvsdEplusMIPs_fit_c[1]->GetXaxis()->SetRangeUser(0.,520.);  
      gr_MIPSvsdEplusMIPs_fit_c[1]->GetYaxis()->SetRangeUser(-10.,10.);  

      c12->Update(); 
      c12->SaveAs("MIPSvsdEplusMIPs_fitresults_constantterm" + particle + ".png"); 
      c12->Write();

 








      //Here for the plot of average energy lost vs radiation length
      c5->cd();
      for (int l=0; l<numberofenergies; l++){
	//dEvsXos[l] = (TGraphErrors*)results[k]->Get(auxnam_tgr.c_str());
	if (l == 0){dEvsXos[l]->SetMarkerColor(4);dEvsXos[l]->SetMarkerStyle(20);}
	if (l == 1){dEvsXos[l]->SetMarkerColor(2);dEvsXos[l]->SetMarkerStyle(21);}
	if (l == 2){dEvsXos[l]->SetMarkerColor(1);dEvsXos[l]->SetMarkerStyle(22);}
	if (l == 3){dEvsXos[l]->SetMarkerColor(3);dEvsXos[l]->SetMarkerStyle(23);}
	if (l == 4){dEvsXos[l]->SetMarkerColor(5);dEvsXos[l]->SetMarkerStyle(28);}
	if (l == 5){dEvsXos[l]->SetMarkerColor(6);dEvsXos[l]->SetMarkerStyle(29);}
	if (l == 6){dEvsXos[l]->SetMarkerColor(7);dEvsXos[l]->SetMarkerStyle(30);}
	if (l == 7){dEvsXos[l]->SetMarkerColor(8);dEvsXos[l]->SetMarkerStyle(31);}
	if (l == 8){dEvsXos[l]->SetMarkerColor(9);dEvsXos[l]->SetMarkerStyle(33);}
	if (l == 9){dEvsXos[l]->SetMarkerColor(12);dEvsXos[l]->SetMarkerStyle(34);}
	if (l == 10){dEvsXos[l]->SetMarkerColor(46);dEvsXos[l]->SetMarkerStyle(25);}
	if (l == 11){dEvsXos[l]->SetMarkerColor(49);dEvsXos[l]->SetMarkerStyle(26);}
 
	dEvsXos[l]->SetTitle( "Energy lost vs. transversed thickness for " + particle); 
	dEvsXos[l]->GetHistogram()->SetXTitle("Transversed thickness (1/X_{0})"); 
	dEvsXos[l]->GetHistogram()->SetYTitle("Average dE (MIPs)");
	dEvsXos[l]->GetXaxis()->SetRangeUser(0.,25.8);  
	l == 0 ? dEvsXos[l]->Draw("APL") : dEvsXos[l]->Draw("PSL");
	l == 0 ? leg_devsxos->Draw() : leg_devsxos->Draw("same");
	dEvsXos[0]->GetXaxis()->SetRangeUser(0.,25.8); 
	c5->Update(); 
	dEvsXos[l]->Write();
      }

      TString cnpng = "dEvsXos.png";
      c5->SaveAs(cnpng);
      c5->Write();

      //Here for the plot of weights vs shower depth
      // c6->cd();
      // TGraphErrors *gr_weights[numberofenergies];

      // for (int l=0; l<numberofenergies; l++){
      // 	std::string auxnam_tgr = "Weights_" + IntToString(l);
      // 	gr_weights[l] = (TGraphErrors *) calib->Clone( auxnam_tgr.c_str()  );
      // 	if (l == 0){gr_weights[l]->SetMarkerColor(4);gr_weights[l]->SetMarkerStyle(20);}
      // 	if (l == 1){gr_weights[l]->SetMarkerColor(2);gr_weights[l]->SetMarkerStyle(21);}
      // 	if (l == 2){gr_weights[l]->SetMarkerColor(1);gr_weights[l]->SetMarkerStyle(22);}
      // 	if (l == 3){gr_weights[l]->SetMarkerColor(3);gr_weights[l]->SetMarkerStyle(23);}
      // 	if (l == 4){gr_weights[l]->SetMarkerColor(5);gr_weights[l]->SetMarkerStyle(28);}
      // 	if (l == 5){gr_weights[l]->SetMarkerColor(6);gr_weights[l]->SetMarkerStyle(29);}
      // 	if (l == 6){gr_weights[l]->SetMarkerColor(7);gr_weights[l]->SetMarkerStyle(30);}
      // 	if (l == 7){gr_weights[l]->SetMarkerColor(8);gr_weights[l]->SetMarkerStyle(31);}
      // 	if (l == 8){gr_weights[l]->SetMarkerColor(9);gr_weights[l]->SetMarkerStyle(33);}
      // 	if (l == 9){gr_weights[l]->SetMarkerColor(12);gr_weights[l]->SetMarkerStyle(34);}
      // 	if (l == 10){gr_weights[l]->SetMarkerColor(46);gr_weights[l]->SetMarkerStyle(25);}
      // 	if (l == 11){gr_weights[l]->SetMarkerColor(49);gr_weights[l]->SetMarkerStyle(26);}
 
      // 	gr_weights[l]->SetTitle( "Weights vs. shower depth for " + particle); 

      // 	for(unsigned iL(1); iL<(*ssvec).size()-1; iL++){
      // 	  Int_t nnp=gr_weights[l]->GetN();
      // 	  gr_weights[l]->SetPoint( nnp ,  h_Showerprofile[(int) iL][l]->GetMean() , wdevsmips[(int) iL][l] );
      // 	  gr_weights[l]->SetPointError(nnp , h_Showerprofile[(int) iL][l]->GetRMS() , wdevsmipserror[(int) iL][l] );
      // 	}
      // 	leg_w->AddEntry(gr_weights[l], procNa[l], "P");

      // 	gr_weights[l]->GetHistogram()->SetXTitle("Shower depth"); 
      // 	gr_weights[l]->GetHistogram()->SetYTitle("Weights");
      // 	// gr_weights[l]->GetXaxis()->SetRangeUser(0.,25.8);  
      // 	l == 0 ? gr_weights[l]->Draw("APL") : gr_weights[l]->Draw("PSL");
      // 	l == 0 ? leg_w->Draw() : leg_w->Draw("same");
      // 	c6->Update(); 
      // }

      // TString cnpngw = "WeightsvsShowerdepth.png";
      // c6->SaveAs(cnpngw);
      // c6->Write();

      // //Comparing different weighting schemes for a 50 GeV electron
      // c7->cd();
      // TGraphErrors *gr_weights_x0;
      // TGraphErrors *gr_weights_dedx;
      // TGraphErrors *gr_weights_lamdba;

      // gr_weights_x0 = (TGraphErrors *) calib->Clone( "Weights_x0"  );
      // gr_weights_dedx = (TGraphErrors *) calib->Clone( "Weights_dedx" );
      // gr_weights_lamdba = (TGraphErrors *) calib->Clone( "Weights_lamdba"  );
    
      // for(unsigned iL(1); iL<(*ssvec).size()-1; iL++){
      // 	Int_t nnpx0=gr_weights_x0->GetN();
      // 	gr_weights_x0->SetPoint( nnpx0 ,  h_Showerprofile[(int) iL][2]->GetMean() , wXos[(int) iL][2] );
      // 	gr_weights_x0->SetPointError(nnpx0 , h_Showerprofile[(int) iL][2]->GetRMS() , 0. );
      // 	Int_t nnpdedx=gr_weights_dedx->GetN();
      // 	gr_weights_dedx->SetPoint( nnpdedx ,  h_Showerprofile[(int) iL][2]->GetMean() , wdEdx[(int) iL][2] );
      // 	gr_weights_dedx->SetPointError(nnpdedx , h_Showerprofile[(int) iL][2]->GetRMS() , 0. );
      // 	Int_t nnplambda=gr_weights_lamdba->GetN();
      // 	gr_weights_lamdba->SetPoint( nnplambda ,  h_Showerprofile[(int) iL][2]->GetMean() , wlamdba[(int) iL][2] );
      // 	gr_weights_lamdba->SetPointError(nnplambda , h_Showerprofile[(int) iL][2]->GetRMS() , 0. );
      // }
      // gr_weights[2]->SetTitle( "Weights vs. shower depth for " + particle); 
      // gr_weights[2]->SetMarkerColor(4);
      // gr_weights[2]->SetMarkerStyle(20);
      // gr_weights[2]->GetHistogram()->SetXTitle("Shower depth"); 
      // gr_weights[2]->GetHistogram()->SetYTitle("Weights");
      // // gr_weights[l]->GetXaxis()->SetRangeUser(0.,25.8);
      // gr_weights[2]->Draw("APL");

      // gr_weights_x0->SetMarkerColor(4);
      // gr_weights_x0 ->SetMarkerStyle(20);
      // gr_weights_x0->Draw("PSL");

      // gr_weights_dedx->SetMarkerColor(2);
      // gr_weights_dedx->SetMarkerStyle(21);
      // gr_weights_dedx->Draw("PSL");

      // gr_weights_lamdba->SetMarkerColor(1);
      // gr_weights_lamdba->SetMarkerStyle(22);
      // gr_weights_lamdba->Draw("PSL");
    
      // leg_allw->AddEntry(gr_weights[2], procNa_weights[0], "P");
      // leg_allw->AddEntry(gr_weights_x0, procNa_weights[1], "P");
      // leg_allw->AddEntry(gr_weights_dedx, procNa_weights[2], "P");
      // leg_allw->AddEntry(gr_weights_lamdba, procNa_weights[3], "P");
      // leg_allw->Draw("PS");

      // c7->Update(); 


      // TString cnpngwe = "WeightsvsShowerdepth_50GeVelectron.png";
      // c7->SaveAs(cnpngwe);
      // c7->Write();
    




















      results_com[k]->Close();
  
    }   //Loop on configs

  }//if to avoid running the second step


  return 0; 

}//main

//===============================================================================
//Definitions here
std::string IntToString (int number1)
{
  std::ostringstream oss1;
  oss1<< number1;
  return oss1.str();
}
