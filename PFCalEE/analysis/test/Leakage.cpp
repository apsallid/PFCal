#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include <boost/algorithm/string.hpp>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TStyle.h"

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TNtuple.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraphErrors.h"
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
  //The root files with the individual plots for all test beam configs 
  const int numberoffiles = 1;
  //======================================================================================
  //Do the files include upstream material?
  bool upstream = false;
  int startlayer = upstream ? 1 : 0;
  //======================================================================================
  //3 layers times 4 Si Pads
  // const int numberofpads = 12; 
  //======================================================================================
  //Particle type and energy
  TString particle = "e-"; //remaining e+, pi+
  const int numberofenergies = 18; // 15 30 50 80 100 120 150 180 200 250 300 500
  //const int numberofenergies = 1; // 15 30 50 80 100 150
 // int particleenergy = 15;
  std::vector<int> particleenergies;
  particleenergies.clear();
  // particleenergies.push_back(15);
  // particleenergies.push_back(30);
  // particleenergies.push_back(50);
  // particleenergies.push_back(80);
  // particleenergies.push_back(100);
  // particleenergies.push_back(120);
  // particleenergies.push_back(150);
  // particleenergies.push_back(180);
  // particleenergies.push_back(200);
  // particleenergies.push_back(250);
  // particleenergies.push_back(300);
  // particleenergies.push_back(500);
  particleenergies.push_back(2);
  particleenergies.push_back(3);
  particleenergies.push_back(5);
  particleenergies.push_back(8);
  particleenergies.push_back(10);
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
  particleenergies.push_back(400);
  particleenergies.push_back(500);

  //======================================================================================
  //Read all files and tree
  TString filename[numberofconfigs][numberofenergies][numberoffiles];
  for (int k=0; k<numberofconfigs; k++){
    for (int l=0; l<numberofenergies; l++){

      // //Here we need to know the exact number of the files 
      // std::string linux_command1 =  "/eos/cms/store/group/phys_b2g/apsallid/SiWEcal/DetectorConfigurations/Config_" + IntToString(configs[k]);
      // std::string linux_command2 =  "/" + particle; 
      // std::string linux_command3 =  "/" + IntToString(particleenergies[l]);
      // std::string linux_command4 = "/results";
      // std::string path =  linux_command1 + linux_command2 + linux_command3 + linux_command4;
      // std::string linux_command5 = " >> pp.txt ";

      // std::string linux_command6 = "eos ls " + path; 
      // std::string linux_command = linux_command6 + linux_command5; 
      // system(linux_command.c_str());
      // std::ifstream pp_file("pp.txt");
      // std::string buffer;
      // getline(pp_file, buffer);
      // std::istringstream is(buffer);
      // is >> numoffilesindataset;
      // pp_file.close();
      // system("rm -rf pp.txt");
 






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
	// files[k]= new TFile(filename[k]);
	//    files[k] = TFile::Open("/afs/cern.ch/work/a/apsallid/CMS/Geant4/SiWEcal/Parallel/testingfile/"+filename[k]+".root");
	// files[k][l][j] = TFile::Open("/afs/cern.ch/work/a/apsallid/CMS/Geant4/SiWEcal/"+filename[k][l][j]+".root");
	// files[k][l] = TFile::Open("/tmp/apsallid/Configs/"+filename[k][l]+".root");
	// lTree[k][l] = (TTree*) files[k][l]->Get("HGCSSTree");
	// lChain[k][l]->Add("/tmp/apsallid/"+filename[k][l][j]+".root");
	// lChain[k][l]->Add("root://eoscms//eos/cms/store/group/phys_b2g/apsallid/PFCal/DetectorConfigurations/Config_12/"+filename[k][l][j]+".root");
	if (upstream){
	  lChain[k][l]->Add("/afs/cern.ch/work/a/apsallid/CMS/Geant4/PFCal/PFCalEE/analysis/data/"+filename[k][l][j]+".root");
	} else {
	  // lChain[k][l]->Add("/afs/cern.ch/work/a/apsallid/CMS/Geant4/PFCal/PFCalEE/analysis/data/NoUpStream"+filename[k][l][j]+".root");
	  lChain[k][l]->Add("/tmp/apsallid/"+filename[k][l][j]+".root");
	}
	// lChain[k][l]->Add("/afs/cern.ch/work/a/apsallid/CMS/Geant4/PFCal/PFCalEE/PFcal.root");
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
  // TCanvas *myc = new TCanvas("myc","myc",500,500);

  //======================================================================================
  // Histos
  //For the loss at the back of the detector
  std::vector<TH1F*> h_Leakagefrombehind;
  h_Leakagefrombehind.clear();
  //For the loss at the front of the detector
  std::vector<TH1F*> h_frontLeakage;
  h_frontLeakage.clear();
  //For the lateral loss
  std::vector<TH1F*> h_lateralLeakage;
  h_lateralLeakage.clear();
  //For the longitudinal leakage
  std::vector<TH1F*> h_longitudinalLeakage;
  h_longitudinalLeakage.clear();
  //For the total leakage
  std::vector<TH1F*> h_totalLeakage;
  h_totalLeakage.clear();
  //For the loss before the calorimeter
  std::vector<TH1F*> h_lossbeforecalo;
  h_lossbeforecalo.clear();
  //Loss in anything but the sensors and the back of the detector
  std::vector<TH1F*> h_lossinpassivelayers;
  h_lossinpassivelayers.clear();
  //Loss in silicon
  std::vector<TH1F*> h_lossinsilicon;
  h_lossinsilicon.clear();
  //Total Loss
  std::vector<TH1F*> h_totalloss;
  h_totalloss.clear();
  //======================================================================================
  // For the fit
  std::vector<TGraphErrors*> ElossUpstreamvsTotalMeasuredMIPs;
  ElossUpstreamvsTotalMeasuredMIPs.clear();

  TGraphErrors *corr = new TGraphErrors();
  for (Int_t k=0; k<numberofenergies; k++){
    h_Leakagefrombehind.push_back(new TH1F(("h_Leakagefrombehind_" + IntToString(particleenergies[k])).c_str(),";E (GeV);SimHits",2500,0.,250.));
    h_frontLeakage.push_back(new TH1F(("h_frontLeakage_" + IntToString(particleenergies[k])).c_str(),";E (GeV);SimHits",2500,0.,250.));
    h_lateralLeakage.push_back(new TH1F(("h_lateralLeakage_" + IntToString(particleenergies[k])).c_str(),";E (GeV);SimHits",2500,0.,250.));
    h_longitudinalLeakage.push_back(new TH1F(("h_longitudinalLeakage_" + IntToString(particleenergies[k])).c_str(),";E (GeV);SimHits",2500,0.,250.));
    h_totalLeakage.push_back(new TH1F(("h_totalLeakage_" + IntToString(particleenergies[k])).c_str(),";E (GeV);SimHits",2500,0.,250.));
    h_lossbeforecalo.push_back(new TH1F(("h_lossbeforecalo_" + IntToString(particleenergies[k])).c_str(),";E (MeV);SimHits",1000,0.,1000.));
    h_lossinpassivelayers.push_back(new TH1F(("h_lossinpassivelayers_" + IntToString(particleenergies[k])).c_str(),";E (GeV);SimHits",520,0.,520.));
    h_lossinsilicon.push_back(new TH1F(("h_lossinsilicon_" + IntToString(particleenergies[k])).c_str(),";E (MeV);SimHits",100,0.,10000.));
    h_totalloss.push_back(new TH1F(("h_totalloss_" + IntToString(particleenergies[k])).c_str(),";E (GeV);SimHits",550,0.,550.));

    ElossUpstreamvsTotalMeasuredMIPs.push_back( (TGraphErrors *) corr->Clone( ("ElossUpstreamvsTotalMeasuredMIPs_" + IntToString(particleenergies[k])).c_str() )  );
    
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
  for (int k=0; k<numberofconfigs; k++){
    results[k]= new TFile(res[k],"recreate");
    std::cout << "Results file " << res[k] << std::endl;

    //======================================================================================
    //Loop on energies
    for (int l=0; l<numberofenergies; l++){
      std::cout << l << std::endl; 
      //======================================================================================
      // Loop on entries
      for (Int_t ievt=0; ievt<lChain[k][l]->GetEntries(); ievt++) {
	// 	// if (ievt != 0) {continue;}
      	lChain[k][l]->GetEntry(ievt);
	if (ievt%(10000)==0) std::cout << "Entry " << ievt << std::endl;
	// if (ievt==1000){break;}
	// std::cout << "Entry " << ievt << std::endl;
	double losseverywherebuttheendandsensors = 0.;
	double lossinsensors = 0.;
	double lossinfirst3sensors = 0.;
	double totalloss = 0.;
	double lossbeforecalo = 0.;
	//Loop on sampling sections
	for(unsigned iL(0); iL<(*ssvec).size(); iL++){
	  // HGCSSSamplingSection lSamSec = (*ssvec)[iL];
	  // double energyabsorber = lSamSec.absorberE();

           //std::cout << "Sampling Section " << iL << " energyabsorber " << (*ssvec)[iL].absorberE()  << std::endl; 
	  
	  //std::cout << " Success " <<(*ssvec)[iL].volNb()  <<  " (*ssvec)[iL].absorberE()  " << (*ssvec)[iL].absorberE() << std::endl; 
	  if ( (iL!=(int)(startlayer-1)) ){
	    h_Leakagefrombehind[l]->Fill( (*ssvec)[iL].rearleakageE()/1000. ); //in GeV
	    h_frontLeakage[l]->Fill( (*ssvec)[iL].frontleakageE()/1000. ); //in GeV
	    h_lateralLeakage[l]->Fill( (*ssvec)[iL].lateralleakageE()/1000. ); //in GeV
	    h_longitudinalLeakage[l]->Fill( ((*ssvec)[iL].rearleakageE()/1000.) + ((*ssvec)[iL].frontleakageE()/1000.) ); //in GeV
	    h_totalLeakage[l]->Fill( (*ssvec)[iL].leakageE()/1000. ); //in GeV
	  }
	  
	  //The loss before calo, that is the 0 sampling section.
	  if ( upstream && (iL == (int)(startlayer-1)) ){
	    lossbeforecalo = (*ssvec)[iL].absorberE();//in MeV
	    h_lossbeforecalo[l]->Fill( (*ssvec)[iL].absorberE() ); //in MeV
	  }
	  // losseverywherebuttheendandsensors = losseverywherebuttheendandsensors + (*ssvec)[iL].passiveE()/1000.;//in GeV
	  lossinsensors = lossinsensors + (*ssvec)[iL].measuredE(); //in MeV
	  if ( ((int)iL)==startlayer || ((int)iL)==(startlayer+1) || ((int)iL)==(startlayer+2) ){  
	    lossinfirst3sensors = lossinfirst3sensors + (*ssvec)[iL].measuredE(); //in MeV
	  }
	  //Total loss including upstream material if any
	  totalloss = totalloss + ((*ssvec)[iL].measuredE()/1000.) + ((*ssvec)[iL].absorberE()/1000.) + (*ssvec)[iL].leakageE()/1000.; 
	  //Below the -1 is in order to fill the histos at the last sampling section of the event
	  if (iL== ((*ssvec).size()-1) ){
	    h_lossinpassivelayers[l]->Fill( losseverywherebuttheendandsensors  );//in GeV
	    h_lossinsilicon[l]->Fill( lossinsensors );//in MeV
	    h_totalloss[l]->Fill( totalloss  );//in GeV
	    // totalloss = 0.;
	    //For the correction in Ereco
	    Int_t np=ElossUpstreamvsTotalMeasuredMIPs[l]->GetN();
	    // std::cout << "np " << np << std::endl;
	    if ( upstream ){
	      ElossUpstreamvsTotalMeasuredMIPs[l]->SetPoint(np, lossinfirst3sensors/(0.081)  , lossbeforecalo  ); //y: ElossUp in MeV, x: TotalMeasuredMIPs infirst 3 sensors in MIPs
	    }
	  }


	} //loop on sampling sections
      
      } //  Loop on entries

      //======================================================================================
      //Write histos 
      h_Leakagefrombehind[l]->Write();
      h_frontLeakage[l]->Write();
      h_lateralLeakage[l]->Write();
      h_longitudinalLeakage[l]->Write();
      h_totalLeakage[l]->Write();
      h_lossinpassivelayers[l]->Write();
      h_lossinsilicon[l]->Write();
      h_totalloss[l]->Write();
      if ( upstream ){
	ElossUpstreamvsTotalMeasuredMIPs[l]->Write();
	h_lossbeforecalo[l]->Write();
      }

      //======================================================================================
      //We should here clear the histograms because we want them empty for the next file. 
      //Reset histos
      h_Leakagefrombehind[l]->Reset();
      h_frontLeakage[l]->Reset();
      h_lateralLeakage[l]->Reset();
      h_longitudinalLeakage[l]->Reset();
      h_totalLeakage[l]->Reset();
      h_lossbeforecalo[l]->Reset();
      h_lossinpassivelayers[l]->Reset();
      h_lossinsilicon[l]->Reset();
      h_totalloss[l]->Reset();
      
      
    }//Loop on energies
   
    results[k]->Close();

  } // Loop on files

  //======================================================================================
  //Make one plot for all different leakage energies
  TCanvas* c1 = new TCanvas("c1", "  ");
  TCanvas* c2 = new TCanvas("c2", "  ");
  TCanvas* c3 = new TCanvas("c3", "  ");
  TCanvas* c4 = new TCanvas("c4", "  ");
  TCanvas* c5 = new TCanvas("c5", "  ");
  TCanvas *c6[numberofenergies];
  for (int l=0; l<numberofenergies; l++){
    c6[l] = new TCanvas( ("c6_"+IntToString(l)) .c_str(), "  ");
  }
  TCanvas* c7 = new TCanvas("c7", "  ");
  TCanvas* c8 = new TCanvas("c8", "  ");
  TCanvas* c9 = new TCanvas("c9", "  ");
  TCanvas* c10 = new TCanvas("c10", "  ");
  TCanvas* c11 = new TCanvas("c11", "  ");
  TCanvas* c12 = new TCanvas("c12", "  ");

  //For the legend
  // TLegend* leg[numberoffiles];
  TLegend* leg1 = new TLegend(0.5, 0.7, 0.8, 0.9);
  leg1->SetHeader("Energy");
  leg1->SetFillColor(17);
  TLegend* leg2 = new TLegend(0.5, 0.7, 0.8, 0.9);
  leg2->SetHeader("Energy");
  leg2->SetFillColor(17);
  TLegend* leg3 = new TLegend(0.5, 0.7, 0.8, 0.9);
  leg3->SetHeader("Energy");
  leg3->SetFillColor(17);
  TLegend* leg4 = new TLegend(0.5, 0.7, 0.8, 0.9);
  leg4->SetHeader("Energy");
  leg4->SetFillColor(17);
  TLegend* leg5 = new TLegend(0.5, 0.7, 0.8, 0.9);
  leg5->SetHeader("Energy");
  leg5->SetFillColor(17);
  TLegend* leg6 = new TLegend(0.5, 0.7, 0.8, 0.9);
  leg6->SetHeader("Energy");
  leg6->SetFillColor(17);
  TLegend* leg7 = new TLegend(0.5, 0.7, 0.8, 0.9);
  leg7->SetHeader("Energy");
  leg7->SetFillColor(17);
  TLegend* leg8 = new TLegend(0.5, 0.7, 0.8, 0.9);
  leg8->SetHeader("Energy");
  leg8->SetFillColor(17);
  TLegend* leg9 = new TLegend(0.5, 0.7, 0.8, 0.9);
  leg9->SetHeader("Energy");
  leg9->SetFillColor(17);

  TString procNa[numberofenergies];
  procNa[0] = "2 GeV";
  procNa[1] = "3 GeV";
  procNa[2] = "5 GeV";
  procNa[3] = "8 GeV";
  procNa[4] = "10 GeV";
  procNa[5] = "15 GeV";
  procNa[6] = "30 GeV";
  procNa[7] = "50 GeV";
  procNa[8] = "80 GeV";
  procNa[9] = "100 GeV";
  procNa[10] = "120 GeV";
  procNa[11] = "150 GeV";
  procNa[12] = "180 GeV";
  procNa[13] = "200 GeV";
  procNa[14] = "250 GeV";
  procNa[15] = "300 GeV";
  procNa[16] = "400 GeV";
  procNa[17] = "500 GeV";


  //======================================================================================
  //The files that we will store the results of this analysis for the combined plot
  TString res_com[numberofconfigs];
  TFile* results_com[numberofconfigs];
  for (int k=0; k<numberofconfigs; k++){
    TString fp = "PFcal_" + IntToString(configs[k]);	
    res_com[k] = fp + "_combinedplots_closure.root";
    //std::cout << res[k] << std::endl;
  }
  TH1F* hist1,* hist2,* hist3,* hist4,* hist5,* hist6,* hist7,* hist8,* hist9;
  TGraphErrors *gr[numberofenergies];
  //For accessing the fit results
  TF1 *fit[numberofenergies];
  TGraphErrors *gr_dEUpstream_totalmeasuredMIPs_fit;
  TGraphErrors *gr_dEUpstream_totalmeasuredMIPs_fit_const;
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  TString titleofplot1,titleofplot2,titleofplot3,titleofplot; 
  for (int k=0; k<numberofconfigs; k++){
    results[k]= TFile::Open(res[k],"read");
    std::cout << "Results file " << res[k] << std::endl;
    
    //======================================================================================
    //Loop on energies
    for (int l=0; l<numberofenergies; l++){
      
      //For the rear leakage
      //------------------------------------------------------------------------------------------------
      hist1 = (TH1F*)results[k]->Get(("h_Leakagefrombehind_" + IntToString(particleenergies[l])).c_str());
      titleofplot1 = "Leakage energy behind calorimeter for "; 
      titleofplot2 =  particle +  " beam particle gun"; 
      titleofplot = titleofplot1 + titleofplot2;
      hist1->SetTitle(titleofplot); 
      hist1->GetXaxis()->SetTitle("Rear Leakage energy (GeV)"); 
      hist1->GetYaxis()->SetTitle("Events/100 MeV");
      c1->SetLogy();
      if ( l == 0){hist1->SetLineColor(4);}
      if ( l == 1){hist1->SetLineColor(2);}
      if ( l == 2){hist1->SetLineColor(1);}
      if ( l == 3){hist1->SetLineColor(3);}
      if ( l == 4){hist1->SetLineColor(5);}
      if ( l == 5){hist1->SetLineColor(6);}
      if ( l == 6){hist1->SetLineColor(7);}
      if ( l == 7){hist1->SetLineColor(8);}
      if ( l == 8){hist1->SetLineColor(9);}
      if ( l == 9){hist1->SetLineColor(12);}
      if ( l == 10){hist1->SetLineColor(46);}
      if ( l == 11){hist1->SetLineColor(41);}
      if ( l == 12){hist1->SetLineColor(38);}
      if ( l == 13){hist1->SetLineColor(40);}
      if ( l == 14){hist1->SetLineColor(30);}
      if ( l == 15){hist1->SetLineColor(24);}
      if ( l == 16){hist1->SetLineColor(29);}
      if ( l == 17){hist1->SetLineColor(49);}
      //hist->SetLineColor(l+1);
      leg1->AddEntry(hist1, procNa[l], "L");
      //hist->GetYaxis()->SetRangeUser(0.,900.);//  21: 100. 22: 900.
      c1->cd();
      c1->Update(); 
      l == 0 ? hist1->Draw("HIST") : hist1->Draw("HISTsame");
      l == 0 ? leg1->Draw() : leg1->Draw("same");
      c1->SaveAs("Leakage_" + particle + ".png"); 
      c1->Update(); 

      //For the front leakage
      //------------------------------------------------------------------------------------------------
      hist6 = (TH1F*)results[k]->Get(("h_frontLeakage_" + IntToString(particleenergies[l])).c_str());
      titleofplot1 = "Leakage energy in front of the calorimeter for "; 
      titleofplot2 =  particle +  " beam particle gun"; 
      titleofplot = titleofplot1 + titleofplot2;
      hist6->SetTitle(titleofplot); 
      hist6->GetXaxis()->SetTitle("Front Leakage energy (GeV)"); 
      hist6->GetYaxis()->SetTitle("Events/100 MeV");
      c9->SetLogy();
      if ( l == 0){hist6->SetLineColor(4);}
      if ( l == 1){hist6->SetLineColor(2);}
      if ( l == 2){hist6->SetLineColor(1);}
      if ( l == 3){hist6->SetLineColor(3);}
      if ( l == 4){hist6->SetLineColor(5);}
      if ( l == 5){hist6->SetLineColor(6);}
      if ( l == 6){hist6->SetLineColor(7);}
      if ( l == 7){hist6->SetLineColor(8);}
      if ( l == 8){hist6->SetLineColor(9);}
      if ( l == 9){hist6->SetLineColor(12);}
      if ( l == 10){hist6->SetLineColor(46);}
      if ( l == 11){hist6->SetLineColor(41);}
      if ( l == 12){hist6->SetLineColor(38);}
      if ( l == 13){hist6->SetLineColor(40);}
      if ( l == 14){hist6->SetLineColor(30);}
      if ( l == 15){hist6->SetLineColor(24);}
      if ( l == 16){hist6->SetLineColor(29);}
      if ( l == 17){hist6->SetLineColor(49);}
     //hist->SetLineColor(l+1);
      leg6->AddEntry(hist6, procNa[l], "L");
      //hist->GetYaxis()->SetRangeUser(0.,900.);//  21: 100. 22: 900.
      c9->cd();
      c9->Update(); 
      l == 0 ? hist6->Draw("HIST") : hist6->Draw("HISTsame");
      l == 0 ? leg6->Draw() : leg6->Draw("same");
      c9->SaveAs("Leakage_" + particle + ".png"); 
      c9->Update(); 

      //For the lateral leakage
      //------------------------------------------------------------------------------------------------
      hist7 = (TH1F*)results[k]->Get(("h_lateralLeakage_" + IntToString(particleenergies[l])).c_str());
      titleofplot1 = "Lateral leakage energy of the calorimeter for "; 
      titleofplot2 =  particle +  " beam particle gun"; 
      titleofplot = titleofplot1 + titleofplot2;
      hist7->SetTitle(titleofplot); 
      hist7->GetXaxis()->SetTitle("Lateral Leakage energy (GeV)"); 
      hist7->GetYaxis()->SetTitle("Events/100 MeV");
      c10->SetLogy();
      if ( l == 0){hist7->SetLineColor(4);}
      if ( l == 1){hist7->SetLineColor(2);}
      if ( l == 2){hist7->SetLineColor(1);}
      if ( l == 3){hist7->SetLineColor(3);}
      if ( l == 4){hist7->SetLineColor(5);}
      if ( l == 5){hist7->SetLineColor(6);}
      if ( l == 6){hist7->SetLineColor(7);}
      if ( l == 7){hist7->SetLineColor(8);}
      if ( l == 8){hist7->SetLineColor(9);}
      if ( l == 9){hist7->SetLineColor(12);}
      if ( l == 10){hist7->SetLineColor(46);}
      if ( l == 11){hist7->SetLineColor(41);}
      if ( l == 12){hist7->SetLineColor(38);}
      if ( l == 13){hist7->SetLineColor(40);}
      if ( l == 14){hist7->SetLineColor(30);}
      if ( l == 15){hist7->SetLineColor(24);}
      if ( l == 16){hist7->SetLineColor(29);}
      if ( l == 17){hist7->SetLineColor(49);}
      //hist->SetLineColor(l+1);
      leg7->AddEntry(hist7, procNa[l], "L");
      //hist->GetYaxis()->SetRangeUser(0.,900.);//  21: 100. 22: 900.
      c10->cd();
      c10->Update(); 
      l == 0 ? hist7->Draw("HIST") : hist7->Draw("HISTsame");
      l == 0 ? leg7->Draw() : leg7->Draw("same");
      c10->SaveAs("Leakage_" + particle + ".png"); 
      c10->Update(); 

      //For the longitudinal leakage
      //------------------------------------------------------------------------------------------------
      hist9 = (TH1F*)results[k]->Get(("h_longitudinalLeakage_" + IntToString(particleenergies[l])).c_str());
      titleofplot1 = "Longitudinal leakage energy of the calorimeter for "; 
      titleofplot2 =  particle +  " beam particle gun"; 
      titleofplot = titleofplot1 + titleofplot2;
      hist9->SetTitle(titleofplot); 
      hist9->GetXaxis()->SetTitle("Longitudinal Leakage energy (GeV)"); 
      hist9->GetYaxis()->SetTitle("Events/100 MeV");
      c12->SetLogy();
      if ( l == 0){hist9->SetLineColor(4);}
      if ( l == 1){hist9->SetLineColor(2);}
      if ( l == 2){hist9->SetLineColor(1);}
      if ( l == 3){hist9->SetLineColor(3);}
      if ( l == 4){hist9->SetLineColor(5);}
      if ( l == 5){hist9->SetLineColor(6);}
      if ( l == 6){hist9->SetLineColor(7);}
      if ( l == 7){hist9->SetLineColor(8);}
      if ( l == 8){hist9->SetLineColor(9);}
      if ( l == 9){hist9->SetLineColor(12);}
      if ( l == 10){hist9->SetLineColor(46);}
      if ( l == 11){hist9->SetLineColor(41);}
      if ( l == 12){hist9->SetLineColor(38);}
      if ( l == 13){hist9->SetLineColor(40);}
      if ( l == 14){hist9->SetLineColor(30);}
      if ( l == 15){hist9->SetLineColor(24);}
      if ( l == 16){hist9->SetLineColor(29);}
      if ( l == 17){hist9->SetLineColor(49);}
      //hist->SetLineColor(l+1);
      leg9->AddEntry(hist9, procNa[l], "L");
      //hist->GetYaxis()->SetRangeUser(0.,900.);//  21: 100. 22: 900.
      c12->cd();
      c12->Update(); 
      l == 0 ? hist9->Draw("HIST") : hist9->Draw("HISTsame");
      l == 0 ? leg9->Draw() : leg9->Draw("same");
      c12->SaveAs("Leakage_" + particle + ".png"); 
      c12->Update(); 

      //For the total leakage
      //------------------------------------------------------------------------------------------------
      hist8 = (TH1F*)results[k]->Get(("h_totalLeakage_" + IntToString(particleenergies[l])).c_str());
      titleofplot1 = "Total leakage energy of the calorimeter for "; 
      titleofplot2 =  particle +  " beam particle gun"; 
      titleofplot = titleofplot1 + titleofplot2;
      hist8->SetTitle(titleofplot); 
      hist8->GetXaxis()->SetTitle("Total Leakage energy (GeV)"); 
      hist8->GetYaxis()->SetTitle("Events/100 MeV");
      c11->SetLogy();
      if ( l == 0){hist8->SetLineColor(4);}
      if ( l == 1){hist8->SetLineColor(2);}
      if ( l == 2){hist8->SetLineColor(1);}
      if ( l == 3){hist8->SetLineColor(3);}
      if ( l == 4){hist8->SetLineColor(5);}
      if ( l == 5){hist8->SetLineColor(6);}
      if ( l == 6){hist8->SetLineColor(7);}
      if ( l == 7){hist8->SetLineColor(8);}
      if ( l == 8){hist8->SetLineColor(9);}
      if ( l == 9){hist8->SetLineColor(12);}
      if ( l == 10){hist8->SetLineColor(46);}
      if ( l == 11){hist8->SetLineColor(41);}
      if ( l == 12){hist8->SetLineColor(38);}
      if ( l == 13){hist8->SetLineColor(40);}
      if ( l == 14){hist8->SetLineColor(30);}
      if ( l == 15){hist8->SetLineColor(24);}
      if ( l == 16){hist8->SetLineColor(29);}
      if ( l == 17){hist8->SetLineColor(49);}
      //hist->SetLineColor(l+1);
      leg8->AddEntry(hist8, procNa[l], "L");
      //hist->GetYaxis()->SetRangeUser(0.,900.);//  21: 100. 22: 900.
      c11->cd();
      c11->Update(); 
      l == 0 ? hist8->Draw("HIST") : hist8->Draw("HISTsame");
      l == 0 ? leg8->Draw() : leg8->Draw("same");
      c11->SaveAs("Leakage_" + particle + ".png"); 
      c11->Update(); 

      //For the loss in sensors
      //------------------------------------------------------------------------------------------------
      c2->cd();
      hist2 = (TH1F*)results[k]->Get(("h_lossinsilicon_" + IntToString(particleenergies[l])).c_str());
      titleofplot1 = "Energy loss in silicon sensors for "; 
      titleofplot2 =  particle +  " beam particle gun"; 
      titleofplot = titleofplot1 + titleofplot2;
      hist2->SetTitle(titleofplot); 
      hist2->GetXaxis()->SetTitle("Energy loss in silicon sensors (MeV)"); 
      hist2->GetYaxis()->SetTitle("Events/100 MeV");
      c2->SetLogy();
      if ( l == 0){hist2->SetLineColor(4);}
      if ( l == 1){hist2->SetLineColor(2);}
      if ( l == 2){hist2->SetLineColor(1);}
      if ( l == 3){hist2->SetLineColor(3);}
      if ( l == 4){hist2->SetLineColor(5);}
      if ( l == 5){hist2->SetLineColor(6);}
      if ( l == 6){hist2->SetLineColor(7);}
      if ( l == 7){hist2->SetLineColor(8);}
      if ( l == 8){hist2->SetLineColor(9);}
      if ( l == 9){hist2->SetLineColor(12);}
      if ( l == 10){hist2->SetLineColor(46);}
      if ( l == 11){hist2->SetLineColor(41);}
      if ( l == 12){hist2->SetLineColor(38);}
      if ( l == 13){hist2->SetLineColor(40);}
      if ( l == 14){hist2->SetLineColor(30);}
      if ( l == 15){hist2->SetLineColor(24);}
      if ( l == 16){hist2->SetLineColor(29);}
      if ( l == 17){hist2->SetLineColor(49);}
      //hist->SetLineColor(l+1);
      leg2->AddEntry(hist2, procNa[l], "L");
      //hist->GetYaxis()->SetRangeUser(0.,900.);//  21: 100. 22: 900.
      c2->Update(); 
      l == 0 ? hist2->Draw("HIST") : hist2->Draw("HISTsame");
      l == 0 ? leg2->Draw() : leg2->Draw("same");
      c2->SaveAs("EnergyLossInSilicon_" + particle + ".png"); 
      c2->Update(); 

      //For the loss before calo
      //------------------------------------------------------------------------------------------------
      if(upstream){
	c3->cd();
	hist3 = (TH1F*)results[k]->Get(("h_lossbeforecalo_" + IntToString(particleenergies[l])).c_str());
	titleofplot1 = "Energy loss before calorimeter for "; 
	titleofplot2 =  particle +  " beam particle gun"; 
	titleofplot = titleofplot1 + titleofplot2;
	hist3->SetTitle(titleofplot); 
	hist3->GetXaxis()->SetTitle("Energy loss before calorimeter (MeV)"); 
	hist3->GetYaxis()->SetTitle("Events/1 MeV");
	c3->SetLogy();
	if ( l == 0){hist3->SetLineColor(4);}
	if ( l == 1){hist3->SetLineColor(2);}
	if ( l == 2){hist3->SetLineColor(1);}
	if ( l == 3){hist3->SetLineColor(3);}
	if ( l == 4){hist3->SetLineColor(5);}
	if ( l == 5){hist3->SetLineColor(6);}
	if ( l == 6){hist3->SetLineColor(7);}
	if ( l == 7){hist3->SetLineColor(8);}
	if ( l == 8){hist3->SetLineColor(9);}
	if ( l == 9){hist3->SetLineColor(12);}
	if ( l == 10){hist3->SetLineColor(46);}
	if ( l == 11){hist3->SetLineColor(41);}
	if ( l == 12){hist3->SetLineColor(38);}
	if ( l == 13){hist3->SetLineColor(40);}
	if ( l == 14){hist3->SetLineColor(30);}
	if ( l == 15){hist3->SetLineColor(24);}
	if ( l == 16){hist3->SetLineColor(29);}
	if ( l == 17){hist3->SetLineColor(49);}
 	//hist->SetLineColor(l+1);
	leg3->AddEntry(hist3, procNa[l], "L");
	//hist->GetYaxis()->SetRangeUser(0.,900.);//  21: 100. 22: 900.
	c3->Update(); 
	l == 0 ? hist3->Draw("HIST") : hist3->Draw("HISTsame");
	l == 0 ? leg3->Draw() : leg3->Draw("same");
	c3->SaveAs("EnergyLossBeforeCalo_" + particle + ".png"); 
	c3->Update(); 
      }
      //For the loss in passive layers
      //------------------------------------------------------------------------------------------------
      c4->cd();
      hist4 = (TH1F*)results[k]->Get(("h_lossinpassivelayers_" + IntToString(particleenergies[l])).c_str());
      titleofplot1 = "Energy loss in passive layers for "; 
      titleofplot2 =  particle +  " beam particle gun"; 
      titleofplot = titleofplot1 + titleofplot2;
      hist4->SetTitle(titleofplot); 
      hist4->GetXaxis()->SetTitle("Energy loss in passive layers (GeV)"); 
      hist4->GetYaxis()->SetTitle("Events/1 GeV");
      c4->SetLogy();
      if ( l == 0){hist4->SetLineColor(4);}
      if ( l == 1){hist4->SetLineColor(2);}
      if ( l == 2){hist4->SetLineColor(1);}
      if ( l == 3){hist4->SetLineColor(3);}
      if ( l == 4){hist4->SetLineColor(5);}
      if ( l == 5){hist4->SetLineColor(6);}
      if ( l == 6){hist4->SetLineColor(7);}
      if ( l == 7){hist4->SetLineColor(8);}
      if ( l == 8){hist4->SetLineColor(9);}
      if ( l == 9){hist4->SetLineColor(12);}
      if ( l == 10){hist4->SetLineColor(46);}
      if ( l == 11){hist4->SetLineColor(41);}
      if ( l == 12){hist4->SetLineColor(38);}
      if ( l == 13){hist4->SetLineColor(40);}
      if ( l == 14){hist4->SetLineColor(30);}
      if ( l == 15){hist4->SetLineColor(24);}
      if ( l == 16){hist4->SetLineColor(29);}
      if ( l == 17){hist4->SetLineColor(49);}
      //hist->SetLineColor(l+1);
      leg4->AddEntry(hist4, procNa[l], "L");
      //hist->GetYaxis()->SetRangeUser(0.,900.);//  21: 100. 22: 900.
      c4->Update(); 
      l == 0 ? hist4->Draw("HIST") : hist4->Draw("HISTsame");
      l == 0 ? leg4->Draw() : leg4->Draw("same");
      c4->SaveAs("EnergyLossInPassiveLayers_" + particle + ".png"); 
      c4->Update(); 
      //For the total loss
      //------------------------------------------------------------------------------------------------
      c5->cd();
      hist5 = (TH1F*)results[k]->Get(("h_totalloss_" + IntToString(particleenergies[l])).c_str());
      titleofplot1 = "Total Energy loss for "; 
      titleofplot2 =  particle +  " beam particle gun"; 
      titleofplot = titleofplot1 + titleofplot2;
      hist5->SetTitle(titleofplot); 
      hist5->GetXaxis()->SetTitle("Total Energy loss (GeV)"); 
      hist5->GetYaxis()->SetTitle("Events/1 GeV");
      c5->SetLogy();
      if ( l == 0){hist5->SetLineColor(4);}
      if ( l == 1){hist5->SetLineColor(2);}
      if ( l == 2){hist5->SetLineColor(1);}
      if ( l == 3){hist5->SetLineColor(3);}
      if ( l == 4){hist5->SetLineColor(5);}
      if ( l == 5){hist5->SetLineColor(6);}
      if ( l == 6){hist5->SetLineColor(7);}
      if ( l == 7){hist5->SetLineColor(8);}
      if ( l == 8){hist5->SetLineColor(9);}
      if ( l == 9){hist5->SetLineColor(12);}
      if ( l == 10){hist5->SetLineColor(46);}
      if ( l == 11){hist5->SetLineColor(41);}
      if ( l == 12){hist5->SetLineColor(38);}
      if ( l == 13){hist5->SetLineColor(40);}
      if ( l == 14){hist5->SetLineColor(30);}
      if ( l == 15){hist5->SetLineColor(24);}
      if ( l == 16){hist5->SetLineColor(29);}
      if ( l == 17){hist5->SetLineColor(49);}
        //hist->SetLineColor(l+1);
      leg5->AddEntry(hist5, procNa[l], "L");
      //hist->GetYaxis()->SetRangeUser(0.,900.);//  21: 100. 22: 900.
      c5->Update(); 
      l == 0 ? hist5->Draw("HIST") : hist5->Draw("HISTsame");
      l == 0 ? leg5->Draw() : leg5->Draw("same");
      c5->SaveAs("TotalEnergyLoss_" + particle + ".png"); 
      c5->Update(); 

      if(upstream){
	//For the fit 
	//------------------------------------------------------------------------------------------------
	gr[l] = (TGraphErrors*)results[k]->Get( ("ElossUpstreamvsTotalMeasuredMIPs_" + IntToString(particleenergies[l])).c_str()  );
	gr[l]->GetXaxis()->SetRangeUser(0.,100000.);
	gr[l]->Fit("pol1");
	
	// prof[(int) iL] = hist[(int) iL]->ProfileX();
	// prof[(int) iL]->Fit("pol1");
	fit[l] = gr[l]->GetFunction("pol1");
	//Results of the fit
	std::cout << "First parameter: " << fit[l]->GetParameter(0) << " Second parameter: " << fit[l]->GetParameter(1) << std::endl;
	
	titleofplot1 = "Energy lost in upstream materials before calorimeter vs measured energy in first 3 layers for "; 
	titleofplot2 =  particle +  " "; 
	titleofplot3 =  IntToString(particleenergies[l]) +  " GeV beam particle gun"; 
	titleofplot = titleofplot1 + titleofplot2 + titleofplot3;
	gr[l]->SetTitle(titleofplot); 
	gr[l]->GetHistogram()->SetXTitle("E_{active} in first 3 layers (MIPs)"); 
	gr[l]->GetHistogram()->SetYTitle("dE_{lossUpStream} (MeV)");
	if ( l == 0){fit[l]->SetLineColor(4);gr[l]->SetMarkerColor(4);gr[l]->SetMarkerStyle(20);}
	if ( l == 1){fit[l]->SetLineColor(2);gr[l]->SetMarkerColor(2);gr[l]->SetMarkerStyle(21);}
	if ( l == 2){fit[l]->SetLineColor(1);gr[l]->SetMarkerColor(1);gr[l]->SetMarkerStyle(22);}
	if ( l == 3){fit[l]->SetLineColor(3);gr[l]->SetMarkerColor(3);gr[l]->SetMarkerStyle(23);}
	if ( l == 4){fit[l]->SetLineColor(5);gr[l]->SetMarkerColor(5);gr[l]->SetMarkerStyle(28);}
	if ( l == 5){fit[l]->SetLineColor(6);gr[l]->SetMarkerColor(6);gr[l]->SetMarkerStyle(29);}
	if ( l == 6){fit[l]->SetLineColor(7);gr[l]->SetMarkerColor(7);gr[l]->SetMarkerStyle(30);}
	if ( l == 7){fit[l]->SetLineColor(8);gr[l]->SetMarkerColor(8);gr[l]->SetMarkerStyle(31);}
	if ( l == 8){fit[l]->SetLineColor(9);gr[l]->SetMarkerColor(9);gr[l]->SetMarkerStyle(33);}
	if ( l == 9){fit[l]->SetLineColor(12);gr[l]->SetMarkerColor(12);gr[l]->SetMarkerStyle(34);}
	if ( l == 10){fit[l]->SetLineColor(46);gr[l]->SetMarkerColor(46);gr[l]->SetMarkerStyle(25);}
	if ( l == 11){fit[l]->SetLineColor(41);gr[l]->SetMarkerColor(41);gr[l]->SetMarkerStyle(26);}
	if ( l == 12){fit[l]->SetLineColor(38);gr[l]->SetMarkerColor(38);gr[l]->SetMarkerStyle(2);}
	if ( l == 13){fit[l]->SetLineColor(40);gr[l]->SetMarkerColor(40);gr[l]->SetMarkerStyle(3);}
	if ( l == 14){fit[l]->SetLineColor(30);gr[l]->SetMarkerColor(30);gr[l]->SetMarkerStyle(5);}
	if ( l == 15){fit[l]->SetLineColor(24);gr[l]->SetMarkerColor(24);gr[l]->SetMarkerStyle(12);}
	if ( l == 16){fit[l]->SetLineColor(29);gr[l]->SetMarkerColor(29);gr[l]->SetMarkerStyle(24);}
	if ( l == 17){fit[l]->SetLineColor(49);gr[l]->SetMarkerColor(49);gr[l]->SetMarkerStyle(32);}
 	//hist[(int) ifit[l]->SetLineColor(4);L]->SetLineColor(l+1);
	// leg_devsmipswithfit[(int) iL]->AddEntry(hist[(int) iL], procNa[l], "P");
	// leg_devsmipswithfit[(int) iL]->AddEntry(gr[(int) iL], procNa[l], "P");
	//hist[(int) iL]->GetYaxis()->SetRangeUser(0.,900.);//  21: 100. 22: 900.
	c6[l]->cd();
	c6[l]->Update(); 
	gr[l]->Draw("AP");
	// l == 0 ? gr[(int) iL]->Draw("AP") : gr[(int) iL]->Draw("APS");
	// leg_devsmipswithfit[(int) iL]->Draw("PS");
	// l == 0 ? leg_devsmipswithfit[(int) iL]->Draw() : leg_devsmipswithfit[(int) iL]->Draw("same");
	// l == 0 ? fit[(int) iL]->Draw() : fit[(int) iL]->Draw("same");
	gr[l]->GetHistogram()->GetXaxis()->SetRangeUser(0.,10000.);

	char buf[500];
	TLatex lat;
	double latx = gr[l]->GetHistogram()->GetXaxis()->GetXmin()+(gr[l]->GetHistogram()->GetXaxis()->GetXmax()-gr[l]->GetHistogram()->GetXaxis()->GetXmin())/20.;
	double laty = gr[l]->GetHistogram()->GetMaximum();
	sprintf(buf,"dE_upstream = p0 + p1 * E_first3layers");
	lat.DrawLatex(latx,laty*0.8,buf);
	sprintf(buf,"p0 = %3.3f +/- %3.3f",fit[l]->GetParameter(0),fit[l]->GetParError(0));
	lat.DrawLatex(latx,laty*0.7,buf);
	sprintf(buf,"p1 = %3.3f +/- %3.3f",fit[l]->GetParameter(1),fit[l]->GetParError(1));
	lat.DrawLatex(latx,laty*0.6,buf);
	sprintf(buf,"#chi^{2}/N = %3.3f/%d = %3.3f",fit[l]->GetChisquare(),fit[l]->GetNDF(),fit[l]->GetChisquare()/fit[l]->GetNDF());
	lat.DrawLatex(latx,laty*0.5,buf);
 

	c6[l]->Update(); 
	TString cantopng1 = "dEUpstreamVStotalmeasuredMIPs_plusfit_"+ particle;
	TString cantopng2 = "_" + IntToString(particleenergies[l]);
	TString cantopng3 = "GeV";
	TString cantopng = cantopng1 + cantopng2 + cantopng3 + ".png";
	
	c6[l]->SaveAs(cantopng);
      }

      
    }//Loop on energies
    
    results_com[k]= new TFile(res_com[k],"recreate");

    //======================================================================================
    //Write canvas with combined plot 
    // c1->Print("Leakage.pdf",".pdf");
    c1->Write();
    c2->Write();
    c4->Write();
    c5->Write();
    if(upstream){
      c3->Write();
      for (int l=0; l<numberofenergies; l++){
	c6[l]->Write();
      }
    }
    c9->Write();
    c10->Write();
    c11->Write();
    c12->Write();


    //Reset histos
    // hist->Reset();

    if(upstream){
      //Here we will make the plot for the fit results
      c7->cd();
      gr_dEUpstream_totalmeasuredMIPs_fit = (TGraphErrors *) corr->Clone( "ElossUpstreamvsTotalMeasuredMIPs_fit");
      gr_dEUpstream_totalmeasuredMIPs_fit->SetMarkerStyle( 20 ); 
      // gr_dEUpstream_totalmeasuredMIPs_fit->SetMarkerSize(0.2); 
      gr_dEUpstream_totalmeasuredMIPs_fit->SetMarkerColor(1);  
      gr_dEUpstream_totalmeasuredMIPs_fit->SetLineColor(1);
      c7->Update(); 

      // Int_t np=gr_dEUpstream_totalmeasuredMIPs_fit->GetN();
      // gr_dEUpstream_totalmeasuredMIPs_fit->SetPoint(numberofenergies, particleenergies[l] , fit[l]->GetParameter(0) );
      // gr_dEUpstream_totalmeasuredMIPs_fit->SetPointError(numberofenergies, 0. , fit[l]->GetParError(0) );
      for (int l=0; l<numberofenergies; l++){
	gr_dEUpstream_totalmeasuredMIPs_fit->SetPoint( l , particleenergies[l] , fit[l]->GetParameter(1) );
	gr_dEUpstream_totalmeasuredMIPs_fit->SetPointError(l , 0. , fit[l]->GetParError(1) );
      }
      gr_dEUpstream_totalmeasuredMIPs_fit->Draw("APL");
      c7->Update(); 

      gr_dEUpstream_totalmeasuredMIPs_fit->SetTitle( "Energy lost Upstream before calorimeter vs. energy measured in first 3 layers fit results for " + particle + " for all beam energies" );
      gr_dEUpstream_totalmeasuredMIPs_fit->GetHistogram()->SetXTitle("Beam Energy (GeV)"); 
      gr_dEUpstream_totalmeasuredMIPs_fit->GetHistogram()->SetYTitle("dE_{lossUpStream} = f(E_{active,3layers}): Slope (MeV/MIP) ");
      gr_dEUpstream_totalmeasuredMIPs_fit->GetXaxis()->SetRangeUser(0.,520.);  
      gr_dEUpstream_totalmeasuredMIPs_fit->GetYaxis()->SetRangeUser(-15.,15.);  

      c7->Update(); 
      TString cnpng = "dEUpstreamVStotalmeasuredMIPs_plusfit_slope"+ particle + ".png";
      c7->SaveAs(cnpng);

      c7->Write();

      //Here we will make the plot for the fit results constant term
      c8->cd();
      gr_dEUpstream_totalmeasuredMIPs_fit_const = (TGraphErrors *) corr->Clone( "ElossUpstreamvsTotalMeasuredMIPs_fit_const");
      gr_dEUpstream_totalmeasuredMIPs_fit_const->SetMarkerStyle( 20 ); 
      // gr_dEUpstream_totalmeasuredMIPs_fit_const->SetMarkerSize(0.2); 
      gr_dEUpstream_totalmeasuredMIPs_fit_const->SetMarkerColor(1);  
      gr_dEUpstream_totalmeasuredMIPs_fit_const->SetLineColor(1);
      c8->Update(); 

      // Int_t np=gr_dEUpstream_totalmeasuredMIPs_fit_const->GetN();
      // gr_dEUpstream_totalmeasuredMIPs_fit_const->SetPoint(numberofenergies, particleenergies[l] , fit[l]->GetParameter(0) );
      // gr_dEUpstream_totalmeasuredMIPs_fit_const->SetPointError(numberofenergies, 0. , fit[l]->GetParError(0) );
      for (int l=0; l<numberofenergies; l++){
	gr_dEUpstream_totalmeasuredMIPs_fit_const->SetPoint( l , particleenergies[l] , fit[l]->GetParameter(0) );
	gr_dEUpstream_totalmeasuredMIPs_fit_const->SetPointError(l , 0. , fit[l]->GetParError(0) );
      }
      gr_dEUpstream_totalmeasuredMIPs_fit_const->Draw("APL");
      c8->Update(); 

      gr_dEUpstream_totalmeasuredMIPs_fit_const->SetTitle( "Energy lost Upstream before calorimeter vs. energy measured in first 3 layers fit results for " + particle + " for all beam energies" );
      gr_dEUpstream_totalmeasuredMIPs_fit_const->GetHistogram()->SetXTitle("Beam Energy (GeV)"); 
      gr_dEUpstream_totalmeasuredMIPs_fit_const->GetHistogram()->SetYTitle("dE_{lossUpStream} = f(E_{active,3layers}): Constant Term (MeV)");
      gr_dEUpstream_totalmeasuredMIPs_fit_const->GetXaxis()->SetRangeUser(0.,520.);  
      gr_dEUpstream_totalmeasuredMIPs_fit_const->GetYaxis()->SetRangeUser(-500.,500.);  

      c8->Update(); 
      TString cnpng1 = "dEUpstreamVStotalmeasuredMIPs_plusfit_const"+ particle + ".png";
      c8->SaveAs(cnpng1);
      c8->Write();

    }

    results_com[k]->Close();

    std::ofstream myfile;
    if(upstream){

      TString titleoffile = "dE_lossUpvsMIPs_" + particle + ".txt";
      myfile.open(titleoffile);
      
      //Here print the w0 and w1 that are going to be used for Ereco_corrected
      for (int l=0; l<numberofenergies; l++){
	std::cout << "---------------------------------------------------------" << std::endl;
	std::cout << "Particle " <<  particle << std::endl;
	std::cout << "Energy of the beam " <<  particleenergies[l] << std::endl;
	std::cout << "LossUpstreamvsmeasuredMIPs slope: " <<  fit[l]->GetParameter(1) << " error slope " << fit[l]->GetParError(1) << std::endl;
	std::cout << "LossUpstreamvsmeasuredMIPs constant: " <<  fit[l]->GetParameter(0) << " error constant " << fit[l]->GetParError(0) << std::endl;
      
	//particle particleenergy constant slope errorconstant errorslope
	myfile << particle << " " << particleenergies[l] << " " << fit[l]->GetParameter(0) << " " << fit[l]->GetParameter(1) << " " << fit[l]->GetParError(0) <<  " " << fit[l]->GetParError(1) << "\n";
      
      }

    }

    myfile.close();

  

    











 
  }//Loop on configs




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
