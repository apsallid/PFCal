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
#include "TGaxis.h"
#include "THStack.h"

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
  //Important boolean parameters related to the running of the code
  bool step1 = true; 
  bool step2 = true; 
  bool dorunfit = true; 
  bool dothe3dplot = true; 
  bool dothedEvsXosplot = true; 
  bool looponlayersforsffit = true; 
  //======================================================================================
  //Do the files include upstream material?
  bool upstream = false;
  int startlayer = upstream ? 1 : 0;
  //======================================================================================
  //Old shower depth: SF = f(X0max/<X0max>)
  //New shower depth: SF = f( Sum(Ei*X0_i)/Sum(E_i)  / <Sum(Ei*X0_i)/Sum(E_i)>)
  bool usenewshowerdepth = true;
  TString shdp = usenewshowerdepth ? " #frac{#sum E X_{0}/#sum E}{<#sum E X_{0}/#sum E>}" : "X_{0,max}/<X_{0,max}>";
  //======================================================================================
  //The root files with the individual plots for all test beam configs 
  const int numberoffiles = 1;
  //======================================================================================
  //3 layers times 4 Si Pads
  // const int numberofpads = 12; 
  //======================================================================================
  //Particle type and energy
  TString particle = "e-"; //mu+
  const int numberofenergies = 18; //15 30 50 80 100 120 150 180 200 250 300 500
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
	if (upstream){
	  lChain[k][l]->Add("/afs/cern.ch/work/a/apsallid/CMS/Geant4/PFCal/PFCalEE/analysis/data/"+filename[k][l][j]+".root");
	} else {
	  // lChain[k][l]->Add("/afs/cern.ch/work/a/apsallid/CMS/Geant4/PFCal/PFCalEE/analysis/data/NoUpStream"+filename[k][l][j]+".root");
	  lChain[k][l]->Add("/tmp/apsallid/"+filename[k][l][j]+".root");
	}
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
  std::map< int, std::vector<TH2F*> > h_MIPSvsdEplusMIPs;
  h_MIPSvsdEplusMIPs.clear();
  //Shower depth
  std::map< int, std::vector<TH1F*> > h_Showerprofile;
  h_Showerprofile.clear();
  //Sampling fractions of layers vs shower depth
  std::map< int, std::vector<TH2F*> > h_sfvsshowerdepthprofile;
  h_sfvsshowerdepthprofile.clear();
  //Shower depth true
  std::vector<TH1F*> h_ShowerDepth_reco;
  h_ShowerDepth_reco.clear();
  //New Shower depth 
  std::vector<TH1F*> h_NewShowerDepth;
  h_NewShowerDepth.clear();
  //Total sampling fraction
  std::vector<TH1F*> h_TotalSamplingFraction;
  h_TotalSamplingFraction.clear();
  //------------------------------------------------------
  //3 sections
  //Total sampling fraction sec 1
  std::vector<TH1F*> h_TotalSamplingFraction_sec1;
  h_TotalSamplingFraction_sec1.clear();
  //Total sampling fraction sec 2
  std::vector<TH1F*> h_TotalSamplingFraction_sec2;
  h_TotalSamplingFraction_sec2.clear();
  //Total sampling fraction sec 3
  std::vector<TH1F*> h_TotalSamplingFraction_sec3;
  h_TotalSamplingFraction_sec3.clear();
  //Total sampling fraction for bins of shower depth
  std::vector<TH2F*> Totalsfvsshowerdepth_forprofile;
  Totalsfvsshowerdepth_forprofile.clear();
  //Total sampling fraction for bins of shower depth sec1
  std::vector<TH2F*> Totalsf_sec1vsshowerdepth_forprofile;
  Totalsf_sec1vsshowerdepth_forprofile.clear();
  //Total sampling fraction for bins of shower depth sec2
  std::vector<TH2F*> Totalsf_sec2vsshowerdepth_forprofile;
  Totalsf_sec2vsshowerdepth_forprofile.clear();
  //Total sampling fraction for bins of shower depth sec3
  std::vector<TH2F*> Totalsf_sec3vsshowerdepth_forprofile;
  Totalsf_sec3vsshowerdepth_forprofile.clear();
  //------------------------------------------------------
  //5 sections plus first layer = 6 sections
  //Total sampling fraction dedxsec 1
  std::vector<TH1F*> h_TotalSamplingFraction_dedxsec1;
  h_TotalSamplingFraction_dedxsec1.clear();
  //Total sampling fraction dedxsec 2
  std::vector<TH1F*> h_TotalSamplingFraction_dedxsec2;
  h_TotalSamplingFraction_dedxsec2.clear();
  //Total sampling fraction dedxsec 3
  std::vector<TH1F*> h_TotalSamplingFraction_dedxsec3;
  h_TotalSamplingFraction_dedxsec3.clear();
  //Total sampling fraction dedxsec 4
  std::vector<TH1F*> h_TotalSamplingFraction_dedxsec4;
  h_TotalSamplingFraction_dedxsec4.clear();
  //Total sampling fraction dedxsec 5
  std::vector<TH1F*> h_TotalSamplingFraction_dedxsec5;
  h_TotalSamplingFraction_dedxsec5.clear();
  //Total sampling fraction dedxsec 6
  std::vector<TH1F*> h_TotalSamplingFraction_dedxsec6;
  h_TotalSamplingFraction_dedxsec6.clear();
  //Total sampling fraction for bins of shower depth dedxsec1
  std::vector<TH2F*> Totalsf_dedxsec1vsshowerdepth_forprofile;
  Totalsf_dedxsec1vsshowerdepth_forprofile.clear();
  //Total sampling fraction for bins of shower depth dedxsec2
  std::vector<TH2F*> Totalsf_dedxsec2vsshowerdepth_forprofile;
  Totalsf_dedxsec2vsshowerdepth_forprofile.clear();
  //Total sampling fraction for bins of shower depth dedxsec3
  std::vector<TH2F*> Totalsf_dedxsec3vsshowerdepth_forprofile;
  Totalsf_dedxsec3vsshowerdepth_forprofile.clear();
  //Total sampling fraction for bins of shower depth dedxsec4
  std::vector<TH2F*> Totalsf_dedxsec4vsshowerdepth_forprofile;
  Totalsf_dedxsec4vsshowerdepth_forprofile.clear();
  //Total sampling fraction for bins of shower depth dedxsec5
  std::vector<TH2F*> Totalsf_dedxsec5vsshowerdepth_forprofile;
  Totalsf_dedxsec5vsshowerdepth_forprofile.clear();
  //Total sampling fraction for bins of shower depth dedxsec6
  std::vector<TH2F*> Totalsf_dedxsec6vsshowerdepth_forprofile;
  Totalsf_dedxsec6vsshowerdepth_forprofile.clear();
  //------------------------------------------------------
  //Total sampling fraction as a function of the shower depth
  std::vector<TGraphErrors*> Totalsfvsshowerdepth;
  Totalsfvsshowerdepth.clear();
  //Total sampling fraction as a function of the shower depth layer unit
  std::vector<TGraphErrors*> Totalsfvsshowerdepthlay;
  Totalsfvsshowerdepthlay.clear();
  //Sampling fraction as a function of shower depth and layer
  std::vector<TH3F*> sfvsshowerdepthprofileandlayer;
  sfvsshowerdepthprofileandlayer.clear();
  //Sampling fraction as a function of shower depth with layer as the weight
  std::vector<TH2F*> sfvsshowerdepthprofile_layerasweight;
  sfvsshowerdepthprofile_layerasweight.clear();
  //Sampling fraction vs normalized shower depth with shower max
  std::vector<TH2F*> h_sfvsnormshowerdepth;
  h_sfvsnormshowerdepth.clear();
  //Energy loss
  std::map< int, std::vector<TH1F*> > h_loss;
  h_loss.clear();
  //Average energy lost as a function of the transversed thickness
  std::vector<TGraphErrors*> dEvsXos;
  dEvsXos.clear();
  //======================================================================================
  // For the fit
  std::map< int, std::vector<TGraphErrors*> > MIPSvsdEplusMIPs_calib;
  MIPSvsdEplusMIPs_calib.clear();

  // int nmbsh = 300; 
  double maxbsh = usenewshowerdepth ? 3. : 10.; 
  // double maxbsh = 10.; 
  // double maxbsh = 10.; 
  lChain[0][0]->GetEntry(0);
  TGraphErrors *calib = new TGraphErrors();
  for(unsigned iL(0); iL<(*ssvec).size(); iL++){
    for (Int_t k=0; k<numberofenergies; k++){
      std::string auxnam1_sf2 = "h_MIPSvsdEplusMIPs_" + IntToString((int) iL);
      std::string auxnam1_tgr2 = "MIPSvsdEplusMIPs_" + IntToString((int) iL);
      std::string auxnam1_tgr1 = "dEvsXos_" + IntToString((int) iL);
      std::string auxnam1_sh = "Showerprofile_" + IntToString((int) iL);
      std::string auxnam1_sf_prof = "h_sfvsshowerdepthprofile_" + IntToString((int) iL);
      std::string auxnam2 = "_" + IntToString(particleenergies[k]);
      std::string auxnam_sf2 = auxnam1_sf2 + auxnam2;
      std::string auxnam_tgr2 = auxnam1_tgr2 + auxnam2;
      std::string auxnam_tgr1 = auxnam1_tgr1 + auxnam2;
      std::string auxnam_sh = auxnam1_sh + auxnam2;
      std::string auxnam_sf_prof = auxnam1_sf_prof + auxnam2;

      h_MIPSvsdEplusMIPs[(int) iL].push_back(new TH2F( auxnam_sf2.c_str(),";E_{active} + E_{passive};E_{active};E (GeV);SimHits",520,0.,520.,1000,0.,1000.));
      h_Showerprofile[(int) iL].push_back(new TH1F( auxnam_sh.c_str(),";E (MIPs);Events",1000,0.,1000.));
      h_loss[(int) iL].push_back(new TH1F( auxnam_tgr1.c_str(),";dE (MIPs);Events",1000,0.,500000.));

      h_sfvsshowerdepthprofile[(int) iL].push_back(new TH2F( auxnam_sf_prof.c_str(),";"+shdp+";Sampling fraction of layer;",100,0.,maxbsh,100,0.006,0.020) );

      MIPSvsdEplusMIPs_calib[(int) iL].push_back( (TGraphErrors *) calib->Clone( auxnam_tgr2.c_str() )  );

    }
  }
  
  for (Int_t k=0; k<numberofenergies; k++){
    std::string auxnam_tgr1 = "Totalsfvsshowerdepth_" + IntToString(particleenergies[k]);
    std::string auxnam_tgr2 = "Totalsfvsshowerdepthlay_" + IntToString(particleenergies[k]);
    std::string auxnam_tgr3 = "dEvsXos_" + IntToString(particleenergies[k]);
    h_ShowerDepth_reco.push_back(new TH1F( ("ShowerDepth_reco_"+ IntToString(particleenergies[k])).c_str(),";Shower depth (1/X_{0});Events",300,0.,30.) );
    h_NewShowerDepth.push_back(new TH1F( ("h_NewShowerDepth_"+ IntToString(particleenergies[k])).c_str(),";#sum_{i=1}^{30}(E_{i}X_{0,i})/#sum_{i=1}^{30}E_{i};Events",300,0.,maxbsh) );
    h_TotalSamplingFraction.push_back(new TH1F( ("TotalSF_"+ IntToString(particleenergies[k])).c_str(),";Total sampling fraction;Events",100,0.006,0.020) );
    h_TotalSamplingFraction_sec1.push_back(new TH1F( ("TotalSF_sec1_"+ IntToString(particleenergies[k])).c_str(),";Total sampling fraction;Events",400,0.,0.04) );
    h_TotalSamplingFraction_sec2.push_back(new TH1F( ("TotalSF_sec2_"+ IntToString(particleenergies[k])).c_str(),";Total sampling fraction;Events",400,0.,0.04) );
    h_TotalSamplingFraction_sec3.push_back(new TH1F( ("TotalSF_sec3_"+ IntToString(particleenergies[k])).c_str(),";Total sampling fraction;Events",400,0.,0.04) );
    Totalsfvsshowerdepth_forprofile.push_back(new TH2F( ("TotalSFvsShowerDepth_"+ IntToString(particleenergies[k])).c_str(),";"+shdp+";Total sampling fraction;",100,0.,maxbsh,100,0.006,0.020) );
    Totalsf_sec1vsshowerdepth_forprofile.push_back(new TH2F( ("TotalSF_sec1vsShowerDepth_"+ IntToString(particleenergies[k])).c_str(),";"+shdp+";Total sampling fraction;",100,0.,maxbsh,400,0.,0.04) );
    Totalsf_sec2vsshowerdepth_forprofile.push_back(new TH2F( ("TotalSF_sec2vsShowerDepth_"+ IntToString(particleenergies[k])).c_str(),";"+shdp+";Total sampling fraction;",100,0.,maxbsh,400,0.,0.04) );
    Totalsf_sec3vsshowerdepth_forprofile.push_back(new TH2F( ("TotalSF_sec3vsShowerDepth_"+ IntToString(particleenergies[k])).c_str(),";"+shdp+";Total sampling fraction;",100,0.,maxbsh,400,0.,0.04) );
    //5 sections plus first layer = 6 sections
    h_TotalSamplingFraction_dedxsec1.push_back(new TH1F( ("TotalSF_dedxsec1_"+ IntToString(particleenergies[k])).c_str(),";Total sampling fraction;Events",400,0.,0.04) );
    h_TotalSamplingFraction_dedxsec2.push_back(new TH1F( ("TotalSF_dedxsec2_"+ IntToString(particleenergies[k])).c_str(),";Total sampling fraction;Events",400,0.,0.04) );
    h_TotalSamplingFraction_dedxsec3.push_back(new TH1F( ("TotalSF_dedxsec3_"+ IntToString(particleenergies[k])).c_str(),";Total sampling fraction;Events",400,0.,0.04) );
    h_TotalSamplingFraction_dedxsec4.push_back(new TH1F( ("TotalSF_dedxsec4_"+ IntToString(particleenergies[k])).c_str(),";Total sampling fraction;Events",400,0.,0.04) );
    h_TotalSamplingFraction_dedxsec5.push_back(new TH1F( ("TotalSF_dedxsec5_"+ IntToString(particleenergies[k])).c_str(),";Total sampling fraction;Events",400,0.,0.04) );
    h_TotalSamplingFraction_dedxsec6.push_back(new TH1F( ("TotalSF_dedxsec6_"+ IntToString(particleenergies[k])).c_str(),";Total sampling fraction;Events",400,0.,0.04) );
    Totalsf_dedxsec1vsshowerdepth_forprofile.push_back(new TH2F( ("TotalSF_dedxsec1vsShowerDepth_"+ IntToString(particleenergies[k])).c_str(),";"+shdp+";Total sampling fraction;",100,0.,maxbsh,400,0.,0.04) );
    Totalsf_dedxsec2vsshowerdepth_forprofile.push_back(new TH2F( ("TotalSF_dedxsec2vsShowerDepth_"+ IntToString(particleenergies[k])).c_str(),";"+shdp+";Total sampling fraction;",100,0.,maxbsh,400,0.,0.04) );
    Totalsf_dedxsec3vsshowerdepth_forprofile.push_back(new TH2F( ("TotalSF_dedxsec3vsShowerDepth_"+ IntToString(particleenergies[k])).c_str(),";"+shdp+";Total sampling fraction;",100,0.,maxbsh,400,0.,0.04) );
    Totalsf_dedxsec4vsshowerdepth_forprofile.push_back(new TH2F( ("TotalSF_dedxsec4vsShowerDepth_"+ IntToString(particleenergies[k])).c_str(),";"+shdp+";Total sampling fraction;",100,0.,maxbsh,400,0.,0.04) );
    Totalsf_dedxsec5vsshowerdepth_forprofile.push_back(new TH2F( ("TotalSF_dedxsec5vsShowerDepth_"+ IntToString(particleenergies[k])).c_str(),";"+shdp+";Total sampling fraction;",100,0.,maxbsh,400,0.,0.04) );
    Totalsf_dedxsec6vsshowerdepth_forprofile.push_back(new TH2F( ("TotalSF_dedxsec6vsShowerDepth_"+ IntToString(particleenergies[k])).c_str(),";"+shdp+";Total sampling fraction;",100,0.,maxbsh,400,0.,0.04) );

    sfvsshowerdepthprofileandlayer.push_back(new TH3F( ("sfvsshowerdepthprofileandlayer_"+ IntToString(particleenergies[k])).c_str(),";"+shdp+";Total sampling fraction;Layer",100,0.,maxbsh,100,0.006,0.020,31,0.,31.) );
    sfvsshowerdepthprofile_layerasweight.push_back(new TH2F( ("sfvsshowerdepthprofile_layerasweight_"+ IntToString(particleenergies[k])).c_str(),";"+shdp+";Sampling fraction of layer;",100,0.,maxbsh,100,0.006,0.020) );
    h_sfvsnormshowerdepth.push_back(new TH2F( ("h_sfvsnormshowerdepth_"+ IntToString(particleenergies[k])).c_str(),";"+shdp+";Sampling fraction of layer;",100,0.,maxbsh,400,0.,0.04) ); 
    Totalsfvsshowerdepth.push_back( (TGraphErrors *) calib->Clone( auxnam_tgr1.c_str() )  );
    Totalsfvsshowerdepthlay.push_back( (TGraphErrors *) calib->Clone( auxnam_tgr2.c_str() )  );
    Totalsfvsshowerdepth[k]->GetXaxis()->SetRangeUser(0.,maxbsh); //27.
    Totalsfvsshowerdepthlay[k]->GetXaxis()->SetRangeUser(0.,30.);
    dEvsXos.push_back( (TGraphErrors *) calib->Clone( auxnam_tgr3.c_str() )  );
  }


  //Now for the shower max
  //RECIPE for the profile of the totsf vs xomax/meanxomax
  //1. Run the Samplingfraction.cpp to make PFcal_12_final.root
  //2. Run the RecoE.cpp to get the meanshowermaxperenergy above. 
  //3. Rerun the Samplingfraction.cpp to get the profile histo
  double meanshowermaxperenergy[numberofenergies];
  meanshowermaxperenergy[0] = 4.95222;
  meanshowermaxperenergy[1] = 5.1993;
  meanshowermaxperenergy[2] = 5.65257;
  meanshowermaxperenergy[3] = 6.07514;
  meanshowermaxperenergy[4] = 6.27705;
  meanshowermaxperenergy[5] = 6.676;
  meanshowermaxperenergy[6] = 7.38304;
  meanshowermaxperenergy[7] = 7.93744;
  meanshowermaxperenergy[8] = 8.41862;
  meanshowermaxperenergy[9] = 8.67068;
  meanshowermaxperenergy[10] = 8.84051;
  meanshowermaxperenergy[11] = 9.09047;
  meanshowermaxperenergy[12] = 9.26026;
  meanshowermaxperenergy[13] = 9.36666;
  meanshowermaxperenergy[14] = 9.57008;
  meanshowermaxperenergy[15] = 9.80897;
  meanshowermaxperenergy[16] = 10.0913;
  meanshowermaxperenergy[17] = 10.3638;
  //Now for the new shower depth variable
  //RECIPE for the new shower depth variable mean value: 
  //Run only step1 and read the values from the terminal
  double meannewshowerdepthperenergy[numberofenergies];
  meannewshowerdepthperenergy[0] = 0.673458;
  meannewshowerdepthperenergy[1] = 0.68317;
  meannewshowerdepthperenergy[2] = 0.698144;
  meannewshowerdepthperenergy[3] = 0.712634;
  meannewshowerdepthperenergy[4] = 0.720433;
  meannewshowerdepthperenergy[5] = 0.733995;
  meannewshowerdepthperenergy[6] = 0.757765;
  meannewshowerdepthperenergy[7] = 0.776155;
  meannewshowerdepthperenergy[8] = 0.792933;
  meannewshowerdepthperenergy[9] = 0.801784;
  meannewshowerdepthperenergy[10] = 0.80831;
  meannewshowerdepthperenergy[11] = 0.816638;
  meannewshowerdepthperenergy[12] = 0.822941;
  meannewshowerdepthperenergy[13] = 0.826855;
  meannewshowerdepthperenergy[14] = 0.834427;
  meannewshowerdepthperenergy[15] = 0.841823;
  meannewshowerdepthperenergy[16] = 0.852045;
  meannewshowerdepthperenergy[17] = 0.860961;

  // Double_t t_showerdepthbinning[(*ssvec).size()-1];//ignore the before calo material
  double currentxotrans = 0.;

  // for(unsigned iL(0); iL<(*ssvec).size(); iL++){
  //   if(iL!=0){
  //     //Count from layer 1 as start of calo and not include before calorimeter.
  //     currentxotrans = currentxotrans + (*ssvec)[iL].volX0trans();
  //     t_showerdepthbinning[iL-1]=currentxotrans;
  //   }

  // }

  //======================================================================================
  //For the dedx definition of the sampling fraction
  TString titleoffile = "data/V5_dedx.txt";
  std::ifstream infile(titleoffile);
  std::string dummyfordedx; 
  //dedx for all layers minus one for the upstream material or zero if no material upstream
  double dedx[(*ssvec).size()-startlayer];

  std::string buffer;
  std::cout << "dedx for layers 1-30 from file " << titleoffile << std::endl;
  std::cout << "dedx   total_dedx    " << titleoffile << std::endl;
  double totdedx = 0.;
  double totdedx_standalone = 0.;

  for (unsigned int k=0; k<(*ssvec).size()-startlayer; k++){
    getline(infile, buffer);
    std::istringstream is(buffer);

    is >> dedx[k];
    totdedx = totdedx + dedx[k];
    totdedx_standalone = totdedx_standalone + (*ssvec)[k+startlayer].voldEdx();
    //Test
    // std::cout << "-------------------------------------" <<std::endl;
    // std::cout << dedx[k] << "   " << totdedx <<  std::endl;
    std::cout << (*ssvec)[k+startlayer].voldEdx() << "   " << totdedx_standalone <<std::endl;
  }
  infile.close();

  //Definition of sampling fraction for dedx
  double SF_dedx = 0;
  //dedx for the silicon per layer
  // double si_dedx = 0.11628; //it is 0.3 mm * 0.3876 MeV/mm
  //total dedx for the silicon in whole calorimeter
  double totalsi_dedx = 3.4884; //30 layers x si_dedx
  //dedx for the absorber
  double abs_dedx = 0.;
  //Do not include the volume for the leakage
  for(unsigned iL(startlayer); iL<(*ssvec).size()-startlayer; iL++){
    abs_dedx = abs_dedx + (*ssvec)[iL].voldEdx();
  }

  SF_dedx = (totalsi_dedx/(abs_dedx+totalsi_dedx));
  
  //So for the e/mip : totalsi_dedx/(totalsi_dedx+) 
  std::cout << "Sampling fraction with dedx for e/mip: " << SF_dedx << std::endl;

  std::vector<double> xotrans;
  std::vector<double> sfpl;
  //======================================================================================
  //The files that we will store the results of this analysis
  TString res[numberofconfigs];
  for (int k=0; k<numberofconfigs; k++){
    TString fp = "PFcal_" + IntToString(configs[k]);	
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
	  // if (ievt==100){break;}
	  // std::cout << "Entry " << ievt << std::endl;
	  // double lossbefcalo = 0.;
	  double totalscalefactor = 0.;
	  double totalscalefactor_sec1 = 0.;
	  double totalscalefactor_sec2 = 0.;
	  double totalscalefactor_sec3 = 0.;
	  double lossinsensors = 0.;
	  double totalloss = 0.;
	  double lossinsensors_sec1 = 0.;
	  double totalloss_sec1 = 0.;
	  double lossinsensors_sec2 = 0.;
	  double totalloss_sec2 = 0.;
	  double lossinsensors_sec3 = 0.;
	  double totalloss_sec3 = 0.;
	  //5 sections plus first layer = 6 sections
	  double totalscalefactor_dedxsec1 = 0.;
	  double totalscalefactor_dedxsec2 = 0.;
	  double totalscalefactor_dedxsec3 = 0.;
	  double totalscalefactor_dedxsec4 = 0.;
	  double totalscalefactor_dedxsec5 = 0.;
	  double totalscalefactor_dedxsec6 = 0.;
	  double lossinsensors_dedxsec1 = 0.;
	  double totalloss_dedxsec1 = 0.;
	  double lossinsensors_dedxsec2 = 0.;
	  double totalloss_dedxsec2 = 0.;
	  double lossinsensors_dedxsec3 = 0.;
	  double totalloss_dedxsec3 = 0.;
	  double lossinsensors_dedxsec4 = 0.;
	  double totalloss_dedxsec4 = 0.;
	  double lossinsensors_dedxsec5 = 0.;
	  double totalloss_dedxsec5 = 0.;
	  double lossinsensors_dedxsec6 = 0.;
	  double totalloss_dedxsec6 = 0.;
	  
	  //For the xotransfered add the before calo xo's
	  xotrans.clear();
	  currentxotrans = 0.;
	  int layermax = -1.;
	  double currmaxenergy = 0.;
	  //Sampling fraction per layer
	  sfpl.clear();
	  //For the sum(Eixoi) in the shower depth
	  double Xoweightedwithenergy = 0.;

	  //Loop on sampling sections
	  for(unsigned iL(0); iL<(*ssvec).size(); iL++){
	    // HGCSSSamplingSection lSamSec = (*ssvec)[iL];
	    // double energyabsorber = lSamSec.absorberE();

	    // std::cout << "Sampling Section " << iL << " energyabsorber " << energyabsorber  << std::endl; 
	  
	    //Esilicon vs (Epassive+Esilicon) for the sampling fraction
	    h_MIPSvsdEplusMIPs[(int) iL][l]->Fill( ( (*ssvec)[iL].measuredE() + (*ssvec)[iL].absorberE() )/(1000.) , (*ssvec)[iL].measuredE()/1000. ); //x is GeV, y is GeV. 


	    Int_t np2=MIPSvsdEplusMIPs_calib[(int) iL][l]->GetN();
	    // std::cout << "np " << np << std::endl;
	    MIPSvsdEplusMIPs_calib[(int) iL][l]->SetPoint(np2, ( (*ssvec)[iL].measuredE() + (*ssvec)[iL].absorberE() )/(1000.)  , (*ssvec)[iL].measuredE()/1000.  );
	    // MIPSvsdEplusMIPs_calib[(int) iL][l]->SetPointError();
	    


	    h_Showerprofile[(int) iL][l]->Fill( ( (*ssvec)[iL].measuredE()/(0.081) ) );//in MIPs
	    h_loss[(int) iL][l]->Fill( (*ssvec)[iL].absorberE()/(0.081) );//in MIPs
	    
	    if( iL!=(int)(startlayer-1) ){
	      //Count from layer 1 as start of calo and not include before calorimeter.
	      //If no material in front starts from zero
	      currentxotrans = currentxotrans + (*ssvec)[iL].volX0trans();
	      xotrans.push_back(currentxotrans);
	      if ( currmaxenergy < (*ssvec)[iL].measuredE() ) { 
		layermax = (int) iL;
		currmaxenergy = (*ssvec)[iL].measuredE();	      
	      }
	      //This is for the new shower depth variable
	      Xoweightedwithenergy = Xoweightedwithenergy + ( (*ssvec)[iL].measuredE() * ( (*ssvec)[iL].volX0trans() ) );
	    }
	    


	    double w = 0.;
	    if ( (iL!=(unsigned)(startlayer-1)) ){
	    // if ( (iL!=(int)(startlayer-1)) && (iL!=31)){
	      totalloss = totalloss + (*ssvec)[iL].measuredE() + (*ssvec)[iL].absorberE(); //in MeV
	      lossinsensors = lossinsensors + (*ssvec)[iL].measuredE(); //in MeV
	      w = (*ssvec)[iL].measuredE()/((*ssvec)[iL].measuredE() + (*ssvec)[iL].absorberE());
	      sfpl.push_back(w);
	    }

	    //The section we are in
	    //!!!!! THIS IS FOR NO UPSTREAM material case
	    //for the total sf per section
	    bool totsf_sec1 = false;
	    bool totsf_sec2 = false;
	    bool totsf_sec3 = false;
	    if ( (iL==0) || (iL==3) || (iL==5) || (iL==7) || (iL==9) ) {totsf_sec1 = true;}
	    if ( (iL==1) || (iL==2) || (iL==4) || (iL==6) || (iL==8) || (iL==10) || (iL==11) || 
		 (iL==12) || (iL==13) || (iL==14) || (iL==15) || (iL==16) || (iL==17) || (iL==18) || 
		 (iL==19) || (iL==20) || (iL==21) || (iL==23) || (iL==25) || (iL==27) || (iL==29) ) {totsf_sec2 = true;}
	    if ( (iL==22) || (iL==24) || (iL==26) || (iL==28) ) {totsf_sec3 = true;}

	    if (totsf_sec1){
	      totalloss_sec1 = totalloss_sec1 + (*ssvec)[iL].measuredE() + (*ssvec)[iL].absorberE(); //in MeV
	      lossinsensors_sec1 = lossinsensors_sec1 + (*ssvec)[iL].measuredE(); //in MeV
	    }
	    if (totsf_sec2){
	      totalloss_sec2 = totalloss_sec2 + (*ssvec)[iL].measuredE() + (*ssvec)[iL].absorberE(); //in MeV
	      lossinsensors_sec2 = lossinsensors_sec2 + (*ssvec)[iL].measuredE(); //in MeV

	    }
	    if (totsf_sec3){
	      totalloss_sec3 = totalloss_sec3 + (*ssvec)[iL].measuredE() + (*ssvec)[iL].absorberE(); //in MeV
	      lossinsensors_sec3 = lossinsensors_sec3 + (*ssvec)[iL].measuredE(); //in MeV
	    }

	    //5 sections plus first layer = 6 sections
	    //!!!!! THIS IS FOR NO UPSTREAM material case
	    //for the total sf per section
	    bool totsf_dedxsec1 = false;
	    bool totsf_dedxsec2 = false;
	    bool totsf_dedxsec3 = false;
	    bool totsf_dedxsec4 = false;
	    bool totsf_dedxsec5 = false;
	    bool totsf_dedxsec6 = false;
	    if (iL==0) {totsf_dedxsec1 = true;}
	    if ( (iL==3) || (iL==5) || (iL==7) || (iL==9) ) {totsf_dedxsec2 = true;}
	    if ( (iL==11) || (iL==13) || (iL==15) || (iL==17) || (iL==19) ) {totsf_dedxsec3 = true;}
	    if ( (iL==1) || (iL==2) || (iL==4) || (iL==6) || (iL==8) || (iL==10)) {totsf_dedxsec4 = true;}
	    if ( (iL==12) || (iL==14) || (iL==16) || (iL==18) || (iL==20) || (iL==21) || (iL==23) || (iL==25) || (iL==27) || (iL==29) ) {totsf_dedxsec5 = true;}
	    if ( (iL==22) || (iL==24) || (iL==26) || (iL==28) ) {totsf_dedxsec6 = true;}

	    if (totsf_dedxsec1){
	      totalloss_dedxsec1 = totalloss_dedxsec1 + (*ssvec)[iL].measuredE() + (*ssvec)[iL].absorberE(); //in MeV
	      lossinsensors_dedxsec1 = lossinsensors_dedxsec1 + (*ssvec)[iL].measuredE(); //in MeV
	    }
	    if (totsf_dedxsec2){
	      totalloss_dedxsec2 = totalloss_dedxsec2 + (*ssvec)[iL].measuredE() + (*ssvec)[iL].absorberE(); //in MeV
	      lossinsensors_dedxsec2 = lossinsensors_dedxsec2 + (*ssvec)[iL].measuredE(); //in MeV

	    }
	    if (totsf_dedxsec3){
	      totalloss_dedxsec3 = totalloss_dedxsec3 + (*ssvec)[iL].measuredE() + (*ssvec)[iL].absorberE(); //in MeV
	      lossinsensors_dedxsec3 = lossinsensors_dedxsec3 + (*ssvec)[iL].measuredE(); //in MeV
	    }
	    if (totsf_dedxsec4){
	      totalloss_dedxsec4 = totalloss_dedxsec4 + (*ssvec)[iL].measuredE() + (*ssvec)[iL].absorberE(); //in MeV
	      lossinsensors_dedxsec4 = lossinsensors_dedxsec4 + (*ssvec)[iL].measuredE(); //in MeV
	    }
	    if (totsf_dedxsec5){
	      totalloss_dedxsec5 = totalloss_dedxsec5 + (*ssvec)[iL].measuredE() + (*ssvec)[iL].absorberE(); //in MeV
	      lossinsensors_dedxsec5 = lossinsensors_dedxsec5 + (*ssvec)[iL].measuredE(); //in MeV

	    }
	    if (totsf_dedxsec6){
	      totalloss_dedxsec6 = totalloss_dedxsec6 + (*ssvec)[iL].measuredE() + (*ssvec)[iL].absorberE(); //in MeV
	      lossinsensors_dedxsec6 = lossinsensors_dedxsec6 + (*ssvec)[iL].measuredE(); //in MeV
	    }

	    
	  } //loop on sampling sections

	  //For the total sampling fraction per event
	  //Sum_1^30(Esilicon_i) / Sum_1^30(Epassive+Esilicon)
	  totalscalefactor = lossinsensors/totalloss;
	  h_TotalSamplingFraction[l]->Fill( totalscalefactor );
	  //per section
	  totalscalefactor_sec1 = lossinsensors_sec1/totalloss_sec1;
	  totalscalefactor_sec2 = lossinsensors_sec2/totalloss_sec2;
	  totalscalefactor_sec3 = lossinsensors_sec3/totalloss_sec3;
	  h_TotalSamplingFraction_sec1[l]->Fill( totalscalefactor_sec1 );
	  h_TotalSamplingFraction_sec2[l]->Fill( totalscalefactor_sec2 );
	  h_TotalSamplingFraction_sec3[l]->Fill( totalscalefactor_sec3 );
	  //for the 5 plus 1 sections
	  totalscalefactor_dedxsec1 = lossinsensors_dedxsec1/totalloss_dedxsec1;
	  totalscalefactor_dedxsec2 = lossinsensors_dedxsec2/totalloss_dedxsec2;
	  totalscalefactor_dedxsec3 = lossinsensors_dedxsec3/totalloss_dedxsec3;
	  totalscalefactor_dedxsec4 = lossinsensors_dedxsec4/totalloss_dedxsec4;
	  totalscalefactor_dedxsec5 = lossinsensors_dedxsec5/totalloss_dedxsec5;
	  totalscalefactor_dedxsec6 = lossinsensors_dedxsec6/totalloss_dedxsec6;
	  h_TotalSamplingFraction_dedxsec1[l]->Fill( totalscalefactor_dedxsec1 );
	  h_TotalSamplingFraction_dedxsec2[l]->Fill( totalscalefactor_dedxsec2 );
	  h_TotalSamplingFraction_dedxsec3[l]->Fill( totalscalefactor_dedxsec3 );
	  h_TotalSamplingFraction_dedxsec4[l]->Fill( totalscalefactor_dedxsec4 );
	  h_TotalSamplingFraction_dedxsec5[l]->Fill( totalscalefactor_dedxsec5 );
	  h_TotalSamplingFraction_dedxsec6[l]->Fill( totalscalefactor_dedxsec6 );
	  

	  //layer max 9 is the 8 element (0,...8) of vector when material in front and 9 when no material upstream
	  h_ShowerDepth_reco[l]->Fill( xotrans[layermax-startlayer]  );
	  // std::cout << "layermax " << layermax << " xo " << xotrans[layermax-startlayer] << std::endl;

	  //Now we know X0max/<X0max> for the current event
	  double samfrac_shomax = xotrans[layermax-startlayer]/meanshowermaxperenergy[l];
	  //This is for the new shower depth variable of the current event defined above see bool
	  double newshowerdepth = Xoweightedwithenergy / lossinsensors;
	  //And of course for the new shower depth variable we should divide with the mean. 
	  double newshowerdepthweight = newshowerdepth / meannewshowerdepthperenergy[l]; 
	  if (usenewshowerdepth){
	    samfrac_shomax = newshowerdepthweight;
	  }

	  h_NewShowerDepth[l]->Fill( newshowerdepth );
	  Totalsfvsshowerdepth_forprofile[l]->Fill(samfrac_shomax,totalscalefactor);
	  Totalsf_sec1vsshowerdepth_forprofile[l]->Fill(samfrac_shomax,totalscalefactor_sec1);
	  Totalsf_sec2vsshowerdepth_forprofile[l]->Fill(samfrac_shomax,totalscalefactor_sec2);
	  Totalsf_sec3vsshowerdepth_forprofile[l]->Fill(samfrac_shomax,totalscalefactor_sec3);
	  //for the 5 plus 1 sections
	  Totalsf_dedxsec1vsshowerdepth_forprofile[l]->Fill(samfrac_shomax,totalscalefactor_dedxsec1);
	  Totalsf_dedxsec2vsshowerdepth_forprofile[l]->Fill(samfrac_shomax,totalscalefactor_dedxsec2);
	  Totalsf_dedxsec3vsshowerdepth_forprofile[l]->Fill(samfrac_shomax,totalscalefactor_dedxsec3);
	  Totalsf_dedxsec4vsshowerdepth_forprofile[l]->Fill(samfrac_shomax,totalscalefactor_dedxsec4);
	  Totalsf_dedxsec5vsshowerdepth_forprofile[l]->Fill(samfrac_shomax,totalscalefactor_dedxsec5);
	  Totalsf_dedxsec6vsshowerdepth_forprofile[l]->Fill(samfrac_shomax,totalscalefactor_dedxsec6);

	  //For the total sampling fraction vs the depth of the shower per event
	  Int_t nnp=Totalsfvsshowerdepth[l]->GetN();
	  // Totalsfvsshowerdepth[l]->SetPoint( nnp, xotrans[layermax-startlayer] , totalscalefactor ); // x is in radiation length, y is the SF
	  Totalsfvsshowerdepth[l]->SetPoint( nnp, samfrac_shomax , totalscalefactor ); // x is in radiation length, y is the SF
	  Totalsfvsshowerdepth[l]->SetPointError( nnp, 0. , 0. ); 

	  if(dothe3dplot){
	    for(unsigned iL(startlayer); iL<(*ssvec).size(); iL++){
	      // sfvsshowerdepthprofileandlayer[l]->Fill(xotrans[layermax-startlayer], sfpl[iL-startlayer] , iL );
	      // sfvsshowerdepthprofile_layerasweight[l]->Fill(xotrans[layermax-startlayer], sfpl[iL-startlayer] , iL );
	      // h_sfvsshowerdepthprofile[iL][l]->Fill(xotrans[layermax-startlayer], sfpl[iL-startlayer] );
	      sfvsshowerdepthprofileandlayer[l]->Fill(samfrac_shomax, sfpl[iL-startlayer] , iL );
	      sfvsshowerdepthprofile_layerasweight[l]->Fill(samfrac_shomax, sfpl[iL-startlayer] , iL );
	      h_sfvsshowerdepthprofile[iL][l]->Fill(samfrac_shomax, sfpl[iL-startlayer] );
	      //SF vs normalized shower depth: sf_i vs t=x0/x0,max
	      // h_sfvsnormshowerdepth[l]->Fill(xotrans[iL-startlayer]/xotrans[layermax-startlayer], sfpl[iL-startlayer] );
	      h_sfvsnormshowerdepth[l]->Fill(samfrac_shomax, sfpl[iL-startlayer] );
	    }
	  }

	  //For the total sampling fraction vs the depth of the shower per event with layer unit
	  Int_t nnp1=Totalsfvsshowerdepthlay[l]->GetN();
	  Totalsfvsshowerdepthlay[l]->SetPoint( nnp1, layermax , totalscalefactor ); // x is layer, y is the SF
	  Totalsfvsshowerdepthlay[l]->SetPointError( nnp1, 0. , 0. ); 

	} //  Loop on entries

	if(dothedEvsXosplot){
	  //Loop again on sampling sections for the average plots
	  for(unsigned iL(startlayer); iL<(*ssvec).size(); iL++){
	    //For the average energy lost as a function of the transversed thickness
	    Int_t nnp=dEvsXos[l]->GetN();
	    dEvsXos[l]->SetPoint( nnp, xotrans[(int) iL] , h_loss[(int) iL][l]->GetMean() ); // y is in MIPs, x is in radiation length
	    dEvsXos[l]->SetPointError( nnp, 0 , h_loss[(int) iL][l]->GetRMS() ); // y is in MIPs, x is in radiation length
	  }
	}

	//======================================================================================
	//Write histos 
	if(dothedEvsXosplot){
	  dEvsXos[l]->Write();
	}
	h_TotalSamplingFraction[l]->Write();
	h_TotalSamplingFraction_sec1[l]->Write();
	h_TotalSamplingFraction_sec2[l]->Write();
	h_TotalSamplingFraction_sec3[l]->Write();
	h_ShowerDepth_reco[l]->Write();
	h_NewShowerDepth[l]->Write();
	Totalsfvsshowerdepth_forprofile[l]->Write();
	Totalsf_sec1vsshowerdepth_forprofile[l]->Write();
	Totalsf_sec2vsshowerdepth_forprofile[l]->Write();
	Totalsf_sec3vsshowerdepth_forprofile[l]->Write();
	//5+1 sections
	h_TotalSamplingFraction_dedxsec1[l]->Write();
	h_TotalSamplingFraction_dedxsec2[l]->Write();
	h_TotalSamplingFraction_dedxsec3[l]->Write();
	h_TotalSamplingFraction_dedxsec4[l]->Write();
	h_TotalSamplingFraction_dedxsec5[l]->Write();
	h_TotalSamplingFraction_dedxsec6[l]->Write();
	Totalsf_dedxsec1vsshowerdepth_forprofile[l]->Write();
	Totalsf_dedxsec2vsshowerdepth_forprofile[l]->Write();
	Totalsf_dedxsec3vsshowerdepth_forprofile[l]->Write();
	Totalsf_dedxsec4vsshowerdepth_forprofile[l]->Write();
	Totalsf_dedxsec5vsshowerdepth_forprofile[l]->Write();
	Totalsf_dedxsec6vsshowerdepth_forprofile[l]->Write();
	sfvsshowerdepthprofileandlayer[l]->Write();
	sfvsshowerdepthprofile_layerasweight[l]->Write();
	h_sfvsnormshowerdepth[l]->Write();
	lChain[0][0]->GetEntry(0);
	for(unsigned iL(0); iL<(*ssvec).size(); iL++){
	  h_MIPSvsdEplusMIPs[(int) iL][l]->Write();
	  MIPSvsdEplusMIPs_calib[(int) iL][l]->Write();
	  h_sfvsshowerdepthprofile[iL][l]->Write();
	  h_loss[iL][l]->Write();
	}

	//We will print the total sampling fraction here 
	std::cout << "Total Sampling fraction for " << particleenergies[l] << " GeV : " << h_TotalSamplingFraction[l]->GetMean() << std::endl;  
	std::cout << "Total Sampling fraction for Sec 1 " << particleenergies[l] << " GeV : " << h_TotalSamplingFraction_sec1[l]->GetMean() << std::endl;  
	std::cout << "Total Sampling fraction for Sec 2 " << particleenergies[l] << " GeV : " << h_TotalSamplingFraction_sec2[l]->GetMean() << std::endl;  
	std::cout << "Total Sampling fraction for Sec 3 " << particleenergies[l] << " GeV : " << h_TotalSamplingFraction_sec3[l]->GetMean() << std::endl;  
	std::cout << "======================================" << std::endl;
	std::cout << "Total Sampling fraction for Dedxsec 1 " << particleenergies[l] << " GeV : " << h_TotalSamplingFraction_dedxsec1[l]->GetMean() << std::endl;  
	std::cout << "Total Sampling fraction for Dedxsec 2 " << particleenergies[l] << " GeV : " << h_TotalSamplingFraction_dedxsec2[l]->GetMean() << std::endl;  
	std::cout << "Total Sampling fraction for Dedxsec 3 " << particleenergies[l] << " GeV : " << h_TotalSamplingFraction_dedxsec3[l]->GetMean() << std::endl;  
	std::cout << "Total Sampling fraction for Dedxsec 4 " << particleenergies[l] << " GeV : " << h_TotalSamplingFraction_dedxsec4[l]->GetMean() << std::endl;  
	std::cout << "Total Sampling fraction for Dedxsec 5 " << particleenergies[l] << " GeV : " << h_TotalSamplingFraction_dedxsec5[l]->GetMean() << std::endl;  
	std::cout << "Total Sampling fraction for Dedxsec 6 " << particleenergies[l] << " GeV : " << h_TotalSamplingFraction_dedxsec6[l]->GetMean() << std::endl;  

	//Here we want to know the average value for the new shower depth variable
	std::cout << "======================================" << std::endl;
	std::cout << "New shower depth:  <Sum(Ei*X0_i)/Sum(E_i)> " << particleenergies[l] << " GeV : " << h_NewShowerDepth[l]->GetMean() << std::endl;  
	


 
	//======================================================================================
	//We should here clear the histograms because we want them empty for the next file. 
	//Reset histos
	for(unsigned iL(0); iL<(*ssvec).size(); iL++){
	  h_MIPSvsdEplusMIPs[(int) iL][l]->Reset();
	}

	h_TotalSamplingFraction[l]->Reset();
	h_TotalSamplingFraction_sec1[l]->Reset();
	h_TotalSamplingFraction_sec2[l]->Reset();
	h_TotalSamplingFraction_sec3[l]->Reset();
	h_TotalSamplingFraction_dedxsec1[l]->Reset();
	h_TotalSamplingFraction_dedxsec2[l]->Reset();
	h_TotalSamplingFraction_dedxsec3[l]->Reset();
	h_TotalSamplingFraction_dedxsec4[l]->Reset();
	h_TotalSamplingFraction_dedxsec5[l]->Reset();
	h_TotalSamplingFraction_dedxsec6[l]->Reset();
	//h_sfvsnormshowerdepth[l]->Reset();
	h_ShowerDepth_reco[l]->Reset();
	h_NewShowerDepth[l]->Reset();
       
      }//Loop on energies
   
      results[k]->Close();

    } //Loop on files

  }// if to avoid running again on input files

  //======================================================================================
  //Make plots
  TCanvas *c1[(*ssvec).size()];
  TCanvas *c2[(*ssvec).size()];
  TCanvas *c4[(*ssvec).size()][numberofenergies];
  TCanvas *c9[(*ssvec).size()];
  TCanvas *c10[(*ssvec).size()][numberofenergies];
  TCanvas *c13[(*ssvec).size()];
  for(unsigned iL(0); iL<(*ssvec).size(); iL++){
    c1[(int) iL] = new TCanvas(("c1_"+ IntToString((int) iL)).c_str(), "  ");
    c9[(int) iL] = new TCanvas(("c9_"+ IntToString((int) iL)).c_str(), "  ");
    c13[(int) iL] = new TCanvas(("c13_"+ IntToString((int) iL)).c_str(), "  ");
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
  TCanvas *c14_1 = new TCanvas("c14_1", "  ");
  TCanvas *c14_2 = new TCanvas("c14_2", "  ");
  TCanvas *c14_3 = new TCanvas("c14_3", "  ");
  TCanvas *c14_4 = new TCanvas("c14_4", "  ");
  TCanvas *c14_5 = new TCanvas("c14_5", "  ");
  TCanvas *c15 = new TCanvas("c15", "  ");
  TCanvas *c16 = new TCanvas("c16", "  ");
  TCanvas *c17_1 = new TCanvas("c17_1", "  ");
  TCanvas *c17_2 = new TCanvas("c17_2", "  ");
  TCanvas *c18 = new TCanvas("c18", "  ");
  TCanvas *c19 = new TCanvas("c19", "  ");
  TCanvas *c20 = new TCanvas("c20", "  ");
  TCanvas *c21 = new TCanvas("c21", "  ");
  TCanvas *c22 = new TCanvas("c22", "  ");
  TCanvas *c23 = new TCanvas("c23", "  ");
  TCanvas *c24 = new TCanvas("c24", "  ");
  TCanvas *c25 = new TCanvas("c25", "  ");
  TCanvas *c26 = new TCanvas("c26", "  ");

  c17_1->Divide(3,3);
  c17_2->Divide(3,3);

  TPad *pad = new TPad("pad","",0,0,1,1);
  //For the legend
  // TLegend* leg[numberoffiles];
  TLegend* leg1 = new TLegend(0.70, 0.2, 0.9, 0.9);
  leg1->SetHeader("Energy");
  leg1->SetFillColor(17);
  TLegend* leg2 = new TLegend(0.5, 0.7, 0.8, 0.9);
  leg2->SetHeader("Energy");
  leg2->SetFillColor(17);
  TLegend* leg[(*ssvec).size()];
  for(unsigned iL(0); iL<(*ssvec).size(); iL++){
    leg[(int) iL] = new TLegend(0.5, 0.7, 0.8, 0.9);
    leg[(int) iL]->SetHeader("Energy");
    leg[(int) iL]->SetFillColor(17);
  }
  TLegend* leg_sf2[(*ssvec).size()];
  for(unsigned iL(0); iL<(*ssvec).size(); iL++){
    leg_sf2[(int) iL] = new TLegend(0.5, 0.7, 0.8, 0.9);
    leg_sf2[(int) iL]->SetHeader("Energy");
    leg_sf2[(int) iL]->SetFillColor(17);
  }
  TLegend* leg_sf[(*ssvec).size()];
  for(unsigned iL(0); iL<(*ssvec).size(); iL++){
    leg_sf[(int) iL] = new TLegend(0.5, 0.7, 0.8, 0.9);
    leg_sf[(int) iL]->SetHeader("Energy");
    leg_sf[(int) iL]->SetFillColor(17);
  }
  TLegend* leg_devsmipswithfit[(*ssvec).size()];
  for(unsigned iL(0); iL<(*ssvec).size(); iL++){
    leg_devsmipswithfit[(int) iL] = new TLegend(0.5, 0.7, 0.8, 0.9);
    leg_devsmipswithfit[(int) iL]->SetHeader("Energy");
    leg_devsmipswithfit[(int) iL]->SetFillColor(17);
  }
  TLegend* leg_sfplvst_prof[(*ssvec).size()];
  for(unsigned iL(0); iL<(*ssvec).size(); iL++){
    leg_sfplvst_prof[(int) iL] = new TLegend(0.83, 0.11, 1.0, 1.0);
    leg_sfplvst_prof[(int) iL]->SetHeader("Energy");
    leg_sfplvst_prof[(int) iL]->SetFillColor(17);
  }

  TLegend* leg_devsxos = new TLegend(0.5, 0.7, 0.8, 0.9);
  leg_devsxos->SetHeader("Energy");
  leg_devsxos->SetFillColor(17);

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
  TLegend* leg_tsfvst = new TLegend(0.5, 0.7, 0.8, 0.9);
  leg_tsfvst->SetHeader("Energy");
  leg_tsfvst->SetFillColor(17);
  TLegend* leg_tsfvst_prof = new TLegend(0.83, 0.11, 1.0, 1.0);
  leg_tsfvst_prof->SetHeader("Energy");
  leg_tsfvst_prof->SetFillColor(17);
  TLegend* leg_tsfsec1vst_prof = new TLegend(0.83, 0.11, 1.0, 1.0);
  leg_tsfsec1vst_prof->SetHeader("Energy");
  leg_tsfsec1vst_prof->SetFillColor(17);
  TLegend* leg_tsfsec2vst_prof = new TLegend(0.83, 0.11, 1.0, 1.0);
  leg_tsfsec2vst_prof->SetHeader("Energy");
  leg_tsfsec2vst_prof->SetFillColor(17);
  TLegend* leg_tsfsec3vst_prof = new TLegend(0.83, 0.11, 1.0, 1.0);
  leg_tsfsec3vst_prof->SetHeader("Energy");
  leg_tsfsec3vst_prof->SetFillColor(17);
  //5+1sections
  TLegend* leg_tsfdedxsec1vst_prof = new TLegend(0.83, 0.11, 1.0, 1.0);
  leg_tsfdedxsec1vst_prof->SetHeader("Energy");
  leg_tsfdedxsec1vst_prof->SetFillColor(17);
  TLegend* leg_tsfdedxsec2vst_prof = new TLegend(0.83, 0.11, 1.0, 1.0);
  leg_tsfdedxsec2vst_prof->SetHeader("Energy");
  leg_tsfdedxsec2vst_prof->SetFillColor(17);
  TLegend* leg_tsfdedxsec3vst_prof = new TLegend(0.83, 0.11, 1.0, 1.0);
  leg_tsfdedxsec3vst_prof->SetHeader("Energy");
  leg_tsfdedxsec3vst_prof->SetFillColor(17);
  TLegend* leg_tsfdedxsec4vst_prof = new TLegend(0.83, 0.11, 1.0, 1.0);
  leg_tsfdedxsec4vst_prof->SetHeader("Energy");
  leg_tsfdedxsec4vst_prof->SetFillColor(17);
  TLegend* leg_tsfdedxsec5vst_prof = new TLegend(0.83, 0.11, 1.0, 1.0);
  leg_tsfdedxsec5vst_prof->SetHeader("Energy");
  leg_tsfdedxsec5vst_prof->SetFillColor(17);
  TLegend* leg_tsfdedxsec6vst_prof = new TLegend(0.83, 0.11, 1.0, 1.0);
  leg_tsfdedxsec6vst_prof->SetHeader("Energy");
  leg_tsfdedxsec6vst_prof->SetFillColor(17);

  TLegend* leg_tsfvstlay = new TLegend(0.5, 0.7, 0.8, 0.9);
  leg_tsfvstlay->SetHeader("Energy");
  leg_tsfvstlay->SetFillColor(17);
  TLegend* leg_w = new TLegend(0.5, 0.7, 0.8, 0.9);
  leg_w->SetHeader("Layers");
  leg_w->SetFillColor(17);
  TLegend* leg_allw = new TLegend(0.5, 0.7, 0.8, 0.9);
  leg_allw->SetHeader("Layers");
  leg_allw->SetFillColor(17);
  TLegend* leg_tsfvsnormdepth_prof = new TLegend(0.83, 0.11, 1.0, 1.0);
  leg_tsfvsnormdepth_prof->SetHeader("Energy");
  leg_tsfvsnormdepth_prof->SetFillColor(17);

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

  TString procNa_weights[4];
  procNa_weights[0] = "devsMIPs";
  procNa_weights[1] = "X_{0}";
  procNa_weights[2] = "dEdx";
  procNa_weights[3] = "#lamdba";


  TString procNa_layers[(*ssvec).size()];
  upstream ? procNa_layers[0] = "Before calo" : procNa_layers[0] = "Layer 1";
  upstream ? procNa_layers[1] = "Layer 1" : procNa_layers[1] = "Layer 2";
  upstream ? procNa_layers[2] = "Layer 2" : procNa_layers[2] = "Layer 3";
  upstream ? procNa_layers[3] = "Layer 3" : procNa_layers[3] = "Layer 4";
  upstream ? procNa_layers[4] = "Layer 4" : procNa_layers[4] = "Layer 5";
  upstream ? procNa_layers[5] = "Layer 5" : procNa_layers[5] = "Layer 6";
  upstream ? procNa_layers[6] = "Layer 6" : procNa_layers[6] = "Layer 7";
  upstream ? procNa_layers[7] = "Layer 7" : procNa_layers[7] = "Layer 8";
  upstream ? procNa_layers[8] = "Layer 8" : procNa_layers[8] = "Layer 9";
  upstream ? procNa_layers[9] = "Layer 9" : procNa_layers[9] = "Layer 10";
  upstream ? procNa_layers[10] = "Layer 10" : procNa_layers[10] = "Layer 11";
  upstream ? procNa_layers[11] = "Layer 11" : procNa_layers[11] = "Layer 12";
  upstream ? procNa_layers[12] = "Layer 12" : procNa_layers[12] = "Layer 13";
  upstream ? procNa_layers[13] = "Layer 13" : procNa_layers[13] = "Layer 14";
  upstream ? procNa_layers[14] = "Layer 14" : procNa_layers[14] = "Layer 15";
  upstream ? procNa_layers[15] = "Layer 15" : procNa_layers[15] = "Layer 16";
  upstream ? procNa_layers[16] = "Layer 16" : procNa_layers[16] = "Layer 17";
  upstream ? procNa_layers[17] = "Layer 17" : procNa_layers[17] = "Layer 18";
  upstream ? procNa_layers[18] = "Layer 18" : procNa_layers[18] = "Layer 19";
  upstream ? procNa_layers[19] = "Layer 19" : procNa_layers[19] = "Layer 20";
  upstream ? procNa_layers[20] = "Layer 20" : procNa_layers[20] = "Layer 21";
  upstream ? procNa_layers[21] = "Layer 21" : procNa_layers[21] = "Layer 22";
  upstream ? procNa_layers[22] = "Layer 22" : procNa_layers[22] = "Layer 23";
  upstream ? procNa_layers[23] = "Layer 23" : procNa_layers[23] = "Layer 24";
  upstream ? procNa_layers[24] = "Layer 24" : procNa_layers[24] = "Layer 25";
  upstream ? procNa_layers[25] = "Layer 25" : procNa_layers[25] = "Layer 26";
  upstream ? procNa_layers[26] = "Layer 26" : procNa_layers[26] = "Layer 27";
  upstream ? procNa_layers[27] = "Layer 27" : procNa_layers[27] = "Layer 28";
  upstream ? procNa_layers[28] = "Layer 28" : procNa_layers[28] = "Layer 29";
  upstream ? procNa_layers[29] = "Layer 29" : procNa_layers[29] = "Layer 30";
  if(upstream){procNa_layers[30] = "Layer 30";}

  //======================================================================================
  //The files that we will store the results of this analysis for the combined plot
  TString res_com[numberofconfigs];
  TFile* results_com[numberofconfigs];
  for (int k=0; k<numberofconfigs; k++){
    TString fp = "PFcal_" + IntToString(configs[k]);	
    res_com[k] = fp + "_combinedplots.root";
    //std::cout << res[k] << std::endl;
  }
  TH1F* hist1;
  TH2F *hist2[(*ssvec).size()];
  TH2F *currenthist1,*currenthist2,*currenthist3,*currenthist4,*currenthistsec1,*currenthistsec2,*currenthistsec3,*currenthistdedxsec1,*currenthistdedxsec2,*currenthistdedxsec3,*currenthistdedxsec4,*currenthistdedxsec5,*currenthistdedxsec6;
  TGraphErrors *gr2[(*ssvec).size()];
  TGraphErrors *gr_totalsf;
  // //use a TProfile to convert the 2-d to 1-d fit problem
  // TProfile *prof[(*ssvec).size()];
  //For accessing the fit results
  TF1 *fit2[(*ssvec).size()][numberofenergies];
  // TF1 *line0 = new TF1("line0","[0]*x+[1]",0.,20);
  // line0->SetParameter(0, 1);
  //line0->SetParameter(1, 0);
  TGraphErrors *gr_MIPSvsdEplusMIPs_fit[(*ssvec).size()];
  TGraphErrors *gr_MIPSvsdEplusMIPs_fit_c[(*ssvec).size()];
  // TMultiGraph* mg = new TMultiGraph();
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  TString titleofplot1,titleofplot2,titleofplot3,titleofplot4,titleofplot; 

  if (step2){

    // TH1F *emp = new TH1F("emp", " ", 400,0.,0.04);
    for (int k=0; k<numberofconfigs; k++){
      results[k]= TFile::Open(res[k],"read");
      std::cout << "Results file " << res[k] << std::endl;
    
      //======================================================================================
      //Loop on energies
      for (int l=0; l<numberofenergies; l++){
	// for (int l=0; l<1; l++){

	// emp->Reset();delete emp;
	// emp = new TH1F("emp", " ", 400,0.,0.04);
	//For the total sampling fraction
	//------------------------------------------------------------------------------------------------
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);
	c3->cd();
	c3->Update(); 
	hist1 = (TH1F*)results[k]->Get(("TotalSF_" + IntToString(particleenergies[l])).c_str());
	titleofplot1 = "Total sampling fraction for "; 
	titleofplot2 =  particle +  " beam particle gun"; 
	titleofplot = titleofplot1 + titleofplot2;
	// emp->SetTitle(titleofplot); 
	// emp->SetXTitle("SF = #sum_{i=1}^{30}Eactive_{i}/#sum_{i=1}^{30}(Epassive_{i}+Eactive_{i}) "); 
	// emp->SetYTitle("Events/0.0001");
	// emp->GetYaxis()->SetTitleOffset(1.4);
	// emp->Draw("axis");
	hist1->SetTitle(titleofplot); 
	hist1->GetXaxis()->SetTitle("SF = #sum_{i=1}^{30}Eactive_{i}/#sum_{i=1}^{30}(Epassive_{i}+Eactive_{i}) "); 
	hist1->GetYaxis()->SetTitle("Events/0.0001");
	hist1->GetXaxis()->SetTitleOffset(1.4);
	gPad->SetBottomMargin(0.14);
	// c3->SetLogy();
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
	hist1->GetYaxis()->SetRangeUser(0.,3000.);//  21: 100. 22: 900.
	// hist1->Draw("HISTsame");
	// leg1->Draw("same");
	l == 0 ? hist1->Draw("HIST") : hist1->Draw("HISTsame");
	l == 0 ? leg1->Draw() : leg1->Draw("same");
	c3->SaveAs("TotalSF_" + particle + ".png"); 
	c3->Update(); 


	if(looponlayersforsffit){

	//======================================================================================
	//Loop on layers
	//Starting from 1 when there is upstream material since in that case 0 is before calorimeter
	for(unsigned iL(startlayer); iL<(*ssvec).size(); iL++){
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
	  if ( l == 11){hist2[(int) iL]->SetMarkerColor(41);hist2[(int) iL]->SetMarkerStyle(26);}
	  if ( l == 12){hist2[(int) iL]->SetMarkerColor(38);hist2[(int) iL]->SetMarkerStyle(2);}
	  if ( l == 13){hist2[(int) iL]->SetMarkerColor(40);hist2[(int) iL]->SetMarkerStyle(3);}
	  if ( l == 14){hist2[(int) iL]->SetMarkerColor(30);hist2[(int) iL]->SetMarkerStyle(5);}
	  if ( l == 15){hist2[(int) iL]->SetMarkerColor(24);hist2[(int) iL]->SetMarkerStyle(21);}
	  if ( l == 16){hist2[(int) iL]->SetMarkerColor(29);hist2[(int) iL]->SetMarkerStyle(22);}
	  if ( l == 17){hist2[(int) iL]->SetMarkerColor(49);hist2[(int) iL]->SetMarkerStyle(32);}
	  //hist2[(int) iL]->SetLineColor(l+1);
	  leg_sf2[(int) iL]->AddEntry(hist2[(int) iL], procNa[l], "P");
	  // leg_sf2[(int) iL]->AddEntry(gr[(int) iL], procNa[l], "P");
	  hist2[(int) iL]->GetYaxis()->SetRangeUser(0.,2.);
	  hist2[(int) iL]->GetXaxis()->SetRangeUser(0.,60.); 
	  c9[(int) iL]->cd();
	  c9[(int) iL]->Update(); 
	  // if (l==5){mg->Draw("ap same");}
	  l == 0 ? hist2[(int) iL]->Draw("HIST") : hist2[(int) iL]->Draw("HISTsame");
	  l == 0 ? leg_sf2[(int) iL]->Draw() : leg_sf2[(int) iL]->Draw("same");
	  // l == 0 ? fit[(int) iL]->Draw() : fit[(int) iL]->Draw("same");
	  c9[(int) iL]->Update(); 
	  c9[(int) iL]->SaveAs("MIPSvsdEplusMIPs_"+particle+"_layer"+ IntToString( ((int) iL) ) + ".png");
	
	  char buf[500];
	  TLatex lat;
	  double latx = 0.;
	  double laty = 0.;

	  TString cantopng1 = "dEvsMIPs_plusfit_"+ particle;
	  TString cantopng2 = "_" + IntToString(particleenergies[l]);
	  TString cantopng3 = "GeV_layer" + IntToString((int) iL);
	  TString cantopng = cantopng1 + cantopng2 + cantopng3 + ".png";
	
	  if (dorunfit){
	  
	    //For the fit 
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
	    if ( l == 11){fit2[(int) iL][l]->SetLineColor(41);gr2[(int) iL]->SetMarkerColor(41);gr2[(int) iL]->SetMarkerStyle(26);}
	    if ( l == 12){fit2[(int) iL][l]->SetLineColor(38);gr2[(int) iL]->SetMarkerColor(38);gr2[(int) iL]->SetMarkerStyle(2);}
	    if ( l == 13){fit2[(int) iL][l]->SetLineColor(40);gr2[(int) iL]->SetMarkerColor(40);gr2[(int) iL]->SetMarkerStyle(3);}
	    if ( l == 14){fit2[(int) iL][l]->SetLineColor(30);gr2[(int) iL]->SetMarkerColor(30);gr2[(int) iL]->SetMarkerStyle(5);}
	    if ( l == 15){fit2[(int) iL][l]->SetLineColor(24);gr2[(int) iL]->SetMarkerColor(24);gr2[(int) iL]->SetMarkerStyle(21);}
	    if ( l == 16){fit2[(int) iL][l]->SetLineColor(29);gr2[(int) iL]->SetMarkerColor(29);gr2[(int) iL]->SetMarkerStyle(22);}
	    if ( l == 17){fit2[(int) iL][l]->SetLineColor(49);gr2[(int) iL]->SetMarkerColor(49);gr2[(int) iL]->SetMarkerStyle(32);}
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

	  }





	}//Loop on layers

	}//closing if 
      
      }//Loop on energies

      results_com[k]= new TFile(res_com[k],"recreate");
    
      //======================================================================================
      //Write canvas with combined plot 
      // c1->Print("Leakage.pdf",".pdf");
      for(unsigned iL(startlayer); iL<(*ssvec).size(); iL++){
	c1[(int) iL]->Write();
	c9[(int) iL]->Write();
	c2[(int) iL]->Write();
	for (int l=0; l<numberofenergies; l++){
	  c4[(int) iL][l]->Write();
	  c10[(int) iL][l]->Write();
	}
	// gr[(int) iL]->Write();
	if(dorunfit){
	  gr2[(int) iL]->Write();
	}
      
      }
      
      c3->Write();

      //Reset histos
      // hist->Reset();

      //======================================================================================
      //Save the fit results for faster running
      std::ofstream myfile;
      TString titleoffile = "MIPSvsdEplusMIPs_" + particle + ".txt";
      myfile.open(titleoffile);

      if(dorunfit){
	//Here print the p0 and p1 that are going to be used for Ereco
	for (int l=0; l<numberofenergies; l++){
	  std::cout << "=========================================================" << std::endl;
	  std::cout << "Particle " <<  particle << std::endl;
	  std::cout << "Energy of the beam " <<  particleenergies[l] << std::endl;
	  for(unsigned iL(startlayer); iL<(*ssvec).size(); iL++){
	    std::cout << "---------------------------------------------------------" << std::endl;
	    std::cout << "Layer " << iL << std::endl;
	    std::cout << "MIPSvsdEplusMIPs slope: " <<  fit2[(int) iL][l]->GetParameter(1) << " error slope " << fit2[(int) iL][l]->GetParError(1) << std::endl;
	    std::cout << "MIPSvsdEplusMIPs constant: " <<  fit2[(int) iL][l]->GetParameter(0) << " error constant " << fit2[(int) iL][l]->GetParError(0) << std::endl;

	    //particle particleenergy layer constant slope errorconstant errorslope
	    myfile << particle << " " << particleenergies[l] << " " << iL << " " << fit2[(int) iL][l]->GetParameter(0) << " " << fit2[(int) iL][l]->GetParameter(1) << " " << fit2[(int) iL][l]->GetParError(0) <<  " " << fit2[(int) iL][l]->GetParError(1) << "\n";


	  }
	}

      }
      //======================================================================================
      //Total sampling fraction vs energy
      c7->cd();
      gr_totalsf = (TGraphErrors *) calib->Clone("Totalsf");
      gr_totalsf->SetMarkerStyle( 20 ); 
      gPad->SetLeftMargin(0.14);
      // gr_totalsf->SetMarkerSize(0.2); 
      gr_totalsf->SetMarkerColor(1);  
      gr_totalsf->SetLineColor(1);
      c7->Update(); 

      // Int_t np=gr_totalsf->GetN();
      // gr_totalsf->SetPoint(numberofenergies, particleenergies[l] , fit[l]->GetParameter(0) );
      // gr_totalsf->SetPointError(numberofenergies, 0. , fit[l]->GetParError(0) );
      for (int l=0; l<numberofenergies; l++){
	hist1 = (TH1F*)results[k]->Get(("TotalSF_" + IntToString(particleenergies[l])).c_str());
	gr_totalsf->SetPoint( l , particleenergies[l] , hist1->GetMean() );
	gr_totalsf->SetPointError(l , 0. , hist1->GetRMS() );
      }
      gr_totalsf->Draw("APL");
      c7->Update(); 

      gr_totalsf->SetTitle( "Total sampling fraction for " + particle + " for all beam energies" );
      gr_totalsf->GetHistogram()->SetXTitle("Beam Energy (GeV)"); 
      gr_totalsf->GetHistogram()->SetYTitle("SF = #sum_{i=1}^{30}Eactive_{i}/#sum_{i=1}^{30}(Epassive_{i}+Eactive_{i}) ");
      gr_totalsf->GetYaxis()->SetTitleOffset(1.7);
      gr_totalsf->GetXaxis()->SetRangeUser(0.,520.);  
      gr_totalsf->GetYaxis()->SetRangeUser(0.006,0.020);  

      //create a transparent pad drawn on top of the main pad
      TPad *overlay = new TPad("overlay","",0,0,1,1);
      overlay->SetFillStyle(4000);
      overlay->SetFillColor(0);
      overlay->SetFrameFillStyle(4000);
      overlay->Draw();
      overlay->cd();
      gPad->SetLeftMargin(0.14);
      TH1F *hframe;
      TGaxis *axis ;
      // overlay->Range(0.,0.006,520.,0.0016);

      // double xmin = pad->GetUxmin();
      // double ymin = pad->GetUymin();
      double ymin = 0.006/SF_dedx;
      // double xmax = pad->GetUxmax();
      // double ymax = pad->GetUymax();
      double ymax = 0.020/SF_dedx;
      
      // Double_t dy = (ymax-ymin)/0.8;

      // hframe =  overlay->DrawFrame(xmin,ymin,xmax,ymax);
      // hframe =  overlay->DrawFrame(0.,ymin,520.,ymax);
      hframe =  overlay->DrawFrame(0.,0.,520.,10.);
      hframe->GetXaxis()->SetLabelOffset(99);
      hframe->GetYaxis()->SetLabelOffset(99);
      hframe->GetYaxis()->SetTickLength(0.001);
      hframe->GetXaxis()->SetTickLength(0.001);

      axis = new TGaxis(520.,0.,520,10.,ymin,ymax,510,"+L");
      axis->SetLineColor(4);
      axis->SetLabelColor(4);
      axis->SetTitle("e/mip");
      axis->SetTitleColor(4);
      axis->Draw();

      c7->Update(); 
      TString cnpng = "TotalSamplingFractionvsEnergy_"+ particle + ".png";
      c7->SaveAs(cnpng);

      c7->Write();

      //======================================================================================
      //Total sampling fraction vs shower depth in radiation length
      c5->cd();
      pad->Draw();
      pad->cd();
      gPad->SetLeftMargin(0.14);
      for (int l=0; l<numberofenergies; l++){
	//Totalsfvsshowerdepth[l] = (TGraphErrors*)results[k]->Get(auxnam_tgr.c_str());
	if (l == 0){Totalsfvsshowerdepth[l]->SetMarkerColor(4);Totalsfvsshowerdepth[l]->SetMarkerStyle(20);}
	if (l == 1){Totalsfvsshowerdepth[l]->SetMarkerColor(2);Totalsfvsshowerdepth[l]->SetMarkerStyle(21);}
	if (l == 2){Totalsfvsshowerdepth[l]->SetMarkerColor(1);Totalsfvsshowerdepth[l]->SetMarkerStyle(22);}
	if (l == 3){Totalsfvsshowerdepth[l]->SetMarkerColor(3);Totalsfvsshowerdepth[l]->SetMarkerStyle(23);}
	if (l == 4){Totalsfvsshowerdepth[l]->SetMarkerColor(5);Totalsfvsshowerdepth[l]->SetMarkerStyle(28);}
	if (l == 5){Totalsfvsshowerdepth[l]->SetMarkerColor(6);Totalsfvsshowerdepth[l]->SetMarkerStyle(29);}
	if (l == 6){Totalsfvsshowerdepth[l]->SetMarkerColor(7);Totalsfvsshowerdepth[l]->SetMarkerStyle(30);}
	if (l == 7){Totalsfvsshowerdepth[l]->SetMarkerColor(8);Totalsfvsshowerdepth[l]->SetMarkerStyle(31);}
	if (l == 8){Totalsfvsshowerdepth[l]->SetMarkerColor(9);Totalsfvsshowerdepth[l]->SetMarkerStyle(33);}
	if (l == 9){Totalsfvsshowerdepth[l]->SetMarkerColor(12);Totalsfvsshowerdepth[l]->SetMarkerStyle(34);}
	if (l == 10){Totalsfvsshowerdepth[l]->SetMarkerColor(46);Totalsfvsshowerdepth[l]->SetMarkerStyle(25);}
	if (l == 11){Totalsfvsshowerdepth[l]->SetMarkerColor(41);Totalsfvsshowerdepth[l]->SetMarkerStyle(26);}
	if (l == 12){Totalsfvsshowerdepth[l]->SetMarkerColor(38);Totalsfvsshowerdepth[l]->SetMarkerStyle(2);}
	if (l == 13){Totalsfvsshowerdepth[l]->SetMarkerColor(40);Totalsfvsshowerdepth[l]->SetMarkerStyle(3);}
	if (l == 14){Totalsfvsshowerdepth[l]->SetMarkerColor(30);Totalsfvsshowerdepth[l]->SetMarkerStyle(5);}
	if (l == 15){Totalsfvsshowerdepth[l]->SetMarkerColor(24);Totalsfvsshowerdepth[l]->SetMarkerStyle(21);}
	if (l == 16){Totalsfvsshowerdepth[l]->SetMarkerColor(29);Totalsfvsshowerdepth[l]->SetMarkerStyle(22);}
	if (l == 17){Totalsfvsshowerdepth[l]->SetMarkerColor(49);Totalsfvsshowerdepth[l]->SetMarkerStyle(32);}
 
	Totalsfvsshowerdepth[l]->SetTitle( "Total sampling fraction vs. "+shdp+" for " + particle); 
	Totalsfvsshowerdepth[l]->GetHistogram()->SetXTitle(shdp); 
	Totalsfvsshowerdepth[l]->GetHistogram()->SetYTitle("SF = #sum_{i=1}^{30}Eactive_{i}/#sum_{i=1}^{30}(Epassive_{i}+Eactive_{i})");
	Totalsfvsshowerdepth[l]->GetYaxis()->SetTitleOffset(1.7);
	leg_tsfvst->AddEntry(Totalsfvsshowerdepth[l], procNa[l], "P");

	l == 0 ? Totalsfvsshowerdepth[l]->Draw("AP") : Totalsfvsshowerdepth[l]->Draw("PS");
	l == 0 ? leg_tsfvst->Draw() : leg_tsfvst->Draw("same");
	// Totalsfvsshowerdepth[0]->GetXaxis()->SetRangeUser(0.,25.); 
	// Totalsfvsshowerdepth[0]->GetYaxis()->SetRangeUser(0.008,0.00135); 
	c5->Update(); 
	Totalsfvsshowerdepth[l]->Write();
      }
      Totalsfvsshowerdepth[0]->GetXaxis()->SetRangeUser(0.,25.); 
      c5->Update(); 

      cnpng = "Totalsfvsshowerdepth.png";
      c5->SaveAs(cnpng);
      c5->Write();

     //======================================================================================
      //Total sampling fraction vs shower depth in radiation length with layer as unit
      c6->cd();
      pad->Draw();
      pad->cd();
      gPad->SetLeftMargin(0.14);
      for (int l=0; l<numberofenergies; l++){
	//Totalsfvsshowerdepthlay[l] = (TGraphErrors*)results[k]->Get(auxnam_tgr.c_str());
	if (l == 0){Totalsfvsshowerdepthlay[l]->SetMarkerColor(4);Totalsfvsshowerdepthlay[l]->SetMarkerStyle(20);}
	if (l == 1){Totalsfvsshowerdepthlay[l]->SetMarkerColor(2);Totalsfvsshowerdepthlay[l]->SetMarkerStyle(21);}
	if (l == 2){Totalsfvsshowerdepthlay[l]->SetMarkerColor(1);Totalsfvsshowerdepthlay[l]->SetMarkerStyle(22);}
	if (l == 3){Totalsfvsshowerdepthlay[l]->SetMarkerColor(3);Totalsfvsshowerdepthlay[l]->SetMarkerStyle(23);}
	if (l == 4){Totalsfvsshowerdepthlay[l]->SetMarkerColor(5);Totalsfvsshowerdepthlay[l]->SetMarkerStyle(28);}
	if (l == 5){Totalsfvsshowerdepthlay[l]->SetMarkerColor(6);Totalsfvsshowerdepthlay[l]->SetMarkerStyle(29);}
	if (l == 6){Totalsfvsshowerdepthlay[l]->SetMarkerColor(7);Totalsfvsshowerdepthlay[l]->SetMarkerStyle(30);}
	if (l == 7){Totalsfvsshowerdepthlay[l]->SetMarkerColor(8);Totalsfvsshowerdepthlay[l]->SetMarkerStyle(31);}
	if (l == 8){Totalsfvsshowerdepthlay[l]->SetMarkerColor(9);Totalsfvsshowerdepthlay[l]->SetMarkerStyle(33);}
	if (l == 9){Totalsfvsshowerdepthlay[l]->SetMarkerColor(12);Totalsfvsshowerdepthlay[l]->SetMarkerStyle(34);}
	if (l == 10){Totalsfvsshowerdepthlay[l]->SetMarkerColor(46);Totalsfvsshowerdepthlay[l]->SetMarkerStyle(25);}
	if (l == 11){Totalsfvsshowerdepthlay[l]->SetMarkerColor(41);Totalsfvsshowerdepthlay[l]->SetMarkerStyle(26);}
	if ( l == 12){Totalsfvsshowerdepthlay[l]->SetMarkerColor(38);Totalsfvsshowerdepthlay[l]->SetMarkerStyle(2);}
	if ( l == 13){Totalsfvsshowerdepthlay[l]->SetMarkerColor(40);Totalsfvsshowerdepthlay[l]->SetMarkerStyle(3);}
	if ( l == 14){Totalsfvsshowerdepthlay[l]->SetMarkerColor(30);Totalsfvsshowerdepthlay[l]->SetMarkerStyle(5);}
	if ( l == 15){Totalsfvsshowerdepthlay[l]->SetMarkerColor(24);Totalsfvsshowerdepthlay[l]->SetMarkerStyle(21);}
	if ( l == 16){Totalsfvsshowerdepthlay[l]->SetMarkerColor(29);Totalsfvsshowerdepthlay[l]->SetMarkerStyle(22);}
	if ( l == 17){Totalsfvsshowerdepthlay[l]->SetMarkerColor(49);Totalsfvsshowerdepthlay[l]->SetMarkerStyle(32);}
 
	Totalsfvsshowerdepthlay[l]->SetTitle( "Total sampling fraction vs. shower depth for " + particle); 
	Totalsfvsshowerdepthlay[l]->GetHistogram()->SetXTitle("Shower depth (layer)"); 
	Totalsfvsshowerdepthlay[l]->GetHistogram()->SetYTitle("SF = #sum_{i=1}^{30}Eactive_{i}/#sum_{i=1}^{30}(Epassive_{i}+Eactive_{i})");
	Totalsfvsshowerdepthlay[l]->GetYaxis()->SetTitleOffset(1.7);
	leg_tsfvstlay->AddEntry(Totalsfvsshowerdepthlay[l], procNa[l], "P");

	l == 0 ? Totalsfvsshowerdepthlay[l]->Draw("AP") : Totalsfvsshowerdepthlay[l]->Draw("PS");
	l == 0 ? leg_tsfvstlay->Draw() : leg_tsfvstlay->Draw("same");
	Totalsfvsshowerdepthlay[0]->GetXaxis()->SetRangeUser(1.,30.); 
	// Totalsfvsshowerdepthlay[0]->GetYaxis()->SetRangeUser(0.008,0.00135); 
	c6->Update(); 
	Totalsfvsshowerdepthlay[l]->Write();
      }

      cnpng = "Totalsfvsshowerdepthlay.png";
      c6->SaveAs(cnpng);
      c6->Write();

      //======================================================================================
      //Total sampling fraction vs X0max/<X0max in radiation length : Profile
      c8->cd();
      pad->Draw();
      pad->cd();
      gPad->SetLeftMargin(0.14);
      // THStack *hs= new THStack("hs","profile");

      // Empty histo
      TH2F * empty = new TH2F("empty", " ", 100,0.,maxbsh,100,0.006,0.020);
      empty->SetTitle( "Profile of total sampling fraction vs. "+shdp+" for " + particle); 
      empty->SetXTitle(shdp); 
      empty->SetYTitle("SF = #sum_{i=1}^{30}Eactive_{i}/#sum_{i=1}^{30}(Epassive_{i}+Eactive_{i})");
      empty->GetYaxis()->SetTitleOffset(1.7);
      empty->Draw("axis");

      TProfile *currentprof;
      for (int l=0; l<numberofenergies; l++){
	//Totalsfvsshowerdepth_forprofile[l]->ProfileX() = (TGraphErrors*)results[k]->Get(auxnam_tgr.c_str());
	currenthist1 = (TH2F*)results[k]->Get(("TotalSFvsShowerDepth_" + IntToString(particleenergies[l])).c_str());
	// currentprof = Totalsfvsshowerdepth_forprofile[l]->ProfileX();
	currentprof = currenthist1->ProfileX();
	if (l == 0){currentprof->SetMarkerColor(4);currentprof->SetMarkerStyle(20);}
	if (l == 1){currentprof->SetMarkerColor(2);currentprof->SetMarkerStyle(21);}
	if (l == 2){currentprof->SetMarkerColor(1);currentprof->SetMarkerStyle(22);}
	if (l == 3){currentprof->SetMarkerColor(3);currentprof->SetMarkerStyle(23);}
	if (l == 4){currentprof->SetMarkerColor(5);currentprof->SetMarkerStyle(28);}
	if (l == 5){currentprof->SetMarkerColor(6);currentprof->SetMarkerStyle(29);}
	if (l == 6){currentprof->SetMarkerColor(7);currentprof->SetMarkerStyle(30);}
	if (l == 7){currentprof->SetMarkerColor(8);currentprof->SetMarkerStyle(31);}
	if (l == 8){currentprof->SetMarkerColor(9);currentprof->SetMarkerStyle(33);}
	if (l == 9){currentprof->SetMarkerColor(12);currentprof->SetMarkerStyle(34);}
	if (l == 10){currentprof->SetMarkerColor(46);currentprof->SetMarkerStyle(25);}
	if (l == 11){currentprof->SetMarkerColor(41);currentprof->SetMarkerStyle(26);}
 	if (l == 12){currentprof->SetMarkerColor(38);currentprof->SetMarkerStyle(2);}
	if (l == 13){currentprof->SetMarkerColor(40);currentprof->SetMarkerStyle(3);}
	if (l == 14){currentprof->SetMarkerColor(30);currentprof->SetMarkerStyle(5);}
	if (l == 15){currentprof->SetMarkerColor(24);currentprof->SetMarkerStyle(21);}
	if (l == 16){currentprof->SetMarkerColor(29);currentprof->SetMarkerStyle(22);}
	if (l == 17){currentprof->SetMarkerColor(49);currentprof->SetMarkerStyle(32);}

	// currentprof->SetTitle( "Profile of total sampling fraction vs. "+shdp+" for " + particle); 
	// currentprof->SetXTitle(shdp); 
	// currentprof->SetYTitle("SF = #sum_{i=1}^{30}Eactive_{i}/#sum_{i=1}^{30}(Epassive_{i}+Eactive_{i})");
	// currentprof->GetYaxis()->SetTitleOffset(1.7);
	leg_tsfvst_prof->AddEntry(currentprof, procNa[l], "P");

	// l == 0 ? currentprof->Draw() : currentprof->Draw("same");
	// l == 0 ? leg_tsfvst_prof->Draw() : leg_tsfvst_prof->Draw("same");
	// currenthist1->SetLineColor(0);
	// currenthist1->Draw();
	// currenthist1->GetYaxis()->SetRangeUser(0.,0.022);

	currentprof->Draw("same");
	// l == 0 ? leg_tsfvst_prof->Draw() : leg_tsfvst_prof->Draw("same");
	leg_tsfvst_prof->Draw("same");


	// currentprof->GetYaxis()->SetRangeUser(0.006,0.0016); 
	// pad->Update();
	// currentprof->GetYaxis()->SetLimits(0.006,0.0016); 
	// Totalsfvsshowerdepth_forprofile[l]->ProfileX()->GetXaxis()->SetRangeUser(0.,.); 
	// Totalsfvsshowerdepth_forprofile[0]->ProfileX()->GetYaxis()->SetRangeUser(0.006,0.0016); 
	// hs->Add(currentprof);
	// if (l==0){
	//   hs.SetHistogram(currentprof);
	// }
	// leg_tsfvst_prof->AddEntry(currentprof, procNa[l], "P");

	currentprof->Write();
      }

      // hs->Draw("nostack");
   
      // THStack *hs_n = (THStack *)0;
      // hs_n = new THStack(*hs,"prof"); // do not modify the "hs_ORIG" THStack, create a copy
      // hs_n->Draw("H NOSTACK");
      // hs->GetHistogram()->SetXTitle("My X Title");
      // hs->GetHistogram()->SetYTitle("My Y Title");
      // hs->GetXaxis()->SetTitleOffset(1.0);
      // hs->GetYaxis()->SetTitleOffset(1.5);
      // hs->GetYaxis()->CenterTitle(kTRUE);
      // hs->GetXaxis()->SetLimits(0, 10.001); // "xmin", "xmax"
      // hs_n->SetMinimum(0.0006); // "ymin"
      // hs_n->SetMaximum(0.0016); // "ymax"
      // leg_tsfvst_prof->Draw();
      // hs.GetYaxis()->SetRangeUser(0.006,0.0010);
      

      c8->Update(); 

      cnpng = "Totalsfvsshowerdepth_profile.png";
      c8->SaveAs(cnpng);
      c8->Write();
     //======================================================================================
      //Total sampling fraction for section 1 vs X0max/<X0max in radiation length : Profile
      c18->cd();
      pad->Draw();
      pad->cd();
      gPad->SetLeftMargin(0.14);

      // Empty histo
      TH2F *emptysec1 = new TH2F("emptysec1", " ", 100,0.,maxbsh,400,0.,0.04);
      emptysec1->SetTitle( "Profile of total sampling fraction for section 1 vs. "+shdp+" for " + particle); 
      emptysec1->SetXTitle(shdp); 
      emptysec1->SetYTitle("SF_{sec1} = #sum_{i=sec1}Eactive_{i}/#sum_{i=sec1}(Epassive_{i}+Eactive_{i})");
      emptysec1->GetYaxis()->SetTitleOffset(1.7);
      emptysec1->Draw("axis");

      TProfile *currentprof_sec1;
      for (int l=0; l<numberofenergies; l++){
	currenthistsec1 = (TH2F*)results[k]->Get(("TotalSF_sec1vsShowerDepth_" + IntToString(particleenergies[l])).c_str());
	currentprof_sec1 = currenthistsec1->ProfileX();
	if (l == 0){currentprof_sec1->SetMarkerColor(4);currentprof_sec1->SetMarkerStyle(20);}
	if (l == 1){currentprof_sec1->SetMarkerColor(2);currentprof_sec1->SetMarkerStyle(21);}
	if (l == 2){currentprof_sec1->SetMarkerColor(1);currentprof_sec1->SetMarkerStyle(22);}
	if (l == 3){currentprof_sec1->SetMarkerColor(3);currentprof_sec1->SetMarkerStyle(23);}
	if (l == 4){currentprof_sec1->SetMarkerColor(5);currentprof_sec1->SetMarkerStyle(28);}
	if (l == 5){currentprof_sec1->SetMarkerColor(6);currentprof_sec1->SetMarkerStyle(29);}
	if (l == 6){currentprof_sec1->SetMarkerColor(7);currentprof_sec1->SetMarkerStyle(30);}
	if (l == 7){currentprof_sec1->SetMarkerColor(8);currentprof_sec1->SetMarkerStyle(31);}
	if (l == 8){currentprof_sec1->SetMarkerColor(9);currentprof_sec1->SetMarkerStyle(33);}
	if (l == 9){currentprof_sec1->SetMarkerColor(12);currentprof_sec1->SetMarkerStyle(34);}
	if (l == 10){currentprof_sec1->SetMarkerColor(46);currentprof_sec1->SetMarkerStyle(25);}
	if (l == 11){currentprof_sec1->SetMarkerColor(41);currentprof_sec1->SetMarkerStyle(26);}
 	if (l == 12){currentprof_sec1->SetMarkerColor(38);currentprof_sec1->SetMarkerStyle(2);}
	if (l == 13){currentprof_sec1->SetMarkerColor(40);currentprof_sec1->SetMarkerStyle(3);}
	if (l == 14){currentprof_sec1->SetMarkerColor(30);currentprof_sec1->SetMarkerStyle(5);}
	if (l == 15){currentprof_sec1->SetMarkerColor(24);currentprof_sec1->SetMarkerStyle(21);}
	if (l == 16){currentprof_sec1->SetMarkerColor(29);currentprof_sec1->SetMarkerStyle(22);}
	if (l == 17){currentprof_sec1->SetMarkerColor(49);currentprof_sec1->SetMarkerStyle(32);}

	leg_tsfsec1vst_prof->AddEntry(currentprof_sec1, procNa[l], "P");
	currentprof_sec1->Draw("same");
	leg_tsfsec1vst_prof->Draw("same");
	currentprof_sec1->Write();
      }

      c18->Update(); 

      cnpng = "Totalsfsec1vsshowerdepth_profile.png";
      c18->SaveAs(cnpng);
      c18->Write();

      

      //======================================================================================
      //Total sampling fraction for section 2 vs X0max/<X0max in radiation length : Profile
      c19->cd();
      // pad->Draw();
      // pad->cd();
      gPad->SetLeftMargin(0.14);

      // Empty histo
      TH2F *emptysec2 = new TH2F("emptysec2", " ", 100,0.,maxbsh,400,0.,0.04);
      emptysec2->SetTitle( "Profile of total sampling fraction for section 2 vs. "+shdp+" for " + particle); 
      emptysec2->SetXTitle(shdp); 
      emptysec2->SetYTitle("SF_{sec2} = #sum_{i=sec2}Eactive_{i}/#sum_{i=sec2}(Epassive_{i}+Eactive_{i})");
      emptysec2->GetYaxis()->SetTitleOffset(1.7);
      emptysec2->Draw("axis");

      TProfile *currentprof_sec2;
      for (int l=0; l<numberofenergies; l++){
	currenthistsec2 = (TH2F*)results[k]->Get(("TotalSF_sec2vsShowerDepth_" + IntToString(particleenergies[l])).c_str());
	currentprof_sec2 = currenthistsec2->ProfileX();
	if (l == 0){currentprof_sec2->SetMarkerColor(4);currentprof_sec2->SetMarkerStyle(20);}
	if (l == 1){currentprof_sec2->SetMarkerColor(2);currentprof_sec2->SetMarkerStyle(21);}
	if (l == 2){currentprof_sec2->SetMarkerColor(1);currentprof_sec2->SetMarkerStyle(22);}
	if (l == 3){currentprof_sec2->SetMarkerColor(3);currentprof_sec2->SetMarkerStyle(23);}
	if (l == 4){currentprof_sec2->SetMarkerColor(5);currentprof_sec2->SetMarkerStyle(28);}
	if (l == 5){currentprof_sec2->SetMarkerColor(6);currentprof_sec2->SetMarkerStyle(29);}
	if (l == 6){currentprof_sec2->SetMarkerColor(7);currentprof_sec2->SetMarkerStyle(30);}
	if (l == 7){currentprof_sec2->SetMarkerColor(8);currentprof_sec2->SetMarkerStyle(31);}
	if (l == 8){currentprof_sec2->SetMarkerColor(9);currentprof_sec2->SetMarkerStyle(33);}
	if (l == 9){currentprof_sec2->SetMarkerColor(12);currentprof_sec2->SetMarkerStyle(34);}
	if (l == 10){currentprof_sec2->SetMarkerColor(46);currentprof_sec2->SetMarkerStyle(25);}
	if (l == 11){currentprof_sec2->SetMarkerColor(41);currentprof_sec2->SetMarkerStyle(26);}
 	if (l == 12){currentprof_sec2->SetMarkerColor(38);currentprof_sec2->SetMarkerStyle(2);}
	if (l == 13){currentprof_sec2->SetMarkerColor(40);currentprof_sec2->SetMarkerStyle(3);}
	if (l == 14){currentprof_sec2->SetMarkerColor(30);currentprof_sec2->SetMarkerStyle(5);}
	if (l == 15){currentprof_sec2->SetMarkerColor(24);currentprof_sec2->SetMarkerStyle(21);}
	if (l == 16){currentprof_sec2->SetMarkerColor(29);currentprof_sec2->SetMarkerStyle(22);}
	if (l == 17){currentprof_sec2->SetMarkerColor(49);currentprof_sec2->SetMarkerStyle(32);}

	leg_tsfsec2vst_prof->AddEntry(currentprof_sec2, procNa[l], "P");
	currentprof_sec2->Draw("same");
	leg_tsfsec2vst_prof->Draw("same");
	currentprof_sec2->Write();
      }

      c19->Update(); 

      cnpng = "Totalsfsec2vsshowerdepth_profile.png";
      c19->SaveAs(cnpng);
      c19->Write();




     //======================================================================================
      //Total sampling fraction for section 3 vs X0max/<X0max in radiation length : Profile
      c20->cd();
      // pad->Draw();
      // pad->cd();
      gPad->SetLeftMargin(0.14);

      // Empty histo
      TH2F *emptysec3 = new TH2F("emptysec3", " ", 100,0.,maxbsh,400,0.,0.04);
      emptysec3->SetTitle( "Profile of total sampling fraction for section 3 vs. "+shdp+" for " + particle); 
      emptysec3->SetXTitle(shdp); 
      emptysec3->SetYTitle("SF_{sec3} = #sum_{i=sec3}Eactive_{i}/#sum_{i=sec3}(Epassive_{i}+Eactive_{i})");
      emptysec3->GetYaxis()->SetTitleOffset(1.7);
      emptysec3->Draw("axis");

      TProfile *currentprof_sec3;
      for (int l=0; l<numberofenergies; l++){
	currenthistsec3 = (TH2F*)results[k]->Get(("TotalSF_sec3vsShowerDepth_" + IntToString(particleenergies[l])).c_str());
	currentprof_sec3 = currenthistsec3->ProfileX();
	if (l == 0){currentprof_sec3->SetMarkerColor(4);currentprof_sec3->SetMarkerStyle(20);}
	if (l == 1){currentprof_sec3->SetMarkerColor(2);currentprof_sec3->SetMarkerStyle(21);}
	if (l == 2){currentprof_sec3->SetMarkerColor(1);currentprof_sec3->SetMarkerStyle(22);}
	if (l == 3){currentprof_sec3->SetMarkerColor(3);currentprof_sec3->SetMarkerStyle(23);}
	if (l == 4){currentprof_sec3->SetMarkerColor(5);currentprof_sec3->SetMarkerStyle(28);}
	if (l == 5){currentprof_sec3->SetMarkerColor(6);currentprof_sec3->SetMarkerStyle(29);}
	if (l == 6){currentprof_sec3->SetMarkerColor(7);currentprof_sec3->SetMarkerStyle(30);}
	if (l == 7){currentprof_sec3->SetMarkerColor(8);currentprof_sec3->SetMarkerStyle(31);}
	if (l == 8){currentprof_sec3->SetMarkerColor(9);currentprof_sec3->SetMarkerStyle(33);}
	if (l == 9){currentprof_sec3->SetMarkerColor(12);currentprof_sec3->SetMarkerStyle(34);}
	if (l == 10){currentprof_sec3->SetMarkerColor(46);currentprof_sec3->SetMarkerStyle(25);}
	if (l == 11){currentprof_sec3->SetMarkerColor(41);currentprof_sec3->SetMarkerStyle(26);}
 	if (l == 12){currentprof_sec3->SetMarkerColor(38);currentprof_sec3->SetMarkerStyle(2);}
	if (l == 13){currentprof_sec3->SetMarkerColor(40);currentprof_sec3->SetMarkerStyle(3);}
	if (l == 14){currentprof_sec3->SetMarkerColor(30);currentprof_sec3->SetMarkerStyle(5);}
	if (l == 15){currentprof_sec3->SetMarkerColor(24);currentprof_sec3->SetMarkerStyle(21);}
	if (l == 16){currentprof_sec3->SetMarkerColor(29);currentprof_sec3->SetMarkerStyle(22);}
	if (l == 17){currentprof_sec3->SetMarkerColor(49);currentprof_sec3->SetMarkerStyle(32);}

	leg_tsfsec3vst_prof->AddEntry(currentprof_sec3, procNa[l], "P");
	currentprof_sec3->Draw("same");
	leg_tsfsec3vst_prof->Draw("same");
	currentprof_sec3->Write();
      }

      c20->Update(); 

      cnpng = "Totalsfsec3vsshowerdepth_profile.png";
      c20->SaveAs(cnpng);
      c20->Write();
     //======================================================================================
     //======================================================================================
      //Total sampling fraction for dedxsection 1 vs X0max/<X0max in radiation length : Profile
      c21->cd();
      pad->Draw();
      pad->cd();
      gPad->SetLeftMargin(0.14);

      // Empty histo
      TH2F *emptydedxsec1 = new TH2F("emptydedxsec1", " ", 100,0.,maxbsh,400,0.,0.04);
      emptydedxsec1->SetTitle( "Profile of total sampling fraction for layer 1 vs. "+shdp+" for " + particle); 
      emptydedxsec1->SetXTitle(shdp); 
      emptydedxsec1->SetYTitle("SF_{lay1} = #sum_{i=lay1}Eactive_{i}/#sum_{i=lay1}(Epassive_{i}+Eactive_{i})");
      emptydedxsec1->GetYaxis()->SetTitleOffset(1.7);
      emptydedxsec1->Draw("axis");

      TProfile *currentprof_dedxsec1;
      for (int l=0; l<numberofenergies; l++){
	currenthistdedxsec1 = (TH2F*)results[k]->Get(("TotalSF_dedxsec1vsShowerDepth_" + IntToString(particleenergies[l])).c_str());
	currentprof_dedxsec1 = currenthistdedxsec1->ProfileX();
	if (l == 0){currentprof_dedxsec1->SetMarkerColor(4);currentprof_dedxsec1->SetMarkerStyle(20);}
	if (l == 1){currentprof_dedxsec1->SetMarkerColor(2);currentprof_dedxsec1->SetMarkerStyle(21);}
	if (l == 2){currentprof_dedxsec1->SetMarkerColor(1);currentprof_dedxsec1->SetMarkerStyle(22);}
	if (l == 3){currentprof_dedxsec1->SetMarkerColor(3);currentprof_dedxsec1->SetMarkerStyle(23);}
	if (l == 4){currentprof_dedxsec1->SetMarkerColor(5);currentprof_dedxsec1->SetMarkerStyle(28);}
	if (l == 5){currentprof_dedxsec1->SetMarkerColor(6);currentprof_dedxsec1->SetMarkerStyle(29);}
	if (l == 6){currentprof_dedxsec1->SetMarkerColor(7);currentprof_dedxsec1->SetMarkerStyle(30);}
	if (l == 7){currentprof_dedxsec1->SetMarkerColor(8);currentprof_dedxsec1->SetMarkerStyle(31);}
	if (l == 8){currentprof_dedxsec1->SetMarkerColor(9);currentprof_dedxsec1->SetMarkerStyle(33);}
	if (l == 9){currentprof_dedxsec1->SetMarkerColor(12);currentprof_dedxsec1->SetMarkerStyle(34);}
	if (l == 10){currentprof_dedxsec1->SetMarkerColor(46);currentprof_dedxsec1->SetMarkerStyle(25);}
	if (l == 11){currentprof_dedxsec1->SetMarkerColor(41);currentprof_dedxsec1->SetMarkerStyle(26);}
 	if (l == 12){currentprof_dedxsec1->SetMarkerColor(38);currentprof_dedxsec1->SetMarkerStyle(2);}
	if (l == 13){currentprof_dedxsec1->SetMarkerColor(40);currentprof_dedxsec1->SetMarkerStyle(3);}
	if (l == 14){currentprof_dedxsec1->SetMarkerColor(30);currentprof_dedxsec1->SetMarkerStyle(5);}
	if (l == 15){currentprof_dedxsec1->SetMarkerColor(24);currentprof_dedxsec1->SetMarkerStyle(21);}
	if (l == 16){currentprof_dedxsec1->SetMarkerColor(29);currentprof_dedxsec1->SetMarkerStyle(22);}
	if (l == 17){currentprof_dedxsec1->SetMarkerColor(49);currentprof_dedxsec1->SetMarkerStyle(32);}

	leg_tsfdedxsec1vst_prof->AddEntry(currentprof_dedxsec1, procNa[l], "P");
	currentprof_dedxsec1->Draw("same");
	leg_tsfdedxsec1vst_prof->Draw("same");
	currentprof_dedxsec1->Write();
      }

      c21->Update(); 

      cnpng = "Totalsfdedxsec1vsshowerdepth_profile.png";
      c21->SaveAs(cnpng);
      c21->Write();

      

      //======================================================================================
      //Total sampling fraction for dedxsection 2 vs X0max/<X0max in radiation length : Profile
      c22->cd();
      // pad->Draw();
      // pad->cd();
      gPad->SetLeftMargin(0.14);

      // Empty histo
      TH2F *emptydedxsec2 = new TH2F("emptydedxsec2", " ", 100,0.,maxbsh,400,0.,0.04);
      emptydedxsec2->SetTitle( "Profile of total sampling fraction for section 1 vs. "+shdp+" for " + particle); 
      emptydedxsec2->SetXTitle(shdp); 
      emptydedxsec2->SetYTitle("SF_{sec1} = #sum_{i=sec1}Eactive_{i}/#sum_{i=sec1}(Epassive_{i}+Eactive_{i})");
      emptydedxsec2->GetYaxis()->SetTitleOffset(1.7);
      emptydedxsec2->Draw("axis");

      TProfile *currentprof_dedxsec2;
      for (int l=0; l<numberofenergies; l++){
	currenthistdedxsec2 = (TH2F*)results[k]->Get(("TotalSF_dedxsec2vsShowerDepth_" + IntToString(particleenergies[l])).c_str());
	currentprof_dedxsec2 = currenthistdedxsec2->ProfileX();
	if (l == 0){currentprof_dedxsec2->SetMarkerColor(4);currentprof_dedxsec2->SetMarkerStyle(20);}
	if (l == 1){currentprof_dedxsec2->SetMarkerColor(2);currentprof_dedxsec2->SetMarkerStyle(21);}
	if (l == 2){currentprof_dedxsec2->SetMarkerColor(1);currentprof_dedxsec2->SetMarkerStyle(22);}
	if (l == 3){currentprof_dedxsec2->SetMarkerColor(3);currentprof_dedxsec2->SetMarkerStyle(23);}
	if (l == 4){currentprof_dedxsec2->SetMarkerColor(5);currentprof_dedxsec2->SetMarkerStyle(28);}
	if (l == 5){currentprof_dedxsec2->SetMarkerColor(6);currentprof_dedxsec2->SetMarkerStyle(29);}
	if (l == 6){currentprof_dedxsec2->SetMarkerColor(7);currentprof_dedxsec2->SetMarkerStyle(30);}
	if (l == 7){currentprof_dedxsec2->SetMarkerColor(8);currentprof_dedxsec2->SetMarkerStyle(31);}
	if (l == 8){currentprof_dedxsec2->SetMarkerColor(9);currentprof_dedxsec2->SetMarkerStyle(33);}
	if (l == 9){currentprof_dedxsec2->SetMarkerColor(12);currentprof_dedxsec2->SetMarkerStyle(34);}
	if (l == 10){currentprof_dedxsec2->SetMarkerColor(46);currentprof_dedxsec2->SetMarkerStyle(25);}
	if (l == 11){currentprof_dedxsec2->SetMarkerColor(41);currentprof_dedxsec2->SetMarkerStyle(26);}
 	if (l == 12){currentprof_dedxsec2->SetMarkerColor(38);currentprof_dedxsec2->SetMarkerStyle(2);}
	if (l == 13){currentprof_dedxsec2->SetMarkerColor(40);currentprof_dedxsec2->SetMarkerStyle(3);}
	if (l == 14){currentprof_dedxsec2->SetMarkerColor(30);currentprof_dedxsec2->SetMarkerStyle(5);}
	if (l == 15){currentprof_dedxsec2->SetMarkerColor(24);currentprof_dedxsec2->SetMarkerStyle(21);}
	if (l == 16){currentprof_dedxsec2->SetMarkerColor(29);currentprof_dedxsec2->SetMarkerStyle(22);}
	if (l == 17){currentprof_dedxsec2->SetMarkerColor(49);currentprof_dedxsec2->SetMarkerStyle(32);}

	leg_tsfdedxsec2vst_prof->AddEntry(currentprof_dedxsec2, procNa[l], "P");
	currentprof_dedxsec2->Draw("same");
	leg_tsfdedxsec2vst_prof->Draw("same");
	currentprof_dedxsec2->Write();
      }

      c22->Update(); 

      cnpng = "Totalsfdedxsec2vsshowerdepth_profile.png";
      c22->SaveAs(cnpng);
      c22->Write();




     //======================================================================================
      //Total sampling fraction for dedxsection 3 vs X0max/<X0max in radiation length : Profile
      c23->cd();
      // pad->Draw();
      // pad->cd();
      gPad->SetLeftMargin(0.14);

      // Empty histo
      TH2F *emptydedxsec3 = new TH2F("emptydedxsec3", " ", 100,0.,maxbsh,400,0.,0.04);
      emptydedxsec3->SetTitle( "Profile of total sampling fraction for section 2 vs. "+shdp+" for " + particle); 
      emptydedxsec3->SetXTitle(shdp); 
      emptydedxsec3->SetYTitle("SF_{sec2} = #sum_{i=sec2}Eactive_{i}/#sum_{i=sec2}(Epassive_{i}+Eactive_{i})");
      emptydedxsec3->GetYaxis()->SetTitleOffset(1.7);
      emptydedxsec3->Draw("axis");

      TProfile *currentprof_dedxsec3;
      for (int l=0; l<numberofenergies; l++){
	currenthistdedxsec3 = (TH2F*)results[k]->Get(("TotalSF_dedxsec3vsShowerDepth_" + IntToString(particleenergies[l])).c_str());
	currentprof_dedxsec3 = currenthistdedxsec3->ProfileX();
	if (l == 0){currentprof_dedxsec3->SetMarkerColor(4);currentprof_dedxsec3->SetMarkerStyle(20);}
	if (l == 1){currentprof_dedxsec3->SetMarkerColor(2);currentprof_dedxsec3->SetMarkerStyle(21);}
	if (l == 2){currentprof_dedxsec3->SetMarkerColor(1);currentprof_dedxsec3->SetMarkerStyle(22);}
	if (l == 3){currentprof_dedxsec3->SetMarkerColor(3);currentprof_dedxsec3->SetMarkerStyle(23);}
	if (l == 4){currentprof_dedxsec3->SetMarkerColor(5);currentprof_dedxsec3->SetMarkerStyle(28);}
	if (l == 5){currentprof_dedxsec3->SetMarkerColor(6);currentprof_dedxsec3->SetMarkerStyle(29);}
	if (l == 6){currentprof_dedxsec3->SetMarkerColor(7);currentprof_dedxsec3->SetMarkerStyle(30);}
	if (l == 7){currentprof_dedxsec3->SetMarkerColor(8);currentprof_dedxsec3->SetMarkerStyle(31);}
	if (l == 8){currentprof_dedxsec3->SetMarkerColor(9);currentprof_dedxsec3->SetMarkerStyle(33);}
	if (l == 9){currentprof_dedxsec3->SetMarkerColor(12);currentprof_dedxsec3->SetMarkerStyle(34);}
	if (l == 10){currentprof_dedxsec3->SetMarkerColor(46);currentprof_dedxsec3->SetMarkerStyle(25);}
	if (l == 11){currentprof_dedxsec3->SetMarkerColor(41);currentprof_dedxsec3->SetMarkerStyle(26);}
 	if (l == 12){currentprof_dedxsec3->SetMarkerColor(38);currentprof_dedxsec3->SetMarkerStyle(2);}
	if (l == 13){currentprof_dedxsec3->SetMarkerColor(40);currentprof_dedxsec3->SetMarkerStyle(3);}
	if (l == 14){currentprof_dedxsec3->SetMarkerColor(30);currentprof_dedxsec3->SetMarkerStyle(5);}
	if (l == 15){currentprof_dedxsec3->SetMarkerColor(24);currentprof_dedxsec3->SetMarkerStyle(21);}
	if (l == 16){currentprof_dedxsec3->SetMarkerColor(29);currentprof_dedxsec3->SetMarkerStyle(22);}
	if (l == 17){currentprof_dedxsec3->SetMarkerColor(49);currentprof_dedxsec3->SetMarkerStyle(32);}

	leg_tsfdedxsec3vst_prof->AddEntry(currentprof_dedxsec3, procNa[l], "P");
	currentprof_dedxsec3->Draw("same");
	leg_tsfdedxsec3vst_prof->Draw("same");
	currentprof_dedxsec3->Write();
      }

      c23->Update(); 

      cnpng = "Totalsfdedxsec3vsshowerdepth_profile.png";
      c23->SaveAs(cnpng);
      c23->Write();
     //======================================================================================
     //======================================================================================
      //Total sampling fraction for dedxsection 4 vs X0max/<X0max in radiation length : Profile
      c24->cd();
      // pad->Draw();
      // pad->cd();
      gPad->SetLeftMargin(0.14);

      // Empty histo
      TH2F *emptydedxsec4 = new TH2F("emptydedxsec4", " ", 100,0.,maxbsh,400,0.,0.04);
      emptydedxsec4->SetTitle( "Profile of total sampling fraction for section 3 vs. "+shdp+" for " + particle); 
      emptydedxsec4->SetXTitle(shdp); 
      emptydedxsec4->SetYTitle("SF_{sec3} = #sum_{i=sec3}Eactive_{i}/#sum_{i=sec3}(Epassive_{i}+Eactive_{i})");
      emptydedxsec4->GetYaxis()->SetTitleOffset(1.7);
      emptydedxsec4->Draw("axis");

      TProfile *currentprof_dedxsec4;
      for (int l=0; l<numberofenergies; l++){
	currenthistdedxsec4 = (TH2F*)results[k]->Get(("TotalSF_dedxsec4vsShowerDepth_" + IntToString(particleenergies[l])).c_str());
	currentprof_dedxsec4 = currenthistdedxsec4->ProfileX();
	if (l == 0){currentprof_dedxsec4->SetMarkerColor(4);currentprof_dedxsec4->SetMarkerStyle(20);}
	if (l == 1){currentprof_dedxsec4->SetMarkerColor(2);currentprof_dedxsec4->SetMarkerStyle(21);}
	if (l == 2){currentprof_dedxsec4->SetMarkerColor(1);currentprof_dedxsec4->SetMarkerStyle(22);}
	if (l == 3){currentprof_dedxsec4->SetMarkerColor(3);currentprof_dedxsec4->SetMarkerStyle(23);}
	if (l == 4){currentprof_dedxsec4->SetMarkerColor(5);currentprof_dedxsec4->SetMarkerStyle(28);}
	if (l == 5){currentprof_dedxsec4->SetMarkerColor(6);currentprof_dedxsec4->SetMarkerStyle(29);}
	if (l == 6){currentprof_dedxsec4->SetMarkerColor(7);currentprof_dedxsec4->SetMarkerStyle(30);}
	if (l == 7){currentprof_dedxsec4->SetMarkerColor(8);currentprof_dedxsec4->SetMarkerStyle(31);}
	if (l == 8){currentprof_dedxsec4->SetMarkerColor(9);currentprof_dedxsec4->SetMarkerStyle(33);}
	if (l == 9){currentprof_dedxsec4->SetMarkerColor(12);currentprof_dedxsec4->SetMarkerStyle(34);}
	if (l == 10){currentprof_dedxsec4->SetMarkerColor(46);currentprof_dedxsec4->SetMarkerStyle(25);}
	if (l == 11){currentprof_dedxsec4->SetMarkerColor(41);currentprof_dedxsec4->SetMarkerStyle(26);}
 	if (l == 12){currentprof_dedxsec4->SetMarkerColor(38);currentprof_dedxsec4->SetMarkerStyle(2);}
	if (l == 13){currentprof_dedxsec4->SetMarkerColor(40);currentprof_dedxsec4->SetMarkerStyle(3);}
	if (l == 14){currentprof_dedxsec4->SetMarkerColor(30);currentprof_dedxsec4->SetMarkerStyle(5);}
	if (l == 15){currentprof_dedxsec4->SetMarkerColor(24);currentprof_dedxsec4->SetMarkerStyle(21);}
	if (l == 16){currentprof_dedxsec4->SetMarkerColor(29);currentprof_dedxsec4->SetMarkerStyle(22);}
	if (l == 17){currentprof_dedxsec4->SetMarkerColor(49);currentprof_dedxsec4->SetMarkerStyle(32);}

	leg_tsfdedxsec4vst_prof->AddEntry(currentprof_dedxsec4, procNa[l], "P");
	currentprof_dedxsec4->Draw("same");
	leg_tsfdedxsec4vst_prof->Draw("same");
	currentprof_dedxsec4->Write();
      }

      c24->Update(); 

      cnpng = "Totalsfdedxsec4vsshowerdepth_profile.png";
      c24->SaveAs(cnpng);
      c24->Write();
    //======================================================================================
      //Total sampling fraction for dedxsection 5 vs X0max/<X0max in radiation length : Profile
      c25->cd();
      // pad->Draw();
      // pad->cd();
      gPad->SetLeftMargin(0.14);

      // Empty histo
      TH2F *emptydedxsec5 = new TH2F("emptydedxsec5", " ", 100,0.,maxbsh,400,0.,0.04);
      emptydedxsec5->SetTitle( "Profile of total sampling fraction for section 4 vs. "+shdp+" for " + particle); 
      emptydedxsec5->SetXTitle(shdp); 
      emptydedxsec5->SetYTitle("SF_{sec4} = #sum_{i=sec4}Eactive_{i}/#sum_{i=sec4}(Epassive_{i}+Eactive_{i})");
      emptydedxsec5->GetYaxis()->SetTitleOffset(1.7);
      emptydedxsec5->Draw("axis");

      TProfile *currentprof_dedxsec5;
      for (int l=0; l<numberofenergies; l++){
	currenthistdedxsec5 = (TH2F*)results[k]->Get(("TotalSF_dedxsec5vsShowerDepth_" + IntToString(particleenergies[l])).c_str());
	currentprof_dedxsec5 = currenthistdedxsec5->ProfileX();
	if (l == 0){currentprof_dedxsec5->SetMarkerColor(4);currentprof_dedxsec5->SetMarkerStyle(20);}
	if (l == 1){currentprof_dedxsec5->SetMarkerColor(2);currentprof_dedxsec5->SetMarkerStyle(21);}
	if (l == 2){currentprof_dedxsec5->SetMarkerColor(1);currentprof_dedxsec5->SetMarkerStyle(22);}
	if (l == 3){currentprof_dedxsec5->SetMarkerColor(3);currentprof_dedxsec5->SetMarkerStyle(23);}
	if (l == 4){currentprof_dedxsec5->SetMarkerColor(5);currentprof_dedxsec5->SetMarkerStyle(28);}
	if (l == 5){currentprof_dedxsec5->SetMarkerColor(6);currentprof_dedxsec5->SetMarkerStyle(29);}
	if (l == 6){currentprof_dedxsec5->SetMarkerColor(7);currentprof_dedxsec5->SetMarkerStyle(30);}
	if (l == 7){currentprof_dedxsec5->SetMarkerColor(8);currentprof_dedxsec5->SetMarkerStyle(31);}
	if (l == 8){currentprof_dedxsec5->SetMarkerColor(9);currentprof_dedxsec5->SetMarkerStyle(33);}
	if (l == 9){currentprof_dedxsec5->SetMarkerColor(12);currentprof_dedxsec5->SetMarkerStyle(34);}
	if (l == 10){currentprof_dedxsec5->SetMarkerColor(46);currentprof_dedxsec5->SetMarkerStyle(25);}
	if (l == 11){currentprof_dedxsec5->SetMarkerColor(41);currentprof_dedxsec5->SetMarkerStyle(26);}
 	if (l == 12){currentprof_dedxsec5->SetMarkerColor(38);currentprof_dedxsec5->SetMarkerStyle(2);}
	if (l == 13){currentprof_dedxsec5->SetMarkerColor(40);currentprof_dedxsec5->SetMarkerStyle(3);}
	if (l == 14){currentprof_dedxsec5->SetMarkerColor(30);currentprof_dedxsec5->SetMarkerStyle(5);}
	if (l == 15){currentprof_dedxsec5->SetMarkerColor(24);currentprof_dedxsec5->SetMarkerStyle(21);}
	if (l == 16){currentprof_dedxsec5->SetMarkerColor(29);currentprof_dedxsec5->SetMarkerStyle(22);}
	if (l == 17){currentprof_dedxsec5->SetMarkerColor(49);currentprof_dedxsec5->SetMarkerStyle(32);}

	leg_tsfdedxsec5vst_prof->AddEntry(currentprof_dedxsec5, procNa[l], "P");
	currentprof_dedxsec5->Draw("same");
	leg_tsfdedxsec5vst_prof->Draw("same");
	currentprof_dedxsec5->Write();
      }

      c25->Update(); 

      cnpng = "Totalsfdedxsec5vsshowerdepth_profile.png";
      c25->SaveAs(cnpng);
      c25->Write();
      //======================================================================================
     //======================================================================================
      //Total sampling fraction for dedxsection 6 vs X0max/<X0max in radiation length : Profile
      c26->cd();
      // pad->Draw();
      // pad->cd();
      gPad->SetLeftMargin(0.14);

      // Empty histo
      TH2F *emptydedxsec6 = new TH2F("emptydedxsec6", " ", 100,0.,maxbsh,400,0.,0.04);
      emptydedxsec6->SetTitle( "Profile of total sampling fraction for section 5 vs. "+shdp+" for " + particle); 
      emptydedxsec6->SetXTitle(shdp); 
      emptydedxsec6->SetYTitle("SF_{sec5} = #sum_{i=sec5}Eactive_{i}/#sum_{i=sec5}(Epassive_{i}+Eactive_{i})");
      emptydedxsec6->GetYaxis()->SetTitleOffset(1.7);
      emptydedxsec6->Draw("axis");

      TProfile *currentprof_dedxsec6;
      for (int l=0; l<numberofenergies; l++){
	currenthistdedxsec6 = (TH2F*)results[k]->Get(("TotalSF_dedxsec6vsShowerDepth_" + IntToString(particleenergies[l])).c_str());
	currentprof_dedxsec6 = currenthistdedxsec6->ProfileX();
	if (l == 0){currentprof_dedxsec6->SetMarkerColor(4);currentprof_dedxsec6->SetMarkerStyle(20);}
	if (l == 1){currentprof_dedxsec6->SetMarkerColor(2);currentprof_dedxsec6->SetMarkerStyle(21);}
	if (l == 2){currentprof_dedxsec6->SetMarkerColor(1);currentprof_dedxsec6->SetMarkerStyle(22);}
	if (l == 3){currentprof_dedxsec6->SetMarkerColor(3);currentprof_dedxsec6->SetMarkerStyle(23);}
	if (l == 4){currentprof_dedxsec6->SetMarkerColor(5);currentprof_dedxsec6->SetMarkerStyle(28);}
	if (l == 5){currentprof_dedxsec6->SetMarkerColor(6);currentprof_dedxsec6->SetMarkerStyle(29);}
	if (l == 6){currentprof_dedxsec6->SetMarkerColor(7);currentprof_dedxsec6->SetMarkerStyle(30);}
	if (l == 7){currentprof_dedxsec6->SetMarkerColor(8);currentprof_dedxsec6->SetMarkerStyle(31);}
	if (l == 8){currentprof_dedxsec6->SetMarkerColor(9);currentprof_dedxsec6->SetMarkerStyle(33);}
	if (l == 9){currentprof_dedxsec6->SetMarkerColor(12);currentprof_dedxsec6->SetMarkerStyle(34);}
	if (l == 10){currentprof_dedxsec6->SetMarkerColor(46);currentprof_dedxsec6->SetMarkerStyle(25);}
	if (l == 11){currentprof_dedxsec6->SetMarkerColor(41);currentprof_dedxsec6->SetMarkerStyle(26);}
 	if (l == 12){currentprof_dedxsec6->SetMarkerColor(38);currentprof_dedxsec6->SetMarkerStyle(2);}
	if (l == 13){currentprof_dedxsec6->SetMarkerColor(40);currentprof_dedxsec6->SetMarkerStyle(3);}
	if (l == 14){currentprof_dedxsec6->SetMarkerColor(30);currentprof_dedxsec6->SetMarkerStyle(5);}
	if (l == 15){currentprof_dedxsec6->SetMarkerColor(24);currentprof_dedxsec6->SetMarkerStyle(21);}
	if (l == 16){currentprof_dedxsec6->SetMarkerColor(29);currentprof_dedxsec6->SetMarkerStyle(22);}
	if (l == 17){currentprof_dedxsec6->SetMarkerColor(49);currentprof_dedxsec6->SetMarkerStyle(32);}

	leg_tsfdedxsec6vst_prof->AddEntry(currentprof_dedxsec6, procNa[l], "P");
	currentprof_dedxsec6->Draw("same");
	leg_tsfdedxsec6vst_prof->Draw("same");
	currentprof_dedxsec6->Write();
      }

      c26->Update(); 

      cnpng = "Totalsfdedxsec6vsshowerdepth_profile.png";
      c26->SaveAs(cnpng);
      c26->Write();
      //======================================================================================
 


























































     
     //======================================================================================
      //Sampling fraction vs normalized shower depth: Profile
      c16->cd();
      // pad->Draw();
      // pad->cd();
      // gPad->SetLeftMargin(0.14);
      // THStack *hs= new THStack("hs","profile");
      empty->Reset();
      delete empty;
      empty = new TH2F("empty", " ", 100,0.,maxbsh,400,0.,0.04);
      empty->SetTitle( "Sampling fraction vs. "+shdp+" for " + particle); 
      empty->SetXTitle(shdp); 
      empty->SetYTitle("SF_{i} = Eactive_{i}/(Epassive_{i}+Eactive_{i})");
      empty->GetYaxis()->SetTitleOffset(1.3);
      empty->Draw("axis");

      TProfile *currentprof1;
      for (int l=0; l<numberofenergies; l++){
	//Totalsfvsshowerdepth_forprofile[l]->ProfileX() = (TGraphErrors*)results[k]->Get(auxnam_tgr.c_str());
	currenthist2 = (TH2F*)results[k]->Get(("h_sfvsnormshowerdepth_" + IntToString(particleenergies[l])).c_str());
	currentprof1 = currenthist2->ProfileX();
	// std::cout << currenthist2->GetMean(1) << std::endl;
	if (l == 0){currentprof1->SetMarkerColor(4);currentprof1->SetMarkerStyle(20);}
	if (l == 1){currentprof1->SetMarkerColor(2);currentprof1->SetMarkerStyle(21);}
	if (l == 2){currentprof1->SetMarkerColor(1);currentprof1->SetMarkerStyle(22);}
	if (l == 3){currentprof1->SetMarkerColor(3);currentprof1->SetMarkerStyle(23);}
	if (l == 4){currentprof1->SetMarkerColor(5);currentprof1->SetMarkerStyle(28);}
	if (l == 5){currentprof1->SetMarkerColor(6);currentprof1->SetMarkerStyle(29);}
	if (l == 6){currentprof1->SetMarkerColor(7);currentprof1->SetMarkerStyle(30);}
	if (l == 7){currentprof1->SetMarkerColor(8);currentprof1->SetMarkerStyle(31);}
	if (l == 8){currentprof1->SetMarkerColor(9);currentprof1->SetMarkerStyle(33);}
	if (l == 9){currentprof1->SetMarkerColor(12);currentprof1->SetMarkerStyle(34);}
	if (l == 10){currentprof1->SetMarkerColor(46);currentprof1->SetMarkerStyle(25);}
	if (l == 11){currentprof1->SetMarkerColor(41);currentprof1->SetMarkerStyle(26);}
 	if (l == 12){currentprof1->SetMarkerColor(38);currentprof1->SetMarkerStyle(2);}
	if (l == 13){currentprof1->SetMarkerColor(40);currentprof1->SetMarkerStyle(3);}
	if (l == 14){currentprof1->SetMarkerColor(30);currentprof1->SetMarkerStyle(5);}
	if (l == 15){currentprof1->SetMarkerColor(24);currentprof1->SetMarkerStyle(21);}
	if (l == 16){currentprof1->SetMarkerColor(29);currentprof1->SetMarkerStyle(22);}
	if (l == 17){currentprof1->SetMarkerColor(49);currentprof1->SetMarkerStyle(32);}

	// currentprof1->SetTitle( "Sampling fraction vs. "+shdp+" for " + particle); 
	// currentprof1->SetXTitle(shdp); 
	// currentprof1->SetYTitle("SF_{i} = Eactive_{i}/(Epassive_{i}+Eactive_{i})");
	// currentprof1->GetYaxis()->SetTitleOffset(1.7);
	// currentprof1->GetYaxis()->SetRangeUser(0.006,0.0016); 
	// currentprof1->GetYaxis()->SetLimits(0.006,0.0016); 
	leg_tsfvsnormdepth_prof->AddEntry(currentprof1, procNa[l], "P");
	currentprof1->Draw("same");
	leg_tsfvsnormdepth_prof->Draw("same");

	// l == 0 ? currentprof1->Draw() : currentprof1->Draw("same");
	// l == 0 ? leg_tsfvsnormdepth_prof->Draw() : leg_tsfvsnormdepth_prof->Draw("same");
	// Totalsfvsshowerdepth_forprofile[l]->ProfileX()->GetXaxis()->SetRangeUser(0.,.); 
	// Totalsfvsshowerdepth_forprofile[0]->ProfileX()->GetYaxis()->SetRangeUser(0.006,0.0016); 
	// hs->Add(currentprof1);
	// if (l==0){
	//   hs.SetHistogram(currentprof1);
	// }
	// leg_tsfvst_prof->AddEntry(currentprof1, procNa[l], "P");

	currentprof1->Write();
	
      }

      c16->Update(); 

      cnpng = "sfvsnormshowerdepth_profile.png";
      c16->SaveAs(cnpng);
      c16->Write();

      TH2F *emp[numberofenergies];
      //Save the above profiles not all together
      for (int l=0; l<numberofenergies; l++){
	double subcan1 = 0.;double subcan2 = 0.;
	switch (particleenergies[l]){ 
	case 2: subcan1 = 1;   break;
	case 3: subcan1 = 2;   break;
	case 5: subcan1 = 3;   break;
	case 8: subcan1 = 4;   break;
	case 10: subcan1 =5 ;   break;
	case 15: subcan1 =6 ;   break;
	case 30: subcan1 =7 ;   break;
	case 50: subcan1 =8 ;   break;
	case 80: subcan1 =9 ;   break;
	case 100: subcan2 = 1 ;   break;
	case 120: subcan2 = 2 ;   break;
	case 150: subcan2 = 3 ;   break;
	case 180: subcan2 = 4 ;   break;
	case 200: subcan2 = 5 ;   break;
	case 250: subcan2 = 6 ;   break;
	case 300: subcan2 = 7 ;   break;
	case 400: subcan2 = 8 ;   break;
	case 500: subcan2 = 9 ;   break;
	}		  
	if ( particleenergies[l] < 90){c17_1->cd(subcan1);}
	if ( particleenergies[l] > 90){c17_2->cd(subcan2);}
	emp[l] = new TH2F( ("Sampling fraction vs. normalized shower depth profile for " + IntToString(particleenergies[l]) + " GeV " +particle), " ", 100,0.,maxbsh,400,0.,0.04);
	emp[l]->SetTitle("Sampling fraction vs. normalized shower depth profile for " + IntToString(particleenergies[l]) + " GeV " +particle); 
	emp[l]->SetXTitle(shdp); 
	emp[l]->SetYTitle("SF_{i} = Eactive_{i}/(Epassive_{i}+Eactive_{i})");
	emp[l]->GetYaxis()->SetTitleOffset(1.3);
	emp[l]->Draw("axis");
	currenthist3 = (TH2F*)results[k]->Get(("h_sfvsnormshowerdepth_" + IntToString(particleenergies[l])).c_str());
	currentprof1 = currenthist3->ProfileX();
	// currentprof1->SetTitle( "Sampling fraction vs. normalized shower depth profile for " + IntToString(particleenergies[l]) + " GeV " +particle); 
	if (l == 0){currentprof1->SetMarkerColor(4);currentprof1->SetMarkerStyle(20);}
	if (l == 1){currentprof1->SetMarkerColor(2);currentprof1->SetMarkerStyle(21);}
	if (l == 2){currentprof1->SetMarkerColor(1);currentprof1->SetMarkerStyle(22);}
	if (l == 3){currentprof1->SetMarkerColor(3);currentprof1->SetMarkerStyle(23);}
	if (l == 4){currentprof1->SetMarkerColor(5);currentprof1->SetMarkerStyle(28);}
	if (l == 5){currentprof1->SetMarkerColor(6);currentprof1->SetMarkerStyle(29);}
	if (l == 6){currentprof1->SetMarkerColor(7);currentprof1->SetMarkerStyle(30);}
	if (l == 7){currentprof1->SetMarkerColor(8);currentprof1->SetMarkerStyle(31);}
	if (l == 8){currentprof1->SetMarkerColor(9);currentprof1->SetMarkerStyle(33);}
	if (l == 9){currentprof1->SetMarkerColor(12);currentprof1->SetMarkerStyle(34);}
	if (l == 10){currentprof1->SetMarkerColor(46);currentprof1->SetMarkerStyle(25);}
	if (l == 11){currentprof1->SetMarkerColor(41);currentprof1->SetMarkerStyle(26);}
 	if (l == 12){currentprof1->SetMarkerColor(38);currentprof1->SetMarkerStyle(2);}
	if (l == 13){currentprof1->SetMarkerColor(40);currentprof1->SetMarkerStyle(3);}
	if (l == 14){currentprof1->SetMarkerColor(30);currentprof1->SetMarkerStyle(5);}
	if (l == 15){currentprof1->SetMarkerColor(24);currentprof1->SetMarkerStyle(21);}
	if (l == 16){currentprof1->SetMarkerColor(29);currentprof1->SetMarkerStyle(22);}
	if (l == 17){currentprof1->SetMarkerColor(49);currentprof1->SetMarkerStyle(32);}
	currentprof1->Draw("same");
	// currentprof1->GetYaxis()->SetTitleOffset(1.3);
	
	if ( particleenergies[l] < 90){c17_1->Update();}
	if ( particleenergies[l] > 90){c17_2->Update();}

	c17_1->SaveAs("sfvsnormshowerdepth_profile_1.png");
	c17_2->SaveAs("sfvsnormshowerdepth_profile_2.png");
	c17_1->SaveAs("sfvsnormshowerdepth_profile_1.root");
	c17_2->SaveAs("sfvsnormshowerdepth_profile_2.root");
	
	c17_1->Write();
	c17_2->Write();
	

      }

     







      if(dothe3dplot){

	TH2F *em[(*ssvec).size()];
	//======================================================================================
	//Sampling fraction per layer vs shower depth in radiation length : Profile 
	TProfile *currprofperlay;
	for(unsigned iL(startlayer); iL<(*ssvec).size(); iL++){

	  c13[iL]->cd();
	  // gPad->SetLeftMargin(0.14);
	  TString auxnam1 = "Profile of sampling fraction in layer " + IntToString((int) iL);
	  TString auxnam2 = " vs. shower depth for " + particle;
	  TString auxnam = auxnam1 + auxnam2;
	  em[iL] = new TH2F( auxnam , " ", 100,0.,maxbsh,400,0.,0.04);
	  em[iL]->SetTitle( auxnam ); 
	  em[iL]->SetXTitle(shdp); 
	  em[iL]->SetYTitle("SF_{i} = Eactive_{i}/(Epassive_{i}+Eactive_{i})");
	  em[iL]->GetYaxis()->SetTitleOffset(1.3);
	  em[iL]->Draw("axis");
	  
	  for (int l=0; l<numberofenergies; l++){
    

	    
	    TString auxnam1_sf_prof = "h_sfvsshowerdepthprofile_" + IntToString((int) iL);
	    TString auxnam2_sf_prof = "_" + IntToString(particleenergies[l]);
	    TString auxnam_sf_prof = auxnam1_sf_prof + auxnam2_sf_prof;

	  
	    //h_sfvsshowerdepthprofile[iL][l]->ProfileX() = (TGraphErrors*)results[k]->Get(auxnam_tgr.c_str());
	    // currprofperlay = h_sfvsshowerdepthprofile[iL][l]->ProfileX();
	    currenthist4 = (TH2F*)results[k]->Get( auxnam_sf_prof );
	    currprofperlay = currenthist4->ProfileX();
	    if (l == 0){currprofperlay->SetMarkerColor(4);currprofperlay->SetMarkerStyle(20);}
	    if (l == 1){currprofperlay->SetMarkerColor(2);currprofperlay->SetMarkerStyle(21);}
	    if (l == 2){currprofperlay->SetMarkerColor(1);currprofperlay->SetMarkerStyle(22);}
	    if (l == 3){currprofperlay->SetMarkerColor(3);currprofperlay->SetMarkerStyle(23);}
	    if (l == 4){currprofperlay->SetMarkerColor(5);currprofperlay->SetMarkerStyle(28);}
	    if (l == 5){currprofperlay->SetMarkerColor(6);currprofperlay->SetMarkerStyle(29);}
	    if (l == 6){currprofperlay->SetMarkerColor(7);currprofperlay->SetMarkerStyle(30);}
	    if (l == 7){currprofperlay->SetMarkerColor(8);currprofperlay->SetMarkerStyle(31);}
	    if (l == 8){currprofperlay->SetMarkerColor(9);currprofperlay->SetMarkerStyle(33);}
	    if (l == 9){currprofperlay->SetMarkerColor(12);currprofperlay->SetMarkerStyle(34);}
	    if (l == 10){currprofperlay->SetMarkerColor(46);currprofperlay->SetMarkerStyle(25);}
	    if (l == 11){currprofperlay->SetMarkerColor(41);currprofperlay->SetMarkerStyle(26);}
	    if (l == 12){currprofperlay->SetMarkerColor(38);currprofperlay->SetMarkerStyle(2);}
	    if (l == 13){currprofperlay->SetMarkerColor(40);currprofperlay->SetMarkerStyle(3);}
	    if (l == 14){currprofperlay->SetMarkerColor(30);currprofperlay->SetMarkerStyle(5);}
	    if (l == 15){currprofperlay->SetMarkerColor(24);currprofperlay->SetMarkerStyle(21);}
	    if (l == 16){currprofperlay->SetMarkerColor(29);currprofperlay->SetMarkerStyle(22);}
	    if (l == 17){currprofperlay->SetMarkerColor(49);currprofperlay->SetMarkerStyle(32);}
 
	    currprofperlay->SetTitle( auxnam ); 
	    currprofperlay->SetXTitle(shdp); 
	    currprofperlay->SetYTitle("SF_{i} = Eactive_{i}/(Epassive_{i}+Eactive_{i})");
	    currprofperlay->GetYaxis()->SetTitleOffset(1.7);
	    // currprofperlay->GetYaxis()->SetRangeUser(0.006,0.0016); 
	    // currprofperlay->GetYaxis()->SetLimits(0.006,0.0016); 
	    leg_sfplvst_prof[iL]->AddEntry(currprofperlay, procNa[l], "P");
	    currprofperlay->Draw("same");
	    leg_sfplvst_prof[iL]->Draw("same");

	    // l == 0 ? currprofperlay->Draw() : currprofperlay->Draw("same");
	    // l == 0 ? leg_sfplvst_prof[iL]->Draw() : leg_sfplvst_prof[iL]->Draw("same");
	    // h_sfvsshowerdepthprofile[iL][l]->ProfileX()->GetXaxis()->SetRangeUser(0.,.); 
	    // h_sfvsshowerdepthprofile[iL][0]->ProfileX()->GetYaxis()->SetRangeUser(0.006,0.0016); 
	    // hs->Add(currprofperlay);
	    // if (l==0){
	    //   hs.SetHistogram(currprofperlay);
	    // }
	    // leg_sfplvst_prof->AddEntry(currprofperlay, procNa[l], "P");

	    currprofperlay->Write();
	    auxnam1 = "sflayer_" + IntToString((int) iL);
	    auxnam2 = "_" + particle;
	    TString auxnam3 = ".png";
	    cnpng = auxnam1 + auxnam2 + auxnam3;
	    c13[iL]->SaveAs(cnpng);
	    c13[iL]->Write();
	  }

	}
	c14_1->Divide(2,3);
	c14_2->Divide(2,3);
	c14_3->Divide(2,3);
	c14_4->Divide(2,3);
	c14_5->Divide(2,3);

	for(unsigned iL(startlayer); iL<(*ssvec).size(); iL++){
	  int thesubcan = iL;
	  if(iL%(5+startlayer)==0){thesubcan = 6;} 

	  if(iL<=(5+startlayer)){c14_1->cd(thesubcan);}
	  if( iL>=(6+startlayer) && iL<=(11+startlayer)){c14_2->cd(thesubcan%6);}
	  if( iL>=(12+startlayer) && iL<=(17+startlayer)){c14_3->cd(thesubcan%6);}
	  if( iL>=(18+startlayer) && iL<=(23+startlayer)){c14_4->cd(thesubcan%6);}
	  if( iL>=(24+startlayer) && iL<=(29+startlayer)){c14_5->cd(thesubcan%6);}
	  // pad->Draw();
	  // pad->cd();
	  c13[iL]->Draw();

	  // c14_1->SaveAs(cnpng);
	  // c14_1->Write();
	  // c14_2->SaveAs(cnpng);
	  // c14_2->Write();
	  // c14_3->SaveAs(cnpng);
	  // c14_3->Write();
	  // c14_4->SaveAs(cnpng);
	  // c14_4->Write();
	  // c14_5->SaveAs(cnpng);
	  // c14_5->Write();

	  if(iL==(5+startlayer)){
	    // c14_1->SaveAs(cnpng);
	    c14_1->Write();
	  }
	  if(iL==(11+startlayer)){
	    // c14_2->SaveAs(cnpng);
	    c14_2->Write();
	  }
	  if(iL==(17+startlayer)){
	    // c14_3->SaveAs(cnpng);
	    c14_3->Write();
	  }
	  if(iL==(23+startlayer)){
	    // c14_4->SaveAs(cnpng);
	    c14_4->Write();
	  }
	  if(iL==(29+startlayer)){
	    // c14_5->SaveAs(cnpng);
	    c14_5->Write();
	  }

	}


	
      }//closing the if for those time consuming plots


     if(dorunfit){


	//=========================================================================================
	//Here we will make the plot for the second fit results slope
	c11->cd();

	for(unsigned iL(startlayer); iL<(*ssvec).size(); iL++){
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
	  ((int) iL) == startlayer ? gr_MIPSvsdEplusMIPs_fit[(int) iL]->Draw("APL") : gr_MIPSvsdEplusMIPs_fit[(int) iL]->Draw("PSL");
	}
	leg_fit2->Draw("PS");
	c11->Update(); 

	gr_MIPSvsdEplusMIPs_fit[startlayer]->SetTitle( "Measured energy vs absorbed plus measured energy fit results for " + particle + " for all beam energies" );
	gr_MIPSvsdEplusMIPs_fit[startlayer]->GetHistogram()->SetXTitle("Beam Energy (GeV)"); 
	gr_MIPSvsdEplusMIPs_fit[startlayer]->GetHistogram()->SetYTitle("E_{active} = f(E_{active} + E_{passive}): Slope ");
	gr_MIPSvsdEplusMIPs_fit[startlayer]->GetXaxis()->SetRangeUser(0.,520.);  
	gr_MIPSvsdEplusMIPs_fit[startlayer]->GetYaxis()->SetRangeUser(-1.,1.);  

	c11->Update(); 
	c11->SaveAs("MIPSvsdEplusMIPs_fitresults_slope" + particle + ".png"); 
	c11->Write();


	//=========================================================================================
	//Here we will make the plot for the second fit results constant term
	c12->cd();

	for(unsigned iL(startlayer); iL<(*ssvec).size(); iL++){
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
	  ((int) iL) == startlayer ? gr_MIPSvsdEplusMIPs_fit_c[(int) iL]->Draw("APL") : gr_MIPSvsdEplusMIPs_fit_c[(int) iL]->Draw("PSL");
	}
	leg_fit2_c->Draw("PS");
	c12->Update(); 

	gr_MIPSvsdEplusMIPs_fit_c[startlayer]->SetTitle( "Measured energy vs absorbed plus measured energy fit results for " + particle + " for all beam energies" );
	gr_MIPSvsdEplusMIPs_fit_c[startlayer]->GetHistogram()->SetXTitle("Beam Energy (GeV)"); 
	gr_MIPSvsdEplusMIPs_fit_c[startlayer]->GetHistogram()->SetYTitle("E_{active} = f(E_{active} + E_{passive}): Constant term (GeV)");
	gr_MIPSvsdEplusMIPs_fit_c[startlayer]->GetXaxis()->SetRangeUser(0.,520.);  
	gr_MIPSvsdEplusMIPs_fit_c[startlayer]->GetYaxis()->SetRangeUser(-10.,10.);  

	c12->Update(); 
	c12->SaveAs("MIPSvsdEplusMIPs_fitresults_constantterm" + particle + ".png"); 
	c12->Write();

 

     }


     if (dothedEvsXosplot){
       //Here for the plot of average energy lost vs radiation length
       c15->cd();
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
	 if (l == 12){dEvsXos[l]->SetMarkerColor(38);dEvsXos[l]->SetMarkerStyle(2);}
	 if (l == 13){dEvsXos[l]->SetMarkerColor(40);dEvsXos[l]->SetMarkerStyle(3);}
	 if (l == 14){dEvsXos[l]->SetMarkerColor(30);dEvsXos[l]->SetMarkerStyle(5);}
	 if (l == 15){dEvsXos[l]->SetMarkerColor(24);dEvsXos[l]->SetMarkerStyle(21);}
	 if (l == 16){dEvsXos[l]->SetMarkerColor(29);dEvsXos[l]->SetMarkerStyle(22);}
	 if (l == 17){dEvsXos[l]->SetMarkerColor(49);dEvsXos[l]->SetMarkerStyle(32);}
  
	 dEvsXos[l]->SetTitle( "Energy lost vs. transversed thickness for " + particle); 
	 dEvsXos[l]->GetHistogram()->SetXTitle("Transversed thickness (1/X_{0})"); 
	 dEvsXos[l]->GetHistogram()->SetYTitle("Average dE (MIPs)");
	 dEvsXos[l]->GetXaxis()->SetRangeUser(0.,25.8);  
	 leg_devsxos->AddEntry(dEvsXos[l], procNa[l], "P");
	 l == 0 ? dEvsXos[l]->Draw("APL") : dEvsXos[l]->Draw("PSL");
	 l == 0 ? leg_devsxos->Draw() : leg_devsxos->Draw("same");
	 c15->Update(); 
	 dEvsXos[l]->Write();
       }

       dEvsXos[0]->GetXaxis()->SetRangeUser(0.,29.8); 
       dEvsXos[0]->GetYaxis()->SetRangeUser(0.,500000.); 
       dEvsXos[0]->GetYaxis()->SetTitleOffset(1.3); 
       c15->Update(); 
       
       TString cnpng = "dEvsXos.png";
       c15->SaveAs(cnpng);
       c15->Write();
     }





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
// Double_t myfunc(Double_t *x, Double_t *) { 
//   return correlateradandlayer->Eval(x[0]);
// }
