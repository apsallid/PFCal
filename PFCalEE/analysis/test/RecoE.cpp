#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include <boost/algorithm/string.hpp>
#include <boost/regex.hpp>

#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TLine.h"

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
std::string DoubleToString (double number1);
struct FitResult{
  double mean;
  double rms;
  double meanerr;
  double rmserr;
  double Fitmean;
  double sigma;
  double Fitmeanerr;
  double sigmaerr;
};
struct FitResultResolution{
  double sigmaStoch;
  double sigmaStochErr;
  double sigmaConst;
  double sigmaConstErr;
  double sigmaNoise;
  double sigmaNoiseErr;
};

FitResult extractMeanEnergy(TCanvas *mycE_1, TCanvas *mycE_2, int color, TH1F *histo, int rebin);
void buildplot(TCanvas *c1, TH1F* hist, TLegend * leg, TString titleofplot, TString xtitle, TString ytitle, TString procNa, int numene);
void fitforlinearity(TCanvas *c1, TGraphErrors * gr, double fitenelow, double fitenemax, std::string yaxisunit);
FitResultResolution fitforresolution(TCanvas *c1, TGraphErrors * res, int color, bool addNoiseTerm);
void buildgraph(TCanvas *c1, TGraphErrors * res, TString titleofplot, TString xtitle, TString ytitle, int point, double energy, FitResult fitres, bool dogaussfit);
void buildresolutiongraph(TCanvas *c1, TGraphErrors * res, TString titleofplot, int point, double energy, FitResult fitres, bool dogaussfit);
void combineresolution(TCanvas *c1, TGraphErrors * res1, int color1, FitResultResolution fitresol1, TGraphErrors * res2, int color2, FitResultResolution fitresol2,  TGraphErrors * res3, int color3, FitResultResolution fitresol3, TGraphErrors * res4, int color4, FitResultResolution fitresol4, TString particle, bool addNoiseTerm);
void combineresolution_noupstream(TCanvas *c1, TGraphErrors * res1, std::string res1str, int color1, FitResultResolution fitresol1, TGraphErrors * res2,  std::string res2str, int color2, FitResultResolution fitresol2,  TGraphErrors * res3, std::string res3str, int color3, FitResultResolution fitresol3, TString particle, bool addNoiseTerm);
void combineresolution_noupstream_sections(TCanvas *c1, TGraphErrors * res1, std::string res1str, int color1, FitResultResolution fitresol1, TGraphErrors * res2,  std::string res2str, int color2, FitResultResolution fitresol2,  TGraphErrors * res3, std::string res3str, int color3, FitResultResolution fitresol3, TGraphErrors * res4, std::string res4str, int color4, FitResultResolution fitresol4, TString particle, bool addNoiseTerm);

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
  //To avoid rerunning the first step : false
  bool step1 = false;
  //======================================================================================
  //Old shower depth: SF = f(X0max/<X0max>)
  //New shower depth: SF = f( Sum(Ei*X0_i)/Sum(E_i)  / <Sum(Ei*X0_i)/Sum(E_i)>)
  bool usenewshowerdepth = true;
  //For shower depth to bin
  double shdepthtobin = usenewshowerdepth ? 0.03 : 0.1;
  // TString shdp = usenewshowerdepth ? " #frac{#sum E X_{0}/#sum E}{<#sum E X_{0}/#sum E>}" : "X_{0,max}/<X_{0,max}>";
  //======================================================================================
  //3 layers times 4 Si Pads
  // const int numberofpads = 12; 
  //======================================================================================
  //Particle type and energy
  TString particle = "e-"; //mu+
  const int numberofenergies = 18; // 15 30 50 80 100 150
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
  //For the rebin
  std::vector<int> rebin_raw, rebin_raw_leak, rebin, rebin_val;
  rebin.clear();rebin_val.clear();rebin_raw.clear();rebin_raw_leak.clear();
  //This is for h_recoEovertrueE
  // rebin.push_back(16);
  // rebin.push_back(16);
  // rebin.push_back(16);
  // rebin.push_back(8);
  // rebin.push_back(16);
  // rebin.push_back(8);
  // rebin.push_back(8);
  // rebin.push_back(4);
  // rebin.push_back(4);
  // rebin.push_back(4);
  // rebin.push_back(2);
  // rebin.push_back(4);
  // rebin.push_back(2);
  // rebin.push_back(2);
  // rebin.push_back(2);
  // rebin.push_back(2);
  // rebin.push_back(2);
  // rebin.push_back(2);
  //This is for h_recoE_raw
  rebin_raw.push_back(1);  rebin_raw_leak.push_back(1);    rebin_val.push_back(1);  rebin.push_back(1); 
  rebin_raw.push_back(1);  rebin_raw_leak.push_back(1);    rebin_val.push_back(1);  rebin.push_back(1); 
  rebin_raw.push_back(1);  rebin_raw_leak.push_back(1);    rebin_val.push_back(1);  rebin.push_back(1); 
  rebin_raw.push_back(2);  rebin_raw_leak.push_back(2);    rebin_val.push_back(1);  rebin.push_back(1); 
  rebin_raw.push_back(2);  rebin_raw_leak.push_back(2);    rebin_val.push_back(1);  rebin.push_back(1); 
  rebin_raw.push_back(4);  rebin_raw_leak.push_back(4);    rebin_val.push_back(1);  rebin.push_back(1); 
  rebin_raw.push_back(4);  rebin_raw_leak.push_back(4);    rebin_val.push_back(1);  rebin.push_back(1); 
  rebin_raw.push_back(8);  rebin_raw_leak.push_back(8);    rebin_val.push_back(1);  rebin.push_back(1); 
  rebin_raw.push_back(8);  rebin_raw_leak.push_back(8);    rebin_val.push_back(1);  rebin.push_back(1); 
  rebin_raw.push_back(16); rebin_raw_leak.push_back(16);   rebin_val.push_back(1);  rebin.push_back(1); 
  rebin_raw.push_back(16); rebin_raw_leak.push_back(16);   rebin_val.push_back(1);  rebin.push_back(1); 
  rebin_raw.push_back(16); rebin_raw_leak.push_back(16);   rebin_val.push_back(1);  rebin.push_back(1); 
  rebin_raw.push_back(25); rebin_raw_leak.push_back(25);   rebin_val.push_back(1);  rebin.push_back(1); 
  rebin_raw.push_back(25); rebin_raw_leak.push_back(25);   rebin_val.push_back(1);  rebin.push_back(1); 
  rebin_raw.push_back(25); rebin_raw_leak.push_back(25);   rebin_val.push_back(1);  rebin.push_back(1); 
  rebin_raw.push_back(25); rebin_raw_leak.push_back(25);   rebin_val.push_back(1);  rebin.push_back(1); 
  rebin_raw.push_back(25); rebin_raw_leak.push_back(25);   rebin_val.push_back(1);  rebin.push_back(1); 
  rebin_raw.push_back(50); rebin_raw_leak.push_back(25);   rebin_val.push_back(1);  rebin.push_back(1);  
  //======================================================================================
  //Factor to recalibrate only for dedx
  std::vector<double> recalibdedx;
  recalibdedx.clear();
  recalibdedx.push_back(0.91514);
  recalibdedx.push_back(0.92068);
  recalibdedx.push_back(0.92425);
  recalibdedx.push_back(0.92741);
  recalibdedx.push_back(0.92937);
  recalibdedx.push_back(0.93171);
  recalibdedx.push_back(0.93593);
  recalibdedx.push_back(0.939594);
  recalibdedx.push_back(0.942661);
  recalibdedx.push_back(0.94405);
  recalibdedx.push_back(0.94488);
  recalibdedx.push_back(0.94630);
  recalibdedx.push_back(0.94714);
  recalibdedx.push_back(0.94781);
  recalibdedx.push_back(0.94891);
  recalibdedx.push_back(0.95001);
  recalibdedx.push_back(0.95161);
  recalibdedx.push_back(0.95275);

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

  //For the total sampling fraction: A single number per energy
  //Also per section
  std::vector<double> totalsfperenergy, totalsfperenergy_sec1, totalsfperenergy_sec2, totalsfperenergy_sec3; 
  totalsfperenergy.clear(); totalsfperenergy_sec1.clear(); totalsfperenergy_sec2.clear(); totalsfperenergy_sec3.clear();
  totalsfperenergy.push_back(0.0113494);   totalsfperenergy_sec1.push_back(0.0135148);  totalsfperenergy_sec2.push_back(0.00996915 ); totalsfperenergy_sec3.push_back(0.00999521 );   
  totalsfperenergy.push_back(0.0111971);   totalsfperenergy_sec1.push_back(0.0139669);  totalsfperenergy_sec2.push_back(0.00987755 ); totalsfperenergy_sec3.push_back(0.0097018  );   
  totalsfperenergy.push_back(0.0109993);   totalsfperenergy_sec1.push_back(0.0144021);  totalsfperenergy_sec2.push_back(0.0097869  ); totalsfperenergy_sec3.push_back(0.00934909 );   
  totalsfperenergy.push_back(0.0108104);   totalsfperenergy_sec1.push_back(0.0147531);  totalsfperenergy_sec2.push_back(0.00971951 ); totalsfperenergy_sec3.push_back(0.00882311 );   
  totalsfperenergy.push_back(0.0107085);   totalsfperenergy_sec1.push_back(0.0149144);  totalsfperenergy_sec2.push_back(0.00967641 ); totalsfperenergy_sec3.push_back(0.00850979 );   
  totalsfperenergy.push_back(0.0105488);   totalsfperenergy_sec1.push_back(0.0150957);  totalsfperenergy_sec2.push_back(0.00964227 ); totalsfperenergy_sec3.push_back(0.00787493 );   
  totalsfperenergy.push_back(0.0102761);   totalsfperenergy_sec1.push_back(0.015342 );  totalsfperenergy_sec2.push_back(0.00958456 ); totalsfperenergy_sec3.push_back(0.00706344 );   
  totalsfperenergy.push_back(0.0100706);   totalsfperenergy_sec1.push_back(0.0155201);  totalsfperenergy_sec2.push_back(0.00954203 ); totalsfperenergy_sec3.push_back(0.00662297 );   
  totalsfperenergy.push_back(0.0098931);   totalsfperenergy_sec1.push_back(0.0155101);  totalsfperenergy_sec2.push_back(0.00950792 ); totalsfperenergy_sec3.push_back(0.00639739 );   
  totalsfperenergy.push_back(0.009804);	   totalsfperenergy_sec1.push_back(0.0155262);  totalsfperenergy_sec2.push_back(0.00949383 ); totalsfperenergy_sec3.push_back(0.00633975 );   
  totalsfperenergy.push_back(0.00974097);  totalsfperenergy_sec1.push_back(0.0156076);  totalsfperenergy_sec2.push_back(0.00948473 ); totalsfperenergy_sec3.push_back(0.00629165 );   
  totalsfperenergy.push_back(0.00965875);  totalsfperenergy_sec1.push_back(0.0153688);  totalsfperenergy_sec2.push_back(0.00946882 ); totalsfperenergy_sec3.push_back(0.0062509  );   
  totalsfperenergy.push_back(0.00959886);  totalsfperenergy_sec1.push_back(0.0153031);  totalsfperenergy_sec2.push_back(0.00945922 ); totalsfperenergy_sec3.push_back(0.00622347 );   
  totalsfperenergy.push_back(0.00956122);  totalsfperenergy_sec1.push_back(0.0152271);  totalsfperenergy_sec2.push_back(0.00945212 ); totalsfperenergy_sec3.push_back(0.00623154 );   
  totalsfperenergy.push_back(0.00949007);  totalsfperenergy_sec1.push_back(0.0155053);  totalsfperenergy_sec2.push_back(0.0094405  ); totalsfperenergy_sec3.push_back(0.00618165 );   
  totalsfperenergy.push_back(0.00942227);  totalsfperenergy_sec1.push_back(0.0138983);  totalsfperenergy_sec2.push_back(0.0094273  ); totalsfperenergy_sec3.push_back(0.00619566 );   
  totalsfperenergy.push_back(0.00932817);  totalsfperenergy_sec1.push_back(0.0156628);  totalsfperenergy_sec2.push_back(0.00940749 ); totalsfperenergy_sec3.push_back(0.00616669 );   
  totalsfperenergy.push_back(0.00925142);  totalsfperenergy_sec1.push_back(0.0154427);  totalsfperenergy_sec2.push_back(0.00939253 ); totalsfperenergy_sec3.push_back(0.00617051 );
  //Also for the 5 sections plus 1 =  6 sections
  std::vector<double> totalsfperenergy_dedxsec1, totalsfperenergy_dedxsec2, totalsfperenergy_dedxsec3, totalsfperenergy_dedxsec4, totalsfperenergy_dedxsec5, totalsfperenergy_dedxsec6; 
  totalsfperenergy_dedxsec1.clear(); totalsfperenergy_dedxsec2.clear(); totalsfperenergy_dedxsec3.clear(); totalsfperenergy_dedxsec4.clear(); totalsfperenergy_dedxsec5.clear(); totalsfperenergy_dedxsec6.clear();
  totalsfperenergy_dedxsec1.push_back(0.0266367 ); totalsfperenergy_dedxsec2.push_back(0.0164444 ); totalsfperenergy_dedxsec3.push_back(0.00892063 );
  totalsfperenergy_dedxsec1.push_back(0.026344  ); totalsfperenergy_dedxsec2.push_back(0.0165866 ); totalsfperenergy_dedxsec3.push_back(0.00903433 );
  totalsfperenergy_dedxsec1.push_back(0.0259283 ); totalsfperenergy_dedxsec2.push_back(0.0168905 ); totalsfperenergy_dedxsec3.push_back(0.00922293 );
  totalsfperenergy_dedxsec1.push_back(0.0255365 ); totalsfperenergy_dedxsec2.push_back(0.0171529 ); totalsfperenergy_dedxsec3.push_back(0.00941933 );
  totalsfperenergy_dedxsec1.push_back(0.0251573 ); totalsfperenergy_dedxsec2.push_back(0.0173019 ); totalsfperenergy_dedxsec3.push_back(0.00948196 );
  totalsfperenergy_dedxsec1.push_back(0.0246383 ); totalsfperenergy_dedxsec2.push_back(0.0175144 ); totalsfperenergy_dedxsec3.push_back(0.00963975 );
  totalsfperenergy_dedxsec1.push_back(0.0235997 ); totalsfperenergy_dedxsec2.push_back(0.0178273 ); totalsfperenergy_dedxsec3.push_back(0.00987294 );
  totalsfperenergy_dedxsec1.push_back(0.0227062 ); totalsfperenergy_dedxsec2.push_back(0.0180658 ); totalsfperenergy_dedxsec3.push_back(0.010035   );
  totalsfperenergy_dedxsec1.push_back(0.0218953 ); totalsfperenergy_dedxsec2.push_back(0.0182689 ); totalsfperenergy_dedxsec3.push_back(0.0101785  );
  totalsfperenergy_dedxsec1.push_back(0.0214512 ); totalsfperenergy_dedxsec2.push_back(0.0183705 ); totalsfperenergy_dedxsec3.push_back(0.0102468  );
  totalsfperenergy_dedxsec1.push_back(0.0213141 ); totalsfperenergy_dedxsec2.push_back(0.0184516 ); totalsfperenergy_dedxsec3.push_back(0.0103043  );
  totalsfperenergy_dedxsec1.push_back(0.0208924 ); totalsfperenergy_dedxsec2.push_back(0.0185342 ); totalsfperenergy_dedxsec3.push_back(0.0103711  );
  totalsfperenergy_dedxsec1.push_back(0.0205356 ); totalsfperenergy_dedxsec2.push_back(0.0186429 ); totalsfperenergy_dedxsec3.push_back(0.010417   );
  totalsfperenergy_dedxsec1.push_back(0.0204838 ); totalsfperenergy_dedxsec2.push_back(0.0186889 ); totalsfperenergy_dedxsec3.push_back(0.0104491  );
  totalsfperenergy_dedxsec1.push_back(0.0202565 ); totalsfperenergy_dedxsec2.push_back(0.0187835 ); totalsfperenergy_dedxsec3.push_back(0.0105103  );
  totalsfperenergy_dedxsec1.push_back(0.0201504 ); totalsfperenergy_dedxsec2.push_back(0.0188686 ); totalsfperenergy_dedxsec3.push_back(0.0105651  );
  totalsfperenergy_dedxsec1.push_back(0.0200637 ); totalsfperenergy_dedxsec2.push_back(0.0189895 ); totalsfperenergy_dedxsec3.push_back(0.0106449  );
  totalsfperenergy_dedxsec1.push_back(0.0201277 ); totalsfperenergy_dedxsec2.push_back(0.0190774 ); totalsfperenergy_dedxsec3.push_back(0.0107149  );

  totalsfperenergy_dedxsec4.push_back(0.0118401 ); totalsfperenergy_dedxsec5.push_back(0.00730203 ); totalsfperenergy_dedxsec6.push_back(0.00326129 );
  totalsfperenergy_dedxsec4.push_back(0.0119534 ); totalsfperenergy_dedxsec5.push_back(0.00731113 ); totalsfperenergy_dedxsec6.push_back(0.0039805  );
  totalsfperenergy_dedxsec4.push_back(0.0120272 ); totalsfperenergy_dedxsec5.push_back(0.0075019  ); totalsfperenergy_dedxsec6.push_back(0.00443284 );
  totalsfperenergy_dedxsec4.push_back(0.0121294 ); totalsfperenergy_dedxsec5.push_back(0.007566   ); totalsfperenergy_dedxsec6.push_back(0.00449145 );
  totalsfperenergy_dedxsec4.push_back(0.0121561 ); totalsfperenergy_dedxsec5.push_back(0.00763044 ); totalsfperenergy_dedxsec6.push_back(0.00459179 );
  totalsfperenergy_dedxsec4.push_back(0.0122721 ); totalsfperenergy_dedxsec5.push_back(0.00770631 ); totalsfperenergy_dedxsec6.push_back(0.00459936 );
  totalsfperenergy_dedxsec4.push_back(0.0124597 ); totalsfperenergy_dedxsec5.push_back(0.00782464 ); totalsfperenergy_dedxsec6.push_back(0.00469215 );
  totalsfperenergy_dedxsec4.push_back(0.0125569 ); totalsfperenergy_dedxsec5.push_back(0.00790313 ); totalsfperenergy_dedxsec6.push_back(0.00474025 );
  totalsfperenergy_dedxsec4.push_back(0.0126658 ); totalsfperenergy_dedxsec5.push_back(0.00795067 ); totalsfperenergy_dedxsec6.push_back(0.00482071 );
  totalsfperenergy_dedxsec4.push_back(0.0127334 ); totalsfperenergy_dedxsec5.push_back(0.00797776 ); totalsfperenergy_dedxsec6.push_back(0.00485709 );
  totalsfperenergy_dedxsec4.push_back(0.0127877 ); totalsfperenergy_dedxsec5.push_back(0.00799488 ); totalsfperenergy_dedxsec6.push_back(0.00488831 );
  totalsfperenergy_dedxsec4.push_back(0.0128429 ); totalsfperenergy_dedxsec5.push_back(0.00801117 ); totalsfperenergy_dedxsec6.push_back(0.00493136 );
  totalsfperenergy_dedxsec4.push_back(0.0128991 ); totalsfperenergy_dedxsec5.push_back(0.00802767 ); totalsfperenergy_dedxsec6.push_back(0.0049548  );
  totalsfperenergy_dedxsec4.push_back(0.0129342 ); totalsfperenergy_dedxsec5.push_back(0.00803325 ); totalsfperenergy_dedxsec6.push_back(0.00497045 );
  totalsfperenergy_dedxsec4.push_back(0.0129835 ); totalsfperenergy_dedxsec5.push_back(0.00804717 ); totalsfperenergy_dedxsec6.push_back(0.00500066 );
  totalsfperenergy_dedxsec4.push_back(0.0130453 ); totalsfperenergy_dedxsec5.push_back(0.00805779 ); totalsfperenergy_dedxsec6.push_back(0.0050345  );
  totalsfperenergy_dedxsec4.push_back(0.0131222 ); totalsfperenergy_dedxsec5.push_back(0.00806754 ); totalsfperenergy_dedxsec6.push_back(0.00508263 );
  totalsfperenergy_dedxsec4.push_back(0.0131873 ); totalsfperenergy_dedxsec5.push_back(0.00807678 ); totalsfperenergy_dedxsec6.push_back(0.00511982 );

  lChain[0][0]->GetEntry(0);
  //For the sampling fraction related to the shower depth
  //Profile histo is 100 + 2 bins (100 bins + Underflow + Overflow) (0-10 range)
  //bin 1:0-0.1, bin2:0.1-0.2 ... 
  const int nb = 102;
  double sfpernormshowerdepth[numberofenergies][nb];
  double totalsfpershmax[numberofenergies][nb];
  double totalsf_sec1pershmax[numberofenergies][nb];
  double totalsf_sec2pershmax[numberofenergies][nb];
  double totalsf_sec3pershmax[numberofenergies][nb];
  //For the 5 plus 1 sections we will later give the first layer fixed sampling fraction
  double totalsf_dedxsec1pershmax[numberofenergies][nb];
  double totalsf_dedxsec2pershmax[numberofenergies][nb];
  double totalsf_dedxsec3pershmax[numberofenergies][nb];
  double totalsf_dedxsec4pershmax[numberofenergies][nb];
  double totalsf_dedxsec5pershmax[numberofenergies][nb];
  double totalsf_dedxsec6pershmax[numberofenergies][nb];

  double sfperlaypershmax[numberofenergies][(*ssvec).size()][nb];

  //Open the file to read the values 
  TFile *fin = new TFile("PFcal_12_combinedplots.root");
  TProfile *currentprof;

  TString proftitle;
  TString auxnam1_sf_prof;
  TString auxnam2_sf_prof;
  TString auxnam_sf_prof;

  for (int l=0; l<numberofenergies; l++){
    proftitle = "h_sfvsnormshowerdepth_" + IntToString(particleenergies[l]) + "_pfx";
    currentprof = (TProfile*) fin->Get(proftitle); 
    for (int i=1; i< currentprof->GetSize()-1; ++i){
      sfpernormshowerdepth[l][i] = currentprof->GetBinContent(i);
      //std::cout << "bin " << i << " low edge " <<currentprof->GetBinLowEdge(i) << std::endl;
    }
    proftitle = "TotalSFvsShowerDepth_" + IntToString(particleenergies[l]) + "_pfx";
    currentprof = (TProfile*) fin->Get(proftitle); 
    for (int i=1; i< currentprof->GetSize()-1; ++i){
      totalsfpershmax[l][i] = currentprof->GetBinContent(i);
      // std::cout << "bin " << i << " low edge " <<currentprof->GetBinLowEdge(i) << std::endl;
    }
    //Total sf section 1
    proftitle = "TotalSF_sec1vsShowerDepth_" + IntToString(particleenergies[l]) + "_pfx";
    currentprof = (TProfile*) fin->Get(proftitle); 
    for (int i=1; i< currentprof->GetSize()-1; ++i){
      totalsf_sec1pershmax[l][i] = currentprof->GetBinContent(i);
      // std::cout << "bin " << i << " low edge " <<currentprof->GetBinLowEdge(i) << std::endl;
    }
    //Total sf section 2
    proftitle = "TotalSF_sec2vsShowerDepth_" + IntToString(particleenergies[l]) + "_pfx";
    currentprof = (TProfile*) fin->Get(proftitle); 
    for (int i=1; i< currentprof->GetSize()-1; ++i){
      totalsf_sec2pershmax[l][i] = currentprof->GetBinContent(i);
      // std::cout << "bin " << i << " low edge " <<currentprof->GetBinLowEdge(i) << std::endl;
    }
    //Total sf section 3
    proftitle = "TotalSF_sec3vsShowerDepth_" + IntToString(particleenergies[l]) + "_pfx";
    currentprof = (TProfile*) fin->Get(proftitle); 
    for (int i=1; i< currentprof->GetSize()-1; ++i){
      totalsf_sec3pershmax[l][i] = currentprof->GetBinContent(i);
      // std::cout << "bin " << i << " low edge " <<currentprof->GetBinLowEdge(i) << std::endl;
    }
    //--------------------------------------
    //Here for the 5 sections ignoring the first one
    //Total sf section 2
    proftitle = "TotalSF_dedxsec2vsShowerDepth_" + IntToString(particleenergies[l]) + "_pfx";
    currentprof = (TProfile*) fin->Get(proftitle); 
    for (int i=1; i< currentprof->GetSize()-1; ++i){
      totalsf_dedxsec2pershmax[l][i] = currentprof->GetBinContent(i);
      // std::cout << "bin " << i << " low edge " <<currentprof->GetBinLowEdge(i) << std::endl;
    }
    //Total sf section 3
    proftitle = "TotalSF_dedxsec3vsShowerDepth_" + IntToString(particleenergies[l]) + "_pfx";
    currentprof = (TProfile*) fin->Get(proftitle); 
    for (int i=1; i< currentprof->GetSize()-1; ++i){
      totalsf_dedxsec3pershmax[l][i] = currentprof->GetBinContent(i);
      // std::cout << "bin " << i << " low edge " <<currentprof->GetBinLowEdge(i) << std::endl;
    }
    //Total sf section 4
    proftitle = "TotalSF_dedxsec4vsShowerDepth_" + IntToString(particleenergies[l]) + "_pfx";
    currentprof = (TProfile*) fin->Get(proftitle); 
    for (int i=1; i< currentprof->GetSize()-1; ++i){
      totalsf_dedxsec4pershmax[l][i] = currentprof->GetBinContent(i);
      // std::cout << "bin " << i << " low edge " <<currentprof->GetBinLowEdge(i) << std::endl;
    }
    //Total sf section 5
    proftitle = "TotalSF_dedxsec5vsShowerDepth_" + IntToString(particleenergies[l]) + "_pfx";
    currentprof = (TProfile*) fin->Get(proftitle); 
    for (int i=1; i< currentprof->GetSize()-1; ++i){
      totalsf_dedxsec5pershmax[l][i] = currentprof->GetBinContent(i);
      // std::cout << "bin " << i << " low edge " <<currentprof->GetBinLowEdge(i) << std::endl;
    }
    //Total sf section 6
    proftitle = "TotalSF_dedxsec6vsShowerDepth_" + IntToString(particleenergies[l]) + "_pfx";
    currentprof = (TProfile*) fin->Get(proftitle); 
    for (int i=1; i< currentprof->GetSize()-1; ++i){
      totalsf_dedxsec6pershmax[l][i] = currentprof->GetBinContent(i);
      // std::cout << "bin " << i << " low edge " <<currentprof->GetBinLowEdge(i) << std::endl;
    }
    //--------------------------------------
    //Here for the per layer per shower
    for(unsigned iL(startlayer); iL<(*ssvec).size(); iL++){
      auxnam1_sf_prof = "h_sfvsshowerdepthprofile_" + IntToString((int) iL);
      auxnam2_sf_prof = "_" + IntToString(particleenergies[l]) + + "_pfx";
      auxnam_sf_prof = auxnam1_sf_prof + auxnam2_sf_prof;
      currentprof = (TProfile*) fin->Get(auxnam_sf_prof); 
      for (int i=1; i< currentprof->GetSize()-1; ++i){
	sfperlaypershmax[l][iL][i] = currentprof->GetBinContent(i);
	// std::cout << "bin " << i << " low edge " <<currentprof->GetBinLowEdge(i) << std::endl;
      }
    }//Close loop on layers 
  }//close loop on energies  
  
  //Some canvases to same them in the same root file
  TCanvas *c_totsfdistrib = (TCanvas*) fin->Get("c3");//c_totsfdistrib->SaveAs();
  TCanvas *c_totsfvsene = (TCanvas*) fin->Get("c7");//c_totsfvsene->SaveAs();
  TCanvas *c_totsfvstprof = (TCanvas*) fin->Get("c8");//c_totsfvstprof->SaveAs();
  TCanvas *c_sfperlayervstprof = (TCanvas*) fin->Get("c16");//c_sfperlayervstprof->SaveAs();
  TCanvas *c_sfperlayervstprof_1 = (TCanvas*) fin->Get("c17_1");//c_sfperlayervstprof_1->SaveAs();
  TCanvas *c_sfperlayervstprof_2 = (TCanvas*) fin->Get("c17_2");//c_sfperlayervstprof_2->SaveAs();
  TCanvas *c_sf_sec1perlayervstprof = (TCanvas*) fin->Get("c18");//c_sfperlayervstprof->SaveAs();
  TCanvas *c_sf_sec2perlayervstprof = (TCanvas*) fin->Get("c19");//c_sfperlayervstprof->SaveAs();
  TCanvas *c_sf_sec3perlayervstprof = (TCanvas*) fin->Get("c20");//c_sfperlayervstprof->SaveAs();
  TCanvas *c_sf_dedxsec1perlayervstprof = (TCanvas*) fin->Get("c21");//c_sfperlayervstprof->SaveAs();
  TCanvas *c_sf_dedxsec2perlayervstprof = (TCanvas*) fin->Get("c22");//c_sfperlayervstprof->SaveAs();
  TCanvas *c_sf_dedxsec3perlayervstprof = (TCanvas*) fin->Get("c23");//c_sfperlayervstprof->SaveAs();
  TCanvas *c_sf_dedxsec4perlayervstprof = (TCanvas*) fin->Get("c24");//c_sfperlayervstprof->SaveAs();
  TCanvas *c_sf_dedxsec5perlayervstprof = (TCanvas*) fin->Get("c25");//c_sfperlayervstprof->SaveAs();
  TCanvas *c_sf_dedxsec6perlayervstprof = (TCanvas*) fin->Get("c26");//c_sfperlayervstprof->SaveAs();

  //fin->Close();

  //Be sure that you run on the correct input root file : PFcal_12_final.root
  TFile *finformax = new TFile("PFcal_12_final.root");
  TH1F *shmax;
  TString shmaxtitle;
  //Now for the shower max
  double meanshowermaxperenergy[numberofenergies];
  for (int l=0; l<numberofenergies; l++){
    shmaxtitle = "ShowerDepth_reco_" + IntToString(particleenergies[l]); 
    shmax = (TH1F*) finformax->Get(shmaxtitle); 
    meanshowermaxperenergy[l] = shmax->GetMean();
    std::cout << "Mean Shower max per energy - Energy:  " << particleenergies[l] << " GeV meanshowermaxperenergy " << meanshowermaxperenergy[l] << std::endl; 
  }
  //Now for the new shower depth variable
  //RECIPE: Go to the Samplingfraction.cpp and see there
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



    
  finformax->Close();

  //RECIPE for the profile of the totsf vs xomax/meanxomax
  //1. Run the Samplingfraction.cpp to make PFcal_12_final.root
  //2. Run the RecoE.cpp to get the meanshowermaxperenergy above. 
  //3. Rerun the Samplingfraction.cpp to get the profile histo

  //Here we save some more canvases for the pdf later
  lChain[0][0]->GetEntry(0);
  // TCanvas *c_sfperlayervst_prof[(*ssvec).size()];
  // for(unsigned iL(startlayer); iL<(*ssvec).size()-startlayer; iL++){
  //   c_sfperlayervst_prof[(int) iL] = (TCanvas*) fin->Get( ("c13_"+ IntToString((int) iL)).c_str() );
  // }

  //======================================================================================
  // Histos
  //Reco energy
  std::vector<TH1F*> h_recoE_raw;
  h_recoE_raw.clear();
  //Reco energy adding leakage
  std::vector<TH1F*> h_recoE_raw_leak;
  h_recoE_raw_leak.clear();
  //Reco energy with dedx
  std::vector<TH1F*> h_recoE_val;
  h_recoE_val.clear();
  //Reco energy using the total sampling fraction
  std::vector<TH1F*> h_recoE_totsf;
  h_recoE_totsf.clear();
  //Reco energy using the total sampling fraction adding leakage
  std::vector<TH1F*> h_recoE_totsf_leak;
  h_recoE_totsf_leak.clear();
  //Reco energy using the total sampling fraction
  std::vector<TH1F*> h_recoE_totsf_dedxweighted;
  h_recoE_totsf_dedxweighted.clear();
  //Reco energy using the total sampling fraction adding leakage
  std::vector<TH1F*> h_recoE_totsf_dedxweighted_leak;
  h_recoE_totsf_dedxweighted_leak.clear();
  //----------------------
  //Reco energy using the total sampling fraction per sections
  std::vector<TH1F*> h_recoE_totsf_sec;
  h_recoE_totsf_sec.clear();
  //Reco energy using the total sampling fraction adding leakage per sections
  std::vector<TH1F*> h_recoE_totsf_sec_leak;
  h_recoE_totsf_sec_leak.clear();
  //Reco energy using the total sampling fraction per sections and shower max
  std::vector<TH1F*> h_recoE_totsfshmax_sec;
  h_recoE_totsfshmax_sec.clear();
  //Reco energy using the total sampling fraction adding leakage per sections and shower max
  std::vector<TH1F*> h_recoE_totsfshmax_sec_leak;
  h_recoE_totsfshmax_sec_leak.clear();
  //----------------------
  //Reco energy using the total sampling fraction per dedx weighted sections 
  std::vector<TH1F*> h_recoE_totsf_sec_dedxweighted;
  h_recoE_totsf_sec_dedxweighted.clear();
  //Reco energy using the total sampling fraction adding leakage per per dedx weighted sections
  std::vector<TH1F*> h_recoE_totsf_sec_dedxweighted_leak;
  h_recoE_totsf_sec_dedxweighted_leak.clear();
  //Reco energy using the total sampling fraction per dedx weighted sections and shower max
  std::vector<TH1F*> h_recoE_totsfshmax_sec_dedxweighted;
  h_recoE_totsfshmax_sec_dedxweighted.clear();
  //Reco energy using the total sampling fraction adding leakage per dedx weighted sections and shower max
  std::vector<TH1F*> h_recoE_totsfshmax_sec_dedxweighted_leak;
  h_recoE_totsfshmax_sec_dedxweighted_leak.clear();
  //----------------------
  //Reco energy using the total sampling fraction per sections
  std::vector<TH1F*> h_recoE_totsf_dedxsec;
  h_recoE_totsf_dedxsec.clear();
  //Reco energy using the total sampling fraction adding leakage per dedxsections
  std::vector<TH1F*> h_recoE_totsf_dedxsec_leak;
  h_recoE_totsf_dedxsec_leak.clear();
  //Reco energy using the total sampling fraction per dedxsections and shower max
  std::vector<TH1F*> h_recoE_totsfshmax_dedxsec;
  h_recoE_totsfshmax_dedxsec.clear();
  //Reco energy using the total sampling fraction adding leakage per dedxsections and shower max
  std::vector<TH1F*> h_recoE_totsfshmax_dedxsec_leak;
  h_recoE_totsfshmax_dedxsec_leak.clear();
  //----------------------
  //----------------------
  //Reco energy using sampling fraction with normalized shower depth info
  std::vector<TH1F*> h_recoE_sfnorm;
  h_recoE_sfnorm.clear();
  //Reco energy using sampling fraction with normalized shower depth info adding leakage
  std::vector<TH1F*> h_recoE_sfnorm_leak;
  h_recoE_sfnorm_leak.clear();
  //Reco energy using sampling fraction with normalized shower depth info and cut
  std::vector<TH1F*> h_recoE_sfnorm_cut;
  h_recoE_sfnorm_cut.clear();
  //Reco energy using sampling fraction with normalized shower depth info adding leakage and cut
  std::vector<TH1F*> h_recoE_sfnorm_cut_leak;
  h_recoE_sfnorm_cut_leak.clear();
  //Reco energy using sampling fraction with normalized shower depth info
  std::vector<TH1F*> h_recoE_totsfnorm;
  h_recoE_totsfnorm.clear();
  //Reco energy using sampling fraction with normalized shower depth info adding leakage
  std::vector<TH1F*> h_recoE_totsfnorm_leak;
  h_recoE_totsfnorm_leak.clear();
  //Reco energy using sampling fraction with normalized shower depth info using dedx weight
  std::vector<TH1F*> h_recoE_totsfnorm_dedxweighted;
  h_recoE_totsfnorm_dedxweighted.clear();
  //Reco energy using sampling fraction with normalized shower depth info using dedx weight adding leakage
  std::vector<TH1F*> h_recoE_totsfnorm_dedxweighted_leak;
  h_recoE_totsfnorm_dedxweighted_leak.clear();
  //Reco energy using sampling fraction with normalized shower depth info and cut
  std::vector<TH1F*> h_recoE_totsfnorm_cut;
  h_recoE_totsfnorm_cut.clear();
  //Reco energy using sampling fraction with normalized shower depth info adding leakage and cut
  std::vector<TH1F*> h_recoE_totsfnorm_cut_leak;
  h_recoE_totsfnorm_cut_leak.clear();
  //Reco energy without upstream material calibration
  std::vector<TH1F*> h_recoE;
  h_recoE.clear();
  //Reco energy without upstream material calibration including leakage
  std::vector<TH1F*> h_recoE_leak;
  h_recoE_leak.clear();
  //Reco energy without upstream material calibration including leakage and offset
  std::vector<TH1F*> h_recoE_leak_offset;
  h_recoE_leak_offset.clear();
  //Reco energy with upstream material calibration
  std::vector<TH1F*> h_recoE_corr;
  h_recoE_corr.clear();
  //Reco energy over true without upstream material calibration
  std::vector<TH1F*> h_recoEovertrueE;
  h_recoEovertrueE.clear();
  //Reco energy over true without upstream material calibration including leakage
  std::vector<TH1F*> h_recoEovertrueE_leak;
  h_recoEovertrueE_leak.clear();
  //Reco energy over true without upstream material calibration including leakage and offset
  std::vector<TH1F*> h_recoEovertrueE_leak_offset;
  h_recoEovertrueE_leak_offset.clear();
  //----------------------
  //Reco energy over true without upstream material calibration with total SF
  std::vector<TH1F*> h_recoEovertrueE_totsf;
  h_recoEovertrueE_totsf.clear();
  //Reco energy over true without upstream material calibration with total SF including leakage
  std::vector<TH1F*> h_recoEovertrueE_totsf_leak;
  h_recoEovertrueE_totsf_leak.clear();
  //Reco energy over true without upstream material calibration with total SF
  std::vector<TH1F*> h_recoEovertrueE_totsf_dedxweighted;
  h_recoEovertrueE_totsf_dedxweighted.clear();
  //Reco energy over true without upstream material calibration with total SF including leakage
  std::vector<TH1F*> h_recoEovertrueE_totsf_dedxweighted_leak;
  h_recoEovertrueE_totsf_dedxweighted_leak.clear();
  //----------------------
  //Reco energy over true without upstream material calibration with total SF
  std::vector<TH1F*> h_recoEovertrueE_totsf_sec;
  h_recoEovertrueE_totsf_sec.clear();
  //Reco energy over true without upstream material calibration with total SF including leakage
  std::vector<TH1F*> h_recoEovertrueE_totsf_sec_leak;
  h_recoEovertrueE_totsf_sec_leak.clear();
  //Reco energy over true without upstream material calibration with total SF and shower max
  std::vector<TH1F*> h_recoEovertrueE_totsfshmax_sec;
  h_recoEovertrueE_totsfshmax_sec.clear();
  //Reco energy over true without upstream material calibration with total SF including leakage and shower max
  std::vector<TH1F*> h_recoEovertrueE_totsfshmax_sec_leak;
  h_recoEovertrueE_totsfshmax_sec_leak.clear();
  //Reco energy over true without upstream material calibration with total SF and shower max
  std::vector<TH1F*> h_recoEovertrueE_totsfshmax_sec_dedxweighted;
  h_recoEovertrueE_totsfshmax_sec_dedxweighted.clear();
  //Reco energy over true without upstream material calibration with total SF including leakage and shower max
  std::vector<TH1F*> h_recoEovertrueE_totsfshmax_sec_dedxweighted_leak;
  h_recoEovertrueE_totsfshmax_sec_dedxweighted_leak.clear();
  //----------------------
  //Reco energy over true without upstream material calibration with total SF weighted with dedx
  std::vector<TH1F*> h_recoEovertrueE_totsf_sec_dedxweighted;
  h_recoEovertrueE_totsf_sec_dedxweighted.clear();
  //Reco energy over true without upstream material calibration with total SF weighted with dedx including leakage
  std::vector<TH1F*> h_recoEovertrueE_totsf_sec_dedxweighted_leak;
  h_recoEovertrueE_totsf_sec_dedxweighted_leak.clear();
  //----------------------
  //Reco energy over true without upstream material calibration with total SF 
  std::vector<TH1F*> h_recoEovertrueE_totsf_dedxsec;
  h_recoEovertrueE_totsf_dedxsec.clear();
  //Reco energy over true without upstream material calibration with total SF including leakage
  std::vector<TH1F*> h_recoEovertrueE_totsf_dedxsec_leak;
  h_recoEovertrueE_totsf_dedxsec_leak.clear();
  //Reco energy over true without upstream material calibration with total SF and shower max
  std::vector<TH1F*> h_recoEovertrueE_totsfshmax_dedxsec;
  h_recoEovertrueE_totsfshmax_dedxsec.clear();
  //Reco energy over true without upstream material calibration with total SF including leakage and shower max
  std::vector<TH1F*> h_recoEovertrueE_totsfshmax_dedxsec_leak;
  h_recoEovertrueE_totsfshmax_dedxsec_leak.clear();
  //----------------------
  //Reco energy over true using sampling fraction with normalized shower depth info
  std::vector<TH1F*> h_recoEovertrueE_sfnorm;
  h_recoEovertrueE_sfnorm.clear();
  //Reco energy over true using sampling fraction with normalized shower depth info adding leakage
  std::vector<TH1F*> h_recoEovertrueE_sfnorm_leak;
  h_recoEovertrueE_sfnorm_leak.clear();
  //Reco energy over true using sampling fraction with normalized shower depth info and cut
  std::vector<TH1F*> h_recoEovertrueE_sfnorm_cut;
  h_recoEovertrueE_sfnorm_cut.clear();
  //Reco energy over true using sampling fraction with normalized shower depth info adding leakage and cut
  std::vector<TH1F*> h_recoEovertrueE_sfnorm_cut_leak;
  h_recoEovertrueE_sfnorm_cut_leak.clear();
  //Reco energy over true without upstream material calibration with total SF and shower max
  std::vector<TH1F*> h_recoEovertrueE_totsfnorm;
  h_recoEovertrueE_totsfnorm.clear();
  //Reco energy over true without upstream material calibration with total SF and shower max including leakage
  std::vector<TH1F*> h_recoEovertrueE_totsfnorm_leak;
  h_recoEovertrueE_totsfnorm_leak.clear();
  //Reco energy over true without upstream material calibration with total SF using dedx weight and shower max
  std::vector<TH1F*> h_recoEovertrueE_totsfnorm_dedxweighted;
  h_recoEovertrueE_totsfnorm_dedxweighted.clear();
  //Reco energy over true without upstream material calibration with total SF using dedx weight and shower max including leakage
  std::vector<TH1F*> h_recoEovertrueE_totsfnorm_dedxweighted_leak;
  h_recoEovertrueE_totsfnorm_dedxweighted_leak.clear();
  //Reco energy over true without upstream material calibration with total SF and shower max cut
  std::vector<TH1F*> h_recoEovertrueE_totsfnorm_cut;
  h_recoEovertrueE_totsfnorm_cut.clear();
  //Reco energy over true without upstream material calibration with total SF and shower max cut including leakage
  std::vector<TH1F*> h_recoEovertrueE_totsfnorm_cut_leak;
  h_recoEovertrueE_totsfnorm_cut_leak.clear();
  //Reco energy over true with upstream material calibration
  std::vector<TH1F*> h_recoE_corrovertrueE;
  h_recoE_corrovertrueE.clear();
  std::vector<TH1F*> h_recoE_valEovertrueE;
  h_recoE_valEovertrueE.clear();
  std::vector<TH1F*> h_recoE_rawEovertrueE;
  h_recoE_rawEovertrueE.clear();
  //Reco energy without upstream material calibration and 4 sections
  std::vector<TH1F*> h_recoE_4sections;
  h_recoE_4sections.clear();
  //Reco energy without upstream material calibration including leakage and 4 sections
  std::vector<TH1F*> h_recoE_4sections_leak;
  h_recoE_4sections_leak.clear();
  //Reco energy without upstream material calibration and 4 sections : Cheat section 1
  std::vector<TH1F*> h_recoE_4sections_ch1;
  h_recoE_4sections_ch1.clear();
  //Reco energy without upstream material calibration including leakage and 4 sections : Cheat section 1
  std::vector<TH1F*> h_recoE_4sections_ch1_leak;
  h_recoE_4sections_ch1_leak.clear();
  //Reco energy without upstream material calibration and 4 sections : Cheat section 2
  std::vector<TH1F*> h_recoE_4sections_ch2;
  h_recoE_4sections_ch2.clear();
  //Reco energy without upstream material calibration including leakage and 4 sections : Cheat section 2
  std::vector<TH1F*> h_recoE_4sections_ch2_leak;
  h_recoE_4sections_ch2_leak.clear();
  //Reco energy without upstream material calibration and 4 sections : Cheat section 3
  std::vector<TH1F*> h_recoE_4sections_ch3;
  h_recoE_4sections_ch3.clear();
  //Reco energy without upstream material calibration including leakage and 4 sections : Cheat section 3
  std::vector<TH1F*> h_recoE_4sections_ch3_leak;
  h_recoE_4sections_ch3_leak.clear();
  //Reco energy without upstream material calibration and 4 sections : Cheat section 4
  std::vector<TH1F*> h_recoE_4sections_ch4;
  h_recoE_4sections_ch4.clear();
  //Reco energy without upstream material calibration including leakage and 4 sections : Cheat section 4
  std::vector<TH1F*> h_recoE_4sections_ch4_leak;
  h_recoE_4sections_ch4_leak.clear();
  //Reco energy without upstream material calibration : Cheat section 1
  std::vector<TH1F*> h_recoE_ch1;
  h_recoE_ch1.clear();
  //Reco energy without upstream material calibration including leakage : Cheat section 1
  std::vector<TH1F*> h_recoE_ch1_leak;
  h_recoE_ch1_leak.clear();
  //Reco energy without upstream material calibration : Cheat section 2
  std::vector<TH1F*> h_recoE_ch2;
  h_recoE_ch2.clear();
  //Reco energy without upstream material calibration including leakage : Cheat section 2
  std::vector<TH1F*> h_recoE_ch2_leak;
  h_recoE_ch2_leak.clear();
  //Reco energy without upstream material calibration : Cheat section 3
  std::vector<TH1F*> h_recoE_ch3;
  h_recoE_ch3.clear();
  //Reco energy without upstream material calibration including leakage : Cheat section 3
  std::vector<TH1F*> h_recoE_ch3_leak;
  h_recoE_ch3_leak.clear();
  //Reco energy without upstream material calibration : Cheat section 4
  std::vector<TH1F*> h_recoE_ch4;
  h_recoE_ch4.clear();
  //Reco energy without upstream material calibration including leakage : Cheat section 4
  std::vector<TH1F*> h_recoE_ch4_leak;
  h_recoE_ch4_leak.clear();
  //----------------------
  //Reco energy over trueE without upstream material calibration : Cheat section 1
  std::vector<TH1F*> h_recoE_ch1_overtrueE;
  h_recoE_ch1_overtrueE.clear();
  //Reco energy over trueE without upstream material calibration including leakage : Cheat section 1
  std::vector<TH1F*> h_recoE_ch1_overtrueE_leak;
  h_recoE_ch1_overtrueE_leak.clear();
  //Reco energy over trueE without upstream material calibration : Cheat section 2
  std::vector<TH1F*> h_recoE_ch2_overtrueE;
  h_recoE_ch2_overtrueE.clear();
  //Reco energy over trueE without upstream material calibration including leakage : Cheat section 2
  std::vector<TH1F*> h_recoE_ch2_overtrueE_leak;
  h_recoE_ch2_overtrueE_leak.clear();
  //Reco energy over trueE without upstream material calibration : Cheat section 3
  std::vector<TH1F*> h_recoE_ch3_overtrueE;
  h_recoE_ch3_overtrueE.clear();
  //Reco energy over trueE without upstream material calibration including leakage : Cheat section 3
  std::vector<TH1F*> h_recoE_ch3_overtrueE_leak;
  h_recoE_ch3_overtrueE_leak.clear();
  //Reco energy over trueE without upstream material calibration : Cheat section 4
  std::vector<TH1F*> h_recoE_ch4_overtrueE;
  h_recoE_ch4_overtrueE.clear();
  //Reco energy over trueE without upstream material calibration including leakage : Cheat section 4
  std::vector<TH1F*> h_recoE_ch4_overtrueE_leak;
  h_recoE_ch4_overtrueE_leak.clear();
  //The shower max distribution
  std::vector<TH1F*> h_Showermaxovermeanshmax;
  h_Showermaxovermeanshmax.clear();
  //Reco energy with the shower max method
  std::vector<TH1F*> h_recoE_shmax;
  h_recoE_shmax.clear();
  //Reco energy with the shower max method including leakage
  std::vector<TH1F*> h_recoE_shmax_leak;
  h_recoE_shmax_leak.clear();
  //Reco energy over trueE with the shower max method
  std::vector<TH1F*> h_recoEovertrueE_shmax;
  h_recoEovertrueE_shmax.clear();
  //Reco energy over trueE with the shower max method including leakage
  std::vector<TH1F*> h_recoEovertrueE_shmax_leak;
  h_recoEovertrueE_shmax_leak.clear();

  for (Int_t k=0; k<numberofenergies; k++){
    // h_recoE.push_back(new TH1F(("h_recoE_" + IntToString(particleenergies[k])).c_str(),";E (MIPs);Events",2048,0.,7000000.));
    h_recoE.push_back(new TH1F(("h_recoE_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",5500,0.,550.));
    h_recoE_leak.push_back(new TH1F(("h_recoE_leak_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",5500,0.,550.));
    h_recoE_leak_offset.push_back(new TH1F(("h_recoE_leak_offset_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",5500,0.,550.));
    h_recoE_corr.push_back(new TH1F(("h_recoE_corr_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",1000.,0.,1000.));
    h_recoEovertrueE.push_back(new TH1F(("h_recoEovertrueE_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",1024.,-1.0,1.0));
    h_recoEovertrueE_leak.push_back(new TH1F(("h_recoEovertrueE_leak_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",1024.,-1.0,1.0));
    h_recoEovertrueE_leak_offset.push_back(new TH1F(("h_recoEovertrueE_leak_offset_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",1024.,-1.0,1.0));
    h_recoE_corrovertrueE.push_back(new TH1F(("h_recoE_corrovertrueE_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",1024.,-1.0,1.0));
    h_recoE_raw.push_back(new TH1F(("h_recoE_raw_" + IntToString(particleenergies[k])).c_str(),";E (MIPs);Events",10000,0.,100000.));
    // h_recoE_raw_leak.push_back(new TH1F(("h_recoE_raw_leak_" + IntToString(particleenergies[k])).c_str(),";E (MIPs);Events",10000,0.,500000.));
    h_recoE_raw_leak.push_back(new TH1F(("h_recoE_raw_leak_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",5500,0.,550.));
    h_recoE_val.push_back(new TH1F(("h_recoE_val_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",5500,0.,550.));
    h_recoE_totsf.push_back(new TH1F(("h_recoE_totsf_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",5500,0.,550.));
    h_recoE_totsf_leak.push_back(new TH1F(("h_recoE_totsf_leak_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",5500,0.,550.));
    h_recoEovertrueE_totsf.push_back(new TH1F(("h_recoEovertrueE_totsf_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",1024.,-1.0,1.0));
    h_recoEovertrueE_totsf_leak.push_back(new TH1F(("h_recoEovertrueE_totsf_leak_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",1024.,-1.0,1.0));
    h_recoE_totsf_dedxweighted.push_back(new TH1F(("h_recoE_totsf_dedxweighted_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",5500,0.,550.));
    h_recoE_totsf_dedxweighted_leak.push_back(new TH1F(("h_recoE_totsf_dedxweighted_leak_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",5500,0.,550.));
    h_recoEovertrueE_totsf_dedxweighted.push_back(new TH1F(("h_recoEovertrueE_totsf_dedxweighted_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",1024.,-1.0,1.0));
    h_recoEovertrueE_totsf_dedxweighted_leak.push_back(new TH1F(("h_recoEovertrueE_totsf_dedxweighted_leak_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",1024.,-1.0,1.0));
    h_recoE_totsf_sec.push_back(new TH1F(("h_recoE_totsf_sec_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",5500,0.,550.));
    h_recoE_totsf_sec_leak.push_back(new TH1F(("h_recoE_totsf_sec_leak_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",5500,0.,550.));
    h_recoEovertrueE_totsf_sec.push_back(new TH1F(("h_recoEovertrueE_totsf_sec_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",1024.,-1.0,1.0));
    h_recoEovertrueE_totsf_sec_leak.push_back(new TH1F(("h_recoEovertrueE_totsf_sec_leak_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",1024.,-1.0,1.0));
    h_recoE_totsf_sec_dedxweighted.push_back(new TH1F(("h_recoE_totsf_sec_dedxweighted_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",5500,0.,550.));
    h_recoE_totsf_sec_dedxweighted_leak.push_back(new TH1F(("h_recoE_totsf_sec_dedxweighted_leak_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",5500,0.,550.));
    h_recoEovertrueE_totsf_sec_dedxweighted.push_back(new TH1F(("h_recoEovertrueE_totsf_sec_dedxweighted_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",1024.,-1.0,1.0));
    h_recoEovertrueE_totsf_sec_dedxweighted_leak.push_back(new TH1F(("h_recoEovertrueE_totsf_sec_dedxweighted_leak_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",1024.,-1.0,1.0));
    h_recoE_totsfshmax_sec.push_back(new TH1F(("h_recoE_totsfshmax_sec_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",5500,0.,550.));
    h_recoE_totsfshmax_sec_leak.push_back(new TH1F(("h_recoE_totsfshmax_sec_leak_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",5500,0.,550.));
    h_recoEovertrueE_totsfshmax_sec.push_back(new TH1F(("h_recoEovertrueE_totsfshmax_sec_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",1024.,-1.0,1.0));
    h_recoEovertrueE_totsfshmax_sec_leak.push_back(new TH1F(("h_recoEovertrueE_totsfshmax_sec_leak_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",1024.,-1.0,1.0));
    h_recoE_totsfshmax_sec_dedxweighted.push_back(new TH1F(("h_recoE_totsfshmax_sec_dedxweighted_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",5500,0.,550.));
    h_recoE_totsfshmax_sec_dedxweighted_leak.push_back(new TH1F(("h_recoE_totsfshmax_sec_dedxweighted_leak_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",5500,0.,550.));
    h_recoEovertrueE_totsfshmax_sec_dedxweighted.push_back(new TH1F(("h_recoEovertrueE_totsfshmax_sec_dedxweighted_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",1024.,-1.0,1.0));
    h_recoEovertrueE_totsfshmax_sec_dedxweighted_leak.push_back(new TH1F(("h_recoEovertrueE_totsfshmax_sec_dedxweighted_leak_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",1024.,-1.0,1.0));
    //dedx more sections
    h_recoE_totsf_dedxsec.push_back(new TH1F(("h_recoE_totsf_dedxsec_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",5500,0.,550.));
    h_recoE_totsf_dedxsec_leak.push_back(new TH1F(("h_recoE_totsf_dedxsec_leak_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",5500,0.,550.));
    h_recoEovertrueE_totsf_dedxsec.push_back(new TH1F(("h_recoEovertrueE_totsf_dedxsec_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",1024.,-1.0,1.0));
    h_recoEovertrueE_totsf_dedxsec_leak.push_back(new TH1F(("h_recoEovertrueE_totsf_dedxsec_leak_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",1024.,-1.0,1.0));
    h_recoE_totsfshmax_dedxsec.push_back(new TH1F(("h_recoE_totsfshmax_dedxsec_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",5500,0.,550.));
    h_recoE_totsfshmax_dedxsec_leak.push_back(new TH1F(("h_recoE_totsfshmax_dedxsec_leak_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",5500,0.,550.));
    h_recoEovertrueE_totsfshmax_dedxsec.push_back(new TH1F(("h_recoEovertrueE_totsfshmax_dedxsec_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",1024.,-1.0,1.0));
    h_recoEovertrueE_totsfshmax_dedxsec_leak.push_back(new TH1F(("h_recoEovertrueE_totsfshmax_dedxsec_leak_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",1024.,-1.0,1.0));
    h_recoEovertrueE_totsfnorm.push_back(new TH1F(("h_recoEovertrueE_totsfnorm_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",1024.,-1.0,1.0));
    h_recoEovertrueE_totsfnorm_leak.push_back(new TH1F(("h_recoEovertrueE_totsfnorm_leak_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",1024.,-1.0,1.0));
    h_recoEovertrueE_totsfnorm_dedxweighted.push_back(new TH1F(("h_recoEovertrueE_totsfnorm_dedxweighted_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",1024.,-1.0,1.0));
    h_recoEovertrueE_totsfnorm_dedxweighted_leak.push_back(new TH1F(("h_recoEovertrueE_totsfnorm_dedxweighted_leak_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",1024.,-1.0,1.0));
    h_recoEovertrueE_totsfnorm_cut.push_back(new TH1F(("h_recoEovertrueE_totsfnorm_cut_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",1024.,-1.0,1.0));
    h_recoEovertrueE_totsfnorm_cut_leak.push_back(new TH1F(("h_recoEovertrueE_totsfnorm_cut_leak_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",1024.,-1.0,1.0));
    h_recoE_sfnorm.push_back(new TH1F(("h_recoE_sfnorm_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",7000,0.,700.));
    h_recoE_sfnorm_leak.push_back(new TH1F(("h_recoE_sfnorm_leak_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",7000,0.,700.));
    h_recoE_sfnorm_cut.push_back(new TH1F(("h_recoE_sfnorm_cut_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",7000,0.,700.));
    h_recoE_sfnorm_cut_leak.push_back(new TH1F(("h_recoE_sfnorm_cut_leak_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",7000,0.,700.));
    h_recoE_totsfnorm.push_back(new TH1F(("h_recoE_totsfnorm_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",7000,0.,700.));
    h_recoE_totsfnorm_leak.push_back(new TH1F(("h_recoE_totsfnorm_leak_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",7000,0.,700.));
    h_recoE_totsfnorm_dedxweighted.push_back(new TH1F(("h_recoE_totsfnorm_dedxweighted_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",7000,0.,700.));
    h_recoE_totsfnorm_dedxweighted_leak.push_back(new TH1F(("h_recoE_totsfnorm_dedxweighted_leak_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",7000,0.,700.));
    h_recoE_totsfnorm_cut.push_back(new TH1F(("h_recoE_totsfnorm_cut_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",7000,0.,700.));
    h_recoE_totsfnorm_cut_leak.push_back(new TH1F(("h_recoE_totsfnorm_cut_leak_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",7000,0.,700.));
    h_recoEovertrueE_sfnorm.push_back(new TH1F(("h_recoEovertrueE_sfnorm_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",1024.,-1.0,1.0));
    h_recoEovertrueE_sfnorm_leak.push_back(new TH1F(("h_recoEovertrueE_sfnorm_leak_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",1024.,-1.0,1.0));
    h_recoEovertrueE_sfnorm_cut.push_back(new TH1F(("h_recoEovertrueE_sfnorm_cut_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",1024.,-1.0,1.0));
    h_recoEovertrueE_sfnorm_cut_leak.push_back(new TH1F(("h_recoEovertrueE_sfnorm_cut_leak_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",1024.,-1.0,1.0));
    h_recoE_valEovertrueE.push_back(new TH1F(("h_recoE_valEovertrueE_" + IntToString(particleenergies[k])).c_str(),";Erec-Etrue/Etrue;Events",1024.,-1.0,1.0));
    h_recoE_rawEovertrueE.push_back(new TH1F(("h_recoE_rawEovertrueE_" + IntToString(particleenergies[k])).c_str(),";Erec-Etrue/Etrue;Events",1024.,-1.03,-0.97));

    h_recoE_4sections.push_back(new TH1F(("h_recoE_4sections_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",5500,0.,550.));
    h_recoE_4sections_leak.push_back(new TH1F(("h_recoE_4sections_leak_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",5500,0.,550.));
    h_recoE_4sections_ch1.push_back(new TH1F(("h_recoE_4sections_ch1_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",5500,0.,550.));
    h_recoE_4sections_ch1_leak.push_back(new TH1F(("h_recoE_4sections_ch1_leak_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",5500,0.,550.));
    h_recoE_4sections_ch2.push_back(new TH1F(("h_recoE_4sections_ch2_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",5500,0.,550.));
    h_recoE_4sections_ch2_leak.push_back(new TH1F(("h_recoE_4sections_ch2_leak_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",5500,0.,550.));
    h_recoE_4sections_ch3.push_back(new TH1F(("h_recoE_4sections_ch3_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",5500,0.,550.));
    h_recoE_4sections_ch3_leak.push_back(new TH1F(("h_recoE_4sections_ch3_leak_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",5500,0.,550.));
    h_recoE_4sections_ch4.push_back(new TH1F(("h_recoE_4sections_ch4_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",5500,0.,550.));
    h_recoE_4sections_ch4_leak.push_back(new TH1F(("h_recoE_4sections_ch4_leak_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",5500,0.,550.));
    h_recoE_ch1.push_back(new TH1F(("h_recoE_ch1_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",5500,0.,550.));
    h_recoE_ch1_leak.push_back(new TH1F(("h_recoE_ch1_leak_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",5500,0.,550.));
    h_recoE_ch2.push_back(new TH1F(("h_recoE_ch2_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",5500,0.,550.));
    h_recoE_ch2_leak.push_back(new TH1F(("h_recoE_ch2_leak_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",5500,0.,550.));
    h_recoE_ch3.push_back(new TH1F(("h_recoE_ch3_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",5500,0.,550.));
    h_recoE_ch3_leak.push_back(new TH1F(("h_recoE_ch3_leak_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",5500,0.,550.));
    h_recoE_ch4.push_back(new TH1F(("h_recoE_ch4_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",5500,0.,550.));
    h_recoE_ch4_leak.push_back(new TH1F(("h_recoE_ch4_leak_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",5500,0.,550.));
    //overtrueE
    h_recoE_ch1_overtrueE.push_back(new TH1F(("h_recoE_ch1_overtrueE_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",1024.,-1.0,1.0));
    h_recoE_ch1_overtrueE_leak.push_back(new TH1F(("h_recoE_ch1_overtrueE_leak_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",1024.,-1.0,1.0));
    h_recoE_ch2_overtrueE.push_back(new TH1F(("h_recoE_ch2_overtrueE_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",1024.,-1.0,1.0));
    h_recoE_ch2_overtrueE_leak.push_back(new TH1F(("h_recoE_ch2_overtrueE_leak_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",1024.,-1.0,1.0));
    h_recoE_ch3_overtrueE.push_back(new TH1F(("h_recoE_ch3_overtrueE_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",1024.,-1.0,1.0));
    h_recoE_ch3_overtrueE_leak.push_back(new TH1F(("h_recoE_ch3_overtrueE_leak_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",1024.,-1.0,1.0));
    h_recoE_ch4_overtrueE.push_back(new TH1F(("h_recoE_ch4_overtrueE_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",1024.,-1.0,1.0));
    h_recoE_ch4_overtrueE_leak.push_back(new TH1F(("h_recoE_ch4_overtrueE_leak_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",1024.,-1.0,1.0));
    h_Showermaxovermeanshmax.push_back(new TH1F( ("h_Showermaxovermeanshmax_"+ IntToString(particleenergies[k])).c_str(),";X_{0,max}/<X_{0,max}>;Events",1000.,0.,10.));
    h_recoE_shmax.push_back(new TH1F(("h_recoE_shmax_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",5500,0.,550.));
    h_recoE_shmax_leak.push_back(new TH1F(("h_recoE_shmax_leak_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",5500,0.,550.));
    h_recoEovertrueE_shmax.push_back(new TH1F(("h_recoEovertrueE_shmax_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",1024.,-1.0,1.0));
    h_recoEovertrueE_shmax_leak.push_back(new TH1F(("h_recoEovertrueE_shmax_leak_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",1024.,-1.0,1.0));
  }

  //Open file and save the scale factors
  TString titleoffile1 = "dE_lossUpvsMIPs_" + particle + ".txt";
  std::ifstream infile1(titleoffile1);
  std::string dummyforparticle; double dummyforenergy;
  double dE_lossUpvsMIPs_fitresults_slope[numberofenergies];
  double dE_lossUpvsMIPs_fitresults_const[numberofenergies];
  double dE_lossUpvsMIPs_fitresults_slope_error[numberofenergies];
  double dE_lossUpvsMIPs_fitresults_const_error[numberofenergies];

  std::string buffer;
  for (Int_t k=0; k<numberofenergies; k++){
    getline(infile1, buffer);
    std::istringstream is(buffer);
    //particle particleenergy constant slope errorconstant errorslope
    is >> dummyforparticle >> dummyforenergy >> dE_lossUpvsMIPs_fitresults_const[k] >> dE_lossUpvsMIPs_fitresults_slope[k] >> dE_lossUpvsMIPs_fitresults_const_error[k] >> dE_lossUpvsMIPs_fitresults_slope_error[k];
    //Test
    std::cout << dummyforparticle << " " << dummyforenergy << " " << dE_lossUpvsMIPs_fitresults_const[k] << " " << dE_lossUpvsMIPs_fitresults_slope[k] << " " << " " << dE_lossUpvsMIPs_fitresults_const_error[k] << " " << dE_lossUpvsMIPs_fitresults_slope_error[k] << std::endl;

  }
  infile1.close();

  TString titleoffile2 = "MIPSvsdEplusMIPs_" + particle + ".txt";
  std::ifstream infile2(titleoffile2);
  const int numoflayers = 31;
  double dEvsMIPs_fitresults_slope[numberofenergies][numoflayers];
  double dEvsMIPs_fitresults_const[numberofenergies][numoflayers];
  double dEvsMIPs_fitresults_slope_error[numberofenergies][numoflayers];
  double dEvsMIPs_fitresults_const_error[numberofenergies][numoflayers];

  while ( getline(infile2,buffer) ) {
    std::istringstream is(buffer);
    //particle particleenergy layer constant slope errorconstant errorslope
    double currfitconst, currfitslope, currfitconsterr, currfitslopeerr;
    int currlayer;
    is >> dummyforparticle >> dummyforenergy >> currlayer >> currfitconst >> currfitslope >> currfitconsterr >> currfitslopeerr;
    int k = 0;
    if (dummyforenergy == 2){k=0;}
    if (dummyforenergy == 3){k=1;}
    if (dummyforenergy == 5){k=2;}
    if (dummyforenergy == 8){k=3;}
    if (dummyforenergy == 10){k=4;}
    if (dummyforenergy == 15){k=5;}
    if (dummyforenergy == 30){k=6;}
    if (dummyforenergy == 50){k=7;}
    if (dummyforenergy == 80){k=8;}
    if (dummyforenergy == 100){k=9;}
    if (dummyforenergy == 120){k=10;}
    if (dummyforenergy == 150){k=11;}
    if (dummyforenergy == 180){k=12;}
    if (dummyforenergy == 200){k=13;}
    if (dummyforenergy == 250){k=14;}
    if (dummyforenergy == 300){k=15;}
    if (dummyforenergy == 400){k=16;}
    if (dummyforenergy == 500){k=17;}

    dEvsMIPs_fitresults_const[k][currlayer] = currfitconst;
    dEvsMIPs_fitresults_slope[k][currlayer] = currfitslope;
    dEvsMIPs_fitresults_const_error[k][currlayer] = currfitconsterr;
    dEvsMIPs_fitresults_slope_error[k][currlayer] = currfitslopeerr;

    std::cout << dummyforparticle << " " << dummyforenergy << " " << currlayer << " " << dEvsMIPs_fitresults_const[k][currlayer] << " " << dEvsMIPs_fitresults_slope[k][currlayer] << " " << " " << dEvsMIPs_fitresults_const_error[k][currlayer] << " " << dEvsMIPs_fitresults_slope_error[k][currlayer] << std::endl;


  }

  infile2.close();
  
  //Grouping in sections 1-6, 7-12, 13-20, 21-30: 4 sections
  double dEvsMIPs_fitresults_slope_4sections[numberofenergies][numoflayers];
  lChain[0][0]->GetEntry(0);
  for (int l=0; l<numberofenergies; l++){
    //1-6
    double dEvsMIPs_sec1 = dEvsMIPs_fitresults_slope[l][0] + dEvsMIPs_fitresults_slope[l][1] + dEvsMIPs_fitresults_slope[l][2] + dEvsMIPs_fitresults_slope[l][3] + dEvsMIPs_fitresults_slope[l][4] + dEvsMIPs_fitresults_slope[l][5]; 
    //7-12
    double dEvsMIPs_sec2 = dEvsMIPs_fitresults_slope[l][6] + dEvsMIPs_fitresults_slope[l][7] + dEvsMIPs_fitresults_slope[l][8] + dEvsMIPs_fitresults_slope[l][9] + dEvsMIPs_fitresults_slope[l][10] + dEvsMIPs_fitresults_slope[l][11]; 
    //13-20
    double dEvsMIPs_sec3 = dEvsMIPs_fitresults_slope[l][12] + dEvsMIPs_fitresults_slope[l][13] + dEvsMIPs_fitresults_slope[l][14] + dEvsMIPs_fitresults_slope[l][15] + dEvsMIPs_fitresults_slope[l][16] + dEvsMIPs_fitresults_slope[l][17] + dEvsMIPs_fitresults_slope[l][18]+ dEvsMIPs_fitresults_slope[l][19];
    //21-30
    double dEvsMIPs_sec4= dEvsMIPs_fitresults_slope[l][20] + dEvsMIPs_fitresults_slope[l][21] + dEvsMIPs_fitresults_slope[l][22] + dEvsMIPs_fitresults_slope[l][23] + dEvsMIPs_fitresults_slope[l][24] + dEvsMIPs_fitresults_slope[l][25] + dEvsMIPs_fitresults_slope[l][26]+ dEvsMIPs_fitresults_slope[l][27]+ dEvsMIPs_fitresults_slope[l][28]+ dEvsMIPs_fitresults_slope[l][29];
    dEvsMIPs_sec1 = dEvsMIPs_sec1/6.;
    dEvsMIPs_sec2 = dEvsMIPs_sec1/6.;
    dEvsMIPs_sec3 = dEvsMIPs_sec1/8.;
    dEvsMIPs_sec4 = dEvsMIPs_sec1/10.;
    for(unsigned iL(0); iL<(*ssvec).size(); iL++){
      if ( (iL!=(int)(startlayer-1)) ){
	if ( ((int) iL)>=0 && iL <=5 ){
	  dEvsMIPs_fitresults_slope_4sections[l][(int) iL] = dEvsMIPs_sec1;
	} else if ( iL>=6 && iL <=11 ){
	  dEvsMIPs_fitresults_slope_4sections[l][(int) iL] = dEvsMIPs_sec2;
	} else if ( iL>=12 && iL <=19 ){
	  dEvsMIPs_fitresults_slope_4sections[l][(int) iL] = dEvsMIPs_sec3;
	} else if ( iL>=20 && iL <=29 ){
	  dEvsMIPs_fitresults_slope_4sections[l][(int) iL] = dEvsMIPs_sec4;
	}
      } //if for the upstream
    } //sampling sections
  } // energies


  //======================================================================================
  //For the dedx
  TString titleoffile = "data/V5_dedx.txt";
  std::ifstream infile(titleoffile);
  std::string dummyfordedx; 
  //dedx for all layers minus one for the upstream material or zero if no material upstream
  double dedx[(*ssvec).size()-startlayer];

  std::cout << "dedx for layers 1-30 from file " << titleoffile << std::endl;
  std::cout << "dedx   total_dedx    " << titleoffile << std::endl;
  double totdedx = 0.;
  double totdedx_standalone = 0.; 
  double dedxweight_sec1 = 0.;
  double dedxweight_sec2 = 0.;
  double dedxweight_sec3 = 0.;
  double dedxweight_totsf = 0.;


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

    unsigned int iL = k+startlayer;
    //----------------------------------------------------
    //This is for the weight based on dedx in the 3 sections case
    if ( (iL==0) || (iL==3) || (iL==5) || (iL==7) || (iL==9) ) {
      dedxweight_sec1 = dedxweight_sec1 + (*ssvec)[k+startlayer].voldEdx();
    }
    if ( (iL==1) || (iL==2) || (iL==4) || (iL==6) || (iL==8) || (iL==10) || (iL==11) || 
	 (iL==12) || (iL==13) || (iL==14) || (iL==15) || (iL==16) || (iL==17) || (iL==18) || 
	 (iL==19) || (iL==20) || (iL==21) || (iL==23) || (iL==25) || (iL==27) || (iL==29) ) 
      {
	dedxweight_sec2 = dedxweight_sec2 + (*ssvec)[k+startlayer].voldEdx();
      }
    if ( (iL==22) || (iL==24) || (iL==26) || (iL==28) ) {
      dedxweight_sec3 = dedxweight_sec3 + (*ssvec)[k+startlayer].voldEdx();
    }
    //----------------------------------------------------
    //This is for the weight based on dedx in total sf case
    dedxweight_totsf = dedxweight_totsf + (*ssvec)[k+startlayer].voldEdx();
  }
  dedxweight_sec1 = dedxweight_sec1 / 5.;
  dedxweight_sec2 = dedxweight_sec2 / 21.;
  dedxweight_sec3 = dedxweight_sec3 / 4.;
  dedxweight_totsf = dedxweight_totsf / 30.;

  infile.close();

  //dedx for the silicon per layer
  double si_dedx = 0.11628; //it is 0.3 mm * 0.3876 MeV/mm
  //total dedx for the silicon in whole calorimeter
  // double totalsi_dedx = 3.4884; //30 layers x si_dedx
  //dedx for the absorber
  double abs_dedx = 0.;
  //Do not include the volume for the leakage
  for(unsigned iL(startlayer); iL<(*ssvec).size()-startlayer; iL++){
    abs_dedx = abs_dedx + (*ssvec)[iL].voldEdx();
  }

  //So dedx for active and absorber
  std::cout << "dedx silicon per layer : " << si_dedx << std::endl;
  std::cout << "dedx absorber total : " << abs_dedx << std::endl;

  std::vector<double> xotrans;
  double currentxotrans = 0.;
  //======================================================================================
  //The files that we will store the results of this analysis
  TString res[numberofconfigs];
  for (int k=0; k<numberofconfigs; k++){
    TString fp = "SiWEcal_" + IntToString(configs[k]);	
    res[k] = fp + "_final.root";
    //std::cout << res[k] << std::endl;
  }
  
  TFile* results[numberofconfigs];

  //Avoid rerun the entries
  if(step1){
    
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
	  if (ievt%(10000)==0) std::cout << "Entry " << ievt << std::endl;
	  // if (ievt==100){break;}
	  // std::cout << "Entry " << ievt << std::endl;
	  double recE = 0.;
	  double recE_leak = 0.;
	  double recE_inGeV = 0.;
	  double recE_leak_inGeV = 0.;
	  double recE_leak_offset_inGeV = 0.;
	  double recE_4sections_inGeV = 0.;
	  double recE_4sections_leak_inGeV = 0.;
	  double recE_4sections_ch1_inGeV = 0.;
	  double recE_4sections_ch1_leak_inGeV = 0.;
	  double recE_4sections_ch2_inGeV = 0.;
	  double recE_4sections_ch2_leak_inGeV = 0.;
	  double recE_4sections_ch3_inGeV = 0.;
	  double recE_4sections_ch3_leak_inGeV = 0.;
	  double recE_4sections_ch4_inGeV = 0.;
	  double recE_4sections_ch4_leak_inGeV = 0.;
	  double recE_ch1_inGeV = 0.;
	  double recE_ch1_leak_inGeV = 0.;
	  double recE_ch2_inGeV = 0.;
	  double recE_ch2_leak_inGeV = 0.;
	  double recE_ch3_inGeV = 0.;
	  double recE_ch3_leak_inGeV = 0.;
	  double recE_ch4_inGeV = 0.;
	  double recE_ch4_leak_inGeV = 0.;
	  double recE_raw = 0.;
	  double recE_raw_leak = 0.;
	  double recE_raw_inGeV = 0.;
	  double recE_raw_leak_inGeV = 0.;
	  double recE_val = 0.;
	  double recE_totsf_inGeV = 0.;
	  double recE_totsf_leak_inGeV = 0.;
	  double recE_totsf_dedxweighted_inGeV = 0.;
	  double recE_totsf_dedxweighted_leak_inGeV = 0.;
	  double recE_totsf_sec_inGeV = 0.;
	  double recE_totsf_sec_leak_inGeV = 0.;
	  double recE_totsf_sec_dedxweighted_inGeV = 0.;
	  double recE_totsf_sec_dedxweighted_leak_inGeV = 0.;
	  double recE_totsf_dedxsec_inGeV = 0.;
	  double recE_totsf_dedxsec_leak_inGeV = 0.;
	  // double recE_sf = 0.;
	  double recE_corr = 0.;
	  double lossinsensors = 0.;
	  double lossinfirst3sensors = 0.;
	  double numofmipsinabs = 0.;
	  //For the xotransfered add the before calo xo's
	  xotrans.clear();
	  currentxotrans = 0.;
	  int layermax = -1.;
	  double currmaxenergy = 0.;
	  //For the sum(Eixoi) in the shower depth
	  double Xoweightedwithenergy = 0.;

	  //Loop on sampling sections
	  for(unsigned iL(0); iL<(*ssvec).size(); iL++){
	    // std::cout << "Sampling Section " << iL << " energyabsorber " << energyabsorber  << std::endl; 
	  
	    //RecoE
	    //!!!! Change this when you do not have the leakage layer
	    if ( (iL!=(int)(startlayer-1)) ){
	      //dEvsMIPs_fitresults_slope[l][(int) iL]: GeV/GeV
	      //So, (MeV/1000.) * (1./GeV/GeV) 
	      recE_raw = recE_raw + ((*ssvec)[iL].measuredE()/(0.081));//in MIPs
	      recE_raw_leak = recE_raw_leak + ((*ssvec)[iL].leakageE()/(0.081)) + ((*ssvec)[iL].measuredE()/(0.081));//in MIPs

	      recE_raw_inGeV = recE_raw_inGeV + (((*ssvec)[iL].measuredE())/1000.); // in GeV
	      recE_raw_leak_inGeV = recE_raw_leak_inGeV + ((*ssvec)[iL].leakageE()/(1000.)) + (((*ssvec)[iL].measuredE())/1000.); // in GeV

	      if (((int)iL)==startlayer){
		numofmipsinabs =  (*ssvec)[iL].measuredE()/(0.081);
	      } else {
		numofmipsinabs =  ( ((*ssvec)[iL].measuredE()/(0.081)) + ( (*ssvec)[iL-1].measuredE()/(0.081) ) )/2.;
	      }
	      recE_val = recE_val + ( ((*ssvec)[iL].measuredE()/(0.081) )*(si_dedx)) + (numofmipsinabs*(*ssvec)[iL].voldEdx()) + (*ssvec)[iL].leakageE();
	      // recE_inGeV = recE_inGeV + ((*ssvec)[iL].measuredE()/0.081) * (dEvsMIPs_fitresults_slope[l][(int) iL]);
	      //Ereco_i = Esilicon_i*(1/slope_i)
	      // recE_inGeV = recE_inGeV + ( ((*ssvec)[iL].measuredE()/0.081) * (dEvsMIPs_fitresults_slope[l][(int) iL]) ) + ((*ssvec)[iL].measuredE()/1000.);

	      recE_inGeV = recE_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./dEvsMIPs_fitresults_slope[l][(int) iL]) )/(1000.) ); // in GeV
	      //Adding back leakage
	      recE_leak_inGeV = recE_leak_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./dEvsMIPs_fitresults_slope[l][(int) iL]) )/(1000.) ) + ((*ssvec)[iL].leakageE()/(1000.)); // in GeV
	      //Adding back leakage plus offset
	      recE_leak_offset_inGeV = recE_leak_offset_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./dEvsMIPs_fitresults_slope[l][(int) iL]) )/(1000.) ) + ((*ssvec)[iL].leakageE()/(1000.)) + dEvsMIPs_fitresults_const[l][(int) iL]; // in GeV

	      //4 sections
	      recE_4sections_inGeV = recE_4sections_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./dEvsMIPs_fitresults_slope_4sections[l][(int) iL]) )/(1000.) ); // in GeV
	      //Adding back leakage
	      recE_4sections_leak_inGeV = recE_4sections_leak_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./dEvsMIPs_fitresults_slope_4sections[l][(int) iL]) )/(1000.) ) + ((*ssvec)[iL].leakageE()/(1000.)); // in GeV


	      if ( ((int) iL)>=0 && iL <=5 ){
		recE_ch1_inGeV = recE_ch1_inGeV + ( ( (*ssvec)[iL].measuredE() + (*ssvec)[iL].absorberE() )/(1000.) ); // in GeV
		recE_ch1_leak_inGeV = recE_ch1_leak_inGeV + ( ( (*ssvec)[iL].measuredE() + (*ssvec)[iL].absorberE() )/(1000.) ) + ((*ssvec)[iL].leakageE()/(1000.)); // in GeV
		//4 sections: Cheat 1
		recE_4sections_ch1_inGeV = recE_4sections_ch1_inGeV + ( ( (*ssvec)[iL].measuredE() + (*ssvec)[iL].absorberE() )/(1000.) ); // in GeV
		//Adding back leakage: Cheat 1
		recE_4sections_ch1_leak_inGeV = recE_4sections_ch1_leak_inGeV + ( ( (*ssvec)[iL].measuredE() + (*ssvec)[iL].absorberE() )/(1000.) ) + ((*ssvec)[iL].leakageE()/(1000.)); // in GeV
	      } else {
		recE_ch1_inGeV = recE_ch1_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./dEvsMIPs_fitresults_slope[l][(int) iL]) )/(1000.) ); // in GeV
		recE_ch1_leak_inGeV = recE_ch1_leak_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./dEvsMIPs_fitresults_slope[l][(int) iL]) )/(1000.) ) + ((*ssvec)[iL].leakageE()/(1000.)); // in GeV
		//4 sections: Cheat 1
		recE_4sections_ch1_inGeV = recE_4sections_ch1_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./dEvsMIPs_fitresults_slope_4sections[l][(int) iL]) )/(1000.) ); // in GeV
		//Adding back leakage: Cheat 1
		recE_4sections_ch1_leak_inGeV = recE_4sections_ch1_leak_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./dEvsMIPs_fitresults_slope_4sections[l][(int) iL]) )/(1000.) ) + ((*ssvec)[iL].leakageE()/(1000.)); // in GeV
	      }   
	      
	      if ( ((int) iL)>=6 && iL <=11 ){
		recE_ch2_inGeV = recE_ch2_inGeV + ( ( (*ssvec)[iL].measuredE() + (*ssvec)[iL].absorberE() )/(1000.) ); // in GeV
		recE_ch2_leak_inGeV = recE_ch2_leak_inGeV + ( ( (*ssvec)[iL].measuredE() + (*ssvec)[iL].absorberE() )/(1000.) ) + ((*ssvec)[iL].leakageE()/(1000.)); // in GeV
		//4 sections: Cheat 2
		recE_4sections_ch2_inGeV = recE_4sections_ch2_inGeV + ( ( (*ssvec)[iL].measuredE() + (*ssvec)[iL].absorberE() )/(1000.) ); // in GeV
		//Adding back leakage: Cheat 2
		recE_4sections_ch2_leak_inGeV = recE_4sections_ch2_leak_inGeV + ( ( (*ssvec)[iL].measuredE() + (*ssvec)[iL].absorberE() )/(1000.) ) + ((*ssvec)[iL].leakageE()/(1000.)); // in GeV
	      } else {
		recE_ch2_inGeV = recE_ch2_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./dEvsMIPs_fitresults_slope[l][(int) iL]) )/(1000.) ); // in GeV
		recE_ch2_leak_inGeV = recE_ch2_leak_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./dEvsMIPs_fitresults_slope[l][(int) iL]) )/(1000.) ) + ((*ssvec)[iL].leakageE()/(1000.)); // in GeV
		//4 sections: Cheat 2
		recE_4sections_ch2_inGeV = recE_4sections_ch2_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./dEvsMIPs_fitresults_slope_4sections[l][(int) iL]) )/(1000.) ); // in GeV
		//Adding back leakage: Cheat 2
		recE_4sections_ch2_leak_inGeV = recE_4sections_ch2_leak_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./dEvsMIPs_fitresults_slope_4sections[l][(int) iL]) )/(1000.) ) + ((*ssvec)[iL].leakageE()/(1000.)); // in GeV
	      }

	      if ( ((int) iL)>=12 && iL <=19 ){
		recE_ch3_inGeV = recE_ch3_inGeV + ( ( (*ssvec)[iL].measuredE() + (*ssvec)[iL].absorberE() )/(1000.) ); // in GeV
		recE_ch3_leak_inGeV = recE_ch3_leak_inGeV + ( ( (*ssvec)[iL].measuredE() + (*ssvec)[iL].absorberE() )/(1000.) ) + ((*ssvec)[iL].leakageE()/(1000.)); // in GeV
		//4 sections: Cheat 3
		recE_4sections_ch3_inGeV = recE_4sections_ch3_inGeV + ( ( (*ssvec)[iL].measuredE() + (*ssvec)[iL].absorberE() )/(1000.) ); // in GeV
		//Adding back leakage: Cheat 3
		recE_4sections_ch3_leak_inGeV = recE_4sections_ch3_leak_inGeV + ( ( (*ssvec)[iL].measuredE() + (*ssvec)[iL].absorberE() )/(1000.) ) + ((*ssvec)[iL].leakageE()/(1000.)); // in GeV
	      } else {
		recE_ch3_inGeV = recE_ch3_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./dEvsMIPs_fitresults_slope[l][(int) iL]) )/(1000.) ); // in GeV
		recE_ch3_leak_inGeV = recE_ch3_leak_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./dEvsMIPs_fitresults_slope[l][(int) iL]) )/(1000.) ) + ((*ssvec)[iL].leakageE()/(1000.)); // in GeV
		//4 sections: Cheat 3
		recE_4sections_ch3_inGeV = recE_4sections_ch3_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./dEvsMIPs_fitresults_slope_4sections[l][(int) iL]) )/(1000.) ); // in GeV
		//Adding back leakage: Cheat 3
		recE_4sections_ch3_leak_inGeV = recE_4sections_ch3_leak_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./dEvsMIPs_fitresults_slope_4sections[l][(int) iL]) )/(1000.) ) + ((*ssvec)[iL].leakageE()/(1000.)); // in GeV
	      }

	      if ( ((int) iL)>=20 && iL <=29 ){
		recE_ch4_inGeV = recE_ch4_inGeV + ( ( (*ssvec)[iL].measuredE() + (*ssvec)[iL].absorberE() )/(1000.) ); // in GeV
		recE_ch4_leak_inGeV = recE_ch4_leak_inGeV + ( ( (*ssvec)[iL].measuredE() + (*ssvec)[iL].absorberE() )/(1000.) ) + ((*ssvec)[iL].leakageE()/(1000.)); // in GeV
		//4 sections: Cheat 4
		recE_4sections_ch4_inGeV = recE_4sections_ch4_inGeV + ( ( (*ssvec)[iL].measuredE() + (*ssvec)[iL].absorberE() )/(1000.) ); // in GeV
		//Adding back leakage: Cheat 4
		recE_4sections_ch4_leak_inGeV = recE_4sections_ch4_leak_inGeV + ( ( (*ssvec)[iL].measuredE() + (*ssvec)[iL].absorberE() )/(1000.) ) + ((*ssvec)[iL].leakageE()/(1000.)); // in GeV
	      } else {
		recE_ch4_inGeV = recE_ch4_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./dEvsMIPs_fitresults_slope[l][(int) iL]) )/(1000.) ); // in GeV
		recE_ch4_leak_inGeV = recE_ch4_leak_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./dEvsMIPs_fitresults_slope[l][(int) iL]) )/(1000.) ) + ((*ssvec)[iL].leakageE()/(1000.)); // in GeV
		//4 sections: Cheat 4
		recE_4sections_ch4_inGeV = recE_4sections_ch4_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./dEvsMIPs_fitresults_slope_4sections[l][(int) iL]) )/(1000.) ); // in GeV
		//Adding back leakage: Cheat 4
		recE_4sections_ch4_leak_inGeV = recE_4sections_ch4_leak_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./dEvsMIPs_fitresults_slope_4sections[l][(int) iL]) )/(1000.) ) + ((*ssvec)[iL].leakageE()/(1000.)); // in GeV
	      }



	      //With the total sampling fraction 
	      recE_totsf_inGeV = recE_totsf_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./totalsfperenergy[l]) )/(1000.) ); // in GeV
	      //Adding back leakage 
	      recE_totsf_leak_inGeV = recE_totsf_leak_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./totalsfperenergy[l]) )/(1000.) ) + ((*ssvec)[iL].leakageE()/(1000.)); // in GeV

	      //With the total sampling fraction 
	      recE_totsf_dedxweighted_inGeV = recE_totsf_dedxweighted_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./totalsfperenergy[l]) * ((*ssvec)[iL].voldEdx()/dedxweight_totsf) )/(1000.) ); // in GeV
	      //Adding back leakage 
	      recE_totsf_dedxweighted_leak_inGeV = recE_totsf_dedxweighted_leak_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./totalsfperenergy[l]) * ((*ssvec)[iL].voldEdx()/dedxweight_totsf) )/(1000.) ) + ((*ssvec)[iL].leakageE()/(1000.)); // in GeV

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
		//With the total sampling fraction 
		recE_totsf_sec_inGeV = recE_totsf_sec_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./totalsfperenergy_sec1[l]) )/(1000.) ); // in GeV
		//Adding back leakage 
		recE_totsf_sec_leak_inGeV = recE_totsf_sec_leak_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./totalsfperenergy_sec1[l]) )/(1000.) ) + ((*ssvec)[iL].leakageE()/(1000.)); // in GeV
	      }
	      if (totsf_sec2){
		//With the total sampling fraction 
		recE_totsf_sec_inGeV = recE_totsf_sec_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./totalsfperenergy_sec2[l]) )/(1000.) ); // in GeV
		//Adding back leakage 
		recE_totsf_sec_leak_inGeV = recE_totsf_sec_leak_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./totalsfperenergy_sec2[l]) )/(1000.) ) + ((*ssvec)[iL].leakageE()/(1000.)); // in GeV
	      }
	      if (totsf_sec3){
		//With the total sampling fraction 
		recE_totsf_sec_inGeV = recE_totsf_sec_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./totalsfperenergy_sec3[l]) )/(1000.) ); // in GeV
		//Adding back leakage 
		recE_totsf_sec_leak_inGeV = recE_totsf_sec_leak_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./totalsfperenergy_sec3[l]) )/(1000.) ) + ((*ssvec)[iL].leakageE()/(1000.)); // in GeV
	      }

	      //----------------------------------------------------------------------------------
	      //For the total sf per section weighted with 
	      if (totsf_sec1){
		//With the total sampling fraction 
		recE_totsf_sec_dedxweighted_inGeV = recE_totsf_sec_dedxweighted_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./totalsfperenergy_sec1[l]) *    ( (*ssvec)[iL].voldEdx() / dedxweight_sec1 )  )/(1000.) ); // in GeV
		//Adding back leakage 
		recE_totsf_sec_dedxweighted_leak_inGeV = recE_totsf_sec_dedxweighted_leak_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./totalsfperenergy_sec1[l]) * ((*ssvec)[iL].voldEdx()/dedxweight_sec1) )/(1000.) ) + ((*ssvec)[iL].leakageE()/(1000.)); // in GeV
	      }
	      if (totsf_sec2){
		//With the total sampling fraction 
		recE_totsf_sec_dedxweighted_inGeV = recE_totsf_sec_dedxweighted_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./totalsfperenergy_sec2[l]) * ((*ssvec)[iL].voldEdx()/dedxweight_sec2) )/(1000.) ); // in GeV
		//Adding back leakage 
		recE_totsf_sec_dedxweighted_leak_inGeV = recE_totsf_sec_dedxweighted_leak_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./totalsfperenergy_sec2[l]) * ((*ssvec)[iL].voldEdx()/dedxweight_sec2) )/(1000.) ) + ((*ssvec)[iL].leakageE()/(1000.)); // in GeV
	      }
	      if (totsf_sec3){
		//With the total sampling fraction 
		recE_totsf_sec_dedxweighted_inGeV = recE_totsf_sec_dedxweighted_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./totalsfperenergy_sec3[l]) * ((*ssvec)[iL].voldEdx()/dedxweight_sec3) )/(1000.) ); // in GeV
		//Adding back leakage 
		recE_totsf_sec_dedxweighted_leak_inGeV = recE_totsf_sec_dedxweighted_leak_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./totalsfperenergy_sec3[l]) * ((*ssvec)[iL].voldEdx()/dedxweight_sec3) )/(1000.) ) + ((*ssvec)[iL].leakageE()/(1000.)); // in GeV
	      }
	      //----------------------------------------------------------------------------------
	      //5 plus 1 sections
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
		//With the total sampling fraction 
		recE_totsf_dedxsec_inGeV = recE_totsf_dedxsec_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./totalsfperenergy_dedxsec1[l]) )/(1000.) ); // in GeV
		//Adding back leakage 
		recE_totsf_dedxsec_leak_inGeV = recE_totsf_dedxsec_leak_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./totalsfperenergy_dedxsec1[l]) )/(1000.) ) + ((*ssvec)[iL].leakageE()/(1000.)); // in GeV
	      }
	      if (totsf_dedxsec2){
		//With the total sampling fraction 
		recE_totsf_dedxsec_inGeV = recE_totsf_dedxsec_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./totalsfperenergy_dedxsec2[l]) )/(1000.) ); // in GeV
		//Adding back leakage 
		recE_totsf_dedxsec_leak_inGeV = recE_totsf_dedxsec_leak_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./totalsfperenergy_dedxsec2[l]) )/(1000.) ) + ((*ssvec)[iL].leakageE()/(1000.)); // in GeV
	      }
	      if (totsf_dedxsec3){
		//With the total sampling fraction 
		recE_totsf_dedxsec_inGeV = recE_totsf_dedxsec_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./totalsfperenergy_dedxsec3[l]) )/(1000.) ); // in GeV
		//Adding back leakage 
		recE_totsf_dedxsec_leak_inGeV = recE_totsf_dedxsec_leak_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./totalsfperenergy_dedxsec3[l]) )/(1000.) ) + ((*ssvec)[iL].leakageE()/(1000.)); // in GeV
	      }
	      if (totsf_dedxsec4){
		//With the total sampling fraction 
		recE_totsf_dedxsec_inGeV = recE_totsf_dedxsec_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./totalsfperenergy_dedxsec4[l]) )/(1000.) ); // in GeV
		//Adding back leakage 
		recE_totsf_dedxsec_leak_inGeV = recE_totsf_dedxsec_leak_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./totalsfperenergy_dedxsec4[l]) )/(1000.) ) + ((*ssvec)[iL].leakageE()/(1000.)); // in GeV
	      }
	      if (totsf_dedxsec5){
		//With the total sampling fraction 
		recE_totsf_dedxsec_inGeV = recE_totsf_dedxsec_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./totalsfperenergy_dedxsec5[l]) )/(1000.) ); // in GeV
		//Adding back leakage 
		recE_totsf_dedxsec_leak_inGeV = recE_totsf_dedxsec_leak_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./totalsfperenergy_dedxsec5[l]) )/(1000.) ) + ((*ssvec)[iL].leakageE()/(1000.)); // in GeV
	      }
	      if (totsf_dedxsec6){
		//With the total sampling fraction 
		recE_totsf_dedxsec_inGeV = recE_totsf_dedxsec_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./totalsfperenergy_dedxsec6[l]) )/(1000.) ); // in GeV
		//Adding back leakage 
		recE_totsf_dedxsec_leak_inGeV = recE_totsf_dedxsec_leak_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./totalsfperenergy_dedxsec6[l]) )/(1000.) ) + ((*ssvec)[iL].leakageE()/(1000.)); // in GeV
	      }


	      recE = recE + ( ( ( (*ssvec)[iL].measuredE() ) * (1./dEvsMIPs_fitresults_slope[l][(int) iL]) )/(0.081) ); //in MIPs
	      //Adding back leakage
	      recE_leak = recE_leak + ( ( ( (*ssvec)[iL].measuredE() ) * (1./dEvsMIPs_fitresults_slope[l][(int) iL]) )/(0.081) )  + ((*ssvec)[iL].leakageE()/(0.081)); //in MIPs
	      lossinsensors = lossinsensors + (*ssvec)[iL].measuredE(); //in MeV
	    }
	    if ( ((int)iL)==startlayer || ((int)iL)==(startlayer+1) || ((int)iL)==(startlayer+2) ){  
	      lossinfirst3sensors = lossinfirst3sensors + (*ssvec)[iL].measuredE(); //in MeV
	    }

	    //This is for the normalized shower depth method. 
	    //We have to know the shower max of the event 
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

	    //Below the -1 is in order to fill the histos at the last sampling section of the event
	    if (iL== ((*ssvec).size()-1) ){
	      // h_recoE[l]->Fill( (recE_inGeV*1000.)/(0.081) );//in MIPs
	      h_recoE[l]->Fill( recE_inGeV );//in GeV
	      h_recoE_leak[l]->Fill( recE_leak_inGeV );//in GeV
	      h_recoE_leak_offset[l]->Fill( recE_leak_offset_inGeV );//in GeV
	      h_recoE_raw[l]->Fill( recE_raw );//in MIPs
	      // h_recoE_raw_leak[l]->Fill( recE_raw_leak );//in MIPs
	      h_recoE_raw_leak[l]->Fill( recE_raw_leak_inGeV );//in GeV
	      //!!!!!! Recalib dedx 
	      h_recoE_val[l]->Fill( recalibdedx[l] * (recE_val/1000.) );//in GeV
	      h_recoE_totsf[l]->Fill( recE_totsf_inGeV );//in GeV
	      h_recoE_totsf_leak[l]->Fill( recE_totsf_leak_inGeV );//in GeV
	      h_recoE_totsf_dedxweighted[l]->Fill( recE_totsf_dedxweighted_inGeV );//in GeV
	      h_recoE_totsf_dedxweighted_leak[l]->Fill( recE_totsf_dedxweighted_leak_inGeV );//in GeV
	      h_recoE_totsf_sec[l]->Fill( recE_totsf_sec_inGeV );//in GeV
	      h_recoE_totsf_sec_leak[l]->Fill( recE_totsf_sec_leak_inGeV );//in GeV
	      h_recoE_totsf_sec_dedxweighted[l]->Fill( recE_totsf_sec_dedxweighted_inGeV );//in GeV
	      h_recoE_totsf_sec_dedxweighted_leak[l]->Fill( recE_totsf_sec_dedxweighted_leak_inGeV );//in GeV
	      h_recoE_totsf_dedxsec[l]->Fill( recE_totsf_dedxsec_inGeV );//in GeV
	      h_recoE_totsf_dedxsec_leak[l]->Fill( recE_totsf_dedxsec_leak_inGeV );//in GeV
	      h_recoE_4sections[l]->Fill( recE_4sections_inGeV );//in GeV
	      // std::cout << recE_4sections_inGeV << std::endl;
	      h_recoE_4sections_leak[l]->Fill( recE_4sections_leak_inGeV );//in GeV
	      h_recoE_4sections_ch1[l]->Fill( recE_4sections_ch1_inGeV );//in GeV
	      h_recoE_4sections_ch1_leak[l]->Fill( recE_4sections_ch1_leak_inGeV );//in GeV
	      h_recoE_4sections_ch2[l]->Fill( recE_4sections_ch2_inGeV );//in GeV
	      h_recoE_4sections_ch2_leak[l]->Fill( recE_4sections_ch2_leak_inGeV );//in GeV
	      h_recoE_4sections_ch3[l]->Fill( recE_4sections_ch3_inGeV );//in GeV
	      h_recoE_4sections_ch3_leak[l]->Fill( recE_4sections_ch3_leak_inGeV );//in GeV
	      h_recoE_4sections_ch4[l]->Fill( recE_4sections_ch4_inGeV );//in GeV
	      h_recoE_4sections_ch4_leak[l]->Fill( recE_4sections_ch4_leak_inGeV );//in GeV
	      h_recoE_ch1[l]->Fill( recE_ch1_inGeV );//in GeV
	      h_recoE_ch1_leak[l]->Fill( recE_ch1_leak_inGeV );//in GeV
	      h_recoE_ch2[l]->Fill( recE_ch2_inGeV );//in GeV
	      h_recoE_ch2_leak[l]->Fill( recE_ch2_leak_inGeV );//in GeV
	      h_recoE_ch3[l]->Fill( recE_ch3_inGeV );//in GeV
	      h_recoE_ch3_leak[l]->Fill( recE_ch3_leak_inGeV );//in GeV
	      h_recoE_ch4[l]->Fill( recE_ch4_inGeV );//in GeV
	      h_recoE_ch4_leak[l]->Fill( recE_ch4_leak_inGeV );//in GeV
	      //overtrueE
	      h_recoE_ch1_overtrueE[l]->Fill( (recE_ch1_inGeV/((double) particleenergies[l])) - 1. );//in GeV
	      h_recoE_ch1_overtrueE_leak[l]->Fill( (recE_ch1_leak_inGeV/((double) particleenergies[l])) - 1. );//in GeV
	      h_recoE_ch2_overtrueE[l]->Fill( (recE_ch2_inGeV/((double) particleenergies[l])) - 1. );//in GeV
	      h_recoE_ch2_overtrueE_leak[l]->Fill( (recE_ch2_leak_inGeV/((double) particleenergies[l])) - 1. );//in GeV
	      h_recoE_ch3_overtrueE[l]->Fill( (recE_ch3_inGeV/((double) particleenergies[l])) - 1. );//in GeV
	      h_recoE_ch3_overtrueE_leak[l]->Fill( (recE_ch3_leak_inGeV/((double) particleenergies[l])) - 1. );//in GeV
	      h_recoE_ch4_overtrueE[l]->Fill( (recE_ch4_inGeV/((double) particleenergies[l])) - 1. );//in GeV
	      h_recoE_ch4_overtrueE_leak[l]->Fill( (recE_ch4_leak_inGeV/((double) particleenergies[l])) - 1. );//in GeV
	      
	      //The correction for the upstream material
	      //Erec_corr = Erec + w0 + w1 * EmeasuredIn3layers
	      //dE_lossUpvsMIPs_fitresults_const[l] : MeV
	      //dE_lossUpvsMIPs_fitresults_slope[l] : MeV/MIP
	      //So, (MeV/0.081) + ( (1./MeV/MIP)* MeV 
	      // recE_corr = recE + (fabs(dE_lossUpvsMIPs_fitresults_const[l])/0.081) + ( fabs(1./dE_lossUpvsMIPs_fitresults_slope[l])*(lossinsensors) ) ;
	      // recE_corr = (recE_inGeV*1000.) + ((dE_lossUpvsMIPs_fitresults_const[l])) + ( (dE_lossUpvsMIPs_fitresults_slope[l])*(lossinsensors/0.081) ) ;
	      if (upstream){
		recE_corr = (recE_inGeV*1000.) + ((dE_lossUpvsMIPs_fitresults_const[l])) + ( (dE_lossUpvsMIPs_fitresults_slope[l])*(lossinfirst3sensors/0.081) ) ;
		// std::cout << "recE " << recE << " recE_corr " << recE_corr << std::endl;
		// h_recoE_corr[l]->Fill( recE_corr );//in MIPs
		h_recoE_corr[l]->Fill( recE_corr/1000. );//in GeV
	      }

	      //Reco energy over beam energy
	      h_recoEovertrueE[l]->Fill( (recE_inGeV/((double) particleenergies[l])) - 1. );//
	      h_recoEovertrueE_leak[l]->Fill( (recE_leak_inGeV/((double) particleenergies[l])) - 1. );//in GeV
	      h_recoEovertrueE_leak_offset[l]->Fill( (recE_leak_offset_inGeV/((double) particleenergies[l])) - 1. );//in GeV
	      h_recoEovertrueE_totsf[l]->Fill( (recE_totsf_inGeV/((double) particleenergies[l])) - 1. );//
	      h_recoEovertrueE_totsf_leak[l]->Fill( (recE_totsf_leak_inGeV/((double) particleenergies[l])) - 1. );//in GeV
	      h_recoEovertrueE_totsf_dedxweighted[l]->Fill( (recE_totsf_dedxweighted_inGeV/((double) particleenergies[l])) - 1. );//
	      h_recoEovertrueE_totsf_dedxweighted_leak[l]->Fill( (recE_totsf_dedxweighted_leak_inGeV/((double) particleenergies[l])) - 1. );//in GeV
	      h_recoEovertrueE_totsf_sec[l]->Fill( (recE_totsf_sec_inGeV/((double) particleenergies[l])) - 1. );//
	      h_recoEovertrueE_totsf_sec_leak[l]->Fill( (recE_totsf_sec_leak_inGeV/((double) particleenergies[l])) - 1. );//in GeV
	      h_recoEovertrueE_totsf_sec_dedxweighted[l]->Fill( (recE_totsf_sec_dedxweighted_inGeV/((double) particleenergies[l])) - 1. );//
	      h_recoEovertrueE_totsf_sec_dedxweighted_leak[l]->Fill( (recE_totsf_sec_dedxweighted_leak_inGeV/((double) particleenergies[l])) - 1. );//in GeV
	      h_recoEovertrueE_totsf_dedxsec[l]->Fill( (recE_totsf_dedxsec_inGeV/((double) particleenergies[l])) - 1. );//
	      h_recoEovertrueE_totsf_dedxsec_leak[l]->Fill( (recE_totsf_dedxsec_leak_inGeV/((double) particleenergies[l])) - 1. );//in GeV
	      if (upstream){
		h_recoE_corrovertrueE[l]->Fill( (recE_corr/1000.)/((double) particleenergies[l]) );//
	      }
	      h_recoE_rawEovertrueE[l]->Fill( ((recE_raw_inGeV)/((double) particleenergies[l])) - 1. );//
	      h_recoE_valEovertrueE[l]->Fill( ((recE_val/1000.)/((double) particleenergies[l])) - 1. );//


	    }
	
	  } //loop on sampling sections

	  double recE_shmax_inGeV = 0.;
	  double recE_shmax_leak_inGeV = 0.;
	  double recE_totsfshmax_sec_inGeV = 0.;
	  double recE_totsfshmax_sec_leak_inGeV = 0.;
	  double recE_totsfshmax_sec_dedxweighted_inGeV = 0.;
	  double recE_totsfshmax_sec_dedxweighted_leak_inGeV = 0.;
	  double recE_totsfshmax_dedxsec_inGeV = 0.;
	  double recE_totsfshmax_dedxsec_leak_inGeV = 0.;
	  //Now we know the x0max/<x0max> for the current event
	  double samfrac_shomax = xotrans[layermax-startlayer]/meanshowermaxperenergy[l];
	  //This is for the new shower depth variable of the current event defined above see bool
	  double newshowerdepth = Xoweightedwithenergy / lossinsensors;
	  //And of course for the new shower depth variable we should divide with the mean. 
	  double newshowerdepthweight = newshowerdepth / meannewshowerdepthperenergy[l]; 
	  if (usenewshowerdepth){
	    samfrac_shomax = newshowerdepthweight;
	  }

	  //Here fill the above sf for writing to root
	  h_Showermaxovermeanshmax[l]->Fill(samfrac_shomax);
	  double sfnormwithshower = 0.;
	  double totsfnormwithshower = 0.;
	  double totsf_sec1normwithshower = 0.;
	  double totsf_sec2normwithshower = 0.;
	  double totsf_sec3normwithshower = 0.;
	  double totsf_dedxsec1normwithshower = 0.;
	  double totsf_dedxsec2normwithshower = 0.;
	  double totsf_dedxsec3normwithshower = 0.;
	  double totsf_dedxsec4normwithshower = 0.;
	  double totsf_dedxsec5normwithshower = 0.;
	  double totsf_dedxsec6normwithshower = 0.;
	  double sfperlaynormwithshower = 0.;
	  double currxooverxomax = 0.;
	  double recE_sfnorm_inGeV = 0.; 
	  double recE_sfnorm_leak_inGeV = 0.; 
	  double recE_sfnorm_cut_inGeV = 0.; 
	  double recE_sfnorm_cut_leak_inGeV = 0.; 
	  double recE_totsfnorm_inGeV = 0.; 
	  double recE_totsfnorm_leak_inGeV = 0.; 
	  double recE_totsfnorm_dedxweighted_inGeV = 0.; 
	  double recE_totsfnorm_dedxweighted_leak_inGeV = 0.; 
	  double recE_totsfnorm_cut_inGeV = 0.; 
	  double recE_totsfnorm_cut_leak_inGeV = 0.; 
	  int binxooverxomax = 0;

	  for(unsigned int iL(startlayer); iL<(*ssvec).size(); iL++){
	    // currxooverxomax = xotrans[iL-startlayer]/xotrans[layermax-startlayer];
	    currxooverxomax = samfrac_shomax;
	    // std::cout <<  currxooverxomax << std::endl;
	    // std::cout <<  (int) ((currxooverxomax*10.) + 1.) << std::endl;
	    //The bin that the current X0/X0,max belongs
	    binxooverxomax = (int) ((currxooverxomax * (1/shdepthtobin)) + 1.);
	    sfnormwithshower = sfpernormshowerdepth[l][binxooverxomax];
	    totsfnormwithshower = totalsfpershmax[l][binxooverxomax];
	    totsf_sec1normwithshower = totalsf_sec1pershmax[l][binxooverxomax];
	    totsf_sec2normwithshower = totalsf_sec2pershmax[l][binxooverxomax];
	    totsf_sec3normwithshower = totalsf_sec3pershmax[l][binxooverxomax];
	    totsf_dedxsec1normwithshower = totalsf_dedxsec1pershmax[l][binxooverxomax];
	    totsf_dedxsec2normwithshower = totalsf_dedxsec2pershmax[l][binxooverxomax];
	    totsf_dedxsec3normwithshower = totalsf_dedxsec3pershmax[l][binxooverxomax];
	    totsf_dedxsec4normwithshower = totalsf_dedxsec4pershmax[l][binxooverxomax];
	    totsf_dedxsec5normwithshower = totalsf_dedxsec5pershmax[l][binxooverxomax];
	    totsf_dedxsec6normwithshower = totalsf_dedxsec6pershmax[l][binxooverxomax];
	    sfperlaynormwithshower = sfperlaypershmax[l][iL][binxooverxomax];
	    //Just checking
	    // std::cout << "Current norm depth " << currxooverxomax << " bin " << binxooverxomax << " sf value " << sfnormwithshower << std::endl;

	    //Total SF per section with shower method
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
	      //With the total sampling fraction 
	      recE_totsfshmax_sec_inGeV = recE_totsfshmax_sec_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./totsf_sec1normwithshower) )/(1000.) ); // in GeV
	      //Adding back leakage 
	      recE_totsfshmax_sec_leak_inGeV = recE_totsfshmax_sec_leak_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./totsf_sec1normwithshower) )/(1000.) ) + ((*ssvec)[iL].leakageE()/(1000.)); // in GeV
	    }
	    if (totsf_sec2){
	      //With the total sampling fraction 
	      recE_totsfshmax_sec_inGeV = recE_totsfshmax_sec_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./totsf_sec2normwithshower) )/(1000.) ); // in GeV
	      //Adding back leakage 
	      recE_totsfshmax_sec_leak_inGeV = recE_totsfshmax_sec_leak_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./totsf_sec2normwithshower) )/(1000.) ) + ((*ssvec)[iL].leakageE()/(1000.)); // in GeV
	    }
	    if (totsf_sec3){
	      //With the total sampling fraction 
	      recE_totsfshmax_sec_inGeV = recE_totsfshmax_sec_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./totsf_sec3normwithshower) )/(1000.) ); // in GeV
	      //Adding back leakage 
	      recE_totsfshmax_sec_leak_inGeV = recE_totsfshmax_sec_leak_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./totsf_sec3normwithshower) )/(1000.) ) + ((*ssvec)[iL].leakageE()/(1000.)); // in GeV
	    }

	    //-------------------------------------------------------------------------------------------
	    if (totsf_sec1){
	      //With the total sampling fraction 
	      recE_totsfshmax_sec_dedxweighted_inGeV = recE_totsfshmax_sec_dedxweighted_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./totsf_sec1normwithshower) * ((*ssvec)[iL].voldEdx()/dedxweight_sec1) )/(1000.) ); // in GeV
	      //Adding back leakage 
	      recE_totsfshmax_sec_dedxweighted_leak_inGeV = recE_totsfshmax_sec_dedxweighted_leak_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./totsf_sec1normwithshower) * ((*ssvec)[iL].voldEdx()/dedxweight_sec1) )/(1000.) ) + ((*ssvec)[iL].leakageE()/(1000.)); // in GeV
	    }
	    if (totsf_sec2){
	      //With the total sampling fraction 
	      recE_totsfshmax_sec_dedxweighted_inGeV = recE_totsfshmax_sec_dedxweighted_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./totsf_sec2normwithshower) * ((*ssvec)[iL].voldEdx()/dedxweight_sec2) )/(1000.) ); // in GeV
	      //Adding back leakage 
	      recE_totsfshmax_sec_dedxweighted_leak_inGeV = recE_totsfshmax_sec_dedxweighted_leak_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./totsf_sec2normwithshower) * ((*ssvec)[iL].voldEdx()/dedxweight_sec2) )/(1000.) ) + ((*ssvec)[iL].leakageE()/(1000.)); // in GeV
	    }
	    if (totsf_sec3){
	      //With the total sampling fraction 
	      recE_totsfshmax_sec_dedxweighted_inGeV = recE_totsfshmax_sec_dedxweighted_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./totsf_sec3normwithshower) * ((*ssvec)[iL].voldEdx()/dedxweight_sec3) )/(1000.) ); // in GeV
	      //Adding back leakage 
	      recE_totsfshmax_sec_dedxweighted_leak_inGeV = recE_totsfshmax_sec_dedxweighted_leak_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./totsf_sec3normwithshower) * ((*ssvec)[iL].voldEdx()/dedxweight_sec3) )/(1000.) ) + ((*ssvec)[iL].leakageE()/(1000.)); // in GeV
	    }



	    //-------------------------------------------------------------------------------------------
	    //For the 5+1 sections here
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
	    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	    //Be careful: First section is layer one which we use the fixed sampling fraction
	    if (totsf_dedxsec1){
	      //With the total sampling fraction 
	      // recE_totsfshmax_dedxsec_inGeV = recE_totsfshmax_dedxsec_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./totsf_dedxsec1normwithshower) )/(1000.) ); // in GeV
	      recE_totsfshmax_dedxsec_inGeV = recE_totsfshmax_dedxsec_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./totalsfperenergy_dedxsec1[l]) )/(1000.) ); // in GeV
	      //Adding back leakage 
	      recE_totsfshmax_dedxsec_leak_inGeV = recE_totsfshmax_dedxsec_leak_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./totalsfperenergy_dedxsec1[l]) )/(1000.) ) + ((*ssvec)[iL].leakageE()/(1000.)); // in GeV
	    }
	    if (totsf_dedxsec2){
	      //With the total sampling fraction 
	      recE_totsfshmax_dedxsec_inGeV = recE_totsfshmax_dedxsec_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./totsf_dedxsec2normwithshower) )/(1000.) ); // in GeV
	      //Adding back leakage 
	      recE_totsfshmax_dedxsec_leak_inGeV = recE_totsfshmax_dedxsec_leak_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./totsf_dedxsec2normwithshower) )/(1000.) ) + ((*ssvec)[iL].leakageE()/(1000.)); // in GeV
	    }
	    if (totsf_dedxsec3){
	      //With the total sampling fraction 
	      recE_totsfshmax_dedxsec_inGeV = recE_totsfshmax_dedxsec_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./totsf_dedxsec3normwithshower) )/(1000.) ); // in GeV
	      //Adding back leakage 
	      recE_totsfshmax_dedxsec_leak_inGeV = recE_totsfshmax_dedxsec_leak_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./totsf_dedxsec3normwithshower) )/(1000.) ) + ((*ssvec)[iL].leakageE()/(1000.)); // in GeV
	    }
	    if (totsf_dedxsec4){
	      //With the total sampling fraction 
	      recE_totsfshmax_dedxsec_inGeV = recE_totsfshmax_dedxsec_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./totsf_dedxsec4normwithshower) )/(1000.) ); // in GeV
	      //Adding back leakage 
	      recE_totsfshmax_dedxsec_leak_inGeV = recE_totsfshmax_dedxsec_leak_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./totsf_dedxsec4normwithshower) )/(1000.) ) + ((*ssvec)[iL].leakageE()/(1000.)); // in GeV
	    }
	    if (totsf_dedxsec5){
	      //With the total sampling fraction 
	      recE_totsfshmax_dedxsec_inGeV = recE_totsfshmax_dedxsec_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./totsf_dedxsec5normwithshower) )/(1000.) ); // in GeV
	      //Adding back leakage 
	      recE_totsfshmax_dedxsec_leak_inGeV = recE_totsfshmax_dedxsec_leak_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./totsf_dedxsec5normwithshower) )/(1000.) ) + ((*ssvec)[iL].leakageE()/(1000.)); // in GeV
	    }
	    if (totsf_dedxsec6){
	      //With the total sampling fraction 
	      recE_totsfshmax_dedxsec_inGeV = recE_totsfshmax_dedxsec_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./totsf_dedxsec6normwithshower) )/(1000.) ); // in GeV
	      //Adding back leakage 
	      recE_totsfshmax_dedxsec_leak_inGeV = recE_totsfshmax_dedxsec_leak_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./totsf_dedxsec6normwithshower) )/(1000.) ) + ((*ssvec)[iL].leakageE()/(1000.)); // in GeV
	    }

	    
	    recE_sfnorm_inGeV = recE_sfnorm_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./sfnormwithshower) )/(1000.) ); // in GeV
	    recE_sfnorm_leak_inGeV =  recE_sfnorm_leak_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./sfnormwithshower) )/(1000.) ) + ((*ssvec)[iL].leakageE()/(1000.)); // in GeV
	    recE_totsfnorm_inGeV = recE_totsfnorm_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./totsfnormwithshower) )/(1000.) ); // in GeV
	    recE_totsfnorm_leak_inGeV =  recE_totsfnorm_leak_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./totsfnormwithshower) )/(1000.) ) + ((*ssvec)[iL].leakageE()/(1000.)); // in GeV
	    recE_totsfnorm_dedxweighted_inGeV = recE_totsfnorm_dedxweighted_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./totsfnormwithshower) * ((*ssvec)[iL].voldEdx()/dedxweight_totsf) )/(1000.) ); // in GeV
	    recE_totsfnorm_dedxweighted_leak_inGeV =  recE_totsfnorm_dedxweighted_leak_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./totsfnormwithshower) * ((*ssvec)[iL].voldEdx()/dedxweight_totsf) )/(1000.) ) + ((*ssvec)[iL].leakageE()/(1000.)); // in GeV
	    //With the shower max method
	    recE_shmax_inGeV = recE_shmax_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./sfperlaynormwithshower) )/(1000.) ); // in GeV
	    //Adding back leakage
	    recE_shmax_leak_inGeV = recE_shmax_leak_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./sfperlaynormwithshower) )/(1000.) ) + ((*ssvec)[iL].leakageE()/(1000.)); // in GeV
	    //Same as above but now cutting at 2. xo/xomax. 
	    if (currxooverxomax >= 2.0){
	      // binxooverxomax = (int) ((2.0*10.) + 1.);
	      binxooverxomax = 21;
	      sfnormwithshower = sfpernormshowerdepth[l][binxooverxomax];
	      recE_sfnorm_cut_inGeV = recE_sfnorm_cut_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./sfnormwithshower) )/(1000.) ); // in GeV
	      recE_sfnorm_cut_leak_inGeV =  recE_sfnorm_cut_leak_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./sfnormwithshower) )/(1000.) ) + ((*ssvec)[iL].leakageE()/(1000.)); // in GeV
	      totsfnormwithshower = totalsfpershmax[l][binxooverxomax];
	      recE_totsfnorm_cut_inGeV = recE_totsfnorm_cut_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./totsfnormwithshower) )/(1000.) ); // in GeV
	      recE_totsfnorm_cut_leak_inGeV =  recE_totsfnorm_cut_leak_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./totsfnormwithshower) )/(1000.) ) + ((*ssvec)[iL].leakageE()/(1000.)); // in GeV
	    } else {
	      binxooverxomax = (int) ((currxooverxomax*10.) + 1.);
	      sfnormwithshower = sfpernormshowerdepth[l][binxooverxomax];
	      recE_sfnorm_cut_inGeV = recE_sfnorm_cut_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./sfnormwithshower) )/(1000.) ); // in GeV
	      recE_sfnorm_cut_leak_inGeV =  recE_sfnorm_cut_leak_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./sfnormwithshower) )/(1000.) ) + ((*ssvec)[iL].leakageE()/(1000.)); // in GeV
	      totsfnormwithshower = sfpernormshowerdepth[l][binxooverxomax];
	      recE_totsfnorm_cut_inGeV = recE_totsfnorm_cut_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./totsfnormwithshower) )/(1000.) ); // in GeV
	      recE_totsfnorm_cut_leak_inGeV =  recE_totsfnorm_cut_leak_inGeV + ( ( ( (*ssvec)[iL].measuredE() ) * (1./totsfnormwithshower) )/(1000.) ) + ((*ssvec)[iL].leakageE()/(1000.)); // in GeV
	    }

	    //Below the -1 is in order to fill the histos at the last sampling section of the event
	    if (iL== ((*ssvec).size()-1) ){
	      h_recoE_totsfshmax_sec[l]->Fill( recE_totsfshmax_sec_inGeV );//in GeV
	      h_recoE_totsfshmax_sec_leak[l]->Fill( recE_totsfshmax_sec_leak_inGeV );//in GeV
	      h_recoE_totsfshmax_sec_dedxweighted[l]->Fill( recE_totsfshmax_sec_dedxweighted_inGeV );//in GeV
	      h_recoE_totsfshmax_sec_dedxweighted_leak[l]->Fill( recE_totsfshmax_sec_dedxweighted_leak_inGeV );//in GeV
	      h_recoE_totsfshmax_dedxsec[l]->Fill( recE_totsfshmax_dedxsec_inGeV );//in GeV
	      h_recoE_totsfshmax_dedxsec_leak[l]->Fill( recE_totsfshmax_dedxsec_leak_inGeV );//in GeV
	      h_recoE_sfnorm[l]->Fill( recE_sfnorm_inGeV );//in GeV
	      h_recoE_sfnorm_leak[l]->Fill( recE_sfnorm_leak_inGeV );//in GeV
	      h_recoE_sfnorm_cut[l]->Fill( recE_sfnorm_cut_inGeV );//in GeV
	      h_recoE_sfnorm_cut_leak[l]->Fill( recE_sfnorm_cut_leak_inGeV );//in GeV
	      h_recoE_totsfnorm[l]->Fill( recE_totsfnorm_inGeV );//in GeV
	      h_recoE_totsfnorm_leak[l]->Fill( recE_totsfnorm_leak_inGeV );//in GeV
	      h_recoE_totsfnorm_dedxweighted[l]->Fill( recE_totsfnorm_dedxweighted_inGeV );//in GeV
	      h_recoE_totsfnorm_dedxweighted_leak[l]->Fill( recE_totsfnorm_dedxweighted_leak_inGeV );//in GeV
	      h_recoE_totsfnorm_cut[l]->Fill( recE_totsfnorm_cut_inGeV );//in GeV
	      h_recoE_totsfnorm_cut_leak[l]->Fill( recE_totsfnorm_cut_leak_inGeV );//in GeV
	      h_recoEovertrueE_sfnorm[l]->Fill( (recE_sfnorm_inGeV/((double) particleenergies[l])) - 1. );//in GeV
	      h_recoEovertrueE_sfnorm_leak[l]->Fill( (recE_sfnorm_leak_inGeV/((double) particleenergies[l])) - 1. );//in GeV
	      h_recoEovertrueE_sfnorm_cut[l]->Fill( (recE_sfnorm_cut_inGeV/((double) particleenergies[l])) - 1. );//in GeV
	      h_recoEovertrueE_sfnorm_cut_leak[l]->Fill( (recE_sfnorm_cut_leak_inGeV/((double) particleenergies[l])) - 1. );//in GeV
	      h_recoEovertrueE_totsfnorm[l]->Fill( (recE_totsfnorm_inGeV/((double) particleenergies[l])) - 1. );//
	      h_recoEovertrueE_totsfnorm_leak[l]->Fill( (recE_totsfnorm_leak_inGeV/((double) particleenergies[l])) - 1. );//in GeV
	      h_recoEovertrueE_totsfnorm_dedxweighted[l]->Fill( (recE_totsfnorm_dedxweighted_inGeV/((double) particleenergies[l])) - 1. );//
	      h_recoEovertrueE_totsfnorm_dedxweighted_leak[l]->Fill( (recE_totsfnorm_dedxweighted_leak_inGeV/((double) particleenergies[l])) - 1. );//in GeV
	      h_recoEovertrueE_totsfnorm_cut[l]->Fill( (recE_totsfnorm_cut_inGeV/((double) particleenergies[l])) - 1. );//
	      h_recoEovertrueE_totsfnorm_cut_leak[l]->Fill( (recE_totsfnorm_cut_leak_inGeV/((double) particleenergies[l])) - 1. );//in GeV
	      
	      //shower max method
	      h_recoE_shmax[l]->Fill( recE_shmax_inGeV );//in GeV
	      h_recoE_shmax_leak[l]->Fill( recE_shmax_leak_inGeV );//in GeV
	      h_recoEovertrueE_shmax[l]->Fill( (recE_shmax_inGeV/((double) particleenergies[l])) - 1. );//in GeV
	      h_recoEovertrueE_shmax_leak[l]->Fill( (recE_shmax_leak_inGeV/((double) particleenergies[l])) - 1. );//in GeV
	      h_recoEovertrueE_totsfshmax_sec[l]->Fill( (recE_totsfshmax_sec_inGeV/((double) particleenergies[l])) - 1. );//
	      h_recoEovertrueE_totsfshmax_sec_leak[l]->Fill( (recE_totsfshmax_sec_leak_inGeV/((double) particleenergies[l])) - 1. );//in GeV
	      h_recoEovertrueE_totsfshmax_sec_dedxweighted[l]->Fill( (recE_totsfshmax_sec_dedxweighted_inGeV/((double) particleenergies[l])) - 1. );//
	      h_recoEovertrueE_totsfshmax_sec_dedxweighted_leak[l]->Fill( (recE_totsfshmax_sec_dedxweighted_leak_inGeV/((double) particleenergies[l])) - 1. );//in GeV
	      h_recoEovertrueE_totsfshmax_dedxsec[l]->Fill( (recE_totsfshmax_dedxsec_inGeV/((double) particleenergies[l])) - 1. );//
	      h_recoEovertrueE_totsfshmax_dedxsec_leak[l]->Fill( (recE_totsfshmax_dedxsec_leak_inGeV/((double) particleenergies[l])) - 1. );//in GeV
	    }

	    
	  }
	  


	} //  Loop on entries

	
	//======================================================================================
	//Write histos 
	h_recoE_raw[l]->Write();
	h_recoE_raw_leak[l]->Write();
	h_recoE_val[l]->Write();
	h_recoE[l]->Write();
	h_recoE_leak[l]->Write();
	h_recoE_leak_offset[l]->Write();
	h_recoEovertrueE[l]->Write();
	h_recoEovertrueE_leak[l]->Write();
	h_recoEovertrueE_leak_offset[l]->Write();
	h_recoE_totsf[l]->Write();
	h_recoE_totsf_leak[l]->Write();
	h_recoE_totsf_dedxweighted[l]->Write();
	h_recoE_totsf_dedxweighted_leak[l]->Write();
	h_recoE_totsf_sec[l]->Write();
	h_recoE_totsf_sec_leak[l]->Write();
	h_recoE_totsf_sec_dedxweighted[l]->Write();
	h_recoE_totsf_sec_dedxweighted_leak[l]->Write();
	h_recoE_totsfshmax_sec[l]->Write();
	h_recoE_totsfshmax_sec_leak[l]->Write();
	h_recoE_totsfshmax_sec_dedxweighted[l]->Write();
	h_recoE_totsfshmax_sec_dedxweighted_leak[l]->Write();
	h_recoE_totsf_dedxsec[l]->Write();
	h_recoE_totsf_dedxsec_leak[l]->Write();
	h_recoE_totsfshmax_dedxsec[l]->Write();
	h_recoE_totsfshmax_dedxsec_leak[l]->Write();
	h_recoE_sfnorm[l]->Write();
	h_recoE_sfnorm_leak[l]->Write();
	h_recoE_sfnorm_cut[l]->Write();
	h_recoE_sfnorm_cut_leak[l]->Write();
	h_recoE_totsfnorm[l]->Write();
	h_recoE_totsfnorm_leak[l]->Write();
	h_recoE_totsfnorm_dedxweighted[l]->Write();
	h_recoE_totsfnorm_dedxweighted_leak[l]->Write();
	h_recoE_totsfnorm_cut[l]->Write();
	h_recoE_totsfnorm_cut_leak[l]->Write();
	h_recoEovertrueE_sfnorm[l]->Write();
	h_recoEovertrueE_sfnorm_leak[l]->Write();
	h_recoEovertrueE_sfnorm_cut[l]->Write();
	h_recoEovertrueE_sfnorm_cut_leak[l]->Write();
	h_recoEovertrueE_totsfnorm[l]->Write();
	h_recoEovertrueE_totsfnorm_leak[l]->Write();
	h_recoEovertrueE_totsfnorm_dedxweighted[l]->Write();
	h_recoEovertrueE_totsfnorm_dedxweighted_leak[l]->Write();
	h_recoEovertrueE_totsfnorm_cut[l]->Write();
	h_recoEovertrueE_totsfnorm_cut_leak[l]->Write();
	h_recoEovertrueE_totsf_leak[l]->Write();
	h_recoEovertrueE_totsf[l]->Write();
	h_recoEovertrueE_totsf_leak[l]->Write();
	h_recoEovertrueE_totsf_dedxweighted[l]->Write();
	h_recoEovertrueE_totsf_dedxweighted_leak[l]->Write();
	h_recoEovertrueE_totsf_sec[l]->Write();
	h_recoEovertrueE_totsf_sec_leak[l]->Write();
	h_recoEovertrueE_totsf_sec_dedxweighted[l]->Write();
	h_recoEovertrueE_totsf_sec_dedxweighted_leak[l]->Write();
	h_recoEovertrueE_totsfshmax_sec[l]->Write();
	h_recoEovertrueE_totsfshmax_sec_leak[l]->Write();
	h_recoEovertrueE_totsfshmax_sec_dedxweighted[l]->Write();
	h_recoEovertrueE_totsfshmax_sec_dedxweighted_leak[l]->Write();
	h_recoEovertrueE_totsf_dedxsec[l]->Write();
	h_recoEovertrueE_totsf_dedxsec_leak[l]->Write();
	h_recoEovertrueE_totsfshmax_dedxsec[l]->Write();
	h_recoEovertrueE_totsfshmax_dedxsec_leak[l]->Write();
	if (upstream){
	  h_recoE_corr[l]->Write();
	  h_recoE_corrovertrueE[l]->Write();
	}
	h_recoE_valEovertrueE[l]->Write();
	h_recoE_rawEovertrueE[l]->Write();
	h_recoE_4sections[l]->Write();
	h_recoE_4sections_leak[l]->Write();
	h_recoE_4sections_ch1[l]->Write();
	h_recoE_4sections_ch1_leak[l]->Write();
	h_recoE_4sections_ch2[l]->Write();
	h_recoE_4sections_ch2_leak[l]->Write();
	h_recoE_4sections_ch3[l]->Write();
	h_recoE_4sections_ch3_leak[l]->Write();
	h_recoE_4sections_ch4[l]->Write();
	h_recoE_4sections_ch4_leak[l]->Write();
	h_recoE_ch1[l]->Write();
	h_recoE_ch1_leak[l]->Write();
	h_recoE_ch2[l]->Write();
	h_recoE_ch2_leak[l]->Write();
	h_recoE_ch3[l]->Write();
	h_recoE_ch3_leak[l]->Write();
	h_recoE_ch4[l]->Write();
	h_recoE_ch4_leak[l]->Write();
	//overtrueE
	h_recoE_ch1_overtrueE[l]->Write();
	h_recoE_ch1_overtrueE_leak[l]->Write();
	h_recoE_ch2_overtrueE[l]->Write();
	h_recoE_ch2_overtrueE_leak[l]->Write();
	h_recoE_ch3_overtrueE[l]->Write();
	h_recoE_ch3_overtrueE_leak[l]->Write();
	h_recoE_ch4_overtrueE[l]->Write();
	h_recoE_ch4_overtrueE_leak[l]->Write();
	h_Showermaxovermeanshmax[l]->Write();
	h_recoE_shmax[l]->Write();
	h_recoE_shmax_leak[l]->Write();
	h_recoEovertrueE_shmax[l]->Write();
	h_recoEovertrueE_shmax_leak[l]->Write();

	//======================================================================================
	//We should here clear the histograms because we want them empty for the next file. 
	//Reset histos
	h_recoE_raw[l]->Reset();
	h_recoE_raw_leak[l]->Reset();
	h_recoE_val[l]->Reset();
	h_recoE[l]->Reset();
	h_recoE_leak[l]->Reset();
	h_recoE_leak_offset[l]->Reset();
	h_recoE_totsf[l]->Reset();
	h_recoE_totsf_leak[l]->Reset();
	h_recoE_totsf_dedxweighted[l]->Reset();
	h_recoE_totsf_dedxweighted_leak[l]->Reset();
	h_recoE_totsf_sec[l]->Reset();
	h_recoE_totsf_sec_leak[l]->Reset();
	h_recoE_totsf_sec_dedxweighted[l]->Reset();
	h_recoE_totsf_sec_dedxweighted_leak[l]->Reset();
	h_recoE_totsfshmax_sec[l]->Reset();
	h_recoE_totsfshmax_sec_leak[l]->Reset();
	h_recoE_totsfshmax_sec_dedxweighted[l]->Reset();
	h_recoE_totsfshmax_sec_dedxweighted_leak[l]->Reset();
	h_recoE_totsf_dedxsec[l]->Reset();
	h_recoE_totsf_dedxsec_leak[l]->Reset();
	h_recoE_totsfshmax_dedxsec[l]->Reset();
	h_recoE_totsfshmax_dedxsec_leak[l]->Reset();
	h_recoE_sfnorm_cut[l]->Reset();
	h_recoE_sfnorm_cut_leak[l]->Reset();
	h_recoE_sfnorm[l]->Reset();
	h_recoE_sfnorm_leak[l]->Reset();
	h_recoEovertrueE_sfnorm[l]->Reset();
	h_recoEovertrueE_sfnorm_leak[l]->Reset();
	h_recoEovertrueE_sfnorm_cut[l]->Reset();
	h_recoEovertrueE_sfnorm_cut_leak[l]->Reset();
	h_recoE_totsfnorm_cut[l]->Reset();
	h_recoE_totsfnorm_cut_leak[l]->Reset();
	h_recoE_totsfnorm[l]->Reset();
	h_recoE_totsfnorm_leak[l]->Reset();
	h_recoE_totsfnorm_dedxweighted[l]->Reset();
	h_recoE_totsfnorm_dedxweighted_leak[l]->Reset();
	h_recoEovertrueE_totsfnorm[l]->Reset();
	h_recoEovertrueE_totsfnorm_leak[l]->Reset();
	h_recoEovertrueE_totsfnorm_dedxweighted[l]->Reset();
	h_recoEovertrueE_totsfnorm_dedxweighted_leak[l]->Reset();
	h_recoEovertrueE_totsfnorm_cut[l]->Reset();
	h_recoEovertrueE_totsfnorm_cut_leak[l]->Reset();
	if (upstream){
	  h_recoE_corr[l]->Reset();
	  h_recoE_corrovertrueE[l]->Reset();
	}     
	h_recoEovertrueE[l]->Reset();
	h_recoEovertrueE_leak[l]->Reset();
	h_recoEovertrueE_leak_offset[l]->Reset();
	h_recoEovertrueE_totsf[l]->Reset();
	h_recoEovertrueE_totsf_leak[l]->Reset();
	h_recoEovertrueE_totsf_dedxweighted[l]->Reset();
	h_recoEovertrueE_totsf_dedxweighted_leak[l]->Reset();
	h_recoEovertrueE_totsf_sec[l]->Reset();
	h_recoEovertrueE_totsf_sec_leak[l]->Reset();
	h_recoEovertrueE_totsf_sec_dedxweighted[l]->Reset();
	h_recoEovertrueE_totsf_sec_dedxweighted_leak[l]->Reset();
	h_recoEovertrueE_totsfshmax_sec[l]->Reset();
	h_recoEovertrueE_totsfshmax_sec_leak[l]->Reset();
	h_recoEovertrueE_totsfshmax_sec_dedxweighted[l]->Reset();
	h_recoEovertrueE_totsfshmax_sec_dedxweighted_leak[l]->Reset();
	h_recoEovertrueE_totsf_dedxsec[l]->Reset();
	h_recoEovertrueE_totsf_dedxsec_leak[l]->Reset();
	h_recoEovertrueE_totsfshmax_dedxsec[l]->Reset();
	h_recoEovertrueE_totsfshmax_dedxsec_leak[l]->Reset();
	h_recoE_valEovertrueE[l]->Reset();
	h_recoE_rawEovertrueE[l]->Reset();
      	h_recoE_4sections[l]->Reset();
	h_recoE_4sections_leak[l]->Reset();
      	h_recoE_4sections_ch1[l]->Reset();
	h_recoE_4sections_ch1_leak[l]->Reset();
     	h_recoE_4sections_ch2[l]->Reset();
	h_recoE_4sections_ch2_leak[l]->Reset();
     	h_recoE_4sections_ch3[l]->Reset();
	h_recoE_4sections_ch3_leak[l]->Reset();
    	h_recoE_4sections_ch4[l]->Reset();
	h_recoE_4sections_ch4_leak[l]->Reset();
      	h_recoE_ch1[l]->Reset();
	h_recoE_ch1_leak[l]->Reset();
     	h_recoE_ch2[l]->Reset();
	h_recoE_ch2_leak[l]->Reset();
     	h_recoE_ch3[l]->Reset();
	h_recoE_ch3_leak[l]->Reset();
    	h_recoE_ch4[l]->Reset();
	h_recoE_ch4_leak[l]->Reset();
	//overtrueE
      	h_recoE_ch1_overtrueE[l]->Reset();
	h_recoE_ch1_overtrueE_leak[l]->Reset();
     	h_recoE_ch2_overtrueE[l]->Reset();
	h_recoE_ch2_overtrueE_leak[l]->Reset();
     	h_recoE_ch3_overtrueE[l]->Reset();
	h_recoE_ch3_overtrueE_leak[l]->Reset();
    	h_recoE_ch4_overtrueE[l]->Reset();
	h_recoE_ch4_overtrueE_leak[l]->Reset();
	h_Showermaxovermeanshmax[l]->Reset();
	h_recoE_shmax[l]->Reset();
	h_recoE_shmax_leak[l]->Reset();
	h_recoEovertrueE_shmax[l]->Reset();
	h_recoEovertrueE_shmax_leak[l]->Reset();

      }//Loop on energies
  
      results[k]->Close();

    } // Loop on files

  }//if to avoid rerun

  //======================================================================================
  //Make the following plots
  //mean(recoE) = f (beam) 
  //mean(recoE_corr) = f (beam)
  
  TCanvas *c1 = new TCanvas("c1", "  ");
  TCanvas *c2 = new TCanvas("c2", "  ");
  TCanvas *c3 = new TCanvas("c3", "  ");
  TCanvas *c4 = new TCanvas("c4", "  ");
  TCanvas *c5 = new TCanvas("c5", "  ");
  TCanvas *c6 = new TCanvas("c6", "  ");
  TCanvas *c7 = new TCanvas("c7", "  ");
  TCanvas *c8 = new TCanvas("c8", "  ");
  TCanvas *c9 = new TCanvas("c9", "  ");
  TCanvas *c10 = new TCanvas("c10", "  ");
  TCanvas *c11 = new TCanvas("c11", "  ");
  TCanvas *c12 = new TCanvas("c12", "  ");
  TCanvas *c13 = new TCanvas("c13", "  ");
  TCanvas *c14 = new TCanvas("c14", "  ");
  TCanvas *c15 = new TCanvas("c15", "  ");
  TCanvas *c16 = new TCanvas("c16", "  ");
  TCanvas *c17 = new TCanvas("c17", "  ");
  TCanvas *c18 = new TCanvas("c18", "  ");
  TCanvas *c19 = new TCanvas("c19", "  ");
  TCanvas *c20 = new TCanvas("c20", "  ");
  TCanvas *c21 = new TCanvas("c21", "  ");
  TCanvas *c22 = new TCanvas("c22", "  ");
  TCanvas *c23 = new TCanvas("c23", "  ");
  TCanvas *c24 = new TCanvas("c24", "  ");
  TCanvas *c25 = new TCanvas("c25", "  ");
  TCanvas *c26 = new TCanvas("c26", "  ");
  TCanvas *c27 = new TCanvas("c27", "  ");
  TCanvas *c28 = new TCanvas("c28", "  ");
  TCanvas *c29 = new TCanvas("c29", "  ");
  TCanvas *c30 = new TCanvas("c30", "  ");
  TCanvas *c31 = new TCanvas("c31", "  ");
  TCanvas *c32 = new TCanvas("c32", "  ");
  TCanvas *c33 = new TCanvas("c33", "  ");
  TCanvas *c34 = new TCanvas("c34", "  ");
  TCanvas *c35 = new TCanvas("c35", "  ");
  TCanvas *c36 = new TCanvas("c36", "  ");
  TCanvas *c37 = new TCanvas("c37", "  ");
  TCanvas *c38 = new TCanvas("c38", "  ");
  TCanvas *c39 = new TCanvas("c39", "  ");
  TCanvas *c40 = new TCanvas("c40", "  ");
  TCanvas *c41 = new TCanvas("c41", "  ");
  TCanvas *c42 = new TCanvas("c42", "  ");
  TCanvas *c43 = new TCanvas("c43", "  ");
  TCanvas *c44 = new TCanvas("c44", "  ");
  TCanvas *c45 = new TCanvas("c45", "  ");
  TCanvas *c46 = new TCanvas("c46", "  ");
  TCanvas *c47 = new TCanvas("c47", "  ");
  TCanvas *c48 = new TCanvas("c48", "  ");
  TCanvas *c49 = new TCanvas("c49", "  ");
  TCanvas *c50 = new TCanvas("c50", "  ");
  TCanvas *c51 = new TCanvas("c51", "  ");
  TCanvas *c52 = new TCanvas("c52", "  ");
  TCanvas *c53 = new TCanvas("c53", "  ");
  TCanvas *c54 = new TCanvas("c54", "  ");
  TCanvas *c55 = new TCanvas("c55", "  ");
  TCanvas *c56 = new TCanvas("c56", "  ");
  TCanvas *c57 = new TCanvas("c57", "  ");
  TCanvas *c58 = new TCanvas("c58", "  ");
  TCanvas *c59 = new TCanvas("c59", "  ");
  TCanvas *c60 = new TCanvas("c60", "  ");
  TCanvas *c61 = new TCanvas("c61", "  ");
  TCanvas *c62 = new TCanvas("c62", "  ");
  TCanvas *c63 = new TCanvas("c63", "  ");
  TCanvas *c64 = new TCanvas("c64", "  ");
  TCanvas *c65 = new TCanvas("c65", "  ");
  TCanvas *c66 = new TCanvas("c66", "  ");
  TCanvas *c67 = new TCanvas("c67", "  ");
  TCanvas *c68 = new TCanvas("c68", "  ");
  TCanvas *c69 = new TCanvas("c69", "  ");
  TCanvas *c70 = new TCanvas("c70", "  ");
  TCanvas *c71 = new TCanvas("c71", "  ");
  TCanvas *c72 = new TCanvas("c72", "  ");
  TCanvas *c73 = new TCanvas("c73", "  ");
  TCanvas *c74 = new TCanvas("c74", "  ");
  TCanvas *c75 = new TCanvas("c75", "  ");
  TCanvas *c76 = new TCanvas("c76", "  ");
  TCanvas *c77 = new TCanvas("c77", "  ");
  TCanvas *c78 = new TCanvas("c78", "  ");
  TCanvas *c79 = new TCanvas("c79", "  ");
  TCanvas *c80 = new TCanvas("c80", "  ");
  TCanvas *c81 = new TCanvas("c81", "  ");
  TCanvas *c82 = new TCanvas("c82", "  ");
  TCanvas *c83 = new TCanvas("c83", "  ");
  TCanvas *c84 = new TCanvas("c84", "  ");
  TCanvas *c85 = new TCanvas("c85", "  ");
  TCanvas *c86 = new TCanvas("c86", "  ");
  TCanvas *c87 = new TCanvas("c87", "  ");
  TCanvas *c88 = new TCanvas("c88", "  ");
  TCanvas *c89 = new TCanvas("c89", "  ");
  TCanvas *c90 = new TCanvas("c90", "  ");
  TCanvas *c91 = new TCanvas("c91", "  ");
  TCanvas *c92 = new TCanvas("c92", "  ");
  TCanvas *c93 = new TCanvas("c93", "  ");
  TCanvas *c94 = new TCanvas("c94", "  ");
  TCanvas *c95 = new TCanvas("c95", "  ");
  TCanvas *c96 = new TCanvas("c96", "  ");
  TCanvas *c97 = new TCanvas("c97", "  ");
  TCanvas *c98 = new TCanvas("c98", "  ");
  TCanvas *c99 = new TCanvas("c99", "  ");
  TCanvas *c100 = new TCanvas("c100", "  ");
  TCanvas *c101 = new TCanvas("c101", "  ");
  TCanvas *c102 = new TCanvas("c102", "  ");
  TCanvas *c103 = new TCanvas("c103", "  ");
  TCanvas *c104 = new TCanvas("c104", "  ");
  TCanvas *c105 = new TCanvas("c105", "  ");
  TCanvas *c106 = new TCanvas("c106", "  ");
  TCanvas *c107 = new TCanvas("c107", "  ");
  TCanvas *c108 = new TCanvas("c108", "  ");
  TCanvas *c109 = new TCanvas("c109", "  ");
  TCanvas *c110 = new TCanvas("c110", "  ");
  TCanvas *c111 = new TCanvas("c111", "  ");
  TCanvas *c112 = new TCanvas("c112", "  ");
  TCanvas *c113 = new TCanvas("c113", "  ");
  TCanvas *c114 = new TCanvas("c114", "  ");
  TCanvas *c115 = new TCanvas("c115", "  ");
  TCanvas *c116 = new TCanvas("c116", "  ");
  TCanvas *c117 = new TCanvas("c117", "  ");
  TCanvas *c118 = new TCanvas("c118", "  ");
  TCanvas *c119 = new TCanvas("c119", "  ");
  TCanvas *c120 = new TCanvas("c120", "  ");
  TCanvas *c121 = new TCanvas("c121", "  ");
  TCanvas *c122 = new TCanvas("c122", "  ");
  TCanvas *c123 = new TCanvas("c123", "  ");
  TCanvas *c124 = new TCanvas("c124", "  ");
  TCanvas *c125 = new TCanvas("c125", "  ");
  TCanvas *c126 = new TCanvas("c126", "  ");
  TCanvas *c127 = new TCanvas("c127", "  ");
  TCanvas *c128 = new TCanvas("c128", "  ");
  TCanvas *c129 = new TCanvas("c129", "  ");
  TCanvas *c130 = new TCanvas("c130", "  ");
  TCanvas *c131 = new TCanvas("c131", "  ");
  TCanvas *c132 = new TCanvas("c132", "  ");
  TCanvas *c133 = new TCanvas("c133", "  ");
  TCanvas *c134 = new TCanvas("c134", "  ");
  TCanvas *c135 = new TCanvas("c135", "  ");
  TCanvas *c136 = new TCanvas("c136", "  ");
  TCanvas *c137 = new TCanvas("c137", "  ");
  TCanvas *c138 = new TCanvas("c138", "  ");
  TCanvas *c139 = new TCanvas("c139", "  ");
  TCanvas *c140 = new TCanvas("c140", "  ");
  TCanvas *c141 = new TCanvas("c141", "  ");
  TCanvas *c142 = new TCanvas("c142", "  ");
  TCanvas *c143 = new TCanvas("c143", "  ");
  TCanvas *c144 = new TCanvas("c144", "  ");
  TCanvas *c145 = new TCanvas("c145", "  ");
  TCanvas *c146 = new TCanvas("c146", "  ");
  TCanvas *c147 = new TCanvas("c147", "  ");
  TCanvas *c148 = new TCanvas("c148", "  ");
  TCanvas *c149 = new TCanvas("c149", "  ");
  TCanvas *c150 = new TCanvas("c150", "  ");
  TCanvas *c151 = new TCanvas("c151", "  ");
  TCanvas *c152 = new TCanvas("c152", "  ");
  TCanvas *c153 = new TCanvas("c153", "  ");
  TCanvas *c154 = new TCanvas("c154", "  ");
  TCanvas *c155 = new TCanvas("c155", "  ");
  TCanvas *c156 = new TCanvas("c156", "  ");
  TCanvas *c157 = new TCanvas("c157", "  ");
  TCanvas *c158 = new TCanvas("c158", "  ");
  TCanvas *c159 = new TCanvas("c159", "  ");
  TCanvas *c160 = new TCanvas("c160", "  ");
  TCanvas *c161 = new TCanvas("c161", "  ");
  TCanvas *c162 = new TCanvas("c162", "  ");
  TCanvas *c163 = new TCanvas("c163", "  ");
  TCanvas *c164 = new TCanvas("c164", "  ");
  TCanvas *c165 = new TCanvas("c165", "  ");
  TCanvas *c166 = new TCanvas("c166", "  ");
  TCanvas *c167 = new TCanvas("c167", "  ");
  TCanvas *c168 = new TCanvas("c168", "  ");
  TCanvas *c169 = new TCanvas("c169", "  ");
  TCanvas *c170 = new TCanvas("c170", "  ");
  TCanvas *c171 = new TCanvas("c171", "  ");
  TCanvas *c172 = new TCanvas("c172", "  ");
  TCanvas *c173 = new TCanvas("c173", "  ");
  TCanvas *c174 = new TCanvas("c174", "  ");
  TCanvas *c175 = new TCanvas("c175", "  ");
  TCanvas *c176 = new TCanvas("c176", "  ");
  TCanvas *c177 = new TCanvas("c177", "  ");
  TCanvas *c178 = new TCanvas("c178", "  ");
  TCanvas *c179 = new TCanvas("c179", "  ");
  TCanvas *c180 = new TCanvas("c180", "  ");
  TCanvas *c181 = new TCanvas("c181", "  ");
  TCanvas *c182 = new TCanvas("c182", "  ");
  TCanvas *c183 = new TCanvas("c183", "  ");
  TCanvas *c184 = new TCanvas("c184", "  ");
  TCanvas *c185 = new TCanvas("c185", "  ");
  TCanvas *c186 = new TCanvas("c186", "  ");
  TCanvas *c187 = new TCanvas("c187", "  ");
  TCanvas *c188 = new TCanvas("c188", "  ");
  TCanvas *c189 = new TCanvas("c189", "  ");
  TCanvas *c190 = new TCanvas("c190", "  ");
  TCanvas *c191 = new TCanvas("c191", "  ");
  TCanvas *c192 = new TCanvas("c192", "  ");
  TCanvas *c193 = new TCanvas("c193", "  ");
  TCanvas *c194 = new TCanvas("c194", "  ");
  TCanvas *c195 = new TCanvas("c195", "  ");
  TCanvas *c196 = new TCanvas("c196", "  ");
  TCanvas *c197 = new TCanvas("c197", "  ");
  TCanvas *c198 = new TCanvas("c198", "  ");
  TCanvas *c199 = new TCanvas("c199", "  ");
  TCanvas *c200 = new TCanvas("c200", "  ");
  TCanvas *c201 = new TCanvas("c201", "  ");
  TCanvas *c202 = new TCanvas("c202", "  ");
  TCanvas *c203 = new TCanvas("c203", "  ");
  TCanvas *c204 = new TCanvas("c204", "  ");
  TCanvas *c205 = new TCanvas("c205", "  ");
  TCanvas *c206 = new TCanvas("c206", "  ");
  TCanvas *c207 = new TCanvas("c207", "  ");
  TCanvas *c208 = new TCanvas("c208", "  ");
  TCanvas *c209 = new TCanvas("c209", "  ");
  TCanvas *c210 = new TCanvas("c210", "  ");
  TCanvas *c211 = new TCanvas("c211", "  ");
  TCanvas *c212 = new TCanvas("c212", "  ");
  TCanvas *c213 = new TCanvas("c213", "  ");
  TCanvas *c214 = new TCanvas("c214", "  ");
  TCanvas *c215 = new TCanvas("c215", "  ");
  TCanvas *c216 = new TCanvas("c216", "  ");
  TCanvas *c217 = new TCanvas("c217", "  ");
  TCanvas *c218 = new TCanvas("c218", "  ");
  TCanvas *c219 = new TCanvas("c219", "  ");
  TCanvas *c220 = new TCanvas("c220", "  ");
  TCanvas *c221 = new TCanvas("c221", "  ");
  TCanvas *c222 = new TCanvas("c222", "  ");
  TCanvas *c223 = new TCanvas("c223", "  ");
  TCanvas *c224 = new TCanvas("c224", "  ");
  TCanvas *c225 = new TCanvas("c225", "  ");
  TCanvas *c226 = new TCanvas("c226", "  ");
  TCanvas *c227 = new TCanvas("c227", "  ");
  TCanvas *c228 = new TCanvas("c228", "  ");
  TCanvas *c229 = new TCanvas("c229", "  ");
  TCanvas *c230 = new TCanvas("c230", "  ");
  TCanvas *c231 = new TCanvas("c231", "  ");
  TCanvas *c232 = new TCanvas("c232", "  ");
  TCanvas *c233 = new TCanvas("c233", "  ");
  TCanvas *c234 = new TCanvas("c234", "  ");
  TCanvas *c235 = new TCanvas("c235", "  ");
  TCanvas *c236 = new TCanvas("c236", "  ");
  TCanvas *c237 = new TCanvas("c237", "  ");
  TCanvas *c238 = new TCanvas("c238", "  ");
  TCanvas *c239 = new TCanvas("c239", "  ");
  TCanvas *c240 = new TCanvas("c240", "  ");
  TCanvas *c241 = new TCanvas("c241", "  ");
  TCanvas *c242 = new TCanvas("c242", "  ");
  TCanvas *c243 = new TCanvas("c243", "  ");
  TCanvas *c244 = new TCanvas("c244", "  ");
  TCanvas *c245 = new TCanvas("c245", "  ");
  TCanvas *c246 = new TCanvas("c246", "  ");
  TCanvas *c247 = new TCanvas("c247", "  ");
  TCanvas *c248 = new TCanvas("c248", "  ");
  TCanvas *c249 = new TCanvas("c249", "  ");
  TCanvas *c250 = new TCanvas("c250", "  ");
  TCanvas *c251 = new TCanvas("c251", "  ");
  TCanvas *c252 = new TCanvas("c252", "  ");
  TCanvas *c253 = new TCanvas("c253", "  ");
  TCanvas *c254 = new TCanvas("c254", "  ");
  TCanvas *c255 = new TCanvas("c255", "  ");
  TCanvas *c256 = new TCanvas("c256", "  ");
  TCanvas *c257 = new TCanvas("c257", "  ");
  TCanvas *c258 = new TCanvas("c258", "  ");
  TCanvas *c259 = new TCanvas("c259", "  ");
  TCanvas *c260 = new TCanvas("c260", "  ");
  TCanvas *c261 = new TCanvas("c261", "  ");
  TCanvas *c262 = new TCanvas("c262", "  ");
  TCanvas *c263 = new TCanvas("c263", "  ");
  TCanvas *c264 = new TCanvas("c264", "  ");
  TCanvas *c265 = new TCanvas("c265", "  ");
  TCanvas *c266 = new TCanvas("c266", "  ");
  TCanvas *c267 = new TCanvas("c267", "  ");
  TCanvas *c268 = new TCanvas("c268", "  ");
  TCanvas *c269 = new TCanvas("c269", "  ");
  TCanvas *c270 = new TCanvas("c270", "  ");
  TCanvas *c271 = new TCanvas("c271", "  ");
  TCanvas *c272 = new TCanvas("c272", "  ");
  TCanvas *c273 = new TCanvas("c273", "  ");
  TCanvas *c274 = new TCanvas("c274", "  ");
  TCanvas *c275 = new TCanvas("c275", "  ");
  TCanvas *c276 = new TCanvas("c276", "  ");
  TCanvas *c277 = new TCanvas("c277", "  ");
  TCanvas *c278 = new TCanvas("c278", "  ");
  TCanvas *c279 = new TCanvas("c279", "  ");
  TCanvas *c280 = new TCanvas("c280", "  ");
  TCanvas *c281 = new TCanvas("c281", "  ");
  TCanvas *c282 = new TCanvas("c282", "  ");
  TCanvas *c283 = new TCanvas("c283", "  ");
  TCanvas *c284 = new TCanvas("c284", "  ");
  TCanvas *c285 = new TCanvas("c285", "  ");
  TCanvas *c286 = new TCanvas("c286", "  ");
  TCanvas *c287 = new TCanvas("c287", "  ");
  TCanvas *c288 = new TCanvas("c288", "  ");
  TCanvas *c289 = new TCanvas("c289", "  ");
  TCanvas *c290 = new TCanvas("c290", "  ");
  TCanvas *c291 = new TCanvas("c291", "  ");
  TCanvas *c292 = new TCanvas("c292", "  ");
  TCanvas *c293 = new TCanvas("c293", "  ");
  TCanvas *c294 = new TCanvas("c294", "  ");
  TCanvas *c295 = new TCanvas("c295", "  ");
  TCanvas *c296 = new TCanvas("c296", "  ");
  TCanvas *c297 = new TCanvas("c297", "  ");
  TCanvas *c298 = new TCanvas("c298", "  ");
  TCanvas *c299 = new TCanvas("c299", "  ");
  TCanvas *c300 = new TCanvas("c300", "  ");
  TCanvas *c301 = new TCanvas("c301", "  ");
  TCanvas *c302 = new TCanvas("c302", "  ");
  TCanvas *c303 = new TCanvas("c303", "  ");
  TCanvas *c304 = new TCanvas("c304", "  ");
  TCanvas *c305 = new TCanvas("c305", "  ");
  TCanvas *c306 = new TCanvas("c306", "  ");
  TCanvas *c307 = new TCanvas("c307", "  ");
  TCanvas *c308 = new TCanvas("c308", "  ");
  TCanvas *c309 = new TCanvas("c309", "  ");
  TCanvas *c310 = new TCanvas("c310", "  ");
  TCanvas *c311 = new TCanvas("c311", "  ");
  TCanvas *c312 = new TCanvas("c312", "  ");
  TCanvas *c313 = new TCanvas("c313", "  ");
  TCanvas *c314 = new TCanvas("c314", "  ");
  TCanvas *c315 = new TCanvas("c315", "  ");
  TCanvas *c316 = new TCanvas("c316", "  ");
  TCanvas *c317 = new TCanvas("c317", "  ");
  TCanvas *c318 = new TCanvas("c318", "  ");
  TCanvas *c319 = new TCanvas("c319", "  ");
  TCanvas *c320 = new TCanvas("c320", "  ");
  TCanvas *c321 = new TCanvas("c321", "  ");
  TCanvas *c322 = new TCanvas("c322", "  ");
  TCanvas *c323 = new TCanvas("c323", "  ");
  TCanvas *c324 = new TCanvas("c324", "  ");
  TCanvas *c325 = new TCanvas("c325", "  ");
  TCanvas *c326 = new TCanvas("c326", "  ");
  TCanvas *c327 = new TCanvas("c327", "  ");
  TCanvas *c328 = new TCanvas("c328", "  ");
  TCanvas *c329 = new TCanvas("c329", "  ");
  TCanvas *c330 = new TCanvas("c330", "  ");
  TCanvas *c331 = new TCanvas("c331", "  ");
  TCanvas *c332 = new TCanvas("c332", "  ");
  TCanvas *c333 = new TCanvas("c333", "  ");
  TCanvas *c334 = new TCanvas("c334", "  ");
  TCanvas *c335 = new TCanvas("c335", "  ");
  TCanvas *c336 = new TCanvas("c336", "  ");
  TCanvas *c337 = new TCanvas("c337", "  ");
  TCanvas *c338 = new TCanvas("c338", "  ");
  TCanvas *c339 = new TCanvas("c339", "  ");
  TCanvas *c340 = new TCanvas("c340", "  ");
  TCanvas *c341 = new TCanvas("c341", "  ");
  TCanvas *c342 = new TCanvas("c342", "  ");

  TCanvas *mycE_1 = new TCanvas("mycE_1","mycE_1",1500,1000);
  TCanvas *mycE_2 = new TCanvas("mycE_2","mycE_2",1500,1000);
  TCanvas *mycE_leak_1 = new TCanvas("mycE_leak_1","mycE_leak_1",1500,1000);
  TCanvas *mycE_leak_2 = new TCanvas("mycE_leak_2","mycE_leak_2",1500,1000);
  TCanvas *mycE_leak_offset_1 = new TCanvas("mycE_leak_offset_1","mycE_leak_offset_1",1500,1000);
  TCanvas *mycE_leak_offset_2 = new TCanvas("mycE_leak_offset_2","mycE_leak_offset_2",1500,1000);
  TCanvas *mycE_4sections_1 = new TCanvas("mycE_4sections_1","mycE_4sections_1",1500,1000);
  TCanvas *mycE_4sections_2 = new TCanvas("mycE_4sections_2","mycE_4sections_2",1500,1000);
  TCanvas *mycE_4sections_leak_1 = new TCanvas("mycE_4sections_leak_1","mycE_4sections_leak_1",1500,1000);
  TCanvas *mycE_4sections_leak_2 = new TCanvas("mycE_4sections_leak_2","mycE_4sections_leak_2",1500,1000);
  TCanvas *mycE_4sections_ch1_1 = new TCanvas("mycE_4sections_ch1_1","mycE_4sections_ch1_1",1500,1000);
  TCanvas *mycE_4sections_ch1_2 = new TCanvas("mycE_4sections_ch1_2","mycE_4sections_ch1_2",1500,1000);
  TCanvas *mycE_4sections_ch1_leak_1 = new TCanvas("mycE_4sections_ch1_leak_1","mycE_4sections_ch1_leak_1",1500,1000);
  TCanvas *mycE_4sections_ch1_leak_2 = new TCanvas("mycE_4sections_ch1_leak_2","mycE_4sections_ch1_leak_2",1500,1000);
  TCanvas *mycE_4sections_ch2_1 = new TCanvas("mycE_4sections_ch2_1","mycE_4sections_ch2_1",1500,1000);
  TCanvas *mycE_4sections_ch2_2 = new TCanvas("mycE_4sections_ch2_2","mycE_4sections_ch2_2",1500,1000);
  TCanvas *mycE_4sections_ch2_leak_1 = new TCanvas("mycE_4sections_ch2_leak_1","mycE_4sections_ch2_leak_1",1500,1000);
  TCanvas *mycE_4sections_ch2_leak_2 = new TCanvas("mycE_4sections_ch2_leak_2","mycE_4sections_ch2_leak_2",1500,1000);
  TCanvas *mycE_4sections_ch3_1 = new TCanvas("mycE_4sections_ch3_1","mycE_4sections_ch3_1",1500,1000);
  TCanvas *mycE_4sections_ch3_2 = new TCanvas("mycE_4sections_ch3_2","mycE_4sections_ch3_2",1500,1000);
  TCanvas *mycE_4sections_ch3_leak_1 = new TCanvas("mycE_4sections_ch3_leak_1","mycE_4sections_ch3_leak_1",1500,1000);
  TCanvas *mycE_4sections_ch3_leak_2 = new TCanvas("mycE_4sections_ch3_leak_2","mycE_4sections_ch3_leak_2",1500,1000);
  TCanvas *mycE_4sections_ch4_1 = new TCanvas("mycE_4sections_ch4_1","mycE_4sections_ch4_1",1500,1000);
  TCanvas *mycE_4sections_ch4_2 = new TCanvas("mycE_4sections_ch4_2","mycE_4sections_ch4_2",1500,1000);
  TCanvas *mycE_4sections_ch4_leak_1 = new TCanvas("mycE_4sections_ch4_leak_1","mycE_4sections_ch4_leak_1",1500,1000);
  TCanvas *mycE_4sections_ch4_leak_2 = new TCanvas("mycE_4sections_ch4_leak_2","mycE_4sections_ch4_leak_2",1500,1000);
  TCanvas *mycE_ch1_1 = new TCanvas("mycE_ch1_1","mycE_ch2_1",1500,1000);
  TCanvas *mycE_ch1_2 = new TCanvas("mycE_ch1_2","mycE_ch2_2",1500,1000);
  TCanvas *mycE_ch1_leak_1 = new TCanvas("mycE_ch1_leak_1","mycE_ch1_leak_1",1500,1000);
  TCanvas *mycE_ch1_leak_2 = new TCanvas("mycE_ch1_leak_2","mycE_ch1_leak_2",1500,1000);
  TCanvas *mycE_ch2_1 = new TCanvas("mycE_ch2_1","mycE_ch2_1",1500,1000);
  TCanvas *mycE_ch2_2 = new TCanvas("mycE_ch2_2","mycE_ch2_2",1500,1000);
  TCanvas *mycE_ch2_leak_1 = new TCanvas("mycE_ch2_leak_1","mycE_ch2_leak_1",1500,1000);
  TCanvas *mycE_ch2_leak_2 = new TCanvas("mycE_ch2_leak_2","mycE_ch2_leak_2",1500,1000);
  TCanvas *mycE_ch3_1 = new TCanvas("mycE_ch3_1","mycE_ch3_1",1500,1000);
  TCanvas *mycE_ch3_2 = new TCanvas("mycE_ch3_2","mycE_ch3_2",1500,1000);
  TCanvas *mycE_ch3_leak_1 = new TCanvas("mycE_ch3_leak_1","mycE_ch3_leak_1",1500,1000);
  TCanvas *mycE_ch3_leak_2 = new TCanvas("mycE_ch3_leak_2","mycE_ch3_leak_2",1500,1000);
  TCanvas *mycE_ch4_1 = new TCanvas("mycE_ch4_1","mycE_ch4_1",1500,1000);
  TCanvas *mycE_ch4_2 = new TCanvas("mycE_ch4_2","mycE_ch4_2",1500,1000);
  TCanvas *mycE_ch4_leak_1 = new TCanvas("mycE_ch4_leak_1","mycE_ch4_leak_1",1500,1000);
  TCanvas *mycE_ch4_leak_2 = new TCanvas("mycE_ch4_leak_2","mycE_ch4_leak_2",1500,1000);
  TCanvas *mycE_totsf_1 = new TCanvas("mycE_totsf_1","mycE_totsf_1",1500,1000);
  TCanvas *mycE_totsf_2 = new TCanvas("mycE_totsf_2","mycE_totsf_2",1500,1000);
  TCanvas *mycE_totsf_leak_1 = new TCanvas("mycE_totsf_leak_1","mycE_totsf_leak_1",1500,1000);
  TCanvas *mycE_totsf_leak_2 = new TCanvas("mycE_totsf_leak_2","mycE_totsf_leak_2",1500,1000);
  TCanvas *mycE_totsf_dedxweighted_1 = new TCanvas("mycE_totsf_dedxweighted_1","mycE_totsf_dedxweighted_1",1500,1000);
  TCanvas *mycE_totsf_dedxweighted_2 = new TCanvas("mycE_totsf_dedxweighted_2","mycE_totsf_dedxweighted_2",1500,1000);
  TCanvas *mycE_totsf_dedxweighted_leak_1 = new TCanvas("mycE_totsf_dedxweighted_leak_1","mycE_totsf_dedxweighted_leak_1",1500,1000);
  TCanvas *mycE_totsf_dedxweighted_leak_2 = new TCanvas("mycE_totsf_dedxweighted_leak_2","mycE_totsf_dedxweighted_leak_2",1500,1000);
  TCanvas *mycE_totsf_sec_1 = new TCanvas("mycE_totsf_sec_1","mycE_totsf_sec_1",1500,1000);
  TCanvas *mycE_totsf_sec_2 = new TCanvas("mycE_totsf_sec_2","mycE_totsf_sec_2",1500,1000);
  TCanvas *mycE_totsf_sec_leak_1 = new TCanvas("mycE_totsf_sec_leak_1","mycE_totsf_sec_leak_1",1500,1000);
  TCanvas *mycE_totsf_sec_leak_2 = new TCanvas("mycE_totsf_sec_leak_2","mycE_totsf_sec_leak_2",1500,1000);
  TCanvas *mycE_totsf_sec_dedxweighted_1 = new TCanvas("mycE_totsf_sec_dedxweighted_1","mycE_totsf_sec_dedxweighted_1",1500,1000);
  TCanvas *mycE_totsf_sec_dedxweighted_2 = new TCanvas("mycE_totsf_sec_dedxweighted_2","mycE_totsf_sec_dedxweighted_2",1500,1000);
  TCanvas *mycE_totsf_sec_dedxweighted_leak_1 = new TCanvas("mycE_totsf_sec_dedxweighted_leak_1","mycE_totsf_sec_dedxweighted_leak_1",1500,1000);
  TCanvas *mycE_totsf_sec_dedxweighted_leak_2 = new TCanvas("mycE_totsf_sec_dedxweighted_leak_2","mycE_totsf_sec_dedxweighted_leak_2",1500,1000);
  TCanvas *mycE_totsfshmax_sec_1 = new TCanvas("mycE_totsfshmax_sec_1","mycE_totsfshmax_sec_1",1500,1000);
  TCanvas *mycE_totsfshmax_sec_2 = new TCanvas("mycE_totsfshmax_sec_2","mycE_totsfshmax_sec_2",1500,1000);
  TCanvas *mycE_totsfshmax_sec_leak_1 = new TCanvas("mycE_totsfshmax_sec_leak_1","mycE_totsfshmax_sec_leak_1",1500,1000);
  TCanvas *mycE_totsfshmax_sec_leak_2 = new TCanvas("mycE_totsfshmax_sec_leak_2","mycE_totsfshmax_sec_leak_2",1500,1000);
  TCanvas *mycE_totsfshmax_sec_dedxweighted_1 = new TCanvas("mycE_totsfshmax_sec_dedxweighted_1","mycE_totsfshmax_sec_dedxweighted_1",1500,1000);
  TCanvas *mycE_totsfshmax_sec_dedxweighted_2 = new TCanvas("mycE_totsfshmax_sec_dedxweighted_2","mycE_totsfshmax_sec_dedxweighted_2",1500,1000);
  TCanvas *mycE_totsfshmax_sec_dedxweighted_leak_1 = new TCanvas("mycE_totsfshmax_sec_dedxweighted_leak_1","mycE_totsfshmax_sec_dedxweighted_leak_1",1500,1000);
  TCanvas *mycE_totsfshmax_sec_dedxweighted_leak_2 = new TCanvas("mycE_totsfshmax_sec_dedxweighted_leak_2","mycE_totsfshmax_sec_dedxweighted_leak_2",1500,1000);
  TCanvas *mycE_totsf_dedxsec_1 = new TCanvas("mycE_totsf_dedxsec_1","mycE_totsf_dedxsec_1",1500,1000);
  TCanvas *mycE_totsf_dedxsec_2 = new TCanvas("mycE_totsf_dedxsec_2","mycE_totsf_dedxsec_2",1500,1000);
  TCanvas *mycE_totsf_dedxsec_leak_1 = new TCanvas("mycE_totsf_dedxsec_leak_1","mycE_totsf_dedxsec_leak_1",1500,1000);
  TCanvas *mycE_totsf_dedxsec_leak_2 = new TCanvas("mycE_totsf_dedxsec_leak_2","mycE_totsf_dedxsec_leak_2",1500,1000);
  TCanvas *mycE_totsfshmax_dedxsec_1 = new TCanvas("mycE_totsfshmax_dedxsec_1","mycE_totsfshmax_dedxsec_1",1500,1000);
  TCanvas *mycE_totsfshmax_dedxsec_2 = new TCanvas("mycE_totsfshmax_dedxsec_2","mycE_totsfshmax_dedxsec_2",1500,1000);
  TCanvas *mycE_totsfshmax_dedxsec_leak_1 = new TCanvas("mycE_totsfshmax_dedxsec_leak_1","mycE_totsfshmax_dedxsec_leak_1",1500,1000);
  TCanvas *mycE_totsfshmax_dedxsec_leak_2 = new TCanvas("mycE_totsfshmax_dedxsec_leak_2","mycE_totsfshmax_dedxsec_leak_2",1500,1000);
  TCanvas *mycE_shmax_1 = new TCanvas("mycE_shmax_1","mycE_shmax_1",1500,1000);
  TCanvas *mycE_shmax_2 = new TCanvas("mycE_shmax_2","mycE_shmax_2",1500,1000);
  TCanvas *mycE_shmax_leak_1 = new TCanvas("mycE_shmax_leak_1","mycE_shmax_leak_1",1500,1000);
  TCanvas *mycE_shmax_leak_2 = new TCanvas("mycE_shmax_leak_2","mycE_shmax_leak_2",1500,1000);
  TCanvas *mycE_sfnorm_1 = new TCanvas("mycE_sfnorm_1","mycE_sfnorm_1",1500,1000);
  TCanvas *mycE_sfnorm_2 = new TCanvas("mycE_sfnorm_2","mycE_sfnorm_2",1500,1000);
  TCanvas *mycE_sfnorm_leak_1 = new TCanvas("mycE_sfnorm_leak_1","mycE_sfnorm_leak_1",1500,1000);
  TCanvas *mycE_sfnorm_leak_2 = new TCanvas("mycE_sfnorm_leak_2","mycE_sfnorm_leak_2",1500,1000);
  TCanvas *mycE_totsfnorm_1 = new TCanvas("mycE_totsfnorm_1","mycE_totsfnorm_1",1500,1000);
  TCanvas *mycE_totsfnorm_2 = new TCanvas("mycE_totsfnorm_2","mycE_totsfnorm_2",1500,1000);
  TCanvas *mycE_totsfnorm_leak_1 = new TCanvas("mycE_totsfnorm_leak_1","mycE_totsfnorm_leak_1",1500,1000);
  TCanvas *mycE_totsfnorm_leak_2 = new TCanvas("mycE_totsfnorm_leak_2","mycE_totsfnorm_leak_2",1500,1000);
  TCanvas *mycE_totsfnorm_dedxweighted_1 = new TCanvas("mycE_totsfnorm_dedxweighted_1","mycE_totsfnorm_dedxweighted_1",1500,1000);
  TCanvas *mycE_totsfnorm_dedxweighted_2 = new TCanvas("mycE_totsfnorm_dedxweighted_2","mycE_totsfnorm_dedxweighted_2",1500,1000);
  TCanvas *mycE_totsfnorm_dedxweighted_leak_1 = new TCanvas("mycE_totsfnorm_dedxweighted_leak_1","mycE_totsfnorm_dedxweighted_leak_1",1500,1000);
  TCanvas *mycE_totsfnorm_dedxweighted_leak_2 = new TCanvas("mycE_totsfnorm_dedxweighted_leak_2","mycE_totsfnorm_dedxweighted_leak_2",1500,1000);
  TCanvas *mycE_raw_1 = new TCanvas("mycE_raw_1","mycE_raw_1",1500,1000);
  TCanvas *mycE_raw_2 = new TCanvas("mycE_raw_2","mycE_raw_2",1500,1000);
  TCanvas *mycE_raw_leak_1 = new TCanvas("mycE_leak_raw_1","mycE_leak_raw_1",1500,1000);
  TCanvas *mycE_raw_leak_2 = new TCanvas("mycE_leak_raw_2","mycE_leak_raw_2",1500,1000);
  TCanvas *mycE_corr_1 = new TCanvas("mycE_corr_1","mycE_corr_1",1500,1000);
  TCanvas *mycE_corr_2 = new TCanvas("mycE_corr_2","mycE_corr_2",1500,1000);
  TCanvas *mycE_val_1 = new TCanvas("mycE_val_1","mycE_val_1",1500,1000);
  TCanvas *mycE_val_2 = new TCanvas("mycE_val_2","mycE_val_2",1500,1000);
  TCanvas *mycE_overtrue_1 = new TCanvas("mycE_overtrue_1","mycE_overtrue_1",1500,1000);
  TCanvas *mycE_overtrue_2 = new TCanvas("mycE_overtrue_2","mycE_overtrue_2",1500,1000);
  TCanvas *mycE_overtrue_leak_1 = new TCanvas("mycE_overtrue_leak_1","mycE_overtrue_leak_1",1500,1000);
  TCanvas *mycE_overtrue_leak_2 = new TCanvas("mycE_overtrue_leak_2","mycE_overtrue_leak_2",1500,1000);
  TCanvas *mycE_overtrue_leak_offset_1 = new TCanvas("mycE_overtrue_leak_offset_1","mycE_overtrue_leak_offset_1",1500,1000);
  TCanvas *mycE_overtrue_leak_offset_2 = new TCanvas("mycE_overtrue_leak_offset_2","mycE_overtrue_leak_offset_2",1500,1000);
  TCanvas *mycE_overtrue_totsf_1 = new TCanvas("mycE_overtrue_totsf_1","mycE_overtrue_totsf_1",1500,1000);
  TCanvas *mycE_overtrue_totsf_2 = new TCanvas("mycE_overtrue_totsf_2","mycE_overtrue_totsf_2",1500,1000);
  TCanvas *mycE_overtrue_totsf_leak_1 = new TCanvas("mycE_overtrue_totsf_leak_1","mycE_overtrue_totsf_leak_1",1500,1000);
  TCanvas *mycE_overtrue_totsf_leak_2 = new TCanvas("mycE_overtrue_totsf_leak_2","mycE_overtrue_totsf_leak_2",1500,1000);
  TCanvas *mycE_overtrue_totsf_dedxweighted_1 = new TCanvas("mycE_overtrue_totsf_dedxweighted_1","mycE_overtrue_totsf_dedxweighted_1",1500,1000);
  TCanvas *mycE_overtrue_totsf_dedxweighted_2 = new TCanvas("mycE_overtrue_totsf_dedxweighted_2","mycE_overtrue_totsf_dedxweighted_2",1500,1000);
  TCanvas *mycE_overtrue_totsf_dedxweighted_leak_1 = new TCanvas("mycE_overtrue_totsf_dedxweighted_leak_1","mycE_overtrue_totsf_dedxweighted_leak_1",1500,1000);
  TCanvas *mycE_overtrue_totsf_dedxweighted_leak_2 = new TCanvas("mycE_overtrue_totsf_dedxweighted_leak_2","mycE_overtrue_totsf_dedxweighted_leak_2",1500,1000);
  TCanvas *mycE_overtrue_totsf_sec_1 = new TCanvas("mycE_overtrue_totsf_sec_1","mycE_overtrue_totsf_sec_1",1500,1000);
  TCanvas *mycE_overtrue_totsf_sec_2 = new TCanvas("mycE_overtrue_totsf_sec_2","mycE_overtrue_totsf_sec_2",1500,1000);
  TCanvas *mycE_overtrue_totsf_sec_leak_1 = new TCanvas("mycE_overtrue_totsf_sec_leak_1","mycE_overtrue_totsf_sec_leak_1",1500,1000);
  TCanvas *mycE_overtrue_totsf_sec_leak_2 = new TCanvas("mycE_overtrue_totsf_sec_leak_2","mycE_overtrue_totsf_sec_leak_2",1500,1000);
  TCanvas *mycE_overtrue_totsf_sec_dedxweighted_1 = new TCanvas("mycE_overtrue_totsf_sec_dedxweighted_1","mycE_overtrue_totsf_sec_dedxweighted_1",1500,1000);
  TCanvas *mycE_overtrue_totsf_sec_dedxweighted_2 = new TCanvas("mycE_overtrue_totsf_sec_dedxweighted_2","mycE_overtrue_totsf_sec_dedxweighted_2",1500,1000);
  TCanvas *mycE_overtrue_totsf_sec_dedxweighted_leak_1 = new TCanvas("mycE_overtrue_totsf_sec_dedxweighted_leak_1","mycE_overtrue_totsf_sec_dedxweighted_leak_1",1500,1000);
  TCanvas *mycE_overtrue_totsf_sec_dedxweighted_leak_2 = new TCanvas("mycE_overtrue_totsf_sec_dedxweighted_leak_2","mycE_overtrue_totsf_sec_dedxweighted_leak_2",1500,1000);
  TCanvas *mycE_overtrue_totsfshmax_sec_1 = new TCanvas("mycE_overtrue_totsfshmax_sec_1","mycE_overtrue_totsfshmax_sec_1",1500,1000);
  TCanvas *mycE_overtrue_totsfshmax_sec_2 = new TCanvas("mycE_overtrue_totsfshmax_sec_2","mycE_overtrue_totsfshmax_sec_2",1500,1000);
  TCanvas *mycE_overtrue_totsfshmax_sec_leak_1 = new TCanvas("mycE_overtrue_totsfshmax_sec_leak_1","mycE_overtrue_totsfshmax_sec_leak_1",1500,1000);
  TCanvas *mycE_overtrue_totsfshmax_sec_leak_2 = new TCanvas("mycE_overtrue_totsfshmax_sec_leak_2","mycE_overtrue_totsfshmax_sec_leak_2",1500,1000);
  TCanvas *mycE_overtrue_totsfshmax_sec_dedxweighted_1 = new TCanvas("mycE_overtrue_totsfshmax_sec_dedxweighted_1","mycE_overtrue_totsfshmax_sec_dedxweighted_1",1500,1000);
  TCanvas *mycE_overtrue_totsfshmax_sec_dedxweighted_2 = new TCanvas("mycE_overtrue_totsfshmax_sec_dedxweighted_2","mycE_overtrue_totsfshmax_sec_dedxweighted_2",1500,1000);
  TCanvas *mycE_overtrue_totsfshmax_sec_dedxweighted_leak_1 = new TCanvas("mycE_overtrue_totsfshmax_sec_dedxweighted_leak_1","mycE_overtrue_totsfshmax_sec_dedxweighted_leak_1",1500,1000);
  TCanvas *mycE_overtrue_totsfshmax_sec_dedxweighted_leak_2 = new TCanvas("mycE_overtrue_totsfshmax_sec_dedxweighted_leak_2","mycE_overtrue_totsfshmax_sec_dedxweighted_leak_2",1500,1000);
  TCanvas *mycE_overtrue_totsf_dedxsec_1 = new TCanvas("mycE_overtrue_totsf_dedxsec_1","mycE_overtrue_totsf_dedxsec_1",1500,1000);
  TCanvas *mycE_overtrue_totsf_dedxsec_2 = new TCanvas("mycE_overtrue_totsf_dedxsec_2","mycE_overtrue_totsf_dedxsec_2",1500,1000);
  TCanvas *mycE_overtrue_totsf_dedxsec_leak_1 = new TCanvas("mycE_overtrue_totsf_dedxsec_leak_1","mycE_overtrue_totsf_dedxsec_leak_1",1500,1000);
  TCanvas *mycE_overtrue_totsf_dedxsec_leak_2 = new TCanvas("mycE_overtrue_totsf_dedxsec_leak_2","mycE_overtrue_totsf_dedxsec_leak_2",1500,1000);
  TCanvas *mycE_overtrue_totsfshmax_dedxsec_1 = new TCanvas("mycE_overtrue_totsfshmax_dedxsec_1","mycE_overtrue_totsfshmax_dedxsec_1",1500,1000);
  TCanvas *mycE_overtrue_totsfshmax_dedxsec_2 = new TCanvas("mycE_overtrue_totsfshmax_dedxsec_2","mycE_overtrue_totsfshmax_dedxsec_2",1500,1000);
  TCanvas *mycE_overtrue_totsfshmax_dedxsec_leak_1 = new TCanvas("mycE_overtrue_totsfshmax_dedxsec_leak_1","mycE_overtrue_totsfshmax_dedxsec_leak_1",1500,1000);
  TCanvas *mycE_overtrue_totsfshmax_dedxsec_leak_2 = new TCanvas("mycE_overtrue_totsfshmax_dedxsec_leak_2","mycE_overtrue_totsfshmax_dedxsec_leak_2",1500,1000);
  TCanvas *mycE_overtrue_totsfnorm_1 = new TCanvas("mycE_overtrue_totsfnorm_1","mycE_overtrue_totsfnorm_1",1500,1000);
  TCanvas *mycE_overtrue_totsfnorm_2 = new TCanvas("mycE_overtrue_totsfnorm_2","mycE_overtrue_totsfnorm_2",1500,1000);
  TCanvas *mycE_overtrue_totsfnorm_leak_1 = new TCanvas("mycE_overtrue_totsfnorm_leak_1","mycE_overtrue_totsfnorm_leak_1",1500,1000);
  TCanvas *mycE_overtrue_totsfnorm_leak_2 = new TCanvas("mycE_overtrue_totsfnorm_leak_2","mycE_overtrue_totsfnorm_leak_2",1500,1000);
  TCanvas *mycE_overtrue_totsfnorm_dedxweighted_1 = new TCanvas("mycE_overtrue_totsfnorm_dedxweighted_1","mycE_overtrue_totsfnorm_dedxweighted_1",1500,1000);
  TCanvas *mycE_overtrue_totsfnorm_dedxweighted_2 = new TCanvas("mycE_overtrue_totsfnorm_dedxweighted_2","mycE_overtrue_totsfnorm_dedxweighted_2",1500,1000);
  TCanvas *mycE_overtrue_totsfnorm_dedxweighted_leak_1 = new TCanvas("mycE_overtrue_totsfnorm_dedxweighted_leak_1","mycE_overtrue_totsfnorm_dedxweighted_leak_1",1500,1000);
  TCanvas *mycE_overtrue_totsfnorm_dedxweighted_leak_2 = new TCanvas("mycE_overtrue_totsfnorm_dedxweighted_leak_2","mycE_overtrue_totsfnorm_dedxweighted_leak_2",1500,1000);
  TCanvas *mycE_overtrue_sfnorm_1 = new TCanvas("mycE_overtrue_sfnorm_1","mycE_overtrue_sfnorm_1",1500,1000);
  TCanvas *mycE_overtrue_sfnorm_2 = new TCanvas("mycE_overtrue_sfnorm_2","mycE_overtrue_sfnorm_2",1500,1000);
  TCanvas *mycE_overtrue_sfnorm_leak_1 = new TCanvas("mycE_overtrue_sfnorm_leak_1","mycE_overtrue_sfnorm_leak_1",1500,1000);
  TCanvas *mycE_overtrue_sfnorm_leak_2 = new TCanvas("mycE_overtrue_sfnorm_leak_2","mycE_overtrue_sfnorm_leak_2",1500,1000);
  TCanvas *mycE_overtrue_corr_1 = new TCanvas("mycE_overtrue_corr_1","mycE_overtrue_corr_1",1500,1000);
  TCanvas *mycE_overtrue_corr_2 = new TCanvas("mycE_overtrue_corr_2","mycE_overtrue_corr_2",1500,1000);
  TCanvas *mycE_overtrue_val_1 = new TCanvas("mycE_overtrue_val_1","mycE_overtrue_val_1",1500,1000);
  TCanvas *mycE_overtrue_val_2 = new TCanvas("mycE_overtrue_val_2","mycE_overtrue_val_2",1500,1000);
  TCanvas *mycE_ch1_overtrueE_1 = new TCanvas("mycE_overtrueE_ch1_1","mycE_ch2_overtrueE_1",1500,1000);
  TCanvas *mycE_ch1_overtrueE_2 = new TCanvas("mycE_overtrueE_ch1_2","mycE_ch2_overtrueE_2",1500,1000);
  TCanvas *mycE_ch1_overtrueE_leak_1 = new TCanvas("mycE_overtrueE_ch1_leak_1","mycE_ch1_overtrueE_leak_1",1500,1000);
  TCanvas *mycE_ch1_overtrueE_leak_2 = new TCanvas("mycE_overtrueE_ch1_leak_2","mycE_ch1_overtrueE_leak_2",1500,1000);
  TCanvas *mycE_ch2_overtrueE_1 = new TCanvas("mycE_ch2_overtrueE_1","mycE_ch2_overtrueE_1",1500,1000);
  TCanvas *mycE_ch2_overtrueE_2 = new TCanvas("mycE_ch2_overtrueE_2","mycE_ch2_overtrueE_2",1500,1000);
  TCanvas *mycE_ch2_overtrueE_leak_1 = new TCanvas("mycE_ch2_overtrueE_leak_1","mycE_ch2_overtrueE_leak_1",1500,1000);
  TCanvas *mycE_ch2_overtrueE_leak_2 = new TCanvas("mycE_ch2_overtrueE_leak_2","mycE_ch2_overtrueE_leak_2",1500,1000);
  TCanvas *mycE_ch3_overtrueE_1 = new TCanvas("mycE_ch3_overtrueE_1","mycE_ch3_overtrueE_1",1500,1000);
  TCanvas *mycE_ch3_overtrueE_2 = new TCanvas("mycE_ch3_overtrueE_2","mycE_ch3_overtrueE_2",1500,1000);
  TCanvas *mycE_ch3_overtrueE_leak_1 = new TCanvas("mycE_ch3_overtrueE_leak_1","mycE_ch3_overtrueE_leak_1",1500,1000);
  TCanvas *mycE_ch3_overtrueE_leak_2 = new TCanvas("mycE_ch3_overtrueE_leak_2","mycE_ch3_overtrueE_leak_2",1500,1000);
  TCanvas *mycE_ch4_overtrueE_1 = new TCanvas("mycE_ch4_overtrueE_1","mycE_ch4_overtrueE_1",1500,1000);
  TCanvas *mycE_ch4_overtrueE_2 = new TCanvas("mycE_ch4_overtrueE_2","mycE_ch4_overtrueE_2",1500,1000);
  TCanvas *mycE_ch4_overtrueE_leak_1 = new TCanvas("mycE_ch4_overtrueE_leak_1","mycE_ch4_overtrueE_leak_1",1500,1000);
  TCanvas *mycE_ch4_overtrueE_leak_2 = new TCanvas("mycE_ch4_overtrueE_leak_2","mycE_ch4_overtrueE_leak_2",1500,1000);
  TCanvas *mycE_overtrue_shmax_1 = new TCanvas("mycE_overtrue_shmax_1","mycE_overtrue_shmax_1",1500,1000);
  TCanvas *mycE_overtrue_shmax_2 = new TCanvas("mycE_overtrue_shmax_2","mycE_overtrue_shmax_2",1500,1000);
  TCanvas *mycE_overtrue_shmax_leak_1 = new TCanvas("mycE_overtrue_shmax_leak_1","mycE_overtrue_shmax_leak_1",1500,1000);
  TCanvas *mycE_overtrue_shmax_leak_2 = new TCanvas("mycE_overtrue_shmax_leak_2","mycE_overtrue_shmax_leak_2",1500,1000);

  mycE_1->Divide(3,3);
  mycE_2->Divide(3,3);
  mycE_leak_1->Divide(3,3);
  mycE_leak_2->Divide(3,3);
  mycE_leak_offset_1->Divide(3,3);
  mycE_leak_offset_2->Divide(3,3);
  mycE_4sections_1->Divide(3,3);
  mycE_4sections_2->Divide(3,3);
  mycE_4sections_leak_1->Divide(3,3);
  mycE_4sections_leak_2->Divide(3,3);
  mycE_4sections_ch1_1->Divide(3,3);
  mycE_4sections_ch1_2->Divide(3,3);
  mycE_4sections_ch1_leak_1->Divide(3,3);
  mycE_4sections_ch1_leak_2->Divide(3,3);
  mycE_4sections_ch2_1->Divide(3,3);
  mycE_4sections_ch2_2->Divide(3,3);
  mycE_4sections_ch2_leak_1->Divide(3,3);
  mycE_4sections_ch2_leak_2->Divide(3,3);
  mycE_4sections_ch3_1->Divide(3,3);
  mycE_4sections_ch3_2->Divide(3,3);
  mycE_4sections_ch3_leak_1->Divide(3,3);
  mycE_4sections_ch3_leak_2->Divide(3,3);
  mycE_4sections_ch4_1->Divide(3,3);
  mycE_4sections_ch4_2->Divide(3,3);
  mycE_4sections_ch4_leak_1->Divide(3,3);
  mycE_4sections_ch4_leak_2->Divide(3,3);
  mycE_ch1_1->Divide(3,3);
  mycE_ch1_2->Divide(3,3);
  mycE_ch1_leak_1->Divide(3,3);
  mycE_ch1_leak_2->Divide(3,3);
  mycE_ch2_1->Divide(3,3);
  mycE_ch2_2->Divide(3,3);
  mycE_ch2_leak_1->Divide(3,3);
  mycE_ch2_leak_2->Divide(3,3);
  mycE_ch3_1->Divide(3,3);
  mycE_ch3_2->Divide(3,3);
  mycE_ch3_leak_1->Divide(3,3);
  mycE_ch3_leak_2->Divide(3,3);
  mycE_ch4_1->Divide(3,3);
  mycE_ch4_2->Divide(3,3);
  mycE_ch4_leak_1->Divide(3,3);
  mycE_ch4_leak_2->Divide(3,3);
  mycE_totsf_1->Divide(3,3);
  mycE_totsf_2->Divide(3,3);
  mycE_totsf_leak_1->Divide(3,3);
  mycE_totsf_leak_2->Divide(3,3);
  mycE_totsf_dedxweighted_1->Divide(3,3);
  mycE_totsf_dedxweighted_2->Divide(3,3);
  mycE_totsf_dedxweighted_leak_1->Divide(3,3);
  mycE_totsf_dedxweighted_leak_2->Divide(3,3);
  mycE_totsf_sec_1->Divide(3,3);
  mycE_totsf_sec_2->Divide(3,3);
  mycE_totsf_sec_leak_1->Divide(3,3);
  mycE_totsf_sec_leak_2->Divide(3,3);
  mycE_totsf_sec_dedxweighted_1->Divide(3,3);
  mycE_totsf_sec_dedxweighted_2->Divide(3,3);
  mycE_totsf_sec_dedxweighted_leak_1->Divide(3,3);
  mycE_totsf_sec_dedxweighted_leak_2->Divide(3,3);
  mycE_totsfshmax_sec_1->Divide(3,3);
  mycE_totsfshmax_sec_2->Divide(3,3);
  mycE_totsfshmax_sec_leak_1->Divide(3,3);
  mycE_totsfshmax_sec_leak_2->Divide(3,3);
  mycE_totsfshmax_sec_dedxweighted_1->Divide(3,3);
  mycE_totsfshmax_sec_dedxweighted_2->Divide(3,3);
  mycE_totsfshmax_sec_dedxweighted_leak_1->Divide(3,3);
  mycE_totsfshmax_sec_dedxweighted_leak_2->Divide(3,3);
  mycE_totsf_dedxsec_1->Divide(3,3);
  mycE_totsf_dedxsec_2->Divide(3,3);
  mycE_totsf_dedxsec_leak_1->Divide(3,3);
  mycE_totsf_dedxsec_leak_2->Divide(3,3);
  mycE_totsfshmax_dedxsec_1->Divide(3,3);
  mycE_totsfshmax_dedxsec_2->Divide(3,3);
  mycE_totsfshmax_dedxsec_leak_1->Divide(3,3);
  mycE_totsfshmax_dedxsec_leak_2->Divide(3,3);
  mycE_shmax_1->Divide(3,3);
  mycE_shmax_2->Divide(3,3);
  mycE_shmax_leak_1->Divide(3,3);
  mycE_shmax_leak_2->Divide(3,3);
  mycE_sfnorm_1->Divide(3,3);
  mycE_sfnorm_2->Divide(3,3);
  mycE_sfnorm_leak_1->Divide(3,3);
  mycE_sfnorm_leak_2->Divide(3,3);
  mycE_totsfnorm_1->Divide(3,3);
  mycE_totsfnorm_2->Divide(3,3);
  mycE_totsfnorm_leak_1->Divide(3,3);
  mycE_totsfnorm_leak_2->Divide(3,3);
  mycE_totsfnorm_dedxweighted_1->Divide(3,3);
  mycE_totsfnorm_dedxweighted_2->Divide(3,3);
  mycE_totsfnorm_dedxweighted_leak_1->Divide(3,3);
  mycE_totsfnorm_dedxweighted_leak_2->Divide(3,3);
  mycE_raw_1->Divide(3,3);
  mycE_raw_2->Divide(3,3);
  mycE_raw_leak_1->Divide(3,3);
  mycE_raw_leak_2->Divide(3,3);
  mycE_corr_1->Divide(3,3);
  mycE_corr_2->Divide(3,3);
  mycE_val_1->Divide(3,3);
  mycE_val_2->Divide(3,3);
  mycE_overtrue_1->Divide(3,3);
  mycE_overtrue_2->Divide(3,3);
  mycE_overtrue_leak_1->Divide(3,3);
  mycE_overtrue_leak_2->Divide(3,3);
  mycE_overtrue_leak_offset_1->Divide(3,3);
  mycE_overtrue_leak_offset_2->Divide(3,3);
  mycE_overtrue_totsf_1->Divide(3,3);
  mycE_overtrue_totsf_2->Divide(3,3);
  mycE_overtrue_totsf_leak_1->Divide(3,3);
  mycE_overtrue_totsf_leak_2->Divide(3,3);
  mycE_overtrue_totsf_dedxweighted_1->Divide(3,3);
  mycE_overtrue_totsf_dedxweighted_2->Divide(3,3);
  mycE_overtrue_totsf_dedxweighted_leak_1->Divide(3,3);
  mycE_overtrue_totsf_dedxweighted_leak_2->Divide(3,3);
  mycE_overtrue_totsf_sec_1->Divide(3,3);
  mycE_overtrue_totsf_sec_2->Divide(3,3);
  mycE_overtrue_totsf_sec_leak_1->Divide(3,3);
  mycE_overtrue_totsf_sec_leak_2->Divide(3,3);
  mycE_overtrue_totsf_sec_dedxweighted_1->Divide(3,3);
  mycE_overtrue_totsf_sec_dedxweighted_2->Divide(3,3);
  mycE_overtrue_totsf_sec_dedxweighted_leak_1->Divide(3,3);
  mycE_overtrue_totsf_sec_dedxweighted_leak_2->Divide(3,3);
  mycE_overtrue_totsfshmax_sec_1->Divide(3,3);
  mycE_overtrue_totsfshmax_sec_2->Divide(3,3);
  mycE_overtrue_totsfshmax_sec_leak_1->Divide(3,3);
  mycE_overtrue_totsfshmax_sec_leak_2->Divide(3,3);
  mycE_overtrue_totsfshmax_sec_dedxweighted_1->Divide(3,3);
  mycE_overtrue_totsfshmax_sec_dedxweighted_2->Divide(3,3);
  mycE_overtrue_totsfshmax_sec_dedxweighted_leak_1->Divide(3,3);
  mycE_overtrue_totsfshmax_sec_dedxweighted_leak_2->Divide(3,3);
  mycE_overtrue_totsf_dedxsec_1->Divide(3,3);
  mycE_overtrue_totsf_dedxsec_2->Divide(3,3);
  mycE_overtrue_totsf_dedxsec_leak_1->Divide(3,3);
  mycE_overtrue_totsf_dedxsec_leak_2->Divide(3,3);
  mycE_overtrue_totsfshmax_dedxsec_1->Divide(3,3);
  mycE_overtrue_totsfshmax_dedxsec_2->Divide(3,3);
  mycE_overtrue_totsfshmax_dedxsec_leak_1->Divide(3,3);
  mycE_overtrue_totsfshmax_dedxsec_leak_2->Divide(3,3);
  mycE_overtrue_totsfnorm_1->Divide(3,3);
  mycE_overtrue_totsfnorm_2->Divide(3,3);
  mycE_overtrue_totsfnorm_leak_1->Divide(3,3);
  mycE_overtrue_totsfnorm_leak_2->Divide(3,3);
  mycE_overtrue_totsfnorm_dedxweighted_1->Divide(3,3);
  mycE_overtrue_totsfnorm_dedxweighted_2->Divide(3,3);
  mycE_overtrue_totsfnorm_dedxweighted_leak_1->Divide(3,3);
  mycE_overtrue_totsfnorm_dedxweighted_leak_2->Divide(3,3);
  mycE_overtrue_sfnorm_1->Divide(3,3);
  mycE_overtrue_sfnorm_2->Divide(3,3);
  mycE_overtrue_sfnorm_leak_1->Divide(3,3);
  mycE_overtrue_sfnorm_leak_2->Divide(3,3);
  mycE_overtrue_shmax_1->Divide(3,3);
  mycE_overtrue_shmax_2->Divide(3,3);
  mycE_overtrue_shmax_leak_1->Divide(3,3);
  mycE_overtrue_shmax_leak_2->Divide(3,3);
  mycE_overtrue_corr_1->Divide(3,3);
  mycE_overtrue_corr_2->Divide(3,3);
  mycE_overtrue_val_1->Divide(3,3);
  mycE_overtrue_val_2->Divide(3,3);
  mycE_ch1_overtrueE_1->Divide(3,3);
  mycE_ch1_overtrueE_2->Divide(3,3);
  mycE_ch1_overtrueE_leak_1->Divide(3,3);
  mycE_ch1_overtrueE_leak_2->Divide(3,3);
  mycE_ch2_overtrueE_1->Divide(3,3);
  mycE_ch2_overtrueE_2->Divide(3,3);
  mycE_ch2_overtrueE_leak_1->Divide(3,3);
  mycE_ch2_overtrueE_leak_2->Divide(3,3);
  mycE_ch3_overtrueE_1->Divide(3,3);
  mycE_ch3_overtrueE_2->Divide(3,3);
  mycE_ch3_overtrueE_leak_1->Divide(3,3);
  mycE_ch3_overtrueE_leak_2->Divide(3,3);
  mycE_ch4_overtrueE_1->Divide(3,3);
  mycE_ch4_overtrueE_2->Divide(3,3);
  mycE_ch4_overtrueE_leak_1->Divide(3,3);
  mycE_ch4_overtrueE_leak_2->Divide(3,3);

  //For the lgend
  TLegend* leg = new TLegend(0.5, 0.7, 0.8, 0.9); leg->SetHeader("Energy"); leg->SetFillColor(17);
  TLegend* leg2 = new TLegend(0.5, 0.7, 0.8, 0.9); leg2->SetHeader("Energy"); leg2->SetFillColor(17);
  TLegend* leg3 = new TLegend(0.5, 0.7, 0.8, 0.9); leg3->SetHeader("Energy"); leg3->SetFillColor(17);
  TLegend* leg4 = new TLegend(0.5, 0.7, 0.8, 0.9); leg4->SetHeader("Energy"); leg4->SetFillColor(17);
  TLegend* leg5 = new TLegend(0.5, 0.7, 0.8, 0.9); leg5->SetHeader("Energy"); leg5->SetFillColor(17);
  TLegend* leg6 = new TLegend(0.5, 0.7, 0.8, 0.9); leg6->SetHeader("Energy"); leg6->SetFillColor(17);
  TLegend* leg7 = new TLegend(0.5, 0.7, 0.8, 0.9); leg7->SetHeader("Energy"); leg7->SetFillColor(17);
  TLegend* leg8 = new TLegend(0.5, 0.7, 0.8, 0.9); leg8->SetHeader("Energy"); leg8->SetFillColor(17);
  TLegend* leg9 = new TLegend(0.5, 0.7, 0.8, 0.9); leg9->SetHeader("Energy"); leg9->SetFillColor(17);
  TLegend* leg10 = new TLegend(0.5, 0.7, 0.8, 0.9); leg10->SetHeader("Energy"); leg10->SetFillColor(17);
  TLegend* leg11 = new TLegend(0.5, 0.7, 0.8, 0.9); leg11->SetHeader("Energy"); leg11->SetFillColor(17);
  TLegend* leg12 = new TLegend(0.5, 0.7, 0.8, 0.9); leg12->SetHeader("Energy"); leg12->SetFillColor(17);
  TLegend* leg13 = new TLegend(0.5, 0.7, 0.8, 0.9); leg13->SetHeader("Energy"); leg13->SetFillColor(17);
  TLegend* leg14 = new TLegend(0.5, 0.7, 0.8, 0.9); leg14->SetHeader("Energy"); leg14->SetFillColor(17);
  TLegend* leg15 = new TLegend(0.5, 0.7, 0.8, 0.9); leg15->SetHeader("Energy"); leg15->SetFillColor(17);
  TLegend* leg16 = new TLegend(0.5, 0.7, 0.8, 0.9); leg16->SetHeader("Energy"); leg16->SetFillColor(17);
  TLegend* leg17 = new TLegend(0.5, 0.7, 0.8, 0.9); leg17->SetHeader("Energy"); leg17->SetFillColor(17);
  TLegend* leg18 = new TLegend(0.5, 0.7, 0.8, 0.9); leg18->SetHeader("Energy"); leg18->SetFillColor(17);
  TLegend* leg19 = new TLegend(0.5, 0.7, 0.8, 0.9); leg19->SetHeader("Energy"); leg19->SetFillColor(17);
  TLegend* leg20 = new TLegend(0.5, 0.7, 0.8, 0.9); leg20->SetHeader("Energy"); leg20->SetFillColor(17);
  TLegend* leg21 = new TLegend(0.5, 0.7, 0.8, 0.9); leg21->SetHeader("Energy"); leg21->SetFillColor(17);
  TLegend* leg22 = new TLegend(0.5, 0.7, 0.8, 0.9); leg22->SetHeader("Energy"); leg22->SetFillColor(17);
  TLegend* leg23 = new TLegend(0.5, 0.7, 0.8, 0.9); leg23->SetHeader("Energy"); leg23->SetFillColor(17);
  TLegend* leg24 = new TLegend(0.5, 0.7, 0.8, 0.9); leg24->SetHeader("Energy"); leg24->SetFillColor(17);
  TLegend* leg25 = new TLegend(0.5, 0.7, 0.8, 0.9); leg25->SetHeader("Energy"); leg25->SetFillColor(17);
  TLegend* leg26 = new TLegend(0.5, 0.7, 0.8, 0.9); leg26->SetHeader("Energy"); leg26->SetFillColor(17);
  TLegend* leg27 = new TLegend(0.5, 0.7, 0.8, 0.9); leg27->SetHeader("Energy"); leg27->SetFillColor(17);
  TLegend* leg28 = new TLegend(0.5, 0.7, 0.8, 0.9); leg28->SetHeader("Energy"); leg28->SetFillColor(17);
  TLegend* leg29 = new TLegend(0.5, 0.7, 0.8, 0.9); leg29->SetHeader("Energy"); leg29->SetFillColor(17);
  TLegend* leg30 = new TLegend(0.5, 0.7, 0.8, 0.9); leg30->SetHeader("Energy"); leg30->SetFillColor(17);
  TLegend* leg31 = new TLegend(0.5, 0.7, 0.8, 0.9); leg31->SetHeader("Energy"); leg31->SetFillColor(17);
  TLegend* leg32 = new TLegend(0.5, 0.7, 0.8, 0.9); leg32->SetHeader("Energy"); leg32->SetFillColor(17);
  TLegend* leg33 = new TLegend(0.5, 0.7, 0.8, 0.9); leg33->SetHeader("Energy"); leg33->SetFillColor(17);
  TLegend* leg34 = new TLegend(0.5, 0.7, 0.8, 0.9); leg34->SetHeader("Energy"); leg34->SetFillColor(17);
  TLegend* leg35 = new TLegend(0.5, 0.7, 0.8, 0.9); leg35->SetHeader("Energy"); leg35->SetFillColor(17);
  TLegend* leg36 = new TLegend(0.5, 0.7, 0.8, 0.9); leg36->SetHeader("Energy"); leg36->SetFillColor(17);
  TLegend* leg37 = new TLegend(0.5, 0.7, 0.8, 0.9); leg37->SetHeader("Energy"); leg37->SetFillColor(17);
  TLegend* leg38 = new TLegend(0.5, 0.7, 0.8, 0.9); leg38->SetHeader("Energy"); leg38->SetFillColor(17);
  TLegend* leg39 = new TLegend(0.5, 0.7, 0.8, 0.9); leg39->SetHeader("Energy"); leg39->SetFillColor(17);
  TLegend* leg40 = new TLegend(0.5, 0.7, 0.8, 0.9); leg40->SetHeader("Energy"); leg40->SetFillColor(17);
  TLegend* leg41 = new TLegend(0.5, 0.7, 0.8, 0.9); leg41->SetHeader("Energy"); leg41->SetFillColor(17);
  TLegend* leg42 = new TLegend(0.5, 0.7, 0.8, 0.9); leg42->SetHeader("Energy"); leg42->SetFillColor(17);
  TLegend* leg43 = new TLegend(0.5, 0.7, 0.8, 0.9); leg43->SetHeader("Energy"); leg43->SetFillColor(17);
  TLegend* leg44 = new TLegend(0.5, 0.7, 0.8, 0.9); leg44->SetHeader("Energy"); leg44->SetFillColor(17);
  TLegend* leg45 = new TLegend(0.5, 0.7, 0.8, 0.9); leg45->SetHeader("Energy"); leg45->SetFillColor(17);
  TLegend* leg46 = new TLegend(0.5, 0.7, 0.8, 0.9); leg46->SetHeader("Energy"); leg46->SetFillColor(17);
  TLegend* leg47 = new TLegend(0.5, 0.7, 0.8, 0.9); leg47->SetHeader("Energy"); leg47->SetFillColor(17);
  TLegend* leg48 = new TLegend(0.5, 0.7, 0.8, 0.9); leg48->SetHeader("Energy"); leg48->SetFillColor(17);
  TLegend* leg49 = new TLegend(0.5, 0.7, 0.8, 0.9); leg49->SetHeader("Energy"); leg49->SetFillColor(17);
  TLegend* leg50 = new TLegend(0.5, 0.7, 0.8, 0.9); leg50->SetHeader("Energy"); leg50->SetFillColor(17);
  TLegend* leg51 = new TLegend(0.5, 0.7, 0.8, 0.9); leg51->SetHeader("Energy"); leg51->SetFillColor(17);
  TLegend* leg52 = new TLegend(0.5, 0.7, 0.8, 0.9); leg52->SetHeader("Energy"); leg52->SetFillColor(17);
  TLegend* leg53 = new TLegend(0.5, 0.7, 0.8, 0.9); leg53->SetHeader("Energy"); leg53->SetFillColor(17);
  TLegend* leg54 = new TLegend(0.5, 0.7, 0.8, 0.9); leg54->SetHeader("Energy"); leg54->SetFillColor(17);
  TLegend* leg55 = new TLegend(0.5, 0.7, 0.8, 0.9); leg55->SetHeader("Energy"); leg55->SetFillColor(17);
  TLegend* leg56 = new TLegend(0.5, 0.7, 0.8, 0.9); leg56->SetHeader("Energy"); leg56->SetFillColor(17);
  TLegend* leg57 = new TLegend(0.5, 0.7, 0.8, 0.9); leg57->SetHeader("Energy"); leg57->SetFillColor(17);
  TLegend* leg58 = new TLegend(0.5, 0.7, 0.8, 0.9); leg58->SetHeader("Energy"); leg58->SetFillColor(17);
  TLegend* leg59 = new TLegend(0.5, 0.7, 0.8, 0.9); leg59->SetHeader("Energy"); leg59->SetFillColor(17);
  TLegend* leg60 = new TLegend(0.5, 0.7, 0.8, 0.9); leg60->SetHeader("Energy"); leg60->SetFillColor(17);
  TLegend* leg61 = new TLegend(0.5, 0.7, 0.8, 0.9); leg61->SetHeader("Energy"); leg61->SetFillColor(17);
  TLegend* leg62 = new TLegend(0.5, 0.7, 0.8, 0.9); leg62->SetHeader("Energy"); leg62->SetFillColor(17);
  TLegend* leg63 = new TLegend(0.5, 0.7, 0.8, 0.9); leg63->SetHeader("Energy"); leg63->SetFillColor(17);
  TLegend* leg64 = new TLegend(0.5, 0.7, 0.8, 0.9); leg64->SetHeader("Energy"); leg64->SetFillColor(17);
  TLegend* leg65 = new TLegend(0.5, 0.7, 0.8, 0.9); leg65->SetHeader("Energy"); leg65->SetFillColor(17);
  TLegend* leg66 = new TLegend(0.5, 0.7, 0.8, 0.9); leg66->SetHeader("Energy"); leg66->SetFillColor(17);
  TLegend* leg67 = new TLegend(0.5, 0.7, 0.8, 0.9); leg67->SetHeader("Energy"); leg67->SetFillColor(17);
  TLegend* leg68 = new TLegend(0.5, 0.7, 0.8, 0.9); leg68->SetHeader("Energy"); leg68->SetFillColor(17);
  TLegend* leg69 = new TLegend(0.5, 0.7, 0.8, 0.9); leg69->SetHeader("Energy"); leg69->SetFillColor(17);
  TLegend* leg70 = new TLegend(0.5, 0.7, 0.8, 0.9); leg70->SetHeader("Energy"); leg70->SetFillColor(17);
  TLegend* leg71 = new TLegend(0.5, 0.7, 0.8, 0.9); leg71->SetHeader("Energy"); leg71->SetFillColor(17);
  TLegend* leg72 = new TLegend(0.5, 0.7, 0.8, 0.9); leg72->SetHeader("Energy"); leg72->SetFillColor(17);
  TLegend* leg73 = new TLegend(0.5, 0.7, 0.8, 0.9); leg73->SetHeader("Energy"); leg73->SetFillColor(17);
  TLegend* leg74 = new TLegend(0.5, 0.7, 0.8, 0.9); leg74->SetHeader("Energy"); leg74->SetFillColor(17);
  TLegend* leg75 = new TLegend(0.5, 0.7, 0.8, 0.9); leg75->SetHeader("Energy"); leg75->SetFillColor(17);
  TLegend* leg76 = new TLegend(0.5, 0.7, 0.8, 0.9); leg76->SetHeader("Energy"); leg76->SetFillColor(17);
  TLegend* leg77 = new TLegend(0.5, 0.7, 0.8, 0.9); leg77->SetHeader("Energy"); leg77->SetFillColor(17);
  TLegend* leg78 = new TLegend(0.5, 0.7, 0.8, 0.9); leg78->SetHeader("Energy"); leg78->SetFillColor(17);
  TLegend* leg79 = new TLegend(0.5, 0.7, 0.8, 0.9); leg79->SetHeader("Energy"); leg79->SetFillColor(17);
  TLegend* leg80 = new TLegend(0.5, 0.7, 0.8, 0.9); leg80->SetHeader("Energy"); leg80->SetFillColor(17);
  TLegend* leg81 = new TLegend(0.5, 0.7, 0.8, 0.9); leg81->SetHeader("Energy"); leg81->SetFillColor(17);
  TLegend* leg82 = new TLegend(0.5, 0.7, 0.8, 0.9); leg82->SetHeader("Energy"); leg82->SetFillColor(17);
  TLegend* leg83 = new TLegend(0.5, 0.7, 0.8, 0.9); leg83->SetHeader("Energy"); leg83->SetFillColor(17);
  TLegend* leg84 = new TLegend(0.5, 0.7, 0.8, 0.9); leg84->SetHeader("Energy"); leg84->SetFillColor(17);
  TLegend* leg85 = new TLegend(0.5, 0.7, 0.8, 0.9); leg85->SetHeader("Energy"); leg85->SetFillColor(17);
  TLegend* leg86 = new TLegend(0.5, 0.7, 0.8, 0.9); leg86->SetHeader("Energy"); leg86->SetFillColor(17);
  TLegend* leg87 = new TLegend(0.5, 0.7, 0.8, 0.9); leg87->SetHeader("Energy"); leg87->SetFillColor(17);
  TLegend* leg88 = new TLegend(0.5, 0.7, 0.8, 0.9); leg88->SetHeader("Energy"); leg88->SetFillColor(17);

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
  TGraphErrors * recoE_rawvsbeamenergy;
  TGraphErrors * recoE_raw_leakvsbeamenergy;
  TGraphErrors * recoE_valvsbeamenergy;
  TGraphErrors * recoEvsbeamenergy;
  TGraphErrors * recoE_totsfvsbeamenergy;
  TGraphErrors * recoE_totsf_leakvsbeamenergy;
  TGraphErrors * recoE_totsf_dedxweightedvsbeamenergy;
  TGraphErrors * recoE_totsf_dedxweighted_leakvsbeamenergy;
  TGraphErrors * recoE_totsf_secvsbeamenergy;
  TGraphErrors * recoE_totsf_sec_leakvsbeamenergy;
  TGraphErrors * recoE_totsf_sec_dedxweightedvsbeamenergy;
  TGraphErrors * recoE_totsf_sec_dedxweighted_leakvsbeamenergy;
  TGraphErrors * recoE_totsfshmax_secvsbeamenergy;
  TGraphErrors * recoE_totsfshmax_sec_leakvsbeamenergy;
  TGraphErrors * recoE_totsfshmax_sec_dedxweightedvsbeamenergy;
  TGraphErrors * recoE_totsfshmax_sec_dedxweighted_leakvsbeamenergy;
  TGraphErrors * recoE_totsf_dedxsecvsbeamenergy;
  TGraphErrors * recoE_totsf_dedxsec_leakvsbeamenergy;
  TGraphErrors * recoE_totsfshmax_dedxsecvsbeamenergy;
  TGraphErrors * recoE_totsfshmax_dedxsec_leakvsbeamenergy;
  TGraphErrors * recoE_shmaxvsbeamenergy;
  TGraphErrors * recoE_shmax_leakvsbeamenergy;
  // TGraphErrors * recoE_totsfovertrueEvsbeamenergy;
  // TGraphErrors * recoE_totsf_leakovertrueEvsbeamenergy;
  TGraphErrors * recoE_sfnormvsbeamenergy;
  TGraphErrors * recoE_sfnorm_leakvsbeamenergy;
  TGraphErrors * recoE_totsfnormvsbeamenergy;
  TGraphErrors * recoE_totsfnorm_leakvsbeamenergy;
  TGraphErrors * recoE_totsfnorm_dedxweightedvsbeamenergy;
  TGraphErrors * recoE_totsfnorm_dedxweighted_leakvsbeamenergy;
  TGraphErrors * recoE_leakvsbeamenergy;
  TGraphErrors * recoE_leak_offsetvsbeamenergy;
  TGraphErrors * recoE_4sectionsvsbeamenergy;
  TGraphErrors * recoE_4sections_leakvsbeamenergy;
  TGraphErrors * recoE_4sections_ch1vsbeamenergy;
  TGraphErrors * recoE_4sections_ch1_leakvsbeamenergy;
  TGraphErrors * recoE_4sections_ch2vsbeamenergy;
  TGraphErrors * recoE_4sections_ch2_leakvsbeamenergy;
  TGraphErrors * recoE_4sections_ch3vsbeamenergy;
  TGraphErrors * recoE_4sections_ch3_leakvsbeamenergy;
  TGraphErrors * recoE_4sections_ch4vsbeamenergy;
  TGraphErrors * recoE_4sections_ch4_leakvsbeamenergy;
  TGraphErrors * recoE_ch1vsbeamenergy;
  TGraphErrors * recoE_ch1_leakvsbeamenergy;
  TGraphErrors * recoE_ch2vsbeamenergy;
  TGraphErrors * recoE_ch2_leakvsbeamenergy;
  TGraphErrors * recoE_ch3vsbeamenergy;
  TGraphErrors * recoE_ch3_leakvsbeamenergy;
  TGraphErrors * recoE_ch4vsbeamenergy;
  TGraphErrors * recoE_ch4_leakvsbeamenergy;
  TGraphErrors * recoE_ch1_overtrueEvsbeamenergy;
  TGraphErrors * recoE_ch1_overtrueE_leakvsbeamenergy;
  TGraphErrors * recoE_ch2_overtrueEvsbeamenergy;
  TGraphErrors * recoE_ch2_overtrueE_leakvsbeamenergy;
  TGraphErrors * recoE_ch3_overtrueEvsbeamenergy;
  TGraphErrors * recoE_ch3_overtrueE_leakvsbeamenergy;
  TGraphErrors * recoE_ch4_overtrueEvsbeamenergy;
  TGraphErrors * recoE_ch4_overtrueE_leakvsbeamenergy;
  TGraphErrors * recoE_corrvsbeamenergy;
  TGraphErrors * recoEovertrueEvsbeamenergy;
  TGraphErrors * recoEovertrueE_leakvsbeamenergy;
  TGraphErrors * recoEovertrueE_leak_offsetvsbeamenergy;
  TGraphErrors * recoEovertrueE_totsfvsbeamenergy;
  TGraphErrors * recoEovertrueE_totsf_leakvsbeamenergy;
  TGraphErrors * recoEovertrueE_totsf_dedxweightedvsbeamenergy;
  TGraphErrors * recoEovertrueE_totsf_dedxweighted_leakvsbeamenergy;
  TGraphErrors * recoEovertrueE_totsf_secvsbeamenergy;
  TGraphErrors * recoEovertrueE_totsf_sec_leakvsbeamenergy;
  TGraphErrors * recoEovertrueE_totsf_sec_dedxweightedvsbeamenergy;
  TGraphErrors * recoEovertrueE_totsf_sec_dedxweighted_leakvsbeamenergy;
  TGraphErrors * recoEovertrueE_totsfshmax_secvsbeamenergy;
  TGraphErrors * recoEovertrueE_totsfshmax_sec_leakvsbeamenergy;
  TGraphErrors * recoEovertrueE_totsfshmax_sec_dedxweightedvsbeamenergy;
  TGraphErrors * recoEovertrueE_totsfshmax_sec_dedxweighted_leakvsbeamenergy;
  TGraphErrors * recoEovertrueE_totsf_dedxsecvsbeamenergy;
  TGraphErrors * recoEovertrueE_totsf_dedxsec_leakvsbeamenergy;
  TGraphErrors * recoEovertrueE_totsfshmax_dedxsecvsbeamenergy;
  TGraphErrors * recoEovertrueE_totsfshmax_dedxsec_leakvsbeamenergy;
  TGraphErrors * recoEovertrueE_totsfnormvsbeamenergy;
  TGraphErrors * recoEovertrueE_totsfnorm_leakvsbeamenergy;
  TGraphErrors * recoEovertrueE_totsfnorm_dedxweightedvsbeamenergy;
  TGraphErrors * recoEovertrueE_totsfnorm_dedxweighted_leakvsbeamenergy;
  TGraphErrors * recoEovertrueE_sfnormvsbeamenergy;
  TGraphErrors * recoEovertrueE_sfnorm_leakvsbeamenergy;
  TGraphErrors * recoEovertrueE_shmaxvsbeamenergy;
  TGraphErrors * recoEovertrueE_shmax_leakvsbeamenergy;
  TGraphErrors * recoE_valovertrueEvsbeamenergy;
  TGraphErrors * recoE_corrovertrueEvsbeamenergy;
  TGraphErrors * resoFit;
  TGraphErrors * resoFit_totsf;
  TGraphErrors * resoFit_totsf_leak;
  TGraphErrors * resoFit_totsf_dedxweighted;
  TGraphErrors * resoFit_totsf_dedxweighted_leak;
  TGraphErrors * resoFit_totsf_sec;
  TGraphErrors * resoFit_totsf_sec_leak;
  TGraphErrors * resoFit_totsf_sec_dedxweighted;
  TGraphErrors * resoFit_totsf_sec_dedxweighted_leak;
  TGraphErrors * resoFit_totsfshmax_sec;
  TGraphErrors * resoFit_totsfshmax_sec_leak;
  TGraphErrors * resoFit_totsfshmax_sec_dedxweighted;
  TGraphErrors * resoFit_totsfshmax_sec_dedxweighted_leak;
  TGraphErrors * resoFit_totsf_dedxsec;
  TGraphErrors * resoFit_totsf_dedxsec_leak;
  TGraphErrors * resoFit_totsfshmax_dedxsec;
  TGraphErrors * resoFit_totsfshmax_dedxsec_leak;
  TGraphErrors * resoFit_shmax;
  TGraphErrors * resoFit_shmax_leak;
  TGraphErrors * resoFit_sfnorm;
  TGraphErrors * resoFit_sfnorm_leak;
  TGraphErrors * resoFit_leak;
  TGraphErrors * resoFit_totsfnorm;
  TGraphErrors * resoFit_totsfnorm_leak;
  TGraphErrors * resoFit_totsfnorm_dedxweighted;
  TGraphErrors * resoFit_totsfnorm_dedxweighted_leak;
  TGraphErrors * resoFit_4sections;
  TGraphErrors * resoFit_4sections_leak;
  TGraphErrors * resoFit_4sections_ch1;
  TGraphErrors * resoFit_4sections_ch1_leak;
  TGraphErrors * resoFit_4sections_ch2;
  TGraphErrors * resoFit_4sections_ch2_leak;
  TGraphErrors * resoFit_4sections_ch3;
  TGraphErrors * resoFit_4sections_ch3_leak;
  TGraphErrors * resoFit_4sections_ch4;
  TGraphErrors * resoFit_4sections_ch4_leak;
  TGraphErrors * resoFit_ch1;
  TGraphErrors * resoFit_ch1_leak;
  TGraphErrors * resoFit_ch2;
  TGraphErrors * resoFit_ch2_leak;
  TGraphErrors * resoFit_ch3;
  TGraphErrors * resoFit_ch3_leak;
  TGraphErrors * resoFit_ch4;
  TGraphErrors * resoFit_ch4_leak;
  TGraphErrors * resoFit_leak_offset;
  TGraphErrors * resoFit_raw;
  TGraphErrors * resoFit_raw_leak;
  TGraphErrors * resoFit_corr;
  TGraphErrors * resoFit_val;
  TGraphErrors *recoEvsbeam = new TGraphErrors();
  recoE_rawvsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoE_rawvsbeamenergy" );
  recoE_raw_leakvsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoE_raw_leakvsbeamenergy" );
  recoE_valvsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoE_valvsbeamenergy" );
  recoE_totsfvsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoE_totsfvsbeamenergy" );
  recoE_totsf_leakvsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoE_totsf_leakvsbeamenergy" );
  recoE_totsf_dedxweightedvsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoE_totsf_dedxweightedvsbeamenergy" );
  recoE_totsf_dedxweighted_leakvsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoE_totsf_dedxweighted_leakvsbeamenergy" );
  recoE_totsf_secvsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoE_totsf_secvsbeamenergy" );
  recoE_totsf_sec_leakvsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoE_totsf_sec_leakvsbeamenergy" );
  recoE_totsf_sec_dedxweightedvsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoE_totsf_sec_dedxweightedvsbeamenergy" );
  recoE_totsf_sec_dedxweighted_leakvsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoE_totsf_sec_dedxweighted_leakvsbeamenergy" );
  recoE_totsfshmax_secvsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoE_totsfshmax_secvsbeamenergy" );
  recoE_totsfshmax_sec_leakvsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoE_totsfshmax_sec_leakvsbeamenergy" );
  recoE_totsfshmax_sec_dedxweightedvsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoE_totsfshmax_sec_dedxweightedvsbeamenergy" );
  recoE_totsfshmax_sec_dedxweighted_leakvsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoE_totsfshmax_sec_dedxweighted_leakvsbeamenergy" );
  recoE_totsf_dedxsecvsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoE_totsf_dedxsecvsbeamenergy" );
  recoE_totsf_dedxsec_leakvsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoE_totsf_dedxsec_leakvsbeamenergy" );
  recoE_totsfshmax_dedxsecvsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoE_totsfshmax_dedxsecvsbeamenergy" );
  recoE_totsfshmax_dedxsec_leakvsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoE_totsfshmax_dedxsec_leakvsbeamenergy" );
  recoEovertrueE_totsfvsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoEovertrueE_totsfvsbeamenergy" );
  recoEovertrueE_totsf_leakvsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoEovertrueE_totsf_leakvsbeamenergy" );
  recoEovertrueE_totsf_dedxweightedvsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoEovertrueE_totsf_dedxweightedvsbeamenergy" );
  recoEovertrueE_totsf_dedxweighted_leakvsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoEovertrueE_totsf_dedxweighted_leakvsbeamenergy" );
  recoEovertrueE_totsf_secvsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoEovertrueE_totsf_secvsbeamenergy" );
  recoEovertrueE_totsf_sec_leakvsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoEovertrueE_totsf_sec_leakvsbeamenergy" );
  recoEovertrueE_totsf_sec_dedxweightedvsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoEovertrueE_totsf_sec_dedxweightedvsbeamenergy" );
  recoEovertrueE_totsf_sec_dedxweighted_leakvsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoEovertrueE_totsf_sec_dedxweighted_leakvsbeamenergy" );
  recoEovertrueE_totsfshmax_secvsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoEovertrueE_totsfshmax_secvsbeamenergy" );
  recoEovertrueE_totsfshmax_sec_leakvsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoEovertrueE_totsfshmax_sec_leakvsbeamenergy" );
  recoEovertrueE_totsfshmax_sec_dedxweightedvsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoEovertrueE_totsfshmax_sec_dedxweightedvsbeamenergy" );
  recoEovertrueE_totsfshmax_sec_dedxweighted_leakvsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoEovertrueE_totsfshmax_sec_dedxweighted_leakvsbeamenergy" );
  recoEovertrueE_totsf_dedxsecvsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoEovertrueE_totsf_dedxsecvsbeamenergy" );
  recoEovertrueE_totsf_dedxsec_leakvsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoEovertrueE_totsf_dedxsec_leakvsbeamenergy" );
  recoEovertrueE_totsfshmax_dedxsecvsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoEovertrueE_totsfshmax_dedxsecvsbeamenergy" );
  recoEovertrueE_totsfshmax_dedxsec_leakvsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoEovertrueE_totsfshmax_dedxsec_leakvsbeamenergy" );
  recoEovertrueE_totsfnormvsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoEovertrueE_totsfnormvsbeamenergy" );
  recoEovertrueE_totsfnorm_leakvsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoEovertrueE_totsfnorm_leakvsbeamenergy" );
  recoEovertrueE_totsfnorm_dedxweightedvsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoEovertrueE_totsfnorm_dedxweightedvsbeamenergy" );
  recoEovertrueE_totsfnorm_dedxweighted_leakvsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoEovertrueE_totsfnorm_dedxweighted_leakvsbeamenergy" );
  recoEovertrueE_sfnormvsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoEovertrueE_sfnormvsbeamenergy" );
  recoEovertrueE_sfnorm_leakvsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoEovertrueE_sfnorm_leakvsbeamenergy" );
  recoE_shmaxvsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoE_shmaxvsbeamenergy" );
  recoE_shmax_leakvsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoE_shmax_leakvsbeamenergy" );
  recoEovertrueE_shmaxvsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoEovertrueE_shmaxvsbeamenergy" );
  recoEovertrueE_shmax_leakvsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoEovertrueE_shmax_leakvsbeamenergy" );
  recoE_sfnormvsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoE_sfnormvsbeamenergy" );
  recoE_sfnorm_leakvsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoE_sfnorm_leakvsbeamenergy" );
  recoE_totsfnormvsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoE_totsfnormvsbeamenergy" );
  recoE_totsfnorm_leakvsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoE_totsfnorm_leakvsbeamenergy" );
  recoE_totsfnorm_dedxweightedvsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoE_totsfnorm_dedxweightedvsbeamenergy" );
  recoE_totsfnorm_dedxweighted_leakvsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoE_totsfnorm_dedxweighted_leakvsbeamenergy" );
  recoEvsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoEvsbeamenergy" );
  recoE_leakvsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoE_leakvsbeamenergy" );
  recoE_leak_offsetvsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoE_leak_offsetvsbeamenergy" );
  recoE_4sectionsvsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoE_4sectionsvsbeamenergy" );
  recoE_4sections_leakvsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoE_4sections_leakvsbeamenergy" );
  recoE_4sections_ch1vsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoE_4sections_ch1vsbeamenergy" );
  recoE_4sections_ch1_leakvsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoE_4sections_ch1_leakvsbeamenergy" );
  recoE_4sections_ch2vsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoE_4sections_ch2vsbeamenergy" );
  recoE_4sections_ch2_leakvsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoE_4sections_ch2_leakvsbeamenergy" );
  recoE_4sections_ch3vsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoE_4sections_ch3vsbeamenergy" );
  recoE_4sections_ch3_leakvsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoE_4sections_ch3_leakvsbeamenergy" );
  recoE_4sections_ch4vsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoE_4sections_ch4vsbeamenergy" );
  recoE_4sections_ch4_leakvsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoE_4sections_ch4_leakvsbeamenergy" );
  recoE_ch1vsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoE_ch1vsbeamenergy" );
  recoE_ch1_leakvsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoE_ch1_leakvsbeamenergy" );
  recoE_ch2vsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoE_ch2vsbeamenergy" );
  recoE_ch2_leakvsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoE_ch2_leakvsbeamenergy" );
  recoE_ch3vsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoE_ch3vsbeamenergy" );
  recoE_ch3_leakvsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoE_ch3_leakvsbeamenergy" );
  recoE_ch4vsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoE_ch4vsbeamenergy" );
  recoE_ch4_leakvsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoE_ch4_leakvsbeamenergy" );
  recoE_ch1_overtrueEvsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoE_ch1_overtrueEvsbeamenergy" );
  recoE_ch1_overtrueE_leakvsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoE_ch1_overtrueE_leakvsbeamenergy" );
  recoE_ch2_overtrueEvsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoE_ch2_overtrueEvsbeamenergy" );
  recoE_ch2_overtrueE_leakvsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoE_ch2_overtrueE_leakvsbeamenergy" );
  recoE_ch3_overtrueEvsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoE_ch3_overtrueEvsbeamenergy" );
  recoE_ch3_overtrueE_leakvsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoE_ch3_overtrueE_leakvsbeamenergy" );
  recoE_ch4_overtrueEvsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoE_ch4_overtrueEvsbeamenergy" );
  recoE_ch4_overtrueE_leakvsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoE_ch4_overtrueE_leakvsbeamenergy" );
  if(upstream){
    recoE_corrvsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoE_corrvsbeamenergy" );
    recoE_corrovertrueEvsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoE_corrovertrueEvsbeamenergy" );
  }
  recoEovertrueEvsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoEovertrueEvsbeamenergy" );
  recoEovertrueE_leakvsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoEovertrueE_leakvsbeamenergy" );
  recoEovertrueE_leak_offsetvsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoEovertrueE_leak_offsetvsbeamenergy" );
  recoE_valovertrueEvsbeamenergy  = (TGraphErrors *) recoEvsbeam->Clone( "recoE_valovertrueEvsbeamenergy" );
  resoFit = (TGraphErrors *) recoEvsbeam->Clone( "resoFit" );
  resoFit_totsf = (TGraphErrors *) recoEvsbeam->Clone( "resoFit_totsf" );
  resoFit_totsf_leak = (TGraphErrors *) recoEvsbeam->Clone( "resoFit_totsf_leak" );
  resoFit_totsf_dedxweighted = (TGraphErrors *) recoEvsbeam->Clone( "resoFit_totsf_dedxweighted" );
  resoFit_totsf_dedxweighted_leak = (TGraphErrors *) recoEvsbeam->Clone( "resoFit_totsf_dedxweighted_leak" );
  resoFit_totsf_sec = (TGraphErrors *) recoEvsbeam->Clone( "resoFit_totsf_sec" );
  resoFit_totsf_sec_leak = (TGraphErrors *) recoEvsbeam->Clone( "resoFit_totsf_sec_leak" );
  resoFit_totsf_sec_dedxweighted = (TGraphErrors *) recoEvsbeam->Clone( "resoFit_totsf_sec_dedxweighted" );
  resoFit_totsf_sec_dedxweighted_leak = (TGraphErrors *) recoEvsbeam->Clone( "resoFit_totsf_sec_dedxweighted_leak" );
  resoFit_totsfshmax_sec = (TGraphErrors *) recoEvsbeam->Clone( "resoFit_totsfshmax_sec" );
  resoFit_totsfshmax_sec_leak = (TGraphErrors *) recoEvsbeam->Clone( "resoFit_totsfshmax_sec_leak" );
  resoFit_totsfshmax_sec_dedxweighted = (TGraphErrors *) recoEvsbeam->Clone( "resoFit_totsfshmax_sec_dedxweighted" );
  resoFit_totsfshmax_sec_dedxweighted_leak = (TGraphErrors *) recoEvsbeam->Clone( "resoFit_totsfshmax_sec_dedxweighted_leak" );
  resoFit_totsf_dedxsec = (TGraphErrors *) recoEvsbeam->Clone( "resoFit_totsf_dedxsec" );
  resoFit_totsf_dedxsec_leak = (TGraphErrors *) recoEvsbeam->Clone( "resoFit_totsf_dedxsec_leak" );
  resoFit_totsfshmax_dedxsec = (TGraphErrors *) recoEvsbeam->Clone( "resoFit_totsfshmax_dedxsec" );
  resoFit_totsfshmax_dedxsec_leak = (TGraphErrors *) recoEvsbeam->Clone( "resoFit_totsfshmax_dedxsec_leak" );
  resoFit_shmax = (TGraphErrors *) recoEvsbeam->Clone( "resoFit_shmax" );
  resoFit_shmax_leak = (TGraphErrors *) recoEvsbeam->Clone( "resoFit_shmax_leak" );
  resoFit_sfnorm = (TGraphErrors *) recoEvsbeam->Clone( "resoFit_sfnorm" );
  resoFit_sfnorm_leak = (TGraphErrors *) recoEvsbeam->Clone( "resoFit_sfnorm_leak" );
  resoFit_totsfnorm = (TGraphErrors *) recoEvsbeam->Clone( "resoFit_totsfnorm" );
  resoFit_totsfnorm_leak = (TGraphErrors *) recoEvsbeam->Clone( "resoFit_totsfnorm_leak" );
  resoFit_totsfnorm_dedxweighted = (TGraphErrors *) recoEvsbeam->Clone( "resoFit_totsfnorm_dedxweighted" );
  resoFit_totsfnorm_dedxweighted_leak = (TGraphErrors *) recoEvsbeam->Clone( "resoFit_totsfnorm_dedxweighted_leak" );
  resoFit_leak = (TGraphErrors *) recoEvsbeam->Clone( "resoFit_leak" );
  resoFit_leak_offset = (TGraphErrors *) recoEvsbeam->Clone( "resoFit_leak_offset" );
  resoFit_4sections = (TGraphErrors *) recoEvsbeam->Clone( "resoFit_4sections" );
  resoFit_4sections_leak = (TGraphErrors *) recoEvsbeam->Clone( "resoFit_4sections_leak" );
  resoFit_4sections_ch1 = (TGraphErrors *) recoEvsbeam->Clone( "resoFit_4sections_ch1" );
  resoFit_4sections_ch1_leak = (TGraphErrors *) recoEvsbeam->Clone( "resoFit_4sections_ch1_leak" );
  resoFit_4sections_ch2 = (TGraphErrors *) recoEvsbeam->Clone( "resoFit_4sections_ch2" );
  resoFit_4sections_ch2_leak = (TGraphErrors *) recoEvsbeam->Clone( "resoFit_4sections_ch2_leak" );
  resoFit_4sections_ch3 = (TGraphErrors *) recoEvsbeam->Clone( "resoFit_4sections_ch3" );
  resoFit_4sections_ch3_leak = (TGraphErrors *) recoEvsbeam->Clone( "resoFit_4sections_ch3_leak" );
  resoFit_4sections_ch4 = (TGraphErrors *) recoEvsbeam->Clone( "resoFit_4sections_ch4" );
  resoFit_4sections_ch4_leak = (TGraphErrors *) recoEvsbeam->Clone( "resoFit_4sections_ch4_leak" );
  resoFit_ch1 = (TGraphErrors *) recoEvsbeam->Clone( "resoFit_ch1" );
  resoFit_ch1_leak = (TGraphErrors *) recoEvsbeam->Clone( "resoFit_ch1_leak" );
  resoFit_ch2 = (TGraphErrors *) recoEvsbeam->Clone( "resoFit_ch2" );
  resoFit_ch2_leak = (TGraphErrors *) recoEvsbeam->Clone( "resoFit_ch2_leak" );
  resoFit_ch3 = (TGraphErrors *) recoEvsbeam->Clone( "resoFit_ch3" );
  resoFit_ch3_leak = (TGraphErrors *) recoEvsbeam->Clone( "resoFit_ch3_leak" );
  resoFit_ch4 = (TGraphErrors *) recoEvsbeam->Clone( "resoFit_ch4" );
  resoFit_ch4_leak = (TGraphErrors *) recoEvsbeam->Clone( "resoFit_ch4_leak" );
  resoFit_raw = (TGraphErrors *) recoEvsbeam->Clone( "resoFit_raw" );
  resoFit_raw_leak = (TGraphErrors *) recoEvsbeam->Clone( "resoFit_raw_leak" );
  if(upstream){
    resoFit_corr = (TGraphErrors *) recoEvsbeam->Clone( "resoFit_corr" );
  }
  resoFit_val = (TGraphErrors *) recoEvsbeam->Clone( "resoFit_val" );
  // TGraphErrors *gr_fit;
  // //For accessing the fit results
  // TF1 *fit[numberofenergies];

  //======================================================================================
  //The files that we will store the results of this analysis for the combined plot
  TString res_com[numberofconfigs];
  TFile* results_com[numberofconfigs];
  for (int k=0; k<numberofconfigs; k++){
    TString fp = "SiWEcal_" + IntToString(configs[k]);	
    res_com[k] = fp + "_combinedplots_recoE.root";
    //std::cout << res[k] << std::endl;
  }
  TH1F *hist,*hist2,*hist3,*hist4,*hist5,*hist6,*hist7,*hist8,*hist9,*hist10,*hist11,*hist12,*hist13,*hist14,*hist15,*hist16,*hist17,*hist18,*hist19,*hist20,*hist21,*hist22,*hist23,*hist24,*hist25,*hist26,*hist27,*hist28,*hist29,*hist30,*hist31,*hist32,*hist33,*hist34,*hist35,*hist36,*hist37,*hist38,*hist39,*hist40,*hist41,*hist42,*hist43,*hist44,*hist45,*hist46,*hist47,*hist48,*hist49,*hist50,*hist51,*hist52,*hist53,*hist54,*hist55,*hist56,*hist57,*hist58,*hist59,*hist60,*hist61,*hist62,*hist63,*hist64,*hist65,*hist66,*hist67,*hist68,*hist69,*hist70,*hist71,*hist72,*hist73,*hist74,*hist75,*hist76,*hist77,*hist78,*hist79,*hist80,*hist81,*hist82,*hist83,*hist84,*hist85,*hist86,*hist87;
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  int color = 1; 
  FitResult fitres;
  FitResultResolution fitresol, fitresol_totsf, fitresol_totsf_leak, fitresol_totsf_dedxweighted, fitresol_totsf_dedxweighted_leak, fitresol_totsf_sec, fitresol_totsf_sec_leak, fitresol_totsf_sec_dedxweighted, fitresol_totsf_sec_dedxweighted_leak,fitresol_totsfshmax_sec, fitresol_totsfshmax_sec_leak,fitresol_totsfshmax_sec_dedxweighted, fitresol_totsfshmax_sec_dedxweighted_leak, fitresol_totsf_dedxsec, fitresol_totsf_dedxsec_leak,fitresol_totsfshmax_dedxsec, fitresol_totsfshmax_dedxsec_leak,fitresol_shmax, fitresol_shmax_leak, fitresol_sfnorm, fitresol_sfnorm_leak, fitresol_totsfnorm, fitresol_totsfnorm_leak, fitresol_totsfnorm_dedxweighted, fitresol_totsfnorm_dedxweighted_leak, fitresol_leak, fitresol_leak_offset, fitresol_raw, fitresol_raw_leak, fitresol_corr, fitresol_val,fitresol_4sections,fitresol_4sections_leak,fitresol_4sections_ch1,fitresol_4sections_ch1_leak,fitresol_4sections_ch2,fitresol_4sections_ch2_leak,fitresol_4sections_ch3,fitresol_4sections_ch3_leak,fitresol_4sections_ch4,fitresol_4sections_ch4_leak,fitresol_ch1,fitresol_ch1_leak,fitresol_ch2,fitresol_ch2_leak,fitresol_ch3,fitresol_ch3_leak,fitresol_ch4,fitresol_ch4_leak;

  TString titleofplot1,titleofplot2,titleofplot; 
  for (int k=0; k<numberofconfigs; k++){
    results[k]= TFile::Open(res[k],"read");
    std::cout << "Results file " << res[k] << std::endl;
    
    //======================================================================================
    //Loop on energies
    for (int l=0; l<numberofenergies; l++){
      
      //-------------------------------------------------------------------------------------------------
      titleofplot1 = "Raw Reconstructed energy for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      hist5 = (TH1F*)results[k]->Get( ("h_recoE_raw_" + IntToString(particleenergies[l])).c_str());
      // histk->GetYaxis()->SetRangeUser(0.0001,100000.);
      buildplot(c9, hist5, leg6, titleofplot, "Raw E_{reco} (MIPs)", "Events/10 MIP", procNa[l], l);
      std::cout << "Working on plot: " << hist5->GetName() << std::endl;
      //-------------------------------------------------------------------------------------------------
      titleofplot1 = "Raw Reconstructed energy adding leakage for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      hist9 = (TH1F*)results[k]->Get( ("h_recoE_raw_leak_" + IntToString(particleenergies[l])).c_str());
      buildplot(c29, hist9, leg10, titleofplot, "Raw E_{reco} (MIPs)", "Events/50 MIP", procNa[l], l);
      std::cout << "Working on plot: " << hist9->GetName() << std::endl;
      //-------------------------------------------------------------------------------------------------
      titleofplot1 = "Reconstructed energy using dEdx for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      hist6 = (TH1F*)results[k]->Get( ("h_recoE_val_" + IntToString(particleenergies[l])).c_str());
      buildplot(c18, hist6, leg7, titleofplot, "E_{reco} (GeV)", "Events/0.1 GeV", procNa[l], l);
      std::cout << "Working on plot: " << hist6->GetName() << std::endl;
      //-------------------------------------------------------------------------------------------------
      titleofplot1 = "Reconstructed energy for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      hist = (TH1F*)results[k]->Get( ("h_recoE_" + IntToString(particleenergies[l])).c_str());
      buildplot(c1, hist, leg, titleofplot, "Reconstructed energy (GeV)", "Events/0.1 GeV", procNa[l], l);
      std::cout << "Working on plot: " << hist->GetName() << std::endl;
      //-------------------------------------------------------------------------------------------------
      //Including leakage
      titleofplot1 = "Reconstructed energy including leakage for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      hist10 = (TH1F*)results[k]->Get( ("h_recoE_leak_" + IntToString(particleenergies[l])).c_str());
      buildplot(c34, hist10, leg11, titleofplot, "Reconstructed energy (GeV)", "Events/0.1 GeV", procNa[l], l);
      std::cout << "Working on plot: " << hist10->GetName() << std::endl;
      //-------------------------------------------------------------------------------------------------
      //Including leakage and offset
      titleofplot1 = "Reconstructed energy including leakage and offset for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      hist11 = (TH1F*)results[k]->Get( ("h_recoE_leak_offset_" + IntToString(particleenergies[l])).c_str());
      buildplot(c39, hist11, leg12, titleofplot, "Reconstructed energy (GeV)", "Events/0.1 GeV", procNa[l], l);
      std::cout << "Working on plot: " << hist11->GetName() << std::endl;
      //-------------------------------------------------------------------------------------------------
      //4 sections is here
      titleofplot1 = "Reconstructed energy using 4 sections for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      hist16 = (TH1F*)results[k]->Get( ("h_recoE_4sections_" + IntToString(particleenergies[l])).c_str());
      buildplot(c66, hist16, leg17, titleofplot, "Reconstructed energy (GeV)", "Events/0.1 GeV", procNa[l], l);
      std::cout << "Working on plot: " << hist16->GetName() << std::endl;
      //-------------------------------------------------------------------------------------------------
      //4 sections is here including leakage
      titleofplot1 = "Reconstructed energy using 4 sections including leakage for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      hist17 = (TH1F*)results[k]->Get( ("h_recoE_4sections_leak_" + IntToString(particleenergies[l])).c_str());
      buildplot(c67, hist17, leg18, titleofplot, "Reconstructed energy (GeV)", "Events/0.1 GeV", procNa[l], l);
      std::cout << "Working on plot: " << hist17->GetName() << std::endl;
      //-------------------------------------------------------------------------------------------------
      //4 sections is here: Cheat 1
      titleofplot1 = "Reconstructed energy using 4 sections but cheat section 1 for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      hist18 = (TH1F*)results[k]->Get( ("h_recoE_4sections_ch1_" + IntToString(particleenergies[l])).c_str());
      buildplot(c111, hist18, leg19, titleofplot, "Reconstructed energy (GeV)", "Events/0.1 GeV", procNa[l], l);
      std::cout << "Working on plot: " << hist18->GetName() << std::endl;
      //-------------------------------------------------------------------------------------------------
      //4 sections is here including leakage: Cheat 1
      titleofplot1 = "Reconstructed energy using 4 sections but cheat section 1 including leakage for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      hist19 = (TH1F*)results[k]->Get( ("h_recoE_4sections_ch1_leak_" + IntToString(particleenergies[l])).c_str());
      buildplot(c112, hist19, leg20, titleofplot, "Reconstructed energy (GeV)", "Events/0.1 GeV", procNa[l], l);
      std::cout << "Working on plot: " << hist19->GetName() << std::endl;
      //-------------------------------------------------------------------------------------------------
      //4 sections is here: Cheat 2
      titleofplot1 = "Reconstructed energy using 4 sections but cheat section 2 for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      hist20 = (TH1F*)results[k]->Get( ("h_recoE_4sections_ch2_" + IntToString(particleenergies[l])).c_str());
      buildplot(c113, hist20, leg21, titleofplot, "Reconstructed energy (GeV)", "Events/0.1 GeV", procNa[l], l);
      std::cout << "Working on plot: " << hist20->GetName() << std::endl;
      //-------------------------------------------------------------------------------------------------
      //4 sections is here including leakage: Cheat 2
      titleofplot1 = "Reconstructed energy using 4 sections but cheat section 2 including leakage for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      hist21 = (TH1F*)results[k]->Get( ("h_recoE_4sections_ch2_leak_" + IntToString(particleenergies[l])).c_str());
      buildplot(c114, hist21, leg22, titleofplot, "Reconstructed energy (GeV)", "Events/0.1 GeV", procNa[l], l);
      std::cout << "Working on plot: " << hist21->GetName() << std::endl;
      //-------------------------------------------------------------------------------------------------
      //4 sections is here: Cheat 3
      titleofplot1 = "Reconstructed energy using 4 sections but cheat section 3 for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      hist22 = (TH1F*)results[k]->Get( ("h_recoE_4sections_ch3_" + IntToString(particleenergies[l])).c_str());
      buildplot(c115, hist22, leg23, titleofplot, "Reconstructed energy (GeV)", "Events/0.1 GeV", procNa[l], l);
      std::cout << "Working on plot: " << hist22->GetName() << std::endl;
      //-------------------------------------------------------------------------------------------------
      //4 sections is here including leakage: Cheat 3
      titleofplot1 = "Reconstructed energy using 4 sections but cheat section 3 including leakage for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      hist23 = (TH1F*)results[k]->Get( ("h_recoE_4sections_ch3_leak_" + IntToString(particleenergies[l])).c_str());
      buildplot(c116, hist23, leg24, titleofplot, "Reconstructed energy (GeV)", "Events/0.1 GeV", procNa[l], l);
      std::cout << "Working on plot: " << hist23->GetName() << std::endl;
      //-------------------------------------------------------------------------------------------------
      //4 sections is here: Cheat 4
      titleofplot1 = "Reconstructed energy using 4 sections but cheat section 4 for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      hist24 = (TH1F*)results[k]->Get( ("h_recoE_4sections_ch4_" + IntToString(particleenergies[l])).c_str());
      buildplot(c117, hist24, leg25, titleofplot, "Reconstructed energy (GeV)", "Events/0.1 GeV", procNa[l], l);
      std::cout << "Working on plot: " << hist24->GetName() << std::endl;
      //-------------------------------------------------------------------------------------------------
      //4 sections is here including leakage: Cheat 4
      titleofplot1 = "Reconstructed energy using 4 sections but cheat section 4 including leakage for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      hist25 = (TH1F*)results[k]->Get( ("h_recoE_4sections_ch4_leak_" + IntToString(particleenergies[l])).c_str());
      buildplot(c118, hist25, leg26, titleofplot, "Reconstructed energy (GeV)", "Events/0.1 GeV", procNa[l], l);
      std::cout << "Working on plot: " << hist25->GetName() << std::endl;
      //-------------------------------------------------------------------------------------------------
      //-------------------------------------------------------------------------------------------------
      //RecoE: Cheat 1
      titleofplot1 = "Reconstructed energy but cheat section 1 for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      hist26 = (TH1F*)results[k]->Get( ("h_recoE_ch1_" + IntToString(particleenergies[l])).c_str());
      buildplot(c119, hist26, leg27, titleofplot, "Reconstructed energy (GeV)", "Events/0.1 GeV", procNa[l], l);
      std::cout << "Working on plot: " << hist26->GetName() << std::endl;
      //-------------------------------------------------------------------------------------------------
      //RecoE including leakage: Cheat 1
      titleofplot1 = "Reconstructed energy but cheat section 1 including leakage for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      hist27 = (TH1F*)results[k]->Get( ("h_recoE_ch1_leak_" + IntToString(particleenergies[l])).c_str());
      buildplot(c120, hist27, leg28, titleofplot, "Reconstructed energy (GeV)", "Events/0.1 GeV", procNa[l], l);
      std::cout << "Working on plot: " << hist27->GetName() << std::endl;
      //-------------------------------------------------------------------------------------------------
      //RecoE: Cheat 2
      titleofplot1 = "Reconstructed energy but cheat section 2 for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      hist28 = (TH1F*)results[k]->Get( ("h_recoE_ch2_" + IntToString(particleenergies[l])).c_str());
      buildplot(c121, hist28, leg29, titleofplot, "Reconstructed energy (GeV)", "Events/0.1 GeV", procNa[l], l);
      std::cout << "Working on plot: " << hist28->GetName() << std::endl;
      //-------------------------------------------------------------------------------------------------
      //RecoE including leakage: Cheat 2
      titleofplot1 = "Reconstructed energy but cheat section 2 including leakage for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      hist29 = (TH1F*)results[k]->Get( ("h_recoE_ch2_leak_" + IntToString(particleenergies[l])).c_str());
      buildplot(c122, hist29, leg30, titleofplot, "Reconstructed energy (GeV)", "Events/0.1 GeV", procNa[l], l);
      std::cout << "Working on plot: " << hist29->GetName() << std::endl;
      //-------------------------------------------------------------------------------------------------
      //RecoE: Cheat 3
      titleofplot1 = "Reconstructed energy but cheat section 3 for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      hist30 = (TH1F*)results[k]->Get( ("h_recoE_ch3_" + IntToString(particleenergies[l])).c_str());
      buildplot(c123, hist30, leg31, titleofplot, "Reconstructed energy (GeV)", "Events/0.1 GeV", procNa[l], l);
      std::cout << "Working on plot: " << hist30->GetName() << std::endl;
      //-------------------------------------------------------------------------------------------------
      //RecoE including leakage: Cheat 3
      titleofplot1 = "Reconstructed energy but cheat section 3 including leakage for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      hist31 = (TH1F*)results[k]->Get( ("h_recoE_ch3_leak_" + IntToString(particleenergies[l])).c_str());
      buildplot(c124, hist31, leg32, titleofplot, "Reconstructed energy (GeV)", "Events/0.1 GeV", procNa[l], l);
      std::cout << "Working on plot: " << hist31->GetName() << std::endl;
      //-------------------------------------------------------------------------------------------------
      //RecoE: Cheat 4
      titleofplot1 = "Reconstructed energy but cheat section 4 for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      hist32 = (TH1F*)results[k]->Get( ("h_recoE_ch4_" + IntToString(particleenergies[l])).c_str());
      buildplot(c125, hist32, leg33, titleofplot, "Reconstructed energy (GeV)", "Events/0.1 GeV", procNa[l], l);
      std::cout << "Working on plot: " << hist32->GetName() << std::endl;
      //-------------------------------------------------------------------------------------------------
      //RecoE including leakage: Cheat 4
      titleofplot1 = "Reconstructed energy but cheat section 4 including leakage for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      hist33 = (TH1F*)results[k]->Get( ("h_recoE_ch4_leak_" + IntToString(particleenergies[l])).c_str());
      buildplot(c126, hist33, leg34, titleofplot, "Reconstructed energy (GeV)", "Events/0.1 GeV", procNa[l], l);
      std::cout << "Working on plot: " << hist33->GetName() << std::endl;
      //-------------------------------------------------------------------------------------------------
      //-------------------------------------------------------------------------------------------------
      //-------------------------------------------------------------------------------------------------
      //-------------------------------------------------------------------------------------------------
      //RecoE overtruE: Cheat 1
      titleofplot1 = "Reconstructed energy minus true energy over true energy but cheat section 1 for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      hist34 = (TH1F*)results[k]->Get( ("h_recoE_ch1_overtrueE_" + IntToString(particleenergies[l])).c_str());
      buildplot(c162, hist34, leg35, titleofplot, "#frac{E_{reco}-E_{true}}{E_{true}}", "Events", procNa[l], l);
      std::cout << "Working on plot: " << hist34->GetName() << std::endl;
      //-------------------------------------------------------------------------------------------------
      //RecoE overtruE including leakage: Cheat 1
      titleofplot1 = "Reconstructed energy minus true energy over true energy but cheat section 1 including leakage for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      hist35 = (TH1F*)results[k]->Get( ("h_recoE_ch1_overtrueE_leak_" + IntToString(particleenergies[l])).c_str());
      buildplot(c163, hist35, leg36, titleofplot, "#frac{E_{reco}-E_{true}}{E_{true}}", "Events", procNa[l], l);
      std::cout << "Working on plot: " << hist35->GetName() << std::endl;
      //-------------------------------------------------------------------------------------------------
      //RecoE overtruE: Cheat 2
      titleofplot1 = "Reconstructed energy minus true energy over true energy but cheat section 2 for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      hist36 = (TH1F*)results[k]->Get( ("h_recoE_ch2_overtrueE_" + IntToString(particleenergies[l])).c_str());
      buildplot(c164, hist36, leg37, titleofplot, "#frac{E_{reco}-E_{true}}{E_{true}}", "Events", procNa[l], l);
      std::cout << "Working on plot: " << hist36->GetName() << std::endl;
      //-------------------------------------------------------------------------------------------------
      //RecoE overtruE including leakage: Cheat 2
      titleofplot1 = "Reconstructed energy minus true energy over true energy but cheat section 2 including leakage for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      hist37 = (TH1F*)results[k]->Get( ("h_recoE_ch2_overtrueE_leak_" + IntToString(particleenergies[l])).c_str());
      buildplot(c165, hist37, leg38, titleofplot, "#frac{E_{reco}-E_{true}}{E_{true}}", "Events", procNa[l], l);
      std::cout << "Working on plot: " << hist37->GetName() << std::endl;
      //-------------------------------------------------------------------------------------------------
      //RecoE overtruE: Cheat 3
      titleofplot1 = "Reconstructed energy minus true energy over true energy but cheat section 3 for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      hist38 = (TH1F*)results[k]->Get( ("h_recoE_ch3_overtrueE_" + IntToString(particleenergies[l])).c_str());
      buildplot(c166, hist38, leg39, titleofplot, "#frac{E_{reco}-E_{true}}{E_{true}}", "Events", procNa[l], l);
      std::cout << "Working on plot: " << hist38->GetName() << std::endl;
      //-------------------------------------------------------------------------------------------------
      //RecoE overtruE including leakage: Cheat 3
      titleofplot1 = "Reconstructed energy minus true energy over true energy but cheat section 3 including leakage for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      hist39 = (TH1F*)results[k]->Get( ("h_recoE_ch3_overtrueE_leak_" + IntToString(particleenergies[l])).c_str());
      buildplot(c167, hist39, leg40, titleofplot, "#frac{E_{reco}-E_{true}}{E_{true}}", "Events", procNa[l], l);
      std::cout << "Working on plot: " << hist39->GetName() << std::endl;
      //-------------------------------------------------------------------------------------------------
      //RecoE overtruE: Cheat 4
      titleofplot1 = "Reconstructed energy minus true energy over true energy but cheat section 4 for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      hist40 = (TH1F*)results[k]->Get( ("h_recoE_ch4_overtrueE_" + IntToString(particleenergies[l])).c_str());
      buildplot(c168, hist40, leg41, titleofplot, "#frac{E_{reco}-E_{true}}{E_{true}}", "Events", procNa[l], l);
      std::cout << "Working on plot: " << hist40->GetName() << std::endl;
      //-------------------------------------------------------------------------------------------------
      //RecoE overtruE including leakage: Cheat 4
      titleofplot1 = "Reconstructed energy minus true energy over true energy but cheat section 4 including leakage for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      hist41 = (TH1F*)results[k]->Get( ("h_recoE_ch4_overtrueE_leak_" + IntToString(particleenergies[l])).c_str());
      buildplot(c169, hist41, leg42, titleofplot, "#frac{E_{reco}-E_{true}}{E_{true}}", "Events", procNa[l], l);
      std::cout << "Working on plot: " << hist41->GetName() << std::endl;
      //-------------------------------------------------------------------------------------------------
      //-------------------------------------------------------------------------------------------------
      //-------------------------------------------------------------------------------------------------
      //-------------------------------------------------------------------------------------------------
      //Using the total sampling fraction per energy here
      titleofplot1 = "Reconstructed energy using total sampling fraction for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      hist12 = (TH1F*)results[k]->Get( ("h_recoE_totsf_" + IntToString(particleenergies[l])).c_str());
      buildplot(c45, hist12, leg13, titleofplot, "Reconstructed energy (GeV)", "Events/0.1 GeV", procNa[l], l);
      std::cout << "Working on plot: " << hist12->GetName() << std::endl;
      //-------------------------------------------------------------------------------------------------
      //Using the total sampling fraction per energy here adding leakage here
      titleofplot1 = "Reconstructed energy using total sampling fraction including leakage for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      hist13 = (TH1F*)results[k]->Get( ("h_recoE_totsf_leak_" + IntToString(particleenergies[l])).c_str());
      buildplot(c53, hist13, leg14, titleofplot, "Reconstructed energy (GeV)", "Events/0.1 GeV", procNa[l], l);
      std::cout << "Working on plot: " << hist13->GetName() << std::endl;

      //-------------------------------------------------------------------------------------------------
      //Using the total sampling fraction per energy here with dedx weight
      titleofplot1 = "Reconstructed energy using total sampling fraction using dEdx weight for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      hist80 = (TH1F*)results[k]->Get( ("h_recoE_totsf_dedxweighted_" + IntToString(particleenergies[l])).c_str());
      buildplot(c314, hist80, leg81, titleofplot, "Reconstructed energy (GeV)", "Events/0.1 GeV", procNa[l], l);
      std::cout << "Working on plot: " << hist80->GetName() << std::endl;
      //-------------------------------------------------------------------------------------------------
      //Using the total sampling fraction per energy here adding leakage here with dedx weight
      titleofplot1 = "Reconstructed energy using total sampling fraction including leakage using dEdx weight for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      hist81 = (TH1F*)results[k]->Get( ("h_recoE_totsf_dedxweighted_leak_" + IntToString(particleenergies[l])).c_str());
      buildplot(c315, hist81, leg82, titleofplot, "Reconstructed energy (GeV)", "Events/0.1 GeV", procNa[l], l);
      std::cout << "Working on plot: " << hist81->GetName() << std::endl;
      //-------------------------------------------------------------------------------------------------



      //-------------------------------------------------------------------------------------------------
      //Using the total sampling fraction per energy here
      titleofplot1 = "Reconstructed energy using total sampling fraction per section for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      hist56 = (TH1F*)results[k]->Get( ("h_recoE_totsf_sec_" + IntToString(particleenergies[l])).c_str());
      buildplot(c222, hist56, leg57, titleofplot, "Reconstructed energy (GeV)", "Events/0.1 GeV", procNa[l], l);
      std::cout << "Working on plot: " << hist56->GetName() << std::endl;
      //-------------------------------------------------------------------------------------------------
      //Using the total sampling fraction per energy here adding leakage here
      titleofplot1 = "Reconstructed energy using total sampling fraction per section including leakage for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      hist57 = (TH1F*)results[k]->Get( ("h_recoE_totsf_sec_leak_" + IntToString(particleenergies[l])).c_str());
      buildplot(c223, hist57, leg58, titleofplot, "Reconstructed energy (GeV)", "Events/0.1 GeV", procNa[l], l);
      std::cout << "Working on plot: " << hist57->GetName() << std::endl;


      //-------------------------------------------------------------------------------------------------
      //Using the total sampling fraction per energy here
      titleofplot1 = "Reconstructed energy using total sampling fraction per section weighted with dedx for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      hist72 = (TH1F*)results[k]->Get( ("h_recoE_totsf_sec_dedxweighted_" + IntToString(particleenergies[l])).c_str());
      buildplot(c284, hist72, leg73, titleofplot, "Reconstructed energy (GeV)", "Events/0.1 GeV", procNa[l], l);
      std::cout << "Working on plot: " << hist72->GetName() << std::endl;
      //-------------------------------------------------------------------------------------------------
      //Using the total sampling fraction per energy here adding leakage here
      titleofplot1 = "Reconstructed energy using total sampling fraction per section weighted with dedx including leakage for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      hist73 = (TH1F*)results[k]->Get( ("h_recoE_totsf_sec_dedxweighted_leak_" + IntToString(particleenergies[l])).c_str());
      buildplot(c285, hist73, leg74, titleofplot, "Reconstructed energy (GeV)", "Events/0.1 GeV", procNa[l], l);
      std::cout << "Working on plot: " << hist73->GetName() << std::endl;


      //-------------------------------------------------------------------------------------------------
      //Using the total sampling fraction per energy per section with shower max method here
      titleofplot1 = "Reconstructed energy using total SF per section with shower max method for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      hist60 = (TH1F*)results[k]->Get( ("h_recoE_totsfshmax_sec_" + IntToString(particleenergies[l])).c_str());
      buildplot(c238, hist60, leg61, titleofplot, "Reconstructed energy (GeV)", "Events/0.1 GeV", procNa[l], l);
      std::cout << "Working on plot: " << hist60->GetName() << std::endl;
      //-------------------------------------------------------------------------------------------------
      //Using the total sampling fraction per energy per section with shower max method here adding leakage here
      titleofplot1 = "Reconstructed energy using total SF per section with shower max method including leakage for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      hist61 = (TH1F*)results[k]->Get( ("h_recoE_totsfshmax_sec_leak_" + IntToString(particleenergies[l])).c_str());
      buildplot(c239, hist61, leg61, titleofplot, "Reconstructed energy (GeV)", "Events/0.1 GeV", procNa[l], l);
      std::cout << "Working on plot: " << hist61->GetName() << std::endl;
      //-------------------------------------------------------------------------------------------------

      //-------------------------------------------------------------------------------------------------
      //Using the total sampling fraction per energy per section with shower max method here
      titleofplot1 = "Reconstructed energy using total SF per dedx weighted section with shower max method for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      hist76 = (TH1F*)results[k]->Get( ("h_recoE_totsfshmax_sec_dedxweighted_" + IntToString(particleenergies[l])).c_str());
      buildplot(c299, hist76, leg77, titleofplot, "Reconstructed energy (GeV)", "Events/0.1 GeV", procNa[l], l);
      std::cout << "Working on plot: " << hist76->GetName() << std::endl;
      //-------------------------------------------------------------------------------------------------
      //Using the total sampling fraction per energy per section with shower max method here adding leakage here
      titleofplot1 = "Reconstructed energy using total SF per dedx weighted section with shower max method including leakage for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      hist77 = (TH1F*)results[k]->Get( ("h_recoE_totsfshmax_sec_dedxweighted_leak_" + IntToString(particleenergies[l])).c_str());
      buildplot(c300, hist77, leg78, titleofplot, "Reconstructed energy (GeV)", "Events/0.1 GeV", procNa[l], l);
      std::cout << "Working on plot: " << hist77->GetName() << std::endl;
      //-------------------------------------------------------------------------------------------------



     //-------------------------------------------------------------------------------------------------
      //Using the total sampling fraction per energy here
      titleofplot1 = "Reconstructed energy using total sampling fraction per section for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      hist64 = (TH1F*)results[k]->Get( ("h_recoE_totsf_dedxsec_" + IntToString(particleenergies[l])).c_str());
      buildplot(c253, hist64, leg65, titleofplot, "Reconstructed energy (GeV)", "Events/0.1 GeV", procNa[l], l);
      std::cout << "Working on plot: " << hist64->GetName() << std::endl;
      //-------------------------------------------------------------------------------------------------
      //Using the total sampling fraction per energy here adding leakage here
      titleofplot1 = "Reconstructed energy using total sampling fraction per section including leakage for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      hist65 = (TH1F*)results[k]->Get( ("h_recoE_totsf_dedxsec_leak_" + IntToString(particleenergies[l])).c_str());
      buildplot(c254, hist65, leg66, titleofplot, "Reconstructed energy (GeV)", "Events/0.1 GeV", procNa[l], l);
      std::cout << "Working on plot: " << hist65->GetName() << std::endl;
      //-------------------------------------------------------------------------------------------------
      //Using the total sampling fraction per energy per section with shower max method here
      titleofplot1 = "Reconstructed energy using total SF per section with shower max method for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      hist66 = (TH1F*)results[k]->Get( ("h_recoE_totsfshmax_dedxsec_" + IntToString(particleenergies[l])).c_str());
      buildplot(c255, hist66, leg67, titleofplot, "Reconstructed energy (GeV)", "Events/0.1 GeV", procNa[l], l);
      std::cout << "Working on plot: " << hist66->GetName() << std::endl;
      //-------------------------------------------------------------------------------------------------
      //Using the total sampling fraction per energy per section with shower max method here adding leakage here
      titleofplot1 = "Reconstructed energy using total SF per section with shower max method including leakage for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      hist67 = (TH1F*)results[k]->Get( ("h_recoE_totsfshmax_dedxsec_leak_" + IntToString(particleenergies[l])).c_str());
      buildplot(c256, hist67, leg68, titleofplot, "Reconstructed energy (GeV)", "Events/0.1 GeV", procNa[l], l);
      std::cout << "Working on plot: " << hist67->GetName() << std::endl;
      //-------------------------------------------------------------------------------------------------













      //-------------------------------------------------------------------------------------------------
      //Using the shower max method here
      titleofplot1 = "Reconstructed energy using shower max method for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      hist46 = (TH1F*)results[k]->Get( ("h_recoE_shmax_" + IntToString(particleenergies[l])).c_str());
      buildplot(c186, hist46, leg47, titleofplot, "Reconstructed energy (GeV)", "Events/0.1 GeV", procNa[l], l);
      std::cout << "Working on plot: " << hist46->GetName() << std::endl;
      //-------------------------------------------------------------------------------------------------
      //Using the shower max method here adding leakage here
      titleofplot1 = "Reconstructed energy using shower max including leakage for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      hist47 = (TH1F*)results[k]->Get( ("h_recoE_shmax_leak_" + IntToString(particleenergies[l])).c_str());
      buildplot(c187, hist47, leg48, titleofplot, "Reconstructed energy (GeV)", "Events/0.1 GeV", procNa[l], l);
      std::cout << "Working on plot: " << hist47->GetName() << std::endl;
      //-------------------------------------------------------------------------------------------------
      //Using the sampling fraction per normalized shower depth
      titleofplot1 = "Reconstructed energy using sampling fraction with normalized shower depth for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      hist14 = (TH1F*)results[k]->Get( ("h_recoE_sfnorm_" + IntToString(particleenergies[l])).c_str());
      buildplot(c56, hist14, leg15, titleofplot, "Reconstructed energy (GeV)", "Events/0.1 GeV", procNa[l], l);
      std::cout << "Working on plot: " << hist14->GetName() << std::endl;
      //-------------------------------------------------------------------------------------------------
      //Using the sampling fraction per normalized shower depth adding leakage
      titleofplot1 = "Reconstructed energy using sampling fraction with normalized shower depth adding leakage for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      hist15 = (TH1F*)results[k]->Get( ("h_recoE_sfnorm_leak_" + IntToString(particleenergies[l])).c_str());
      buildplot(c61, hist15, leg16, titleofplot, "Reconstructed energy (GeV)", "Events/0.1 GeV", procNa[l], l);
      std::cout << "Working on plot: " << hist15->GetName() << std::endl;
      
      //-------------------------------------------------------------------------------------------------
      //Using the sampling fraction per normalized shower depth
      titleofplot1 = "Reconstructed energy using total SF with X_{0,max}/<X_{0,max}> for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      hist50 = (TH1F*)results[k]->Get( ("h_recoE_totsfnorm_" + IntToString(particleenergies[l])).c_str());
      buildplot(c201, hist50, leg51, titleofplot, "Reconstructed energy (GeV)", "Events/0.1 GeV", procNa[l], l);
      std::cout << "Working on plot: " << hist50->GetName() << std::endl;
      //-------------------------------------------------------------------------------------------------
      //Using the sampling fraction per normalized shower depth adding leakage
      titleofplot1 = "Reconstructed energy using total SF with X_{0,max}/<X_{0,max}> adding leakage for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      hist51 = (TH1F*)results[k]->Get( ("h_recoE_totsfnorm_leak_" + IntToString(particleenergies[l])).c_str());
      buildplot(c202, hist51, leg52, titleofplot, "Reconstructed energy (GeV)", "Events/0.1 GeV", procNa[l], l);
      std::cout << "Working on plot: " << hist51->GetName() << std::endl;
      


      //-------------------------------------------------------------------------------------------------
      //Using the sampling fraction per normalized shower depth using dedx weight 
      titleofplot1 = "Reconstructed energy using total SF with X_{0,max}/<X_{0,max}> using dedx weight for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      hist84 = (TH1F*)results[k]->Get( ("h_recoE_totsfnorm_dedxweighted_" + IntToString(particleenergies[l])).c_str());
      buildplot(c329, hist84, leg85, titleofplot, "Reconstructed energy (GeV)", "Events/0.1 GeV", procNa[l], l);
      std::cout << "Working on plot: " << hist84->GetName() << std::endl;
      //-------------------------------------------------------------------------------------------------
      //Using the sampling fraction per normalized shower depth using dedx weight adding leakage
      titleofplot1 = "Reconstructed energy using total SF with X_{0,max}/<X_{0,max}> using dedx weight adding leakage for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      hist85 = (TH1F*)results[k]->Get( ("h_recoE_totsfnorm_dedxweighted_leak_" + IntToString(particleenergies[l])).c_str());
      buildplot(c330, hist85, leg86, titleofplot, "Reconstructed energy (GeV)", "Events/0.1 GeV", procNa[l], l);
      std::cout << "Working on plot: " << hist85->GetName() << std::endl;

      if(upstream){
	//-------------------------------------------------------------------------------------------------
	titleofplot1 = "Reconstructed energy with upstream material correction for "; 
	titleofplot2 =  particle +  " and all beam energies"; 
	titleofplot = titleofplot1 + titleofplot2;
	hist2 = (TH1F*)results[k]->Get( ("h_recoE_corr_" + IntToString(particleenergies[l])).c_str());
	buildplot(c2, hist2, leg2, titleofplot, "Reconstructed energy corrected (GeV)", "Events/GeV", procNa[l], l);
	std::cout << "Working on plot: " << hist2->GetName() << std::endl;
      }

      //-------------------------------------------------------------------------------------------------
      upstream ?  titleofplot1 = "Reconstructed energy without upstream material correction over beam energy for ": titleofplot1 = "Reconstructed energy minus true energy over true energy for "; 
      titleofplot2 =  particle; 
      titleofplot = titleofplot1 + titleofplot2;
      hist3 = (TH1F*)results[k]->Get( ("h_recoEovertrueE_" + IntToString(particleenergies[l])).c_str());
      buildplot(c5, hist3, leg4, titleofplot, "#frac{E_{reco}-E_{true}}{E_{true}}", "Events", procNa[l], l);
      std::cout << "Working on plot: " << hist3->GetName() << std::endl;

      //-------------------------------------------------------------------------------------------------
      upstream ?  titleofplot1 = "Reconstructed energy without upstream material correction over beam energy adding leakage for ": titleofplot1 = "Reconstructed energy minus true energy over true energy adding leakage for "; 
      titleofplot2 =  particle; 
      titleofplot = titleofplot1 + titleofplot2;
      hist42 = (TH1F*)results[k]->Get( ("h_recoEovertrueE_leak_" + IntToString(particleenergies[l])).c_str());
      buildplot(c178, hist42, leg43, titleofplot, "#frac{E_{reco,leak}-E_{true}}{E_{true}}", "Events", procNa[l], l);
      std::cout << "Working on plot: " << hist42->GetName() << std::endl;

      //-------------------------------------------------------------------------------------------------
      upstream ?  titleofplot1 = "Reconstructed energy without upstream material correction over beam energy adding leakage and offset for ": titleofplot1 = "Reconstructed energy minus true energy over true energy adding leakage and offset for "; 
      titleofplot2 =  particle; 
      titleofplot = titleofplot1 + titleofplot2;
      hist43 = (TH1F*)results[k]->Get( ("h_recoEovertrueE_leak_offset_" + IntToString(particleenergies[l])).c_str());
      buildplot(c179, hist43, leg44, titleofplot, "#frac{E_{reco,leak,offset}-E_{true}}{E_{true}}", "Events", procNa[l], l);
      std::cout << "Working on plot: " << hist43->GetName() << std::endl;

      //-------------------------------------------------------------------------------------------------
      titleofplot1 = "Reconstructed energy minus true energy over true energy for total SF for "; 
      titleofplot2 =  particle; 
      titleofplot = titleofplot1 + titleofplot2;
      hist44 = (TH1F*)results[k]->Get( ("h_recoEovertrueE_totsf_" + IntToString(particleenergies[l])).c_str());
      buildplot(c184, hist44, leg45, titleofplot, "#frac{E_{reco,totsf}-E_{true}}{E_{true}}", "Events", procNa[l], l);
      std::cout << "Working on plot: " << hist44->GetName() << std::endl;

      //-------------------------------------------------------------------------------------------------
      titleofplot1 = "Reconstructed energy minus true energy over true energy for total SF adding leakage for "; 
      titleofplot2 =  particle; 
      titleofplot = titleofplot1 + titleofplot2;
      hist45 = (TH1F*)results[k]->Get( ("h_recoEovertrueE_totsf_leak_" + IntToString(particleenergies[l])).c_str());
      buildplot(c185, hist45, leg46, titleofplot, "#frac{E_{reco,totsf,leak}-E_{true}}{E_{true}}", "Events", procNa[l], l);
      std::cout << "Working on plot: " << hist45->GetName() << std::endl;

      //-------------------------------------------------------------------------------------------------
      titleofplot1 = "Reconstructed energy minus true energy over true energy for total SF using dEdx weight for "; 
      titleofplot2 =  particle; 
      titleofplot = titleofplot1 + titleofplot2;
      hist82 = (TH1F*)results[k]->Get( ("h_recoEovertrueE_totsf_dedxweighted_" + IntToString(particleenergies[l])).c_str());
      buildplot(c316, hist82, leg83, titleofplot, "#frac{E_{reco,totsfdedx}-E_{true}}{E_{true}}", "Events", procNa[l], l);
      std::cout << "Working on plot: " << hist82->GetName() << std::endl;

      //-------------------------------------------------------------------------------------------------
      titleofplot1 = "Reconstructed energy minus true energy over true energy for total SF using dEdx weight adding leakage for "; 
      titleofplot2 =  particle; 
      titleofplot = titleofplot1 + titleofplot2;
      hist83 = (TH1F*)results[k]->Get( ("h_recoEovertrueE_totsf_dedxweighted_leak_" + IntToString(particleenergies[l])).c_str());
      buildplot(c317, hist83, leg84, titleofplot, "#frac{E_{reco,totsfdedx,leak}-E_{true}}{E_{true}}", "Events", procNa[l], l);
      std::cout << "Working on plot: " << hist83->GetName() << std::endl;

      //-------------------------------------------------------------------------------------------------

      //-------------------------------------------------------------------------------------------------
      titleofplot1 = "Reconstructed energy minus true energy over true energy for total SF per section for "; 
      titleofplot2 =  particle; 
      titleofplot = titleofplot1 + titleofplot2;
      hist58 = (TH1F*)results[k]->Get( ("h_recoEovertrueE_totsf_sec_" + IntToString(particleenergies[l])).c_str());
      buildplot(c224, hist58, leg59, titleofplot, "#frac{E_{reco,totsfsec}-E_{true}}{E_{true}}", "Events", procNa[l], l);
      std::cout << "Working on plot: " << hist58->GetName() << std::endl;

      //-------------------------------------------------------------------------------------------------
      titleofplot1 = "Reconstructed energy minus true energy over true energy for total SF per section adding leakage for "; 
      titleofplot2 =  particle; 
      titleofplot = titleofplot1 + titleofplot2;
      hist59 = (TH1F*)results[k]->Get( ("h_recoEovertrueE_totsf_sec_leak_" + IntToString(particleenergies[l])).c_str());
      buildplot(c225, hist59, leg60, titleofplot, "#frac{E_{reco,totsfsec,leak}-E_{true}}{E_{true}}", "Events", procNa[l], l);
      std::cout << "Working on plot: " << hist59->GetName() << std::endl;

      //-------------------------------------------------------------------------------------------------
      titleofplot1 = "Reconstructed energy minus true energy over true energy for total SF per section weighted with dedx for "; 
      titleofplot2 =  particle; 
      titleofplot = titleofplot1 + titleofplot2;
      hist74 = (TH1F*)results[k]->Get( ("h_recoEovertrueE_totsf_sec_dedxweighted_" + IntToString(particleenergies[l])).c_str());
      buildplot(c286, hist74, leg75, titleofplot, "#frac{E_{reco,totsfsecdedx}-E_{true}}{E_{true}}", "Events", procNa[l], l);
      std::cout << "Working on plot: " << hist74->GetName() << std::endl;

      //-------------------------------------------------------------------------------------------------
      titleofplot1 = "Reconstructed energy minus true energy over true energy for total SF per section weighted with dedx adding leakage for "; 
      titleofplot2 =  particle; 
      titleofplot = titleofplot1 + titleofplot2;
      hist75 = (TH1F*)results[k]->Get( ("h_recoEovertrueE_totsf_sec_dedxweighted_leak_" + IntToString(particleenergies[l])).c_str());
      buildplot(c287, hist75, leg76, titleofplot, "#frac{E_{reco,totsfsecdedx,leak}-E_{true}}{E_{true}}", "Events", procNa[l], l);
      std::cout << "Working on plot: " << hist75->GetName() << std::endl;

      //-------------------------------------------------------------------------------------------------
      titleofplot1 = "Reconstructed energy minus true energy over true energy for total SF per section with shower max method for "; 
      titleofplot2 =  particle; 
      titleofplot = titleofplot1 + titleofplot2;
      hist62 = (TH1F*)results[k]->Get( ("h_recoEovertrueE_totsfshmax_sec_" + IntToString(particleenergies[l])).c_str());
      buildplot(c240, hist62, leg63, titleofplot, "#frac{E_{reco,totsfsec}-E_{true}}{E_{true}}", "Events", procNa[l], l);
      std::cout << "Working on plot: " << hist62->GetName() << std::endl;

      //-------------------------------------------------------------------------------------------------
      titleofplot1 = "Reconstructed energy minus true energy over true energy for total SF per section with shower max method adding leakage for "; 
      titleofplot2 =  particle; 
      titleofplot = titleofplot1 + titleofplot2;
      hist63 = (TH1F*)results[k]->Get( ("h_recoEovertrueE_totsfshmax_sec_leak_" + IntToString(particleenergies[l])).c_str());
      buildplot(c241, hist63, leg64, titleofplot, "#frac{E_{reco,totsfsec,leak}-E_{true}}{E_{true}}", "Events", procNa[l], l);
      std::cout << "Working on plot: " << hist63->GetName() << std::endl;

      //-------------------------------------------------------------------------------------------------
      titleofplot1 = "Reconstructed energy minus true energy over true energy for total SF per dedx weighted section with shower max method for "; 
      titleofplot2 =  particle; 
      titleofplot = titleofplot1 + titleofplot2;
      hist78 = (TH1F*)results[k]->Get( ("h_recoEovertrueE_totsfshmax_sec_dedxweighted_" + IntToString(particleenergies[l])).c_str());
      buildplot(c301, hist78, leg79, titleofplot, "#frac{E_{reco,totsfsec}-E_{true}}{E_{true}}", "Events", procNa[l], l);
      std::cout << "Working on plot: " << hist78->GetName() << std::endl;

      //-------------------------------------------------------------------------------------------------
      titleofplot1 = "Reconstructed energy minus true energy over true energy for total SF per dedx weighted section with shower max method adding leakage for "; 
      titleofplot2 =  particle; 
      titleofplot = titleofplot1 + titleofplot2;
      hist79 = (TH1F*)results[k]->Get( ("h_recoEovertrueE_totsfshmax_sec_dedxweighted_leak_" + IntToString(particleenergies[l])).c_str());
      buildplot(c302, hist79, leg80, titleofplot, "#frac{E_{reco,totsfsec,leak}-E_{true}}{E_{true}}", "Events", procNa[l], l);
      std::cout << "Working on plot: " << hist79->GetName() << std::endl;





      //-------------------------------------------------------------------------------------------------
      titleofplot1 = "Reconstructed energy minus true energy over true energy for total SF per section for "; 
      titleofplot2 =  particle; 
      titleofplot = titleofplot1 + titleofplot2;
      hist68 = (TH1F*)results[k]->Get( ("h_recoEovertrueE_totsf_dedxsec_" + IntToString(particleenergies[l])).c_str());
      buildplot(c257, hist68, leg69, titleofplot, "#frac{E_{reco,totsfdedxsec}-E_{true}}{E_{true}}", "Events", procNa[l], l);
      std::cout << "Working on plot: " << hist68->GetName() << std::endl;

      //-------------------------------------------------------------------------------------------------
      titleofplot1 = "Reconstructed energy minus true energy over true energy for total SF per section adding leakage for "; 
      titleofplot2 =  particle; 
      titleofplot = titleofplot1 + titleofplot2;
      hist69 = (TH1F*)results[k]->Get( ("h_recoEovertrueE_totsf_dedxsec_leak_" + IntToString(particleenergies[l])).c_str());
      buildplot(c258, hist69, leg70, titleofplot, "#frac{E_{reco,totsfdedxsec,leak}-E_{true}}{E_{true}}", "Events", procNa[l], l);
      std::cout << "Working on plot: " << hist69->GetName() << std::endl;

      //-------------------------------------------------------------------------------------------------
      titleofplot1 = "Reconstructed energy minus true energy over true energy for total SF per section with shower max method for "; 
      titleofplot2 =  particle; 
      titleofplot = titleofplot1 + titleofplot2;
      hist70 = (TH1F*)results[k]->Get( ("h_recoEovertrueE_totsfshmax_dedxsec_" + IntToString(particleenergies[l])).c_str());
      buildplot(c259, hist70, leg71, titleofplot, "#frac{E_{reco,totsfdedxsec}-E_{true}}{E_{true}}", "Events", procNa[l], l);
      std::cout << "Working on plot: " << hist70->GetName() << std::endl;

      //-------------------------------------------------------------------------------------------------
      titleofplot1 = "Reconstructed energy minus true energy over true energy for total SF per section with shower max method adding leakage for "; 
      titleofplot2 =  particle; 
      titleofplot = titleofplot1 + titleofplot2;
      hist71 = (TH1F*)results[k]->Get( ("h_recoEovertrueE_totsfshmax_dedxsec_leak_" + IntToString(particleenergies[l])).c_str());
      buildplot(c260, hist71, leg72, titleofplot, "#frac{E_{reco,totsfdedxsec,leak}-E_{true}}{E_{true}}", "Events", procNa[l], l);
      std::cout << "Working on plot: " << hist71->GetName() << std::endl;




      //-------------------------------------------------------------------------------------------------
      titleofplot1 = "Reconstructed energy minus true energy over true energy for total SF with X_{0,max}/X_{0,max} for "; 
      titleofplot2 =  particle; 
      titleofplot = titleofplot1 + titleofplot2;
      hist52 = (TH1F*)results[k]->Get( ("h_recoEovertrueE_totsfnorm_" + IntToString(particleenergies[l])).c_str());
      buildplot(c203, hist52, leg53, titleofplot, "#frac{E_{reco,totsfshmax}-E_{true}}{E_{true}}", "Events", procNa[l], l);
      std::cout << "Working on plot: " << hist52->GetName() << std::endl;

      //-------------------------------------------------------------------------------------------------
      titleofplot1 = "Reconstructed energy minus true energy over true energy for total SF with X0,max/X_{0,max} adding leakage for "; 
      titleofplot2 =  particle; 
      titleofplot = titleofplot1 + titleofplot2;
      hist53 = (TH1F*)results[k]->Get( ("h_recoEovertrueE_totsfnorm_leak_" + IntToString(particleenergies[l])).c_str());
      buildplot(c204, hist53, leg54, titleofplot, "#frac{E_{reco,totsfshmax,leak}-E_{true}}{E_{true}}", "Events", procNa[l], l);
      std::cout << "Working on plot: " << hist53->GetName() << std::endl;


      //-------------------------------------------------------------------------------------------------
      titleofplot1 = "Reconstructed energy minus true energy over true energy for total SF with X_{0,max}/X_{0,max} using dedx weight for "; 
      titleofplot2 =  particle; 
      titleofplot = titleofplot1 + titleofplot2;
      hist86 = (TH1F*)results[k]->Get( ("h_recoEovertrueE_totsfnorm_dedxweighted_" + IntToString(particleenergies[l])).c_str());
      buildplot(c331, hist86, leg87, titleofplot, "#frac{E_{reco,totsfshmaxdedx}-E_{true}}{E_{true}}", "Events", procNa[l], l);
      std::cout << "Working on plot: " << hist86->GetName() << std::endl;

      //-------------------------------------------------------------------------------------------------
      titleofplot1 = "Reconstructed energy minus true energy over true energy for total SF with X0,max/X_{0,max} using dedx weight adding leakage for "; 
      titleofplot2 =  particle; 
      titleofplot = titleofplot1 + titleofplot2;
      hist87 = (TH1F*)results[k]->Get( ("h_recoEovertrueE_totsfnorm_dedxweighted_leak_" + IntToString(particleenergies[l])).c_str());
      buildplot(c332, hist87, leg88, titleofplot, "#frac{E_{reco,totsfshmaxdedx,leak}-E_{true}}{E_{true}}", "Events", procNa[l], l);
      std::cout << "Working on plot: " << hist87->GetName() << std::endl;




      //-------------------------------------------------------------------------------------------------
      titleofplot1 = "Reconstructed energy minus true energy over true energy for total SF with X_{0,max}/X_{0,max} for "; 
      titleofplot2 =  particle; 
      titleofplot = titleofplot1 + titleofplot2;
      hist54 = (TH1F*)results[k]->Get( ("h_recoEovertrueE_sfnorm_" + IntToString(particleenergies[l])).c_str());
      buildplot(c216, hist54, leg55, titleofplot, "#frac{E_{reco,sfshmax}-E_{true}}{E_{true}}", "Events", procNa[l], l);
      std::cout << "Working on plot: " << hist54->GetName() << std::endl;

      //-------------------------------------------------------------------------------------------------
      titleofplot1 = "Reconstructed energy minus true energy over true energy for total SF with X_{0,max}/<X_{0,max}> adding leakage for "; 
      titleofplot2 =  particle; 
      titleofplot = titleofplot1 + titleofplot2;
      hist55 = (TH1F*)results[k]->Get( ("h_recoEovertrueE_sfnorm_leak_" + IntToString(particleenergies[l])).c_str());
      buildplot(c217, hist55, leg56, titleofplot, "#frac{E_{reco,sfshmax,leak}-E_{true}}{E_{true}}", "Events", procNa[l], l);
      std::cout << "Working on plot: " << hist55->GetName() << std::endl;

      //-------------------------------------------------------------------------------------------------
      titleofplot1 = "Reconstructed energy minus true energy over true energy using SF_{i} and shower max method for "; 
      titleofplot2 =  particle; 
      titleofplot = titleofplot1 + titleofplot2;
      hist48 = (TH1F*)results[k]->Get( ("h_recoEovertrueE_shmax_" + IntToString(particleenergies[l])).c_str());
      buildplot(c188, hist48, leg49, titleofplot, "#frac{E_{reco,shmax}-E_{true}}{E_{true}}", "Events", procNa[l], l);
      std::cout << "Working on plot: " << hist48->GetName() << std::endl;

      //-------------------------------------------------------------------------------------------------
      titleofplot1 = "Reconstructed energy minus true energy over true energy using SF_{i} and shower max method adding leakage for "; 
      titleofplot2 =  particle; 
      titleofplot = titleofplot1 + titleofplot2;
      hist49 = (TH1F*)results[k]->Get( ("h_recoEovertrueE_shmax_leak_" + IntToString(particleenergies[l])).c_str());
      buildplot(c189, hist49, leg50, titleofplot, "#frac{E_{reco,shmax,leak}-E_{true}}{E_{true}}", "Events", procNa[l], l);
      std::cout << "Working on plot: " << hist49->GetName() << std::endl;
      //-------------------------------------------------------------------------------------------------

      if(upstream){
	titleofplot1 = "Reconstructed energy with upstream material correction minus true energy over true energy for "; 
	titleofplot2 =  particle; 
	titleofplot = titleofplot1 + titleofplot2;
	hist4 = (TH1F*)results[k]->Get( ("h_recoE_corrovertrueE_" + IntToString(particleenergies[l])).c_str());
	buildplot(c6, hist4, leg5, titleofplot, "#frac{E_{reco}-E_{true}}{E_{true}}", "Events", procNa[l], l);
	std::cout << "Working on plot: " << hist4->GetName() << std::endl;
      }

      //-------------------------------------------------------------------------------------------------
      titleofplot1 = "Reconstructed energy using dEdx minus true energy over true energy for "; 
      titleofplot2 =  particle; 
      titleofplot = titleofplot1 + titleofplot2;
      hist7 = (TH1F*)results[k]->Get( ("h_recoE_valEovertrueE_" + IntToString(particleenergies[l])).c_str());
      buildplot(c26, hist7, leg8, titleofplot, "#frac{E_{reco}-E_{true}}{E_{true}}", "Events", procNa[l], l);
      std::cout << "Working on plot: " << hist7->GetName() << std::endl;
 
      //-------------------------------------------------------------------------------------------------
      titleofplot1 = "Raw reconstructed energy minus true energy over true energy for "; 
      titleofplot2 =  particle; 
      titleofplot = titleofplot1 + titleofplot2;
      hist8 = (TH1F*)results[k]->Get( ("h_recoE_rawEovertrueE_" + IntToString(particleenergies[l])).c_str());
      buildplot(c28, hist8, leg9, titleofplot, "#frac{E_{sim}-E_{true}}{E_{true}}", "Events", procNa[l], l);
      std::cout << "Working on plot: " << hist8->GetName() << std::endl;

      //-------------------------------------------------------------------------------------------------
      bool dogaussfit = true;
      fitres = extractMeanEnergy(mycE_1, mycE_2, color, hist, rebin[l]);
      //For the linearity plot
      titleofplot1 = "Reconstructed energy for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c3,recoEvsbeamenergy,titleofplot,"Beam Energy (GeV)","E_{reco} (GeV)", l, particleenergies[l], fitres, dogaussfit);

      //Resolution
      titleofplot1 = "Relative energy resolution for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildresolutiongraph(c14, resoFit, titleofplot, l, particleenergies[l], fitres, dogaussfit);

      //Including leakage
      fitres = extractMeanEnergy(mycE_leak_1, mycE_leak_2, color, hist10, rebin[l]);
      //For the linearity plot
      titleofplot1 = "Reconstructed energy including leakage for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c35,recoE_leakvsbeamenergy,titleofplot,"Beam Energy (GeV)","E_{reco} (GeV)", l, particleenergies[l], fitres, dogaussfit);

      //Resolution
      titleofplot1 = "Relative energy resolution including leakage for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildresolutiongraph(c36, resoFit_leak, titleofplot, l, particleenergies[l], fitres, dogaussfit);
      
      //Including leakage plus the offset
      fitres = extractMeanEnergy(mycE_leak_offset_1, mycE_leak_offset_2, color, hist11, rebin[l]);
      //For the linearity plot
      titleofplot1 = "Reconstructed energy including leakage and offset for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c40,recoE_leak_offsetvsbeamenergy,titleofplot,"Beam Energy (GeV)","E_{reco} (GeV)", l, particleenergies[l], fitres, dogaussfit);

      //Resolution
      titleofplot1 = "Relative energy resolution including leakage and offset for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildresolutiongraph(c41, resoFit_leak_offset, titleofplot, l, particleenergies[l], fitres, dogaussfit);
      
      //For the 4 sections
      fitres = extractMeanEnergy(mycE_4sections_1, mycE_4sections_2, color, hist16, rebin[l]);
      //For the linearity plot
      titleofplot1 = "Reconstructed energy using 4 sections for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c68,recoE_4sectionsvsbeamenergy,titleofplot,"Beam Energy (GeV)","E_{reco,4 sections} (GeV)", l, particleenergies[l], fitres, dogaussfit);

      //Resolution
      titleofplot1 = "Relative energy resolution using 4 sections for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildresolutiongraph(c69, resoFit_4sections, titleofplot, l, particleenergies[l], fitres, dogaussfit);

      //For the 4 sections including leakage
      fitres = extractMeanEnergy(mycE_4sections_leak_1, mycE_4sections_leak_2, color, hist17, rebin[l]);
      //For the linearity plot
      titleofplot1 = "Reconstructed energy using 4 sections including leakage for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c70,recoE_4sections_leakvsbeamenergy,titleofplot,"Beam Energy (GeV)","E_{reco,4 sections, leak} (GeV)", l, particleenergies[l], fitres, dogaussfit);

      //Resolution
      titleofplot1 = "Relative energy resolution using 4 sections including leakage for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildresolutiongraph(c71, resoFit_4sections_leak, titleofplot, l, particleenergies[l], fitres, dogaussfit);
      
      //For the 4 sections: Cheat 1
      fitres = extractMeanEnergy(mycE_4sections_ch1_1, mycE_4sections_ch1_2, color, hist18, rebin[l]);
      //For the linearity plot
      titleofplot1 = "Reconstructed energy using 4 sections but cheat section 1 for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c78,recoE_4sections_ch1vsbeamenergy,titleofplot,"Beam Energy (GeV)","E_{reco,4 sections,ch1} (GeV)", l, particleenergies[l], fitres, dogaussfit);

      //Resolution
      titleofplot1 = "Relative energy resolution using 4 sections but cheat section 1 for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildresolutiongraph(c79, resoFit_4sections_ch1, titleofplot, l, particleenergies[l], fitres, dogaussfit);

      //For the 4 sections including leakage: Cheat 1
      fitres = extractMeanEnergy(mycE_4sections_ch1_leak_1, mycE_4sections_ch1_leak_2, color, hist19, rebin[l]);
      //For the linearity plot
      titleofplot1 = "Reconstructed energy using 4 sections but cheat section 1 including leakage for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c80,recoE_4sections_ch1_leakvsbeamenergy,titleofplot,"Beam Energy (GeV)","E_{reco,4 sections,ch1,leak} (GeV)", l, particleenergies[l], fitres, dogaussfit);

      //Resolution
      titleofplot1 = "Relative energy resolution using 4 sections but cheat section 1 including leakage for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildresolutiongraph(c81, resoFit_4sections_ch1_leak, titleofplot, l, particleenergies[l], fitres, dogaussfit);

      //For the 4 sections: Cheat 2
      fitres = extractMeanEnergy(mycE_4sections_ch2_1, mycE_4sections_ch2_2, color, hist20, rebin[l]);
      //For the linearity plot
      titleofplot1 = "Reconstructed energy using 4 sections but cheat section 2 for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c82,recoE_4sections_ch2vsbeamenergy,titleofplot,"Beam Energy (GeV)","E_{reco,4 sections,ch2} (GeV)", l, particleenergies[l], fitres, dogaussfit);

      //Resolution
      titleofplot1 = "Relative energy resolution using 4 sections but cheat section 2 for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildresolutiongraph(c83, resoFit_4sections_ch2, titleofplot, l, particleenergies[l], fitres, dogaussfit);

      //For the 4 sections including leakage: Cheat 2
      fitres = extractMeanEnergy(mycE_4sections_ch2_leak_1, mycE_4sections_ch2_leak_2, color, hist21, rebin[l]);
      //For the linearity plot
      titleofplot1 = "Reconstructed energy using 4 sections but cheat section 2 including leakage for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c84,recoE_4sections_ch2_leakvsbeamenergy,titleofplot,"Beam Energy (GeV)","E_{reco,4 sections,ch2,leak} (GeV)", l, particleenergies[l], fitres, dogaussfit);

      //Resolution
      titleofplot1 = "Relative energy resolution using 4 sections but cheat section 2 including leakage for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildresolutiongraph(c85, resoFit_4sections_ch2_leak, titleofplot, l, particleenergies[l], fitres, dogaussfit);

      //For the 4 sections: Cheat 3
      fitres = extractMeanEnergy(mycE_4sections_ch3_1, mycE_4sections_ch3_2, color, hist22, rebin[l]);
      //For the linearity plot
      titleofplot1 = "Reconstructed energy using 4 sections but cheat section 3 for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c86,recoE_4sections_ch3vsbeamenergy,titleofplot,"Beam Energy (GeV)","E_{reco,4 sections,ch3} (GeV)", l, particleenergies[l], fitres, dogaussfit);

      //Resolution
      titleofplot1 = "Relative energy resolution using 4 sections but cheat section 3 for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildresolutiongraph(c87, resoFit_4sections_ch3, titleofplot, l, particleenergies[l], fitres, dogaussfit);

      //For the 4 sections including leakage: Cheat 3
      fitres = extractMeanEnergy(mycE_4sections_ch3_leak_1, mycE_4sections_ch3_leak_2, color, hist23, rebin[l]);
      //For the linearity plot
      titleofplot1 = "Reconstructed energy using 4 sections but cheat section 3 including leakage for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c88,recoE_4sections_ch3_leakvsbeamenergy,titleofplot,"Beam Energy (GeV)","E_{reco,4 sections,ch3,leak} (GeV)", l, particleenergies[l], fitres, dogaussfit);

      //Resolution
      titleofplot1 = "Relative energy resolution using 4 sections but cheat section 3 including leakage for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildresolutiongraph(c89, resoFit_4sections_ch3_leak, titleofplot, l, particleenergies[l], fitres, dogaussfit);

      //For the 4 sections: Cheat 4
      fitres = extractMeanEnergy(mycE_4sections_ch4_1, mycE_4sections_ch4_2, color, hist24, rebin[l]);
      //For the linearity plot
      titleofplot1 = "Reconstructed energy using 4 sections but cheat section 4 for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c90,recoE_4sections_ch4vsbeamenergy,titleofplot,"Beam Energy (GeV)","E_{reco,4 sections,ch4} (GeV)", l, particleenergies[l], fitres, dogaussfit);

      //Resolution
      titleofplot1 = "Relative energy resolution using 4 sections but cheat section 4 for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildresolutiongraph(c91, resoFit_4sections_ch4, titleofplot, l, particleenergies[l], fitres, dogaussfit);

      //For the 4 sections including leakage: Cheat 4
      fitres = extractMeanEnergy(mycE_4sections_ch4_leak_1, mycE_4sections_ch4_leak_2, color, hist25, rebin[l]);
      //For the linearity plot
      titleofplot1 = "Reconstructed energy using 4 sections but cheat section 4 including leakage for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c92,recoE_4sections_ch4_leakvsbeamenergy,titleofplot,"Beam Energy (GeV)","E_{reco,4 sections,ch4,leak} (GeV)", l, particleenergies[l], fitres, dogaussfit);

      //Resolution
      titleofplot1 = "Relative energy resolution using 4 sections but cheat section 4 including leakage for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildresolutiongraph(c93, resoFit_4sections_ch4_leak, titleofplot, l, particleenergies[l], fitres, dogaussfit);

      //-------------------------------------------------------------------------------------------------
      //-------------------------------------------------------------------------------------------------
      //RecoE: Cheat 1
      fitres = extractMeanEnergy(mycE_ch1_1, mycE_ch1_2, color, hist26, rebin[l]);
      //For the linearity plot
      titleofplot1 = "Reconstructed energy but cheat section 1 for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c127,recoE_ch1vsbeamenergy,titleofplot,"Beam Energy (GeV)","E_{reco,ch1} (GeV)", l, particleenergies[l], fitres, dogaussfit);

      //Resolution
      titleofplot1 = "Relative energy resolution but cheat section 1 for ";  
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildresolutiongraph(c128, resoFit_ch1, titleofplot, l, particleenergies[l], fitres, dogaussfit);
      //Hack here
      titleofplot1 = "Relative energy resolution but cheat section 1 for ";  
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildresolutiongraph(c160, resoFit_ch1, titleofplot, l, particleenergies[l], fitres, dogaussfit);

      //RecoE including leakage: Cheat 1
      fitres = extractMeanEnergy(mycE_ch1_leak_1, mycE_ch1_leak_2, color, hist27, rebin[l]);
      //For the linearity plot
      titleofplot1 = "Reconstructed energy but cheat section 1 including leakage for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c129,recoE_ch1_leakvsbeamenergy,titleofplot,"Beam Energy (GeV)","E_{reco,ch1,leak} (GeV)", l, particleenergies[l], fitres, dogaussfit);

      //Resolution
      titleofplot1 = "Relative energy resolution but cheat section 1 including leakage for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildresolutiongraph(c130, resoFit_ch1_leak, titleofplot, l, particleenergies[l], fitres, dogaussfit);

      //RecoE: Cheat 2
      fitres = extractMeanEnergy(mycE_ch2_1, mycE_ch2_2, color, hist28, rebin[l]);
      //For the linearity plot
      titleofplot1 = "Reconstructed energy but cheat section 2 for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c131,recoE_ch2vsbeamenergy,titleofplot,"Beam Energy (GeV)","E_{reco,ch2} (GeV)", l, particleenergies[l], fitres, dogaussfit);

      //Resolution
      titleofplot1 = "Relative energy resolution but cheat section 2 for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildresolutiongraph(c132, resoFit_ch2, titleofplot, l, particleenergies[l], fitres, dogaussfit);

      //RecoE including leakage: Cheat 2
      fitres = extractMeanEnergy(mycE_ch2_leak_1, mycE_ch2_leak_2, color, hist29, rebin[l]);
      //For the linearity plot
      titleofplot1 = "Reconstructed energy but cheat section 2 including leakage for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c133,recoE_ch2_leakvsbeamenergy,titleofplot,"Beam Energy (GeV)","E_{reco,ch2,leak} (GeV)", l, particleenergies[l], fitres, dogaussfit);

      //Resolution
      titleofplot1 = "Relative energy resolution but cheat section 2 including leakage for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildresolutiongraph(c134, resoFit_ch2_leak, titleofplot, l, particleenergies[l], fitres, dogaussfit);

      //RecoE: Cheat 3
      fitres = extractMeanEnergy(mycE_ch3_1, mycE_ch3_2, color, hist30, rebin[l]);
      //For the linearity plot
      titleofplot1 = "Reconstructed energy but cheat section 3 for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c135,recoE_ch3vsbeamenergy,titleofplot,"Beam Energy (GeV)","E_{reco,ch3} (GeV)", l, particleenergies[l], fitres, dogaussfit);

      //Resolution
      titleofplot1 = "Relative energy resolution but cheat section 3 for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildresolutiongraph(c136, resoFit_ch3, titleofplot, l, particleenergies[l], fitres, dogaussfit);

      //RecoE including leakage: Cheat 3
      fitres = extractMeanEnergy(mycE_ch3_leak_1, mycE_ch3_leak_2, color, hist31, rebin[l]);
      //For the linearity plot
      titleofplot1 = "Reconstructed energy but cheat section 3 including leakage for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c137,recoE_ch3_leakvsbeamenergy,titleofplot,"Beam Energy (GeV)","E_{reco,ch3,leak} (GeV)", l, particleenergies[l], fitres, dogaussfit);

      //Resolution
      titleofplot1 = "Relative energy resolution but cheat section 3 including leakage for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildresolutiongraph(c138, resoFit_ch3_leak, titleofplot, l, particleenergies[l], fitres, dogaussfit);

      //RecoE: Cheat 4
      fitres = extractMeanEnergy(mycE_ch4_1, mycE_ch4_2, color, hist32, rebin[l]);
      //For the linearity plot
      titleofplot1 = "Reconstructed energy but cheat section 4 for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c139,recoE_ch4vsbeamenergy,titleofplot,"Beam Energy (GeV)","E_{reco,ch4} (GeV)", l, particleenergies[l], fitres, dogaussfit);

      //Resolution
      titleofplot1 = "Relative energy resolution but cheat section 4 for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildresolutiongraph(c140, resoFit_ch4, titleofplot, l, particleenergies[l], fitres, dogaussfit);

      //RecoE including leakage: Cheat 4
      fitres = extractMeanEnergy(mycE_ch4_leak_1, mycE_ch4_leak_2, color, hist33, rebin[l]);
      //For the linearity plot
      titleofplot1 = "Reconstructed energy but cheat section 4 including leakage for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c141,recoE_ch4_leakvsbeamenergy,titleofplot,"Beam Energy (GeV)","E_{reco,ch4,leak} (GeV)", l, particleenergies[l], fitres, dogaussfit);

      //Resolution
      titleofplot1 = "Relative energy resolution but cheat section 4 including leakage for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildresolutiongraph(c142, resoFit_ch4_leak, titleofplot, l, particleenergies[l], fitres, dogaussfit);
      //-------------------------------------------------------------------------------------------------
      //-------------------------------------------------------------------------------------------------
      //-------------------------------------------------------------------------------------------------
      //-------------------------------------------------------------------------------------------------
      //RecoE overtruE: Cheat 1
      fitres = extractMeanEnergy(mycE_ch1_overtrueE_1, mycE_ch1_overtrueE_2, color, hist34, rebin[l]);
      //For the linearity plot
      titleofplot1 = "Reconstructed energy minus true energy over true energy but cheat section 1 for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c170,recoE_ch1_overtrueEvsbeamenergy,titleofplot,"Beam Energy (GeV)","#frac{E_{reco,ch1}-E_{true}}{E_{true}}", l, particleenergies[l], fitres, dogaussfit);

      //RecoE overtruE including leakage: Cheat 1
      fitres = extractMeanEnergy(mycE_ch1_overtrueE_leak_1, mycE_ch1_overtrueE_leak_2, color, hist35, rebin[l]);
      //For the linearity plot
      titleofplot1 = "Reconstructed energy minus true energy over true energy but cheat section 1 including leakage for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c171,recoE_ch1_overtrueE_leakvsbeamenergy,titleofplot,"Beam Energy (GeV)","#frac{E_{reco,ch1,leak}-E_{true}}{E_{true}}", l, particleenergies[l], fitres, dogaussfit);

      //RecoE overtruE: Cheat 2
      fitres = extractMeanEnergy(mycE_ch2_overtrueE_1, mycE_ch2_overtrueE_2, color, hist36, rebin[l]);
      //For the linearity plot
      titleofplot1 = "Reconstructed energy minus true energy over true energy but cheat section 2 for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c172,recoE_ch2_overtrueEvsbeamenergy,titleofplot,"Beam Energy (GeV)","#frac{E_{reco,ch2}-E_{true}}{E_{true}}", l, particleenergies[l], fitres, dogaussfit);

      //RecoE overtruE including leakage: Cheat 2
      fitres = extractMeanEnergy(mycE_ch2_overtrueE_leak_1, mycE_ch2_overtrueE_leak_2, color, hist37, rebin[l]);
      //For the linearity plot
      titleofplot1 = "Reconstructed energy minus true energy over true energy but cheat section 2 including leakage for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c173,recoE_ch2_overtrueE_leakvsbeamenergy,titleofplot,"Beam Energy (GeV)","#frac{E_{reco,ch2,leak}-E_{true}}{E_{true}}", l, particleenergies[l], fitres, dogaussfit);

      //RecoE overtruE: Cheat 3
      fitres = extractMeanEnergy(mycE_ch3_overtrueE_1, mycE_ch3_overtrueE_2, color, hist38, rebin[l]);
      //For the linearity plot
      titleofplot1 = "Reconstructed energy minus true energy over true energy but cheat section 3 for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c174,recoE_ch3_overtrueEvsbeamenergy,titleofplot,"Beam Energy (GeV)","#frac{E_{reco,ch3}-E_{true}}{E_{true}}", l, particleenergies[l], fitres, dogaussfit);

      //RecoE overtruE including leakage: Cheat 3
      fitres = extractMeanEnergy(mycE_ch3_overtrueE_leak_1, mycE_ch3_overtrueE_leak_2, color, hist39, rebin[l]);
      //For the linearity plot
      titleofplot1 = "Reconstructed energy minus true energy over true energy but cheat section 3 including leakage for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c175,recoE_ch3_overtrueE_leakvsbeamenergy,titleofplot,"Beam Energy (GeV)","#frac{E_{reco,ch3,leak}-E_{true}}{E_{true}}", l, particleenergies[l], fitres, dogaussfit);

      //RecoE overtruE: Cheat 4
      fitres = extractMeanEnergy(mycE_ch4_overtrueE_1, mycE_ch4_overtrueE_2, color, hist40, rebin[l]);
      //For the linearity plot
      titleofplot1 = "Reconstructed energy minus true energy over true energy but cheat section 4 for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c176,recoE_ch4_overtrueEvsbeamenergy,titleofplot,"Beam Energy (GeV)","#frac{E_{reco,ch4}-E_{true}}{E_{true}}", l, particleenergies[l], fitres, dogaussfit);

      //RecoE overtruE including leakage: Cheat 4
      fitres = extractMeanEnergy(mycE_ch4_overtrueE_leak_1, mycE_ch4_overtrueE_leak_2, color, hist41, rebin[l]);
      //For the linearity plot
      titleofplot1 = "Reconstructed energy minus true energy over true energy but cheat section 4 including leakage for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c177,recoE_ch4_overtrueE_leakvsbeamenergy,titleofplot,"Beam Energy (GeV)","#frac{E_{reco,ch4,leak}-E_{true}}{E_{true}}", l, particleenergies[l], fitres, dogaussfit);

      //-------------------------------------------------------------------------------------------------
      //-------------------------------------------------------------------------------------------------
      //-------------------------------------------------------------------------------------------------
      //-------------------------------------------------------------------------------------------------

      //Using the total sampling fraction
      fitres = extractMeanEnergy(mycE_totsf_1, mycE_totsf_2, color, hist12, rebin[l]);
      //For the linearity plot
      titleofplot1 = "Reconstructed energy using the total sampling fraction for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c46,recoE_totsfvsbeamenergy,titleofplot,"Beam Energy (GeV)","E_{reco} (GeV)", l, particleenergies[l], fitres, dogaussfit);

      //Resolution
      titleofplot1 = "Relative energy resolution using the total sampling fraction for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildresolutiongraph(c47, resoFit_totsf, titleofplot, l, particleenergies[l], fitres, dogaussfit);

      //Using the total sampling fraction adding leakage
      fitres = extractMeanEnergy(mycE_totsf_leak_1, mycE_totsf_leak_2, color, hist13, rebin[l]);
      //For the linearity plot
      titleofplot1 = "Reconstructed energy using the total sampling fraction including leakage for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c54,recoE_totsf_leakvsbeamenergy,titleofplot,"Beam Energy (GeV)","E_{reco} (GeV)", l, particleenergies[l], fitres, dogaussfit);

      //Resolution
      titleofplot1 = "Relative energy resolution using the total sampling fraction including leakage for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildresolutiongraph(c55, resoFit_totsf_leak, titleofplot, l, particleenergies[l], fitres, dogaussfit);




      //-------------------------------------------------------------------------------------------------

      //Using the total sampling fraction
      fitres = extractMeanEnergy(mycE_totsf_dedxweighted_1, mycE_totsf_dedxweighted_2, color, hist80, rebin[l]);
      //For the linearity plot
      titleofplot1 = "Reconstructed energy using the total sampling fraction using dEdx weight for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c318,recoE_totsf_dedxweightedvsbeamenergy,titleofplot,"Beam Energy (GeV)","E_{reco} (GeV)", l, particleenergies[l], fitres, dogaussfit);

      //Resolution
      titleofplot1 = "Relative energy resolution using the total sampling fraction using dEdx weight for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildresolutiongraph(c319, resoFit_totsf_dedxweighted, titleofplot, l, particleenergies[l], fitres, dogaussfit);

      //Using the total sampling fraction adding leakage
      fitres = extractMeanEnergy(mycE_totsf_dedxweighted_leak_1, mycE_totsf_dedxweighted_leak_2, color, hist81, rebin[l]);
      //For the linearity plot
      titleofplot1 = "Reconstructed energy using the total sampling fraction including leakage using dEdx weight for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c320,recoE_totsf_dedxweighted_leakvsbeamenergy,titleofplot,"Beam Energy (GeV)","E_{reco} (GeV)", l, particleenergies[l], fitres, dogaussfit);

      //Resolution
      titleofplot1 = "Relative energy resolution using the total sampling fraction using dEdx weight including leakage for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildresolutiongraph(c321, resoFit_totsf_dedxweighted_leak, titleofplot, l, particleenergies[l], fitres, dogaussfit);





      //-------------------------------------------------------------------------------------------------
      //Using the total sampling fraction per section 
      fitres = extractMeanEnergy(mycE_totsf_sec_1, mycE_totsf_sec_2, color, hist56, rebin[l]);
      //For the linearity plot
      titleofplot1 = "Reconstructed energy using the total sampling fraction per section for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c226,recoE_totsf_secvsbeamenergy,titleofplot,"Beam Energy (GeV)","E_{reco} (GeV)", l, particleenergies[l], fitres, dogaussfit);

      //Resolution
      titleofplot1 = "Relative energy resolution using the total sampling fraction per section for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildresolutiongraph(c227, resoFit_totsf_sec, titleofplot, l, particleenergies[l], fitres, dogaussfit);

      //Using the total sampling fraction per section adding leakage
      fitres = extractMeanEnergy(mycE_totsf_sec_leak_1, mycE_totsf_sec_leak_2, color, hist57, rebin[l]);
      //For the linearity plot
      titleofplot1 = "Reconstructed energy using the total sampling fraction per section including leakage for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c228,recoE_totsf_sec_leakvsbeamenergy,titleofplot,"Beam Energy (GeV)","E_{reco} (GeV)", l, particleenergies[l], fitres, dogaussfit);

      //Resolution
      titleofplot1 = "Relative energy resolution using the total sampling fraction per section including leakage for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildresolutiongraph(c229, resoFit_totsf_sec_leak, titleofplot, l, particleenergies[l], fitres, dogaussfit);
      //-------------------------------------------------------------------------------------------------

      //-------------------------------------------------------------------------------------------------
      //Using the total sampling fraction per section weighted with dedx
      fitres = extractMeanEnergy(mycE_totsf_sec_dedxweighted_1, mycE_totsf_sec_dedxweighted_2, color, hist72, rebin[l]);
      //For the linearity plot
      titleofplot1 = "Reconstructed energy using the total sampling fraction per section weighted with dedx for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c288,recoE_totsf_sec_dedxweightedvsbeamenergy,titleofplot,"Beam Energy (GeV)","E_{reco} (GeV)", l, particleenergies[l], fitres, dogaussfit);

      //Resolution
      titleofplot1 = "Relative energy resolution using the total sampling fraction per section weighted with dedx for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildresolutiongraph(c289, resoFit_totsf_sec_dedxweighted, titleofplot, l, particleenergies[l], fitres, dogaussfit);

      //Using the total sampling fraction per section weighted with dedx adding leakage
      fitres = extractMeanEnergy(mycE_totsf_sec_dedxweighted_leak_1, mycE_totsf_sec_dedxweighted_leak_2, color, hist73, rebin[l]);
      //For the linearity plot
      titleofplot1 = "Reconstructed energy using the total sampling fraction per section weighted with dedx including leakage for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c290,recoE_totsf_sec_dedxweighted_leakvsbeamenergy,titleofplot,"Beam Energy (GeV)","E_{reco} (GeV)", l, particleenergies[l], fitres, dogaussfit);

      //Resolution
      titleofplot1 = "Relative energy resolution using the total sampling fraction per section weighted with dedx including leakage for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildresolutiongraph(c291, resoFit_totsf_sec_dedxweighted_leak, titleofplot, l, particleenergies[l], fitres, dogaussfit);





      //-------------------------------------------------------------------------------------------------
      //Using the total sampling fraction per section and shower max method 
      fitres = extractMeanEnergy(mycE_totsfshmax_sec_1, mycE_totsfshmax_sec_2, color, hist60, rebin[l]);
      //For the linearity plot
      titleofplot1 = "Reconstructed energy using the total SF per section and shower max method for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c242,recoE_totsfshmax_secvsbeamenergy,titleofplot,"Beam Energy (GeV)","E_{reco} (GeV)", l, particleenergies[l], fitres, dogaussfit);

      //Resolution
      titleofplot1 = "Relative energy resolution using the total SF per section and shower max method for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildresolutiongraph(c243, resoFit_totsfshmax_sec, titleofplot, l, particleenergies[l], fitres, dogaussfit);

      //Using the total sampling fraction per section and shower max method adding leakage
      fitres = extractMeanEnergy(mycE_totsfshmax_sec_leak_1, mycE_totsfshmax_sec_leak_2, color, hist61, rebin[l]);
      //For the linearity plot
      titleofplot1 = "Reconstructed energy using the total SF per section and shower max method including leakage for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c244,recoE_totsfshmax_sec_leakvsbeamenergy,titleofplot,"Beam Energy (GeV)","E_{reco} (GeV)", l, particleenergies[l], fitres, dogaussfit);

      //Resolution
      titleofplot1 = "Relative energy resolution using the total SF per section and shower max method including leakage for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildresolutiongraph(c245, resoFit_totsfshmax_sec_leak, titleofplot, l, particleenergies[l], fitres, dogaussfit);
      //-------------------------------------------------------------------------------------------------




      //-------------------------------------------------------------------------------------------------
      //Using the total sampling fraction per section and shower max method 
      fitres = extractMeanEnergy(mycE_totsfshmax_sec_dedxweighted_1, mycE_totsfshmax_sec_dedxweighted_2, color, hist76, rebin[l]);
      //For the linearity plot
      titleofplot1 = "Reconstructed energy using the total SF per section and shower max method for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c303,recoE_totsfshmax_sec_dedxweightedvsbeamenergy,titleofplot,"Beam Energy (GeV)","E_{reco} (GeV)", l, particleenergies[l], fitres, dogaussfit);

      //Resolution
      titleofplot1 = "Relative energy resolution using the total SF per section and shower max method for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildresolutiongraph(c304, resoFit_totsfshmax_sec_dedxweighted, titleofplot, l, particleenergies[l], fitres, dogaussfit);

      //Using the total sampling fraction per section and shower max method adding leakage
      fitres = extractMeanEnergy(mycE_totsfshmax_sec_dedxweighted_leak_1, mycE_totsfshmax_sec_dedxweighted_leak_2, color, hist77, rebin[l]);
      //For the linearity plot
      titleofplot1 = "Reconstructed energy using the total SF per section and shower max method including leakage for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c305,recoE_totsfshmax_sec_dedxweighted_leakvsbeamenergy,titleofplot,"Beam Energy (GeV)","E_{reco} (GeV)", l, particleenergies[l], fitres, dogaussfit);

      //Resolution
      titleofplot1 = "Relative energy resolution using the total SF per section and shower max method including leakage for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildresolutiongraph(c306, resoFit_totsfshmax_sec_dedxweighted_leak, titleofplot, l, particleenergies[l], fitres, dogaussfit);
      //-------------------------------------------------------------------------------------------------





      //-------------------------------------------------------------------------------------------------
      //Using the total sampling fraction per section 
      fitres = extractMeanEnergy(mycE_totsf_dedxsec_1, mycE_totsf_dedxsec_2, color, hist64, rebin[l]);
      //For the linearity plot
      titleofplot1 = "Reconstructed energy using the total sampling fraction per section for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c261,recoE_totsf_dedxsecvsbeamenergy,titleofplot,"Beam Energy (GeV)","E_{reco} (GeV)", l, particleenergies[l], fitres, dogaussfit);

      //Resolution
      titleofplot1 = "Relative energy resolution using the total sampling fraction per section for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildresolutiongraph(c262, resoFit_totsf_dedxsec, titleofplot, l, particleenergies[l], fitres, dogaussfit);

      //Using the total sampling fraction per section adding leakage
      fitres = extractMeanEnergy(mycE_totsf_dedxsec_leak_1, mycE_totsf_dedxsec_leak_2, color, hist65, rebin[l]);
      //For the linearity plot
      titleofplot1 = "Reconstructed energy using the total sampling fraction per section including leakage for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c263,recoE_totsf_dedxsec_leakvsbeamenergy,titleofplot,"Beam Energy (GeV)","E_{reco} (GeV)", l, particleenergies[l], fitres, dogaussfit);

      //Resolution
      titleofplot1 = "Relative energy resolution using the total sampling fraction per section including leakage for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildresolutiongraph(c264, resoFit_totsf_dedxsec_leak, titleofplot, l, particleenergies[l], fitres, dogaussfit);
      //-------------------------------------------------------------------------------------------------
      //Using the total sampling fraction per section and shower max method 
      fitres = extractMeanEnergy(mycE_totsfshmax_dedxsec_1, mycE_totsfshmax_dedxsec_2, color, hist66, rebin[l]);
      //For the linearity plot
      titleofplot1 = "Reconstructed energy using the total SF per section and shower max method for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c265,recoE_totsfshmax_dedxsecvsbeamenergy,titleofplot,"Beam Energy (GeV)","E_{reco} (GeV)", l, particleenergies[l], fitres, dogaussfit);

      //Resolution
      titleofplot1 = "Relative energy resolution using the total SF per section and shower max method for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildresolutiongraph(c266, resoFit_totsfshmax_dedxsec, titleofplot, l, particleenergies[l], fitres, dogaussfit);

      //Using the total sampling fraction per section and shower max method adding leakage
      fitres = extractMeanEnergy(mycE_totsfshmax_dedxsec_leak_1, mycE_totsfshmax_dedxsec_leak_2, color, hist67, rebin[l]);
      //For the linearity plot
      titleofplot1 = "Reconstructed energy using the total SF per section and shower max method including leakage for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c267,recoE_totsfshmax_dedxsec_leakvsbeamenergy,titleofplot,"Beam Energy (GeV)","E_{reco} (GeV)", l, particleenergies[l], fitres, dogaussfit);

      //Resolution
      titleofplot1 = "Relative energy resolution using the total SF per section and shower max method including leakage for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildresolutiongraph(c268, resoFit_totsfshmax_dedxsec_leak, titleofplot, l, particleenergies[l], fitres, dogaussfit);
      //-------------------------------------------------------------------------------------------------






















      //-------------------------------------------------------------------------------------------------
      //Using the shower max method
      fitres = extractMeanEnergy(mycE_shmax_1, mycE_shmax_2, color, hist46, rebin[l]);
      //For the linearity plot
      titleofplot1 = "Reconstructed energy using SF_{i} and the shower max method for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c190,recoE_shmaxvsbeamenergy,titleofplot,"Beam Energy (GeV)","E_{reco} (GeV)", l, particleenergies[l], fitres, dogaussfit);

      //Resolution
      titleofplot1 = "Relative energy resolution using SF_{i} and the shower max method for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildresolutiongraph(c192, resoFit_shmax, titleofplot, l, particleenergies[l], fitres, dogaussfit);

      //Using the shower max method adding leakage
      fitres = extractMeanEnergy(mycE_shmax_leak_1, mycE_shmax_leak_2, color, hist47, rebin[l]);
      //For the linearity plot
      titleofplot1 = "Reconstructed energy using SF_{i} and and the shower max method including leakage for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c191,recoE_shmax_leakvsbeamenergy,titleofplot,"Beam Energy (GeV)","E_{reco} (GeV)", l, particleenergies[l], fitres, dogaussfit);

      //Resolution
      titleofplot1 = "Relative energy resolution using SF_{i} and and the shower max method including leakage for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildresolutiongraph(c193, resoFit_shmax_leak, titleofplot, l, particleenergies[l], fitres, dogaussfit);
      //-------------------------------------------------------------------------------------------------
      //Using the sampling fraction with normalized shower depth
      fitres = extractMeanEnergy(mycE_sfnorm_1, mycE_sfnorm_2, color, hist14, rebin[l]);
      //For the linearity plot
      titleofplot1 = "Reconstructed energy using the sampling fraction with X_{0,max}/<X_{0,max}> for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c57,recoE_sfnormvsbeamenergy,titleofplot,"Beam Energy (GeV)","E_{reco} (GeV)", l, particleenergies[l], fitres, dogaussfit);

      //Resolution
      titleofplot1 = "Relative energy resolution using the sampling fraction with X_{0,max}/<X_{0,max}> for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildresolutiongraph(c58, resoFit_sfnorm, titleofplot, l, particleenergies[l], fitres, dogaussfit);

      //Using the sampling fraction with normalized shower depth adding leakage
      fitres = extractMeanEnergy(mycE_sfnorm_leak_1, mycE_sfnorm_leak_2, color, hist15, rebin[l]);
      //For the linearity plot
      titleofplot1 = "Reconstructed energy using the sampling fraction with X_{0,max}/<X_{0,max}> adding leakage for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c62,recoE_sfnorm_leakvsbeamenergy,titleofplot,"Beam Energy (GeV)","E_{reco} (GeV)", l, particleenergies[l], fitres, dogaussfit);

      //Resolution
      titleofplot1 = "Relative energy resolution using the sampling fraction with X_{0,max}/<X_{0,max}> adding leakage for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildresolutiongraph(c63, resoFit_sfnorm_leak, titleofplot, l, particleenergies[l], fitres, dogaussfit);

      //-------------------------------------------------------------------------------------------------
      //Using the total SF with normalized shower depth
      fitres = extractMeanEnergy(mycE_totsfnorm_1, mycE_totsfnorm_2, color, hist50, rebin[l]);
      //For the linearity plot
      titleofplot1 = "Reconstructed energy using total SF with X_{0,max}/<X_{0,max}> for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c205,recoE_totsfnormvsbeamenergy,titleofplot,"Beam Energy (GeV)","E_{reco,totsfshmax} (GeV)", l, particleenergies[l], fitres, dogaussfit);

      //Resolution
      titleofplot1 = "Relative energy resolution using total SF with X_{0,max}/<X_{0,max}> for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildresolutiongraph(c206, resoFit_totsfnorm, titleofplot, l, particleenergies[l], fitres, dogaussfit);

      //Using the total SF with normalized shower depth adding leakage
      fitres = extractMeanEnergy(mycE_totsfnorm_leak_1, mycE_totsfnorm_leak_2, color, hist51, rebin[l]);
      //For the linearity plot
      titleofplot1 = "Reconstructed energy using total SF with X_{0,max}/<X_{0,max}> adding leakage for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c207,recoE_totsfnorm_leakvsbeamenergy,titleofplot,"Beam Energy (GeV)","E_{reco,totsfshmax,leak} (GeV)", l, particleenergies[l], fitres, dogaussfit);

      //Resolution
      titleofplot1 = "Relative energy resolution using total SF with X_{0,max}/<X_{0,max}> adding leakage for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildresolutiongraph(c208, resoFit_totsfnorm_leak, titleofplot, l, particleenergies[l], fitres, dogaussfit);
      //-------------------------------------------------------------------------------------------------
     


      //-------------------------------------------------------------------------------------------------
      //Using the total SF with normalized shower depth
      fitres = extractMeanEnergy(mycE_totsfnorm_dedxweighted_1, mycE_totsfnorm_dedxweighted_2, color, hist84, rebin[l]);
      //For the linearity plot
      titleofplot1 = "Reconstructed energy using total SF with X_{0,max}/<X_{0,max}> using dedx weight for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c333,recoE_totsfnorm_dedxweightedvsbeamenergy,titleofplot,"Beam Energy (GeV)","E_{reco,totsfshmaxdedx} (GeV)", l, particleenergies[l], fitres, dogaussfit);

      //Resolution
      titleofplot1 = "Relative energy resolution using total SF with X_{0,max}/<X_{0,max}> using dedx weight for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildresolutiongraph(c334, resoFit_totsfnorm_dedxweighted, titleofplot, l, particleenergies[l], fitres, dogaussfit);

      //Using the total SF with normalized shower depth adding leakage
      fitres = extractMeanEnergy(mycE_totsfnorm_dedxweighted_leak_1, mycE_totsfnorm_dedxweighted_leak_2, color, hist85, rebin[l]);
      //For the linearity plot
      titleofplot1 = "Reconstructed energy using total SF with X_{0,max}/<X_{0,max}> using dedx weight adding leakage for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c335,recoE_totsfnorm_dedxweighted_leakvsbeamenergy,titleofplot,"Beam Energy (GeV)","E_{reco,totsfshmaxdedx,leak} (GeV)", l, particleenergies[l], fitres, dogaussfit);

      //Resolution
      titleofplot1 = "Relative energy resolution using total SF with X_{0,max}/<X_{0,max}> using dedx weight adding leakage for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildresolutiongraph(c336, resoFit_totsfnorm_dedxweighted_leak, titleofplot, l, particleenergies[l], fitres, dogaussfit);
      //-------------------------------------------------------------------------------------------------





      //For the linearity plot
      fitres = extractMeanEnergy(mycE_raw_1, mycE_raw_2, color, hist5, rebin_raw[l]);
      titleofplot1 = "Raw Reconstructed energy for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c10,recoE_rawvsbeamenergy,titleofplot,"Beam Energy (GeV)","Raw E_{reco} (MIPs)", l, particleenergies[l], fitres, dogaussfit);
 
      //Here we fit the (Erec-Etrue)/Etrue
      // fitres = extractMeanEnergy(mycE_raw_1, mycE_raw_2, color, hist8, rebin[l]);
      //Resolution
      titleofplot1 = "Relative energy resolution for raw energy for "; 
      titleofplot2 =  particle; 
      titleofplot = titleofplot1 + titleofplot2;
      buildresolutiongraph(c16, resoFit_raw, titleofplot, l, particleenergies[l], fitres, dogaussfit);

      //For the linearity plot
      fitres = extractMeanEnergy(mycE_raw_leak_1, mycE_raw_leak_2, color, hist9, rebin_raw_leak[l]);
      titleofplot1 = "Raw Reconstructed energy adding leakage for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c30,recoE_raw_leakvsbeamenergy,titleofplot,"Beam Energy (GeV)","Raw E_{reco} (MIPs)", l, particleenergies[l], fitres, dogaussfit);
 
      //Resolution
      titleofplot1 = "Relative energy resolution for raw energy adding leakage for "; 
      titleofplot2 =  particle; 
      titleofplot = titleofplot1 + titleofplot2;
      buildresolutiongraph(c31, resoFit_raw_leak, titleofplot, l, particleenergies[l], fitres, dogaussfit);

      if(upstream){
	//For the linearity plot
	fitres = extractMeanEnergy(mycE_corr_1, mycE_corr_2, color, hist2, rebin[l]);
	titleofplot1 = "Reconstructed energy with upstream material correction for "; 
	titleofplot2 =  particle +  " and all beam energies"; 
	titleofplot = titleofplot1 + titleofplot2;
	buildgraph(c4,recoE_corrvsbeamenergy,titleofplot,"Beam Energy (GeV)","E_{reco} (GeV)", l, particleenergies[l], fitres, dogaussfit);
 
	//Resolution
	titleofplot1 = "Relative energy resolution for energy corrected for upstream material for "; 
	titleofplot2 =  particle; 
	titleofplot = titleofplot1 + titleofplot2;
	buildresolutiongraph(c22, resoFit_corr, titleofplot, l, particleenergies[l], fitres, dogaussfit);
      }

      //-------------------------------------------------------------------------------------------------
      //Here we fit the (Erec-Etrue)/Etrue
      fitres = extractMeanEnergy(mycE_overtrue_1, mycE_overtrue_2, color, hist3, rebin[l]);
      upstream ? titleofplot1 = "Reconstructed and true energy without upstream material correction for " : titleofplot1 = "Reconstructed and true energy for "; 
      titleofplot2 =  particle; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c7,recoEovertrueEvsbeamenergy,titleofplot,"Beam Energy (GeV)","#frac{E_{reco}-E_{true}}{E_{true}}", l, particleenergies[l], fitres, dogaussfit);

      //Here we fit the (Erec-Etrue)/Etrue with leakage
      fitres = extractMeanEnergy(mycE_overtrue_leak_1, mycE_overtrue_leak_2, color, hist42, rebin[l]);
      upstream ? titleofplot1 = "Reconstructed and true energy without upstream material correction adding leakage for " : titleofplot1 = "Reconstructed and true energy adding leakage for "; 
      titleofplot2 =  particle; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c180,recoEovertrueE_leakvsbeamenergy,titleofplot,"Beam Energy (GeV)","#frac{E_{reco,leak}-E_{true}}{E_{true}}", l, particleenergies[l], fitres, dogaussfit);

      //Here we fit the (Erec-Etrue)/Etrue with leakage plus offset
      fitres = extractMeanEnergy(mycE_overtrue_leak_offset_1, mycE_overtrue_leak_offset_2, color, hist43, rebin[l]);
      upstream ? titleofplot1 = "Reconstructed and true energy without upstream material correction adding leakage and offset for " : titleofplot1 = "Reconstructed and true energy adding leakage and offset for "; 
      titleofplot2 =  particle; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c181,recoEovertrueE_leak_offsetvsbeamenergy,titleofplot,"Beam Energy (GeV)","#frac{E_{reco,leak,offset}-E_{true}}{E_{true}}", l, particleenergies[l], fitres, dogaussfit);
      //-------------------------------------------------------------------------------------------------
      //Here we fit the (Erec-Etrue)/Etrue with total SF
      fitres = extractMeanEnergy(mycE_overtrue_totsf_1, mycE_overtrue_totsf_2, color, hist44, rebin[l]);
      titleofplot1 = "Reconstructed and true energy with total sampling fraction for "; 
      titleofplot2 =  particle; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c182,recoEovertrueE_totsfvsbeamenergy,titleofplot,"Beam Energy (GeV)","#frac{E_{reco,totsf}-E_{true}}{E_{true}}", l, particleenergies[l], fitres, dogaussfit);

      //Here we fit the (Erec-Etrue)/Etrue  with total SF with leakage
      fitres = extractMeanEnergy(mycE_overtrue_totsf_leak_1, mycE_overtrue_totsf_leak_2, color, hist45, rebin[l]);
      titleofplot1 = "Reconstructed and true energy  with total sampling fraction adding leakage for "; 
      titleofplot2 =  particle; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c183,recoEovertrueE_totsf_leakvsbeamenergy,titleofplot,"Beam Energy (GeV)","#frac{E_{reco,totsf,leak}-E_{true}}{E_{true}}", l, particleenergies[l], fitres, dogaussfit);

      //-------------------------------------------------------------------------------------------------
      //Here we fit the (Erec-Etrue)/Etrue with total SF weighted with dedx
      fitres = extractMeanEnergy(mycE_overtrue_totsf_dedxweighted_1, mycE_overtrue_totsf_dedxweighted_2, color, hist82, rebin[l]);
      titleofplot1 = "Reconstructed and true energy with total sampling fraction using dEdx weight for "; 
      titleofplot2 =  particle; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c322,recoEovertrueE_totsf_dedxweightedvsbeamenergy,titleofplot,"Beam Energy (GeV)","#frac{E_{reco,totsf_dedxweighteddedx}-E_{true}}{E_{true}}", l, particleenergies[l], fitres, dogaussfit);

      //Here we fit the (Erec-Etrue)/Etrue  with total SF with leakage
      fitres = extractMeanEnergy(mycE_overtrue_totsf_dedxweighted_leak_1, mycE_overtrue_totsf_dedxweighted_leak_2, color, hist83, rebin[l]);
      titleofplot1 = "Reconstructed and true energy  with total sampling fraction using dEdx weight adding leakage for "; 
      titleofplot2 =  particle; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c323,recoEovertrueE_totsf_dedxweighted_leakvsbeamenergy,titleofplot,"Beam Energy (GeV)","#frac{E_{reco,totsf_dedxweighted,leak}-E_{true}}{E_{true}}", l, particleenergies[l], fitres, dogaussfit);
      //-------------------------------------------------------------------------------------------------


      //-------------------------------------------------------------------------------------------------
      //Here we fit the (Erec-Etrue)/Etrue with total SF per section 
      fitres = extractMeanEnergy(mycE_overtrue_totsf_sec_1, mycE_overtrue_totsf_sec_2, color, hist58, rebin[l]);
      titleofplot1 = "Reconstructed and true energy with total sampling fraction per section for "; 
      titleofplot2 =  particle; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c230,recoEovertrueE_totsf_secvsbeamenergy,titleofplot,"Beam Energy (GeV)","#frac{E_{reco,totsfsec}-E_{true}}{E_{true}}", l, particleenergies[l], fitres, dogaussfit);

      //Here we fit the (Erec-Etrue)/Etrue  with total SF with leakage
      fitres = extractMeanEnergy(mycE_overtrue_totsf_sec_leak_1, mycE_overtrue_totsf_sec_leak_2, color, hist59, rebin[l]);
      titleofplot1 = "Reconstructed and true energy  with total sampling fraction per section adding leakage for "; 
      titleofplot2 =  particle; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c231,recoEovertrueE_totsf_sec_leakvsbeamenergy,titleofplot,"Beam Energy (GeV)","#frac{E_{reco,totsfsec,leak}-E_{true}}{E_{true}}", l, particleenergies[l], fitres, dogaussfit);

      //-------------------------------------------------------------------------------------------------
      //Here we fit the (Erec-Etrue)/Etrue with total SF per section weighted with dedx
      fitres = extractMeanEnergy(mycE_overtrue_totsf_sec_dedxweighted_1, mycE_overtrue_totsf_sec_dedxweighted_2, color, hist74, rebin[l]);
      titleofplot1 = "Reconstructed and true energy with total sampling fraction per section for "; 
      titleofplot2 =  particle; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c292,recoEovertrueE_totsf_sec_dedxweightedvsbeamenergy,titleofplot,"Beam Energy (GeV)","#frac{E_{reco,totsfsecdedx}-E_{true}}{E_{true}}", l, particleenergies[l], fitres, dogaussfit);

      //Here we fit the (Erec-Etrue)/Etrue  with total SF weighted with dedx with leakage
      fitres = extractMeanEnergy(mycE_overtrue_totsf_sec_dedxweighted_leak_1, mycE_overtrue_totsf_sec_dedxweighted_leak_2, color, hist75, rebin[l]);
      titleofplot1 = "Reconstructed and true energy  with total sampling fraction per section adding leakage for "; 
      titleofplot2 =  particle; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c293,recoEovertrueE_totsf_sec_dedxweighted_leakvsbeamenergy,titleofplot,"Beam Energy (GeV)","#frac{E_{reco,totsfsecdedx,leak}-E_{true}}{E_{true}}", l, particleenergies[l], fitres, dogaussfit);

       //-------------------------------------------------------------------------------------------------
      //Here we fit the (Erec-Etrue)/Etrue with total SF per section and shower max method 
      fitres = extractMeanEnergy(mycE_overtrue_totsfshmax_sec_1, mycE_overtrue_totsfshmax_sec_2, color, hist62, rebin[l]);
      titleofplot1 = "Reconstructed and true energy with total SF per section and shower max method for "; 
      titleofplot2 =  particle; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c246,recoEovertrueE_totsfshmax_secvsbeamenergy,titleofplot,"Beam Energy (GeV)","#frac{E_{reco,totsfshmaxsec}-E_{true}}{E_{true}}", l, particleenergies[l], fitres, dogaussfit);

      //Here we fit the (Erec-Etrue)/Etrue  with total SF with and shower max method leakage
      fitres = extractMeanEnergy(mycE_overtrue_totsfshmax_sec_leak_1, mycE_overtrue_totsfshmax_sec_leak_2, color, hist63, rebin[l]);
      titleofplot1 = "Reconstructed and true energy  with total SF per section and shower max method adding leakage for "; 
      titleofplot2 =  particle; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c247,recoEovertrueE_totsfshmax_sec_leakvsbeamenergy,titleofplot,"Beam Energy (GeV)","#frac{E_{reco,totsfshmaxsec,leak}-E_{true}}{E_{true}}", l, particleenergies[l], fitres, dogaussfit);

      //-------------------------------------------------------------------------------------------------

      //-------------------------------------------------------------------------------------------------
      //Here we fit the (Erec-Etrue)/Etrue with total SF per section and shower max method 
      fitres = extractMeanEnergy(mycE_overtrue_totsfshmax_sec_dedxweighted_1, mycE_overtrue_totsfshmax_sec_dedxweighted_2, color, hist78, rebin[l]);
      titleofplot1 = "Reconstructed and true energy with total SF per section and shower max method for "; 
      titleofplot2 =  particle; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c307,recoEovertrueE_totsfshmax_sec_dedxweightedvsbeamenergy,titleofplot,"Beam Energy (GeV)","#frac{E_{reco,totsfshmaxsecdedx}-E_{true}}{E_{true}}", l, particleenergies[l], fitres, dogaussfit);

      //Here we fit the (Erec-Etrue)/Etrue  with total SF with and shower max method leakage
      fitres = extractMeanEnergy(mycE_overtrue_totsfshmax_sec_dedxweighted_leak_1, mycE_overtrue_totsfshmax_sec_dedxweighted_leak_2, color, hist79, rebin[l]);
      titleofplot1 = "Reconstructed and true energy  with total SF per section and shower max method adding leakage for "; 
      titleofplot2 =  particle; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c308,recoEovertrueE_totsfshmax_sec_dedxweighted_leakvsbeamenergy,titleofplot,"Beam Energy (GeV)","#frac{E_{reco,totsfshmaxsecdedx,leak}-E_{true}}{E_{true}}", l, particleenergies[l], fitres, dogaussfit);

      //-------------------------------------------------------------------------------------------------





      //-------------------------------------------------------------------------------------------------
      //Here we fit the (Erec-Etrue)/Etrue with total SF per section 
      fitres = extractMeanEnergy(mycE_overtrue_totsf_dedxsec_1, mycE_overtrue_totsf_dedxsec_2, color, hist68, rebin[l]);
      titleofplot1 = "Reconstructed and true energy with total sampling fraction per section for "; 
      titleofplot2 =  particle; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c269,recoEovertrueE_totsf_dedxsecvsbeamenergy,titleofplot,"Beam Energy (GeV)","#frac{E_{reco,totsfdedxsec}-E_{true}}{E_{true}}", l, particleenergies[l], fitres, dogaussfit);

      //Here we fit the (Erec-Etrue)/Etrue  with total SF with leakage
      fitres = extractMeanEnergy(mycE_overtrue_totsf_dedxsec_leak_1, mycE_overtrue_totsf_dedxsec_leak_2, color, hist69, rebin[l]);
      titleofplot1 = "Reconstructed and true energy  with total sampling fraction per section adding leakage for "; 
      titleofplot2 =  particle; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c270,recoEovertrueE_totsf_dedxsec_leakvsbeamenergy,titleofplot,"Beam Energy (GeV)","#frac{E_{reco,totsfdedxsec,leak}-E_{true}}{E_{true}}", l, particleenergies[l], fitres, dogaussfit);

      //-------------------------------------------------------------------------------------------------
      //Here we fit the (Erec-Etrue)/Etrue with total SF per section and shower max method 
      fitres = extractMeanEnergy(mycE_overtrue_totsfshmax_dedxsec_1, mycE_overtrue_totsfshmax_dedxsec_2, color, hist70, rebin[l]);
      titleofplot1 = "Reconstructed and true energy with total SF per section and shower max method for "; 
      titleofplot2 =  particle; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c271,recoEovertrueE_totsfshmax_dedxsecvsbeamenergy,titleofplot,"Beam Energy (GeV)","#frac{E_{reco,totsfshmaxdedxsec}-E_{true}}{E_{true}}", l, particleenergies[l], fitres, dogaussfit);

      //Here we fit the (Erec-Etrue)/Etrue  with total SF with and shower max method leakage
      fitres = extractMeanEnergy(mycE_overtrue_totsfshmax_dedxsec_leak_1, mycE_overtrue_totsfshmax_dedxsec_leak_2, color, hist71, rebin[l]);
      titleofplot1 = "Reconstructed and true energy  with total SF per section and shower max method adding leakage for "; 
      titleofplot2 =  particle; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c272,recoEovertrueE_totsfshmax_dedxsec_leakvsbeamenergy,titleofplot,"Beam Energy (GeV)","#frac{E_{reco,totsfshmaxdedxsec,leak}-E_{true}}{E_{true}}", l, particleenergies[l], fitres, dogaussfit);

      //-------------------------------------------------------------------------------------------------











      //-------------------------------------------------------------------------------------------------
      //Here we fit the (Erec-Etrue)/Etrue with total SF with X0max/<X0max>
      fitres = extractMeanEnergy(mycE_overtrue_totsfnorm_1, mycE_overtrue_totsfnorm_2, color, hist52, rebin[l]);
      titleofplot1 = "Reconstructed and true energy with total SF using X_{0,max}/<X_{0,max}> for "; 
      titleofplot2 =  particle; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c209,recoEovertrueE_totsfnormvsbeamenergy,titleofplot,"Beam Energy (GeV)","#frac{E_{reco,totsfshmax}-E_{true}}{E_{true}}", l, particleenergies[l], fitres, dogaussfit);

      //Here we fit the (Erec-Etrue)/Etrue  with total SF with leakage
      fitres = extractMeanEnergy(mycE_overtrue_totsfnorm_leak_1, mycE_overtrue_totsfnorm_leak_2, color, hist53, rebin[l]);
      titleofplot1 = "Reconstructed and true energy  with total SF using X_{0,max}/<X_{0,max}> adding leakage for "; 
      titleofplot2 =  particle; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c210,recoEovertrueE_totsfnorm_leakvsbeamenergy,titleofplot,"Beam Energy (GeV)","#frac{E_{reco,totsfshmax,leak}-E_{true}}{E_{true}}", l, particleenergies[l], fitres, dogaussfit);



      //-------------------------------------------------------------------------------------------------
      //Here we fit the (Erec-Etrue)/Etrue with total SF with X0max/<X0max> using dedx weight
      fitres = extractMeanEnergy(mycE_overtrue_totsfnorm_dedxweighted_1, mycE_overtrue_totsfnorm_dedxweighted_2, color, hist86, rebin[l]);
      titleofplot1 = "Reconstructed and true energy with total SF using X_{0,max}/<X_{0,max}> and dedx weight for "; 
      titleofplot2 =  particle; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c337,recoEovertrueE_totsfnorm_dedxweightedvsbeamenergy,titleofplot,"Beam Energy (GeV)","#frac{E_{reco,totsfshmaxdedx}-E_{true}}{E_{true}}", l, particleenergies[l], fitres, dogaussfit);

      //Here we fit the (Erec-Etrue)/Etrue  with total SF using dedx weight with leakage
      fitres = extractMeanEnergy(mycE_overtrue_totsfnorm_dedxweighted_leak_1, mycE_overtrue_totsfnorm_dedxweighted_leak_2, color, hist87, rebin[l]);
      titleofplot1 = "Reconstructed and true energy  with total SF using X_{0,max}/<X_{0,max}> and dedx weight adding leakage for "; 
      titleofplot2 =  particle; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c338,recoEovertrueE_totsfnorm_dedxweighted_leakvsbeamenergy,titleofplot,"Beam Energy (GeV)","#frac{E_{reco,totsfshmaxdedx,leak}-E_{true}}{E_{true}}", l, particleenergies[l], fitres, dogaussfit);

 


      //-------------------------------------------------------------------------------------------------
      //Here we fit the (Erec-Etrue)/Etrue with  SF with X0max/<X0max>
      fitres = extractMeanEnergy(mycE_overtrue_sfnorm_1, mycE_overtrue_sfnorm_2, color, hist54, rebin[l]);
      titleofplot1 = "Reconstructed and true energy with SF per layer and energy using X_{0,max}/<X_{0,max}> for "; 
      titleofplot2 =  particle; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c218,recoEovertrueE_sfnormvsbeamenergy,titleofplot,"Beam Energy (GeV)","#frac{E_{reco,sfshmax}-E_{true}}{E_{true}}", l, particleenergies[l], fitres, dogaussfit);

      //Here we fit the (Erec-Etrue)/Etrue  with  SF with leakage
      fitres = extractMeanEnergy(mycE_overtrue_sfnorm_leak_1, mycE_overtrue_sfnorm_leak_2, color, hist55, rebin[l]);
      titleofplot1 = "Reconstructed and true energy  with SF per layer and energy using X_{0,max}/<X_{0,max}> adding leakage for "; 
      titleofplot2 =  particle; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c219,recoEovertrueE_sfnorm_leakvsbeamenergy,titleofplot,"Beam Energy (GeV)","#frac{E_{reco,sfshmax,leak}-E_{true}}{E_{true}}", l, particleenergies[l], fitres, dogaussfit);

      //-------------------------------------------------------------------------------------------------
      //Here we fit the (Erec-Etrue)/Etrue with shower max method
      fitres = extractMeanEnergy(mycE_overtrue_shmax_1, mycE_overtrue_shmax_2, color, hist48, rebin[l]);
      titleofplot1 = "Reconstructed and true energy with SF_{i} and shower max method for "; 
      titleofplot2 =  particle; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c194,recoEovertrueE_shmaxvsbeamenergy,titleofplot,"Beam Energy (GeV)","#frac{E_{reco,shmax}-E_{true}}{E_{true}}", l, particleenergies[l], fitres, dogaussfit);

      //Here we fit the (Erec-Etrue)/Etrue shower max method with leakage
      fitres = extractMeanEnergy(mycE_overtrue_shmax_leak_1, mycE_overtrue_shmax_leak_2, color, hist49, rebin[l]);
      titleofplot1 = "Reconstructed and true energy  with SF_{i} and shower max method adding leakage for "; 
      titleofplot2 =  particle; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c195,recoEovertrueE_shmax_leakvsbeamenergy,titleofplot,"Beam Energy (GeV)","#frac{E_{reco,shmax,leak}-E_{true}}{E_{true}}", l, particleenergies[l], fitres, dogaussfit);
      //-------------------------------------------------------------------------------------------------

        if(upstream){
	//Here we fit the (Erec-Etrue)/Etrue 
	fitres = extractMeanEnergy(mycE_overtrue_corr_1, mycE_overtrue_corr_2, color, hist4, rebin[l]);
	titleofplot1 = "Reconstructed energy with upstream material correction for "; 
	titleofplot2 =  particle; 
	titleofplot = titleofplot1 + titleofplot2;
	buildgraph(c8,recoE_corrovertrueEvsbeamenergy,titleofplot,"Beam Energy (GeV)","#frac{E_{reco}-E_{true}}{E_{true}}", l, particleenergies[l], fitres, dogaussfit);
      }


      //For the linearity plot 
      fitres = extractMeanEnergy(mycE_val_1, mycE_val_2, color, hist6, rebin_val[l]);
      titleofplot1 = "Reconstructed energy using dEdx correction for "; 
      titleofplot2 =  particle; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c21,recoE_valvsbeamenergy,titleofplot,"Beam Energy (GeV)","E_{reco} (GeV)", l, particleenergies[l], fitres, dogaussfit);

      //Resolution
      titleofplot1 = "Relative energy resolution for dedx weighting energy for "; 
      titleofplot2 =  particle; 
      titleofplot = titleofplot1 + titleofplot2;
      buildresolutiongraph(c19, resoFit_val, titleofplot, l, particleenergies[l], fitres, dogaussfit);

      //Here we fit the (Erec-Etrue)/Etrue
      fitres = extractMeanEnergy(mycE_overtrue_val_1, mycE_overtrue_val_2, color, hist7, rebin[l]);
      titleofplot1 = "Reconstructed energy using dEdx correction for "; 
      titleofplot2 =  particle; 
      titleofplot = titleofplot1 + titleofplot2;
      buildgraph(c27,recoE_valovertrueEvsbeamenergy,titleofplot,"Beam Energy (GeV)","#frac{E_{reco,dEdx}-E_{true}}{E_{true}}", l, particleenergies[l], fitres, dogaussfit);

 


    }//Loop on energies

    //Fit for linearity
    fitforlinearity(c11,recoEvsbeamenergy,0.,520.,"GeV");
    fitforlinearity(c48,recoE_totsfvsbeamenergy,0.,520.,"GeV");
    fitforlinearity(c51,recoE_totsf_leakvsbeamenergy,0.,520.,"GeV");
    fitforlinearity(c324,recoE_totsf_dedxweightedvsbeamenergy,0.,520.,"GeV");
    fitforlinearity(c325,recoE_totsf_dedxweighted_leakvsbeamenergy,0.,520.,"GeV");
    fitforlinearity(c232,recoE_totsf_secvsbeamenergy,0.,520.,"GeV");
    fitforlinearity(c233,recoE_totsf_sec_leakvsbeamenergy,0.,520.,"GeV");
    fitforlinearity(c294,recoE_totsf_sec_dedxweightedvsbeamenergy,0.,520.,"GeV");
    fitforlinearity(c295,recoE_totsf_sec_dedxweighted_leakvsbeamenergy,0.,520.,"GeV");
    fitforlinearity(c248,recoE_totsfshmax_secvsbeamenergy,0.,520.,"GeV");
    fitforlinearity(c249,recoE_totsfshmax_sec_leakvsbeamenergy,0.,520.,"GeV");
    fitforlinearity(c309,recoE_totsfshmax_sec_dedxweightedvsbeamenergy,0.,520.,"GeV");
    fitforlinearity(c310,recoE_totsfshmax_sec_dedxweighted_leakvsbeamenergy,0.,520.,"GeV");
    fitforlinearity(c273,recoE_totsf_dedxsecvsbeamenergy,0.,520.,"GeV");
    fitforlinearity(c274,recoE_totsf_dedxsec_leakvsbeamenergy,0.,520.,"GeV");
    fitforlinearity(c275,recoE_totsfshmax_dedxsecvsbeamenergy,0.,520.,"GeV");
    fitforlinearity(c276,recoE_totsfshmax_dedxsec_leakvsbeamenergy,0.,520.,"GeV");
    // fitforlinearity(c196,recoE_shmaxvsbeamenergy,0.,520.,"GeV");
    // fitforlinearity(c197,recoE_shmax_leakvsbeamenergy,0.,520.,"GeV");
    fitforlinearity(c59,recoE_sfnormvsbeamenergy,0.,520.,"GeV");
    fitforlinearity(c64,recoE_sfnorm_leakvsbeamenergy,0.,520.,"GeV");
    fitforlinearity(c211,recoE_totsfnormvsbeamenergy,0.,520.,"GeV");
    fitforlinearity(c212,recoE_totsfnorm_leakvsbeamenergy,0.,520.,"GeV");
    fitforlinearity(c339,recoE_totsfnorm_dedxweightedvsbeamenergy,0.,520.,"GeV");
    fitforlinearity(c340,recoE_totsfnorm_dedxweighted_leakvsbeamenergy,0.,520.,"GeV");
    fitforlinearity(c37,recoE_leakvsbeamenergy,0.,520.,"GeV");
    fitforlinearity(c42,recoE_leak_offsetvsbeamenergy,0.,520.,"GeV");
    fitforlinearity(c72,recoE_4sectionsvsbeamenergy,0.,520.,"GeV");
    fitforlinearity(c73,recoE_4sections_leakvsbeamenergy,0.,520.,"GeV");
    fitforlinearity(c94,recoE_4sections_ch1vsbeamenergy,0.,520.,"GeV");
    fitforlinearity(c95,recoE_4sections_ch1_leakvsbeamenergy,0.,520.,"GeV");
    fitforlinearity(c96,recoE_4sections_ch2vsbeamenergy,0.,520.,"GeV");
    fitforlinearity(c97,recoE_4sections_ch2_leakvsbeamenergy,0.,520.,"GeV");
    fitforlinearity(c98,recoE_4sections_ch3vsbeamenergy,0.,520.,"GeV");
    fitforlinearity(c99,recoE_4sections_ch3_leakvsbeamenergy,0.,520.,"GeV");
    fitforlinearity(c100,recoE_4sections_ch4vsbeamenergy,0.,520.,"GeV");
    fitforlinearity(c101,recoE_4sections_ch4_leakvsbeamenergy,0.,520.,"GeV");
    fitforlinearity(c143,recoE_ch1vsbeamenergy,0.,520.,"GeV");
    fitforlinearity(c144,recoE_ch1_leakvsbeamenergy,0.,520.,"GeV");
    fitforlinearity(c145,recoE_ch2vsbeamenergy,0.,520.,"GeV");
    fitforlinearity(c146,recoE_ch2_leakvsbeamenergy,0.,520.,"GeV");
    fitforlinearity(c147,recoE_ch3vsbeamenergy,0.,520.,"GeV");
    fitforlinearity(c148,recoE_ch3_leakvsbeamenergy,0.,520.,"GeV");
    fitforlinearity(c149,recoE_ch4vsbeamenergy,0.,520.,"GeV");
    fitforlinearity(c150,recoE_ch4_leakvsbeamenergy,0.,520.,"GeV");
    fitforlinearity(c12,recoE_rawvsbeamenergy,0.,520.,"MIP");
    fitforlinearity(c32,recoE_raw_leakvsbeamenergy,0.,520.,"MIP");
    if(upstream){
      fitforlinearity(c13,recoE_corrvsbeamenergy,0.,520.,"GeV");
    }
    fitforlinearity(c24,recoE_valvsbeamenergy,0.,520.,"GeV");

    //Fit for resolution
    bool addNoiseTerm = false;
    int color_reco = 1;
    int color_sim = 2;
    int color_reco_val = 3;
    int color_reco_corr = 4;
    fitresol = fitforresolution(c15,resoFit,color_reco, addNoiseTerm);
    fitresol_totsf = fitforresolution(c49,resoFit_totsf,color_reco, addNoiseTerm);
    fitresol_totsf_leak = fitforresolution(c50,resoFit_totsf_leak,color_reco, addNoiseTerm);
    fitresol_totsf_dedxweighted = fitforresolution(c326,resoFit_totsf_dedxweighted,color_reco, addNoiseTerm);
    fitresol_totsf_dedxweighted_leak = fitforresolution(c327,resoFit_totsf_dedxweighted_leak,color_reco, addNoiseTerm);
    fitresol_totsf_sec = fitforresolution(c234,resoFit_totsf_sec,color_reco, addNoiseTerm);
    fitresol_totsf_sec_leak = fitforresolution(c235,resoFit_totsf_sec_leak,color_reco, addNoiseTerm);
    fitresol_totsf_sec_dedxweighted = fitforresolution(c296,resoFit_totsf_sec_dedxweighted,color_reco, addNoiseTerm);
    fitresol_totsf_sec_dedxweighted_leak = fitforresolution(c297,resoFit_totsf_sec_dedxweighted_leak,color_reco, addNoiseTerm);
    fitresol_totsfshmax_sec = fitforresolution(c250,resoFit_totsfshmax_sec,color_reco, addNoiseTerm);
    fitresol_totsfshmax_sec_leak = fitforresolution(c251,resoFit_totsfshmax_sec_leak,color_reco, addNoiseTerm);
    fitresol_totsfshmax_sec_dedxweighted = fitforresolution(c311,resoFit_totsfshmax_sec_dedxweighted,color_reco, addNoiseTerm);
    fitresol_totsfshmax_sec_dedxweighted_leak = fitforresolution(c312,resoFit_totsfshmax_sec_dedxweighted_leak,color_reco, addNoiseTerm);
    fitresol_totsf_dedxsec = fitforresolution(c277,resoFit_totsf_dedxsec,color_reco, addNoiseTerm);
    fitresol_totsf_dedxsec_leak = fitforresolution(c278,resoFit_totsf_dedxsec_leak,color_reco, addNoiseTerm);
    fitresol_totsfshmax_dedxsec = fitforresolution(c279,resoFit_totsfshmax_dedxsec,color_reco, addNoiseTerm);
    fitresol_totsfshmax_dedxsec_leak = fitforresolution(c280,resoFit_totsfshmax_dedxsec_leak,color_reco, addNoiseTerm);
    fitresol_shmax = fitforresolution(c198,resoFit_shmax,color_reco, addNoiseTerm);
    fitresol_shmax_leak = fitforresolution(c199,resoFit_shmax_leak,color_reco, addNoiseTerm);
    fitresol_sfnorm = fitforresolution(c60,resoFit_sfnorm,color_reco, addNoiseTerm);
    fitresol_sfnorm_leak = fitforresolution(c65,resoFit_sfnorm_leak,color_reco, addNoiseTerm);
    fitresol_totsfnorm = fitforresolution(c213,resoFit_totsfnorm,color_reco, addNoiseTerm);
    fitresol_totsfnorm_leak = fitforresolution(c214,resoFit_totsfnorm_leak,color_reco, addNoiseTerm);
    fitresol_totsfnorm_dedxweighted = fitforresolution(c341,resoFit_totsfnorm_dedxweighted,color_reco, addNoiseTerm);
    fitresol_totsfnorm_dedxweighted_leak = fitforresolution(c342,resoFit_totsfnorm_dedxweighted_leak,color_reco, addNoiseTerm);
    fitresol_leak = fitforresolution(c38,resoFit_leak,color_reco, addNoiseTerm);
    fitresol_leak_offset = fitforresolution(c43,resoFit_leak_offset,color_reco, addNoiseTerm);
    fitresol_4sections = fitforresolution(c74,resoFit_4sections,color_reco, addNoiseTerm);
    fitresol_4sections_leak = fitforresolution(c75,resoFit_4sections_leak,color_reco, addNoiseTerm);
    fitresol_4sections_ch1 = fitforresolution(c102,resoFit_4sections_ch1,color_reco, addNoiseTerm);
    fitresol_4sections_ch1_leak = fitforresolution(c103,resoFit_4sections_ch1_leak,color_reco, addNoiseTerm);
    fitresol_4sections_ch2 = fitforresolution(c104,resoFit_4sections_ch2,color_reco, addNoiseTerm);
    fitresol_4sections_ch2_leak = fitforresolution(c105,resoFit_4sections_ch2_leak,color_reco, addNoiseTerm);
    fitresol_4sections_ch3 = fitforresolution(c106,resoFit_4sections_ch3,color_reco, addNoiseTerm);
    fitresol_4sections_ch3_leak = fitforresolution(c107,resoFit_4sections_ch3_leak,color_reco, addNoiseTerm);
    fitresol_4sections_ch4 = fitforresolution(c108,resoFit_4sections_ch4,color_reco, addNoiseTerm);
    fitresol_4sections_ch4_leak = fitforresolution(c109,resoFit_4sections_ch4_leak,color_reco, addNoiseTerm);
    fitresol_ch1 = fitforresolution(c151,resoFit_ch1,color_reco, addNoiseTerm);
    fitresol_ch1_leak = fitforresolution(c152,resoFit_ch1_leak,color_reco, addNoiseTerm);
    fitresol_ch2 = fitforresolution(c153,resoFit_ch2,color_reco, addNoiseTerm);
    fitresol_ch2_leak = fitforresolution(c154,resoFit_ch2_leak,color_reco, addNoiseTerm);
    fitresol_ch3 = fitforresolution(c155,resoFit_ch3,color_reco, addNoiseTerm);
    fitresol_ch3_leak = fitforresolution(c156,resoFit_ch3_leak,color_reco, addNoiseTerm);
    fitresol_ch4 = fitforresolution(c157,resoFit_ch4,color_reco, addNoiseTerm);
    fitresol_ch4_leak = fitforresolution(c158,resoFit_ch4_leak,color_reco, addNoiseTerm);
    fitresol_raw = fitforresolution(c17,resoFit_raw, color_sim, addNoiseTerm);
    fitresol_raw_leak = fitforresolution(c33,resoFit_raw_leak, color_sim, addNoiseTerm);
    if(upstream){
      fitresol_corr = fitforresolution(c23,resoFit_corr,color_reco_corr, addNoiseTerm);
    }
    fitresol_val = fitforresolution(c20,resoFit_val,color_reco_val, addNoiseTerm);

    //Combined resolution plot
    if(upstream){
      combineresolution(c25, resoFit, color_reco, fitresol, resoFit_raw, color_sim, fitresol_raw, resoFit_corr, color_reco_val, fitresol_corr, resoFit_val, color_reco_val, fitresol_val, particle, addNoiseTerm);
    } else {
      combineresolution_noupstream(c25, resoFit_totsf, "E_{totsf}: ", color_sim, fitresol_totsf, resoFit, "E_{reco}: ", color_reco, fitresol, resoFit_val, "E_{reco,dEdx corr}: ", color_reco_val, fitresol_val, particle, addNoiseTerm);
    }
    combineresolution_noupstream(c44, resoFit, "E_{reco}: ", color_reco, fitresol, resoFit_leak, "E_{reco,leak}: ", color_sim, fitresol_leak, resoFit_leak_offset, "E_{reco,leak,offset}: ", color_reco_val, fitresol_leak_offset, particle, addNoiseTerm);
    combineresolution_noupstream(c52, resoFit_totsf, "E_{reco,totsf}: ", color_reco, fitresol_totsf, resoFit_totsf_leak, "E_{reco,totsf,leak}: ", color_sim, fitresol_totsf_leak, resoFit, "E_{reco}: ", color_reco_val, fitresol, particle, addNoiseTerm);

    combineresolution_noupstream(c76, resoFit_totsf, "E_{reco,totsf}: ", color_reco, fitresol_totsf, resoFit_4sections, "E_{reco,4 sections}: ", color_sim, fitresol_4sections, resoFit, "E_{reco}: ", color_reco_val, fitresol, particle, addNoiseTerm);

    combineresolution_noupstream(c77, resoFit_totsf_leak, "E_{reco,totsf,leak}: ", color_reco, fitresol_totsf_leak, resoFit_4sections_leak, "E_{reco,4 sections,leak}: ", color_sim, fitresol_4sections_leak, resoFit_leak, "E_{reco,leak}: ", color_reco_val, fitresol_leak, particle, addNoiseTerm);

    // combineresolution_noupstream(c110, resoFit_4sections_ch1, "E_{reco,4 sections,ch1}: ", color_reco, fitresol_4sections_ch1, resoFit_4sections_ch2, "E_{reco,4 sections,ch2}: ", color_sim, fitresol_4sections_ch2, resoFit_4sections_ch3, "E_{reco,4 sections,ch3}: ", color_reco_val, fitresol_4sections_ch3, particle, addNoiseTerm);
    combineresolution_noupstream_sections(c110, resoFit_4sections_ch1, "E_{reco,4 sections,ch1}: ", color_reco, fitresol_4sections_ch1, resoFit_4sections_ch2, "E_{reco,4 sections,ch2}: ", color_sim, fitresol_4sections_ch2, resoFit_4sections_ch3, "E_{reco,4 sections,ch3}: ", color_reco_val, fitresol_4sections_ch3, resoFit_4sections_ch4, "E_{reco,4 sections,ch4}: ", color_reco_corr, fitresol_4sections_ch4, particle, addNoiseTerm);

    // combineresolution_noupstream(c159, resoFit_ch1, "E_{reco,ch1}: ", color_reco, fitresol_ch1, resoFit_ch2, "E_{reco,ch2}: ", color_sim, fitresol_ch2, resoFit_ch3, "E_{reco,ch3}: ", color_reco_val, fitresol_ch3, particle, addNoiseTerm);

    combineresolution_noupstream_sections(c159, resoFit_ch1, "E_{reco,ch1}: ", color_reco, fitresol_ch1, resoFit_ch2, "E_{reco,ch2}: ", color_sim, fitresol_ch2, resoFit_ch3, "E_{reco,ch3}: ", color_reco_val, fitresol_ch3, resoFit_ch4, "E_{reco,ch4}: ", color_reco_corr, fitresol_ch4, particle, addNoiseTerm);

    combineresolution_noupstream_sections(c161, resoFit_ch1_leak, "E_{reco,ch1,leak}: ", color_reco, fitresol_ch1_leak, resoFit_ch2_leak, "E_{reco,ch2,leak}: ", color_sim, fitresol_ch2_leak, resoFit_ch3_leak, "E_{reco,ch3,leak}: ", color_reco_val, fitresol_ch3_leak, resoFit_ch4_leak, "E_{reco,ch4,leak}: ", color_reco_corr, fitresol_ch4_leak, particle, addNoiseTerm);

    combineresolution_noupstream(c200, resoFit_shmax, "E_{reco,shmax}: ", color_reco, fitresol_shmax, resoFit_shmax_leak, "E_{reco,shmax,leak}: ", color_sim, fitresol_shmax_leak, resoFit, "E_{reco}: ", color_reco_val, fitresol, particle, addNoiseTerm);

    combineresolution_noupstream(c215, resoFit_totsfnorm, "E_{reco,totsfshmax}: ", color_reco, fitresol_totsfnorm, resoFit_totsfnorm_leak, "E_{reco,totsfshmax,leak}: ", color_sim, fitresol_totsfnorm_leak, resoFit, "E_{reco}: ", color_reco_val, fitresol, particle, addNoiseTerm);

    combineresolution_noupstream_sections(c221, resoFit_totsfnorm, "E_{reco,totsfshmax}: ", color_reco, fitresol_totsfnorm, resoFit_totsfnorm_leak, "E_{reco,totsfshmax,leak}: ", color_sim, fitresol_totsfnorm_leak, resoFit_totsf, "E_{reco,totsf}: ", color_reco_val, fitresol_totsf, resoFit_totsf_leak, "E_{reco,totsf,leak}: ", color_reco_corr, fitresol_totsf_leak, particle, addNoiseTerm);

    combineresolution_noupstream(c220, resoFit_sfnorm, "E_{reco,sfshmax}: ", color_reco, fitresol_sfnorm, resoFit_sfnorm_leak, "E_{reco,sfshmax,leak}: ", color_sim, fitresol_sfnorm_leak, resoFit, "E_{reco}: ", color_reco_val, fitresol, particle, addNoiseTerm);

   combineresolution_noupstream(c236, resoFit_totsf_sec, "E_{reco,totsf_sec}: ", color_reco, fitresol_totsf_sec, resoFit_totsf_sec_leak, "E_{reco,totsf_sec,leak}: ", color_sim, fitresol_totsf_sec_leak, resoFit, "E_{reco}: ", color_reco_val, fitresol, particle, addNoiseTerm);

   combineresolution_noupstream(c237, resoFit_totsf_leak, "E_{reco,totsf,leak}: ", color_reco, fitresol_totsf_leak, resoFit_totsf_sec_leak, "E_{reco,totsf_sec,leak}: ", color_sim, fitresol_totsf_sec_leak, resoFit_val, "E_{reco,dEdx}: ", color_reco_val, fitresol_val, particle, addNoiseTerm);

   combineresolution_noupstream(c298, resoFit_totsf_leak, "E_{reco,totsf,leak}: ", color_reco, fitresol_totsf_leak, resoFit_totsf_sec_dedxweighted_leak, "E_{reco,totsf_sec_dedxweighted,leak}: ", color_sim, fitresol_totsf_sec_dedxweighted_leak, resoFit_val, "E_{reco,dEdx}: ", color_reco_val, fitresol_val, particle, addNoiseTerm);

   combineresolution_noupstream(c328, resoFit_totsf_dedxweighted_leak, "E_{reco,totsf_dedxweighted,leak}: ", color_reco, fitresol_totsf_dedxweighted_leak, resoFit_totsf_sec_dedxweighted_leak, "E_{reco,totsf_sec_dedxweighted,leak}: ", color_sim, fitresol_totsf_sec_dedxweighted_leak, resoFit_val, "E_{reco,dEdx}: ", color_reco_val, fitresol_val, particle, addNoiseTerm);

   combineresolution_noupstream(c252, resoFit_totsfnorm_leak, "E_{reco,totsfshmax,leak}: ", color_reco, fitresol_totsfnorm_leak, resoFit_totsfshmax_sec_leak, "E_{reco,totsfshmaxsec,leak}: ", color_sim, fitresol_totsfshmax_sec_leak, resoFit_val, "E_{reco,dEdx,leak}: ", color_reco_val, fitresol_val, particle, addNoiseTerm);

   combineresolution_noupstream(c313, resoFit_totsfnorm_dedxweighted_leak, "E_{reco,totsfshmaxdedx,leak}: ", color_reco, fitresol_totsfnorm_dedxweighted_leak, resoFit_totsfshmax_sec_dedxweighted_leak, "E_{reco,totsfshmaxsecdedx,leak}: ", color_sim, fitresol_totsfshmax_sec_dedxweighted_leak, resoFit_val, "E_{reco,dEdx,leak}: ", color_reco_val, fitresol_val, particle, addNoiseTerm);

   //------------------
   combineresolution_noupstream(c281, resoFit_totsf_dedxsec, "E_{reco,totsf_dedxsec}: ", color_reco, fitresol_totsf_dedxsec, resoFit_totsf_dedxsec_leak, "E_{reco,totsf_dedxsec,leak}: ", color_sim, fitresol_totsf_dedxsec_leak, resoFit, "E_{reco}: ", color_reco_val, fitresol, particle, addNoiseTerm);

   combineresolution_noupstream(c282, resoFit_totsf_leak, "E_{reco,totsf,leak}: ", color_reco, fitresol_totsf_leak, resoFit_totsf_dedxsec_leak, "E_{reco,totsf_dedxsec,leak}: ", color_sim, fitresol_totsf_dedxsec_leak, resoFit_val, "E_{reco,dEdx}: ", color_reco_val, fitresol_val, particle, addNoiseTerm);

   combineresolution_noupstream(c283, resoFit_totsfnorm_leak, "E_{reco,totsfshmax,leak}: ", color_reco, fitresol_totsfnorm_leak, resoFit_totsfshmax_dedxsec_leak, "E_{reco,totsfshmaxdedxsec,leak}: ", color_sim, fitresol_totsfshmax_dedxsec_leak, resoFit_val, "E_{reco,dEdx,leak}: ", color_reco_val, fitresol_val, particle, addNoiseTerm);
   






    results_com[k]= new TFile(res_com[k],"recreate");
    
    //======================================================================================
    //Write canvas with combined plot 
    c1->Write();
    c2->Write();
    c3->Write();
    c3->SaveAs("RecoEvsBeamEne.png");
    c4->Write();
    c4->SaveAs("RecoE_corrvsBeamEne.png");
    c5->Write();
    // c5->SaveAs("RecoEovertrueEvsBeamEne.png");
    c6->Write();
    // c6->SaveAs("RecoE_corrovertrueEvsBeamEne.png");
    c7->Write();
    c7->SaveAs("RecoEovertrueEvsBeamEne.png");
    c8->Write();
    c8->SaveAs("RecoE_corrovertrueEvsBeamEne.png");
    c9->Write();
    c10->Write();
    c10->SaveAs("RecoE_rawvsBeamEne.png");
    c11->Write();
    c11->SaveAs("RecoEvsBeamEne_plusfit.png");
    c12->Write();
    c12->SaveAs("RecoE_rawvsBeamEne_plusfit.png");
    c13->Write();
    c13->SaveAs("RecoE_corrvsBeamEne_plusfit.png");
    c14->Write();
    c14->SaveAs("Resolution.png");
    c15->Write();
    c15->SaveAs("Resolution_plusfit.png");
    c16->Write();
    c16->SaveAs("Resolution_RecoRaw.png");
    c17->Write();
    c17->SaveAs("Resolution_RecoRaw_plusfit.png");
    c18->Write();
    c18->SaveAs("RecoE_valvsBeamEne.png");
    c19->Write();
    c19->SaveAs("Resolution_RecoVal.png");
    c20->Write();
    c20->SaveAs("Resolution_RecoVal_plusfit.png");
    c21->Write();
    c21->SaveAs("RecoE_valvsBeamEne.png");
    c22->Write();
    c22->SaveAs("Resolution_Recocorr.png");
    c23->Write();
    c23->SaveAs("Resolution_Recocorr_plusfit.png");
    c24->Write();
    c24->SaveAs("RecoE_valvsBeamEne_plusfit.png");
    c25->Write();
    c25->SaveAs("Combinedresolution.png");
    c26->Write();
    c27->Write();
    c27->SaveAs("RecoE_valovertrueEvsBeamEne.png");
    c28->Write();
    c29->Write();
    c30->Write();
    c31->Write();
    c32->Write();
    c33->Write();
    c34->Write();
    c35->Write();
    c36->Write();
    c37->Write();
    c38->Write();
    c39->Write();
    c40->Write();
    c41->Write();
    c42->Write();
    c43->Write();
    c44->Write();
    c45->Write();
    c46->Write();
    c47->Write();
    c48->Write();
    c49->Write();
    c50->Write();
    c51->Write();
    c52->Write();
    c53->Write();
    c54->Write();
    c55->Write();
    c56->Write();
    c57->Write();
    c58->Write();
    c59->Write();
    c60->Write();
    c61->Write();
    c62->Write();
    c63->Write();
    c64->Write();
    c65->Write();
    c66->Write();
    c67->Write();
    c68->Write();
    c69->Write();
    c70->Write();
    c71->Write();
    c72->Write();
    c73->Write();
    c74->Write();
    c75->Write();
    c76->Write();
    c77->Write();
    c78->Write();
    c79->Write();
    c80->Write();
    c81->Write();
    c82->Write();
    c83->Write();
    c84->Write();
    c85->Write();
    c86->Write();
    c87->Write();
    c88->Write();
    c89->Write();
    c90->Write();
    c91->Write();
    c92->Write();
    c93->Write();
    c94->Write();
    c95->Write();
    c96->Write();
    c97->Write();
    c98->Write();
    c99->Write();
    c100->Write();
    c101->Write();
    c102->Write();
    c103->Write();
    c104->Write();
    c105->Write();
    c106->Write();
    c107->Write();
    c108->Write();
    c109->Write();
    c110->Write();
    c111->Write();
    c112->Write();
    c113->Write();
    c114->Write();
    c115->Write();
    c116->Write();
    c117->Write();
    c118->Write();
    c119->Write();
    c120->Write();
    c121->Write();
    c122->Write();
    c123->Write();
    c124->Write();
    c125->Write();
    c126->Write();
    c127->Write();
    c128->Write();
    c129->Write();
    c130->Write();
    c131->Write();
    c132->Write();
    c133->Write();
    c134->Write();
    c135->Write();
    c136->Write();
    c137->Write();
    c138->Write();
    c139->Write();
    c140->Write();
    c141->Write();
    c142->Write();
    c143->Write();
    c144->Write();
    c145->Write();
    c146->Write();
    c147->Write();
    c148->Write();
    c149->Write();
    c150->Write();
    c151->Write();
    c152->Write();
    c153->Write();
    c154->Write();
    c155->Write();
    c156->Write();
    c157->Write();
    c158->Write();
    c159->Write();
    c161->Write();
    c162->Write();
    c163->Write();
    c164->Write();
    c165->Write();
    c166->Write();
    c167->Write();
    c168->Write();
    c169->Write();
    c170->Write();
    c171->Write();
    c172->Write();
    c173->Write();
    c174->Write();
    c175->Write();
    c176->Write();
    c177->Write();
    c178->Write();
    c179->Write();
    c180->Write();
    c181->Write();
    c182->Write();
    c183->Write();
    c184->Write();
    c185->Write();
    c186->Write();
    c187->Write();
    c188->Write();
    c189->Write();
    c190->Write();
    c191->Write();
    c192->Write();
    c193->Write();
    c194->Write();
    c195->Write();
    c196->Write();
    c197->Write();
    c198->Write();
    c199->Write();
    c200->Write();
    c201->Write();
    c202->Write();
    c203->Write();
    c204->Write();
    c205->Write();
    c206->Write();
    c207->Write();
    c208->Write();
    c209->Write();
    c210->Write();
    c211->Write();
    c212->Write();
    c213->Write();
    c214->Write();
    c215->Write();
    c216->Write();
    c217->Write();
    c218->Write();
    c219->Write();
    c220->Write();
    c221->Write();
    c222->Write();
    c223->Write();
    c224->Write();
    c225->Write();
    c226->Write();
    c227->Write();
    c228->Write();
    c229->Write();
    c230->Write();
    c231->Write();
    c232->Write();
    c233->Write();
    c234->Write();
    c235->Write();
    c236->Write();
    c237->Write();
    c238->Write();
    c239->Write();
    c240->Write();
    c241->Write();
    c242->Write();
    c243->Write();
    c244->Write();
    c245->Write();
    c246->Write();
    c247->Write();
    c248->Write();
    c249->Write();
    c250->Write();
    c251->Write();
    c252->Write();
    c253->Write();
    c254->Write();
    c255->Write();
    c256->Write();
    c257->Write();
    c258->Write();
    c259->Write();
    c260->Write();
    c261->Write();
    c262->Write();
    c263->Write();
    c264->Write();
    c265->Write();
    c266->Write();
    c267->Write();
    c268->Write();
    c269->Write();
    c270->Write();
    c271->Write();
    c272->Write();
    c273->Write();
    c274->Write();
    c275->Write();
    c276->Write();
    c277->Write();
    c278->Write();
    c279->Write();
    c280->Write();
    c281->Write();
    c282->Write();
    c283->Write();
    c284->Write();
    c285->Write();
    c286->Write();
    c287->Write();
    c288->Write();
    c289->Write();
    c290->Write();
    c291->Write();
    c292->Write();
    c293->Write();
    c294->Write();
    c295->Write();
    c296->Write();
    c297->Write();
    c298->Write();

    mycE_1->Write();
    mycE_2->Write();
    mycE_leak_1->Write();
    mycE_leak_2->Write();
    mycE_totsf_1->Write();
    mycE_totsf_2->Write();
    mycE_totsf_leak_1->Write();
    mycE_totsf_leak_2->Write();
    mycE_shmax_1->Write();
    mycE_shmax_2->Write();
    mycE_shmax_leak_1->Write();
    mycE_shmax_leak_2->Write();
    mycE_sfnorm_1->Write();
    mycE_sfnorm_2->Write();
    mycE_sfnorm_leak_1->Write();
    mycE_sfnorm_leak_2->Write();
    mycE_totsfnorm_1->Write();
    mycE_totsfnorm_2->Write();
    mycE_totsfnorm_leak_1->Write();
    mycE_totsfnorm_leak_2->Write();
    mycE_leak_offset_1->Write();
    mycE_leak_offset_2->Write();
    mycE_4sections_1->Write();
    mycE_4sections_2->Write();
    mycE_4sections_leak_1->Write();
    mycE_4sections_leak_2->Write();
    mycE_4sections_ch1_1->Write();
    mycE_4sections_ch1_2->Write();
    mycE_4sections_ch1_leak_1->Write();
    mycE_4sections_ch1_leak_2->Write();
    mycE_4sections_ch2_1->Write();
    mycE_4sections_ch2_2->Write();
    mycE_4sections_ch2_leak_1->Write();
    mycE_4sections_ch2_leak_2->Write();
    mycE_4sections_ch3_1->Write();
    mycE_4sections_ch3_2->Write();
    mycE_4sections_ch3_leak_1->Write();
    mycE_4sections_ch3_leak_2->Write();
    mycE_4sections_ch4_1->Write();
    mycE_4sections_ch4_2->Write();
    mycE_4sections_ch4_leak_1->Write();
    mycE_4sections_ch4_leak_2->Write();
    mycE_ch1_1->Write();
    mycE_ch1_2->Write();
    mycE_ch1_leak_1->Write();
    mycE_ch1_leak_2->Write();
    mycE_ch2_1->Write();
    mycE_ch2_2->Write();
    mycE_ch2_leak_1->Write();
    mycE_ch2_leak_2->Write();
    mycE_ch3_1->Write();
    mycE_ch3_2->Write();
    mycE_ch3_leak_1->Write();
    mycE_ch3_leak_2->Write();
    mycE_ch4_1->Write();
    mycE_ch4_2->Write();
    mycE_ch4_leak_1->Write();
    mycE_ch4_leak_2->Write();
    mycE_ch1_overtrueE_1->Write();
    mycE_ch1_overtrueE_2->Write();
    mycE_ch1_overtrueE_leak_1->Write();
    mycE_ch1_overtrueE_leak_2->Write();
    mycE_ch2_overtrueE_1->Write();
    mycE_ch2_overtrueE_2->Write();
    mycE_ch2_overtrueE_leak_1->Write();
    mycE_ch2_overtrueE_leak_2->Write();
    mycE_ch3_overtrueE_1->Write();
    mycE_ch3_overtrueE_2->Write();
    mycE_ch3_overtrueE_leak_1->Write();
    mycE_ch3_overtrueE_leak_2->Write();
    mycE_ch4_overtrueE_1->Write();
    mycE_ch4_overtrueE_2->Write();
    mycE_ch4_overtrueE_leak_1->Write();
    mycE_ch4_overtrueE_leak_2->Write();
    mycE_raw_1->Write();
    mycE_raw_2->Write();
    mycE_raw_leak_1->Write();
    mycE_raw_leak_2->Write();
    mycE_corr_1->Write();
    mycE_corr_2->Write();
    mycE_val_1->Write();
    mycE_val_2->Write();
    mycE_overtrue_1->Write();
    mycE_overtrue_2->Write();
    mycE_overtrue_corr_1->Write();
    mycE_overtrue_corr_2->Write();
    mycE_overtrue_val_1->Write();
    mycE_overtrue_val_2->Write();
    
    //The name of the pdf we want to produce
    TString pdfname = "Plots_spread.pdf"; //"Plots_nospread.pdf"
    //Auxillary canvas to write some titles inside pdf
    TCanvas * tpdf = new TCanvas("tpdf", " ");
    tpdf->cd();
    TLatex lat; lat.SetTextSize(0.08);
    lat.DrawLatex(0.1,0.5,"Scenario: SF per layer per energy");
    lat.DrawLatex(0.1,0.4,"but cheating sections");
    lat.DrawLatex(0.1,0.3,"1. Linearity plots");
    tpdf->Update();
    tpdf->Print(pdfname+"(");
    //This is what I thought it was the linearity plots. Uncomment if needed. 
    // c143->Print(pdfname);
    // c144->Print(pdfname);
    // c145->Print(pdfname);
    // c146->Print(pdfname);
    // c147->Print(pdfname);
    // c148->Print(pdfname);
    // c149->Print(pdfname);
    // c150->Print(pdfname);
    c170->Print(pdfname);
    c171->Print(pdfname);
    c172->Print(pdfname);
    c173->Print(pdfname);
    c174->Print(pdfname);
    c175->Print(pdfname);
    c176->Print(pdfname);
    c177->Print(pdfname);
    tpdf->Clear();
    tpdf->cd();
    lat.DrawLatex(0.1,0.5,"Scenario: SF per layer per energy");
    lat.DrawLatex(0.1,0.4,"but cheating sections");
    lat.DrawLatex(0.1,0.3,"2. Resolution plots");
    tpdf->Update();
    tpdf->Print(pdfname);
    //This is a hack because we have re-used resoFit_ch1
    resoFit_ch1->SetTitle("Relative energy resolution but cheat section 1 for e- and all beam energies");
    c151->Update();
    c151->Print(pdfname);
    //This is a hack because we have re-used resoFit_ch1_leak
    resoFit_ch1_leak->SetTitle("Relative energy resolution but cheat section 1 including leakage for e- and all beam energies");
    c152->Update();
    c152->Print(pdfname);
    c153->Print(pdfname);
    c154->Print(pdfname);
    c155->Print(pdfname);
    c156->Print(pdfname);
    c157->Print(pdfname);
    c158->Print(pdfname);
    tpdf->Clear();
    tpdf->cd();
    lat.DrawLatex(0.1,0.5,"Scenario: SF per layer per energy");
    lat.DrawLatex(0.1,0.4,"but cheating sections");
    lat.DrawLatex(0.1,0.3,"3. Combine Resolution plots");
    tpdf->Update();
    tpdf->Print(pdfname);
    //This is a hack because we have re-used resoFit_ch1
    resoFit_ch1->SetTitle("Relative energy resolution for e-");
    resoFit_ch1->GetFunction("reso")->SetLineColor(1);
    resoFit_ch2->GetFunction("reso")->SetLineColor(2);
    resoFit_ch3->GetFunction("reso")->SetLineColor(3);
    resoFit_ch4->GetFunction("reso")->SetLineColor(4);
    c159->Update();
    c159->Print(pdfname);
    resoFit_ch1_leak->SetTitle("Relative energy resolution for e-");
    resoFit_ch1_leak->GetFunction("reso")->SetLineColor(1);
    resoFit_ch2_leak->GetFunction("reso")->SetLineColor(2);
    resoFit_ch3_leak->GetFunction("reso")->SetLineColor(3);
    resoFit_ch4_leak->GetFunction("reso")->SetLineColor(4);
     c161->Update();
    c161->Print(pdfname);
    tpdf->Clear();
    tpdf->cd();
    lat.DrawLatex(0.1,0.5,"Scenario: SF per layer per energy");
    lat.DrawLatex(0.1,0.4,"but cheating sections");
    lat.DrawLatex(0.1,0.3,"4. The fits used for previous plots");
    tpdf->Update();
    tpdf->Print(pdfname);
    mycE_ch1_1->Print(pdfname);
    mycE_ch1_2->Print(pdfname);
    mycE_ch1_leak_1->Print(pdfname);
    mycE_ch1_leak_2->Print(pdfname);
    mycE_ch2_1->Print(pdfname);
    mycE_ch2_2->Print(pdfname);
    mycE_ch2_leak_1->Print(pdfname);
    mycE_ch2_leak_2->Print(pdfname);
    mycE_ch3_1->Print(pdfname);
    mycE_ch3_2->Print(pdfname);
    mycE_ch3_leak_1->Print(pdfname);
    mycE_ch3_leak_2->Print(pdfname);
    mycE_ch4_1->Print(pdfname);
    mycE_ch4_2->Print(pdfname);
    mycE_ch4_leak_1->Print(pdfname);
    mycE_ch4_leak_2->Print(pdfname);
    tpdf->Clear();
    tpdf->cd();
    lat.DrawLatex(0.1,0.5,"Scenario: SF per layer per energy");
    lat.DrawLatex(0.1,0.4,"but cheating sections");
    lat.DrawLatex(0.1,0.3,"4. No fits for  #frac{E_{reco}-E_{true}}{E_{true}} plots");
    tpdf->Update();
    tpdf->Print(pdfname);
    mycE_ch1_overtrueE_1->Print(pdfname);
    mycE_ch1_overtrueE_2->Print(pdfname);
    mycE_ch1_overtrueE_leak_1->Print(pdfname);
    mycE_ch1_overtrueE_leak_2->Print(pdfname);
    mycE_ch2_overtrueE_1->Print(pdfname);
    mycE_ch2_overtrueE_2->Print(pdfname);
    mycE_ch2_overtrueE_leak_1->Print(pdfname);
    mycE_ch2_overtrueE_leak_2->Print(pdfname);
    mycE_ch3_overtrueE_1->Print(pdfname);
    mycE_ch3_overtrueE_2->Print(pdfname);
    mycE_ch3_overtrueE_leak_1->Print(pdfname);
    mycE_ch3_overtrueE_leak_2->Print(pdfname);
    mycE_ch4_overtrueE_1->Print(pdfname);
    mycE_ch4_overtrueE_2->Print(pdfname);
    mycE_ch4_overtrueE_leak_1->Print(pdfname);
    mycE_ch4_overtrueE_leak_2->Print(pdfname);
    tpdf->Clear();
    tpdf->cd();
    lat.DrawLatex(0.1,0.5,"Scenario: SF per layer per energy");
    lat.DrawLatex(0.1,0.4,"1. Linearity plots");
    lat.SetTextSize(0.08);
    tpdf->Update();
    tpdf->Print(pdfname);
    //This is what I thought it was the linearity plots. Uncomment if needed. 
    // c11->Print(pdfname);
    // c37->Print(pdfname);
    // c42->Print(pdfname);
    c7->Print(pdfname);
    c180->Print(pdfname);
    c181->Print(pdfname);
    tpdf->Clear();
    tpdf->cd();
    lat.DrawLatex(0.1,0.5,"Scenario: SF per layer per energy");
    lat.DrawLatex(0.1,0.4,"2. Resolution plots");
    tpdf->Update();
    tpdf->Print(pdfname);
    c15->Print(pdfname);
    c38->Print(pdfname);
    c43->Print(pdfname);
    tpdf->Clear();
    tpdf->cd();
    lat.DrawLatex(0.1,0.5,"Scenario: SF per layer per energy");
    lat.DrawLatex(0.1,0.4,"3. Combine resolution plots");
    tpdf->Update();
    tpdf->Print(pdfname);
    resoFit->SetMarkerColor(1);resoFit->SetLineColor(1);resoFit->SetMarkerStyle(20);
    resoFit_leak->SetMarkerColor(2);resoFit_leak->SetLineColor(2);resoFit_leak->SetMarkerStyle(21);
    resoFit_leak_offset->SetMarkerColor(3);resoFit_leak_offset->SetLineColor(3);resoFit_leak_offset->SetMarkerStyle(22);
    resoFit->GetFunction("reso")->SetLineColor(1);
    resoFit_leak->GetFunction("reso")->SetLineColor(2);
    resoFit_leak_offset->GetFunction("reso")->SetLineColor(3);
    c44->Update();
    c44->Print(pdfname);
    // tpdf->Clear();
    // tpdf->cd();
    // lat.DrawLatex(0.1,0.5,"Scenario: SF per layer per energy");
    // lat.DrawLatex(0.1,0.3,"4. #frac{E_{reco}-E_{true}}{E_{true}}");
    // tpdf->Update();
    // tpdf->Print(pdfname);
    // c7->Print(pdfname);
    tpdf->Clear();
    tpdf->cd();
    lat.DrawLatex(0.1,0.5,"Scenario: SF per layer per energy");
    lat.DrawLatex(0.1,0.4,"5. Fits for the previous pages");
    tpdf->Update();
    tpdf->Print(pdfname);
    mycE_1->Print(pdfname);
    mycE_2->Print(pdfname);
    mycE_leak_1->Print(pdfname);
    mycE_leak_2->Print(pdfname);
    mycE_leak_offset_1->Print(pdfname);
    mycE_leak_offset_2->Print(pdfname);
    mycE_overtrue_1->Print(pdfname);
    mycE_overtrue_2->Print(pdfname);
    tpdf->Clear();
    tpdf->cd();
    lat.SetTextSize(0.06);
    lat.DrawLatex(0.1,0.5,"Scenario: One total SF for each energy");
    lat.DrawLatex(0.1,0.4,"1. Linearity plots");
    tpdf->Update();
    tpdf->Print(pdfname);
    //This is what I thought it was the linearity plots. Uncomment if needed. 
    // c48->Print(pdfname);
    // c51->Print(pdfname);
    c182->Print(pdfname);
    c183->Print(pdfname);
    tpdf->Clear();
    tpdf->cd();
    lat.DrawLatex(0.1,0.5,"Scenario: One total SF for each energy");
    lat.DrawLatex(0.1,0.4,"2. Resolution plots");
    tpdf->Update();
    tpdf->Print(pdfname);
    c49->Print(pdfname);
    c50->Print(pdfname);
    tpdf->Clear();
    tpdf->cd();
    lat.DrawLatex(0.1,0.5,"Scenario: One total SF for each energy");
    lat.DrawLatex(0.1,0.4,"3. Combine resolution plots");
    tpdf->Update();
    tpdf->Print(pdfname);
    resoFit_totsf->SetMarkerColor(1);resoFit_totsf->SetLineColor(1);resoFit_totsf->SetMarkerStyle(20);
    resoFit_totsf_leak->SetMarkerColor(2);resoFit_totsf_leak->SetLineColor(2);resoFit_totsf_leak->SetMarkerStyle(21);
    resoFit->SetMarkerColor(3);resoFit->SetLineColor(3);resoFit->SetMarkerStyle(22);
    resoFit_totsf->GetFunction("reso")->SetLineColor(1);
    resoFit_totsf_leak->GetFunction("reso")->SetLineColor(2);
    resoFit->GetFunction("reso")->SetLineColor(3);
    c52->Update();
    c52->Print(pdfname);
    tpdf->Clear();
    tpdf->cd();
    lat.DrawLatex(0.1,0.5,"Scenario: One total SF for each energy");
    lat.DrawLatex(0.1,0.4,"4. Fits for the previous plots");
    tpdf->Update();
    tpdf->Print(pdfname);
    mycE_totsf_1->Print(pdfname);
    mycE_totsf_2->Print(pdfname);
    mycE_totsf_leak_1->Print(pdfname);
    mycE_totsf_leak_2->Print(pdfname);
    tpdf->Clear();
    tpdf->cd();
    lat.DrawLatex(0.1,0.5,"Scenario: One total SF for each energy");
    lat.DrawLatex(0.1,0.4,"5. No fits for  #frac{E_{reco}-E_{true}}{E_{true}} plots");
    tpdf->Update();
    tpdf->Print(pdfname);
    mycE_overtrue_totsf_1->Print(pdfname);
    mycE_overtrue_totsf_2->Print(pdfname);
    mycE_overtrue_totsf_leak_1->Print(pdfname);
    mycE_overtrue_totsf_leak_2->Print(pdfname);
    tpdf->Clear();
    tpdf->cd();
    lat.DrawLatex(0.1,0.5,"Scenario: dEdx");
    lat.DrawLatex(0.1,0.4,"1. Linearity plots");
    tpdf->Update();
    tpdf->Print(pdfname);
    //This is what I thought it was the linearity plots. Uncomment if needed. 
    // c24->Print(pdfname);
    c27->Print(pdfname);
    tpdf->Clear();
    tpdf->cd();
    lat.DrawLatex(0.1,0.5,"Scenario: dEdx");
    lat.DrawLatex(0.1,0.4,"2. Resolution plots");
    tpdf->Update();
    tpdf->Print(pdfname);
    c20->Print(pdfname);
    tpdf->Clear();
    tpdf->cd();
    lat.SetTextSize(0.08);
    lat.DrawLatex(0.1,0.5,"Scenario: dEdx");
    lat.DrawLatex(0.1,0.4,"3. Combine resolution plots");
    tpdf->Update();
    tpdf->Print(pdfname);
    resoFit_raw->SetMarkerColor(1);resoFit_raw->SetLineColor(1);resoFit_raw->SetMarkerStyle(20);
    resoFit->SetMarkerColor(2);resoFit->SetLineColor(2);resoFit->SetMarkerStyle(21);
    resoFit_val->SetMarkerColor(3);resoFit_val->SetLineColor(3);resoFit_val->SetMarkerStyle(22);
    resoFit_raw->GetFunction("reso")->SetLineColor(1);
    resoFit->GetFunction("reso")->SetLineColor(2);
    resoFit_val->GetFunction("reso")->SetLineColor(3);
    c25->Update();
    c25->Print(pdfname);
    // tpdf->Clear();
    // tpdf->cd();
    // lat.DrawLatex(0.1,0.5,"Scenario: dEdx");
    // lat.DrawLatex(0.1,0.3,"4. #frac{E_{reco}-E_{true}}{E_{true}}");
    // tpdf->Update();
    // tpdf->Print(pdfname);
    // c27->Print(pdfname);
    tpdf->Clear();
    tpdf->cd();
    lat.DrawLatex(0.1,0.5,"Scenario: dEdx");
    lat.DrawLatex(0.1,0.4,"5. Fits used for the previous pages");
    tpdf->Update();
    tpdf->Print(pdfname);
    mycE_val_1->Print(pdfname);
    mycE_val_2->Print(pdfname);
    tpdf->Clear();
    tpdf->cd();
    lat.DrawLatex(0.1,0.5,"Scenario: dEdx");
    lat.DrawLatex(0.1,0.4,"4. No fits for  #frac{E_{reco}-E_{true}}{E_{true}} plots");
    tpdf->Update();
    tpdf->Print(pdfname);
    mycE_overtrue_val_1->Print(pdfname);
    mycE_overtrue_val_2->Print(pdfname);
    // tpdf->Clear();
    // tpdf->cd();
    // lat.SetTextSize(0.06);
    // lat.DrawLatex(0.1,0.5,"Scenario: Total SF and shower max method");
    // lat.DrawLatex(0.1,0.4,"1. Sampling fraction dependency");
    // tpdf->Update();
    // tpdf->Print(pdfname);
    // c_totsfdistrib->Print(pdfname);
    // c_totsfvsene->Print(pdfname);
    // c_totsfvstprof->Print(pdfname);
    tpdf->Clear();
    tpdf->cd();
    lat.SetTextSize(0.06);
    lat.DrawLatex(0.1,0.5,"Scenario: Total SF and shower max method");
    lat.DrawLatex(0.1,0.4,"1. Linearity plots");
    tpdf->Update();
    tpdf->Print(pdfname);
    //This is what I thought it was the linearity plots. Uncomment if needed. 
    // c211->Print(pdfname);
    // c212->Print(pdfname);
    c209->Print(pdfname);
    c210->Print(pdfname);
    tpdf->Clear();
    tpdf->cd();
    lat.DrawLatex(0.1,0.5,"Scenario: Total SF and shower max method");
    lat.DrawLatex(0.1,0.4,"2. Resolution plots");
    tpdf->Update();
    tpdf->Print(pdfname);
    resoFit_totsfnorm->SetTitle("Relative energy resolution using total SF with X_{0,max}/<X_{0,max}> for e- and all beam energies");
    c213->Update();
    c213->Print(pdfname);
    c214->Print(pdfname);
    tpdf->Clear();
    tpdf->cd();
    lat.DrawLatex(0.1,0.5,"Scenario: Total SF and shower max method");
    lat.DrawLatex(0.1,0.4,"3. Combine resolution plots");
    tpdf->Update();
    tpdf->Print(pdfname);
    resoFit_totsfnorm->SetMarkerColor(1);resoFit_totsfnorm->SetLineColor(1);resoFit_totsfnorm->SetMarkerStyle(20);
    resoFit_totsfnorm_leak->SetMarkerColor(2);resoFit_totsfnorm_leak->SetLineColor(2);resoFit_totsfnorm_leak->SetMarkerStyle(21);
    resoFit->SetMarkerColor(3);resoFit->SetLineColor(3);resoFit->SetMarkerStyle(22);
    resoFit_totsfnorm->GetFunction("reso")->SetLineColor(1);
    resoFit_totsfnorm_leak->GetFunction("reso")->SetLineColor(2);
    resoFit->GetFunction("reso")->SetLineColor(3);
    c215->Update();
    c215->Print(pdfname);
    resoFit_totsfnorm->SetMarkerColor(1);resoFit_totsfnorm->SetLineColor(1);resoFit_totsfnorm->SetMarkerStyle(20);
    resoFit_totsfnorm_leak->SetMarkerColor(2);resoFit_totsfnorm_leak->SetLineColor(2);resoFit_totsfnorm_leak->SetMarkerStyle(21);
    resoFit_totsf->SetMarkerColor(3);resoFit_totsf->SetLineColor(3);resoFit_totsf->SetMarkerStyle(22);
    resoFit_totsf_leak->SetMarkerColor(4);resoFit_totsf_leak->SetLineColor(4);resoFit_totsf_leak->SetMarkerStyle(23);
    resoFit_totsfnorm->GetFunction("reso")->SetLineColor(1);
    resoFit_totsfnorm_leak->GetFunction("reso")->SetLineColor(2);
    resoFit_totsf->GetFunction("reso")->SetLineColor(3);
    resoFit_totsf_leak->GetFunction("reso")->SetLineColor(4);
    c221->Update();
    c221->Print(pdfname);
    tpdf->Clear();
    tpdf->cd();
    lat.DrawLatex(0.1,0.5,"Scenario: Total SF and shower max method");
    lat.DrawLatex(0.1,0.4,"4. Fits for the previous plots");
    tpdf->Update();
    tpdf->Print(pdfname);
    mycE_totsfnorm_1->Print(pdfname);
    mycE_totsfnorm_2->Print(pdfname);
    mycE_totsfnorm_leak_1->Print(pdfname);
    mycE_totsfnorm_leak_2->Print(pdfname);
    tpdf->Clear();
    tpdf->cd();
    lat.DrawLatex(0.1,0.5,"Scenario: Total SF and shower max method");
    lat.DrawLatex(0.1,0.4,"5. No fits for  #frac{E_{reco,totsfnorm}-E_{true}}{E_{true}} plots");
    tpdf->Update();
    tpdf->Print(pdfname);
    mycE_overtrue_totsfnorm_1->Print(pdfname);
    mycE_overtrue_totsfnorm_2->Print(pdfname);
    mycE_overtrue_totsfnorm_leak_1->Print(pdfname);
    mycE_overtrue_totsfnorm_leak_2->Print(pdfname);
    tpdf->Clear();
    tpdf->cd();
    lat.SetTextSize(0.04);
    lat.DrawLatex(0.1,0.5,"Scenario: SF per layer and per energy and shower max method");
    lat.DrawLatex(0.1,0.4,"1. Linearity plots");
    tpdf->Update();
    tpdf->Print(pdfname);
    // c_sfperlayervstprof->Print(pdfname);
    // c_sfperlayervstprof_1->Print(pdfname);
    // c_sfperlayervstprof_2->Print(pdfname);
    c218->Print(pdfname);
    c219->Print(pdfname);
    tpdf->Clear();
    tpdf->cd();
    lat.DrawLatex(0.1,0.5,"Scenario: SF per layer and per energy and shower max method");
    lat.DrawLatex(0.1,0.4,"2. Resolution plots");
    tpdf->Update();
    tpdf->Print(pdfname);
    resoFit_sfnorm->SetTitle("Relative energy resolution using the sampling fraction with X_{0,max}/<X_{0,max}> for e- and all beam energies"); 
    c60->Update();
    c60->Print(pdfname);
    c65->Print(pdfname);
    tpdf->Clear();
    tpdf->cd();
    lat.DrawLatex(0.1,0.5,"Scenario: SF per layer and per energy and shower max method");
    lat.DrawLatex(0.1,0.4,"3. Combine resolution plots");
    tpdf->Update();
    tpdf->Print(pdfname);
    resoFit_sfnorm->SetMarkerColor(1);resoFit_sfnorm->SetLineColor(1);resoFit_sfnorm->SetMarkerStyle(20);
    resoFit_sfnorm_leak->SetMarkerColor(2);resoFit_sfnorm_leak->SetLineColor(2);resoFit_sfnorm_leak->SetMarkerStyle(21);
    resoFit->SetMarkerColor(3);resoFit->SetLineColor(3);resoFit->SetMarkerStyle(22);
    resoFit_sfnorm->GetFunction("reso")->SetLineColor(1);
    resoFit_sfnorm_leak->GetFunction("reso")->SetLineColor(2);
    resoFit->GetFunction("reso")->SetLineColor(3);
    c220->Update();
    c220->Print(pdfname);
    tpdf->Clear();
    tpdf->cd();
    lat.DrawLatex(0.1,0.5,"Scenario: SF per layer and per energy and shower max method");
    lat.DrawLatex(0.1,0.4,"4. Fits for the previous plots");
    tpdf->Update();
    tpdf->Print(pdfname);
    mycE_sfnorm_1->Print(pdfname);
    mycE_sfnorm_2->Print(pdfname);
    mycE_sfnorm_leak_1->Print(pdfname);
    mycE_sfnorm_leak_2->Print(pdfname);
    tpdf->Clear();
    tpdf->cd();
    lat.DrawLatex(0.1,0.5,"Scenario: SF per layer and per energy and shower max method");
    lat.DrawLatex(0.1,0.4,"5. No fits for  #frac{E_{reco,sfnorm}-E_{true}}{E_{true}} plots");
    tpdf->Update();
    tpdf->Print(pdfname);
    mycE_overtrue_sfnorm_1->Print(pdfname);
    mycE_overtrue_sfnorm_2->Print(pdfname);
    mycE_overtrue_sfnorm_leak_1->Print(pdfname);
    mycE_overtrue_sfnorm_leak_2->Print(pdfname);
    tpdf->Clear();
    tpdf->cd();
    lat.SetTextSize(0.04);
    lat.DrawLatex(0.1,0.5,"Scenario: SF_{i} and shower max method");
    lat.DrawLatex(0.1,0.4,"1. Linearity plots");
    tpdf->Update();
    tpdf->Print(pdfname);
    // c_sfperlayervstprof->Print(pdfname);
    // c_sfperlayervstprof_1->Print(pdfname);
    // c_sfperlayervstprof_2->Print(pdfname);
    c194->Print(pdfname);
    c195->Print(pdfname);
    tpdf->Clear();
    tpdf->cd();
    lat.DrawLatex(0.1,0.5,"Scenario: SF_{i} and shower max method");
    lat.DrawLatex(0.1,0.4,"2. Resolution plots");
    tpdf->Update();
    tpdf->Print(pdfname);
    resoFit_shmax->SetTitle("Relative energy resolution using SF_{i} and the shower max method for e- and all beam energies"); 
    c198->Update();
    c198->Print(pdfname);
    c199->Print(pdfname);
    tpdf->Clear();
    tpdf->cd();
    lat.DrawLatex(0.1,0.5,"Scenario: SF_{i} and shower max method");
    lat.DrawLatex(0.1,0.4,"3. Combine resolution plots");
    tpdf->Update();
    tpdf->Print(pdfname);
    resoFit_shmax->SetMarkerColor(1);resoFit_shmax->SetLineColor(1);resoFit_shmax->SetMarkerStyle(20);
    resoFit_shmax_leak->SetMarkerColor(2);resoFit_shmax_leak->SetLineColor(2);resoFit_shmax_leak->SetMarkerStyle(21);
    resoFit->SetMarkerColor(3);resoFit->SetLineColor(3);resoFit->SetMarkerStyle(22);
    resoFit_shmax->GetFunction("reso")->SetLineColor(1);
    resoFit_shmax_leak->GetFunction("reso")->SetLineColor(2);
    resoFit->GetFunction("reso")->SetLineColor(3);
    c200->Update();
    c200->Print(pdfname);
    tpdf->Clear();
    tpdf->cd();
    lat.DrawLatex(0.1,0.5,"Scenario: SF_{i} and shower max method");
    lat.DrawLatex(0.1,0.4,"4. Fits for the previous plots");
    tpdf->Update();
    tpdf->Print(pdfname);
    mycE_shmax_1->Print(pdfname);
    mycE_shmax_2->Print(pdfname);
    mycE_shmax_leak_1->Print(pdfname);
    mycE_shmax_leak_2->Print(pdfname);
    tpdf->Clear();
    tpdf->cd();
    lat.DrawLatex(0.1,0.5,"Scenario: SF_{i} and shower max method");
    lat.DrawLatex(0.1,0.4,"4. No fits for  #frac{E_{reco,shmax}-E_{true}}{E_{true}} plots");
    tpdf->Update();
    tpdf->Print(pdfname);
    mycE_overtrue_shmax_1->Print(pdfname);
    mycE_overtrue_shmax_2->Print(pdfname);
    mycE_overtrue_shmax_leak_1->Print(pdfname);
    mycE_overtrue_shmax_leak_2->Print(pdfname);
    tpdf->Clear();
    tpdf->cd();
    lat.SetTextSize(0.06);
    lat.DrawLatex(0.1,0.5,"Scenario: Total SF per section for each energy");
    lat.DrawLatex(0.1,0.4,"1. Linearity plots");
    tpdf->Update();
    tpdf->Print(pdfname);
    //This is what I thought it was the linearity plots. Uncomment if needed. 
    // c232->Print(pdfname);
    // c233->Print(pdfname);
    c230->Print(pdfname);
    c231->Print(pdfname);
    tpdf->Clear();
    tpdf->cd();
    lat.DrawLatex(0.1,0.5,"Scenario: Total SF per section for each energy");
    lat.DrawLatex(0.1,0.4,"2. Resolution plots");
    tpdf->Update();
    tpdf->Print(pdfname);
    resoFit_totsf_sec->SetTitle("Relative energy resolution using the total sampling fraction per section for e- and all beam energies");
    c234->Update();
    c234->Print(pdfname);
    resoFit_totsf_sec_leak->SetTitle("Relative energy resolution using the total sampling fraction per section adding leakage for e- and all beam energies");
    c235->Update();
    c235->Print(pdfname);
    tpdf->Clear();
    tpdf->cd();
    lat.DrawLatex(0.1,0.5,"Scenario: Total SF per section for each energy");
    lat.DrawLatex(0.1,0.4,"3. Combine resolution plots");
    tpdf->Update();
    tpdf->Print(pdfname);
    resoFit_totsf_sec->SetMarkerColor(1);resoFit_totsf_sec->SetLineColor(1);resoFit_totsf_sec->SetMarkerStyle(20);
    resoFit_totsf_sec_leak->SetMarkerColor(2);resoFit_totsf_sec_leak->SetLineColor(2);resoFit_totsf_sec_leak->SetMarkerStyle(21);
    resoFit->SetMarkerColor(3);resoFit->SetLineColor(3);resoFit->SetMarkerStyle(22);
    resoFit_totsf_sec->GetFunction("reso")->SetLineColor(1);
    resoFit_totsf_sec_leak->GetFunction("reso")->SetLineColor(2);
    resoFit->GetFunction("reso")->SetLineColor(3);
    c236->Update();
    c236->Print(pdfname);
    tpdf->Clear();
    tpdf->cd();
    lat.DrawLatex(0.1,0.5,"Scenario: Total SF per section for each energy");
    lat.DrawLatex(0.1,0.4,"4. Fits for the previous plots");
    tpdf->Update();
    tpdf->Print(pdfname);
    mycE_totsf_sec_1->Print(pdfname);
    mycE_totsf_sec_2->Print(pdfname);
    mycE_totsf_sec_leak_1->Print(pdfname);
    mycE_totsf_sec_leak_2->Print(pdfname);
    tpdf->Clear();
    tpdf->cd();
    lat.DrawLatex(0.1,0.5,"Scenario: Total SF per section for each energy");
    lat.DrawLatex(0.1,0.4,"5. No fits for  #frac{E_{reco}-E_{true}}{E_{true}} plots");
    tpdf->Update();
    tpdf->Print(pdfname);
    mycE_overtrue_totsf_sec_1->Print(pdfname);
    mycE_overtrue_totsf_sec_2->Print(pdfname);
    mycE_overtrue_totsf_sec_leak_1->Print(pdfname);
    mycE_overtrue_totsf_sec_leak_2->Print(pdfname+")");

    //============================================================================================
    //The name of the pdf we want to produce
    pdfname = usenewshowerdepth ? "Plots_spread_1SF_3SF_dEdx_newshvar.pdf" : "Plots_spread_1SF_3SF_dEdx_oldshvar.pdf"; //"Plots_nospread.pdf"
    tpdf->Clear();
    tpdf->cd();
    lat.SetTextSize(0.04);
    lat.DrawLatex(0.1,0.5,"Scenario: One total SF for each energy");
    lat.DrawLatex(0.1,0.4,"1. Linearity plots");
    tpdf->Update();
    tpdf->Print(pdfname+"(");
    //This is what I thought it was the linearity plots. Uncomment if needed. 
    // c48->Print(pdfname);
    // c51->Print(pdfname);
    // c182->Print(pdfname); //without leakage
    c183->Print(pdfname);
    tpdf->Clear();
    tpdf->cd();
    lat.DrawLatex(0.1,0.5,"Scenario: One total SF for each energy");
    lat.DrawLatex(0.1,0.4,"2. Resolution plots");
    tpdf->Update();
    tpdf->Print(pdfname);
    // c49->Print(pdfname); //without leakage
    c50->Print(pdfname);
    tpdf->Clear();
    tpdf->cd();
    lat.DrawLatex(0.1,0.5,"Scenario: One total SF for each energy");
    lat.DrawLatex(0.1,0.4,"3. Fits for the previous plots");
    tpdf->Update();
    tpdf->Print(pdfname);
    // mycE_totsf_1->Print(pdfname);
    // mycE_totsf_2->Print(pdfname);
    mycE_totsf_leak_1->Print(pdfname);
    mycE_totsf_leak_2->Print(pdfname);
    tpdf->Clear();
    tpdf->cd();
    lat.DrawLatex(0.1,0.5,"Scenario: One total SF for each energy");
    lat.DrawLatex(0.1,0.4,"4. No fits for  #frac{E_{reco,totsf}-E_{true}}{E_{true}} plots");
    tpdf->Update();
    tpdf->Print(pdfname);
    // mycE_overtrue_totsf_1->Print(pdfname);
    // mycE_overtrue_totsf_2->Print(pdfname);
    mycE_overtrue_totsf_leak_1->Print(pdfname);
    mycE_overtrue_totsf_leak_2->Print(pdfname);
    
    //Total sf using dedx weight 
    tpdf->Clear();
    tpdf->cd();
    lat.SetTextSize(0.04);
    lat.DrawLatex(0.1,0.5,"Scenario: One total SF for each energy using dedx weight");
    lat.DrawLatex(0.1,0.4,"1. Linearity plots");
    tpdf->Update();
    tpdf->Print(pdfname+"(");
    //This is what I thought it was the linearity plots. Uncomment if needed. 
    // c->Print(pdfname);
    // c->Print(pdfname);
    // c322->Print(pdfname); //without leakage
    c323->Print(pdfname);
    tpdf->Clear();
    tpdf->cd();
    lat.DrawLatex(0.1,0.5,"Scenario: One total SF for each energy using dedx weight");
    lat.DrawLatex(0.1,0.4,"2. Resolution plots");
    tpdf->Update();
    tpdf->Print(pdfname);
    // c326->Print(pdfname); //without leakage
    c327->Print(pdfname);
    tpdf->Clear();
    tpdf->cd();
    lat.DrawLatex(0.1,0.5,"Scenario: One total SF for each energy using dedx weight");
    lat.DrawLatex(0.1,0.4,"3. Fits for the previous plots");
    tpdf->Update();
    tpdf->Print(pdfname);
    // mycE_totsf_dedxweighted_1->Print(pdfname);
    // mycE_totsf_dedxweighted_2->Print(pdfname);
    mycE_totsf_dedxweighted_leak_1->Print(pdfname);
    mycE_totsf_dedxweighted_leak_2->Print(pdfname);
    tpdf->Clear();
    tpdf->cd();
    lat.DrawLatex(0.1,0.5,"Scenario: One total SF for each energy using dedx weight");
    lat.DrawLatex(0.1,0.4,"4. No fits for  #frac{E_{reco,totsf_dedxweighted}-E_{true}}{E_{true}} plots");
    tpdf->Update();
    tpdf->Print(pdfname);
    // mycE_overtrue_totsf_dedxweighted_1->Print(pdfname);
    // mycE_overtrue_totsf_dedxweighted_2->Print(pdfname);
    mycE_overtrue_totsf_dedxweighted_leak_1->Print(pdfname);
    mycE_overtrue_totsf_dedxweighted_leak_2->Print(pdfname);










    tpdf->Clear();
    tpdf->cd();
    lat.SetTextSize(0.06);
    lat.DrawLatex(0.1,0.5,"Scenario: Total SF per section for each energy");
    lat.DrawLatex(0.1,0.4,"1. Linearity plots");
    tpdf->Update();
    tpdf->Print(pdfname);
    //This is what I thought it was the linearity plots. Uncomment if needed. 
    // c232->Print(pdfname);
    // c233->Print(pdfname);
    // c230->Print(pdfname);
    c231->Print(pdfname);
    tpdf->Clear();
    tpdf->cd();
    lat.DrawLatex(0.1,0.5,"Scenario: Total SF per section for each energy");
    lat.DrawLatex(0.1,0.4,"2. Resolution plots");
    tpdf->Update();
    tpdf->Print(pdfname);
    // resoFit_totsf_sec->SetTitle("Relative energy resolution using the total sampling fraction per section for e- and all beam energies");
    // c234->Update();
    // c234->Print(pdfname);
    resoFit_totsf_sec_leak->SetTitle("Relative energy resolution using the total sampling fraction per section adding leakage for e- and all beam energies");
    c235->Update();
    c235->Print(pdfname);
    tpdf->Clear();
    tpdf->cd();
    lat.DrawLatex(0.1,0.5,"Scenario: Total SF per section for each energy");
    lat.DrawLatex(0.1,0.4,"3. Fits for the previous plots");
    tpdf->Update();
    tpdf->Print(pdfname);
    // mycE_totsf_sec_1->Print(pdfname);
    // mycE_totsf_sec_2->Print(pdfname);
    mycE_totsf_sec_leak_1->Print(pdfname);
    mycE_totsf_sec_leak_2->Print(pdfname);
    tpdf->Clear();
    tpdf->cd();
    lat.DrawLatex(0.1,0.5,"Scenario: Total SF per section for each energy");
    lat.DrawLatex(0.1,0.4,"4. No fits for  #frac{E_{reco,totsfsec}-E_{true}}{E_{true}} plots");
    tpdf->Update();
    tpdf->Print(pdfname);
    // mycE_overtrue_totsf_sec_1->Print(pdfname);
    // mycE_overtrue_totsf_sec_2->Print(pdfname);
    mycE_overtrue_totsf_sec_leak_1->Print(pdfname);
    mycE_overtrue_totsf_sec_leak_2->Print(pdfname);


    tpdf->Clear();
    tpdf->cd();
    lat.SetTextSize(0.04);
    lat.DrawLatex(0.1,0.5,"Scenario: Total SF per section weigthed with dedx for each energy");
    lat.DrawLatex(0.1,0.4,"1. Linearity plots");
    tpdf->Update();
    tpdf->Print(pdfname);
    //This is what I thought it was the linearity plots. Uncomment if needed. 
    // c->Print(pdfname);
    // c->Print(pdfname);
    // c292->Print(pdfname);
    c293->Print(pdfname);
    tpdf->Clear();
    tpdf->cd();
    lat.DrawLatex(0.1,0.5,"Scenario: Total SF per section weigthed with dedx for each energy");
    lat.DrawLatex(0.1,0.4,"2. Resolution plots");
    tpdf->Update();
    tpdf->Print(pdfname);
    // resoFit_totsf_sec_dedxweighted->SetTitle("Relative energy resolution using the total sampling fraction per section for e- and all beam energies");
    // c297->Update();
    // c297->Print(pdfname);
    resoFit_totsf_sec_dedxweighted_leak->SetTitle("Relative energy resolution using the total sampling fraction per section adding leakage for e- and all beam energies");
    c297->Update();
    c297->Print(pdfname);
    tpdf->Clear();
    tpdf->cd();
    lat.DrawLatex(0.1,0.5,"Scenario: Total SF per section weigthed with dedx for each energy");
    lat.DrawLatex(0.1,0.4,"3. Fits for the previous plots");
    tpdf->Update();
    tpdf->Print(pdfname);
    // mycE_totsf_sec_dedxweighted_1->Print(pdfname);
    // mycE_totsf_sec_dedxweighted_2->Print(pdfname);
    mycE_totsf_sec_dedxweighted_leak_1->Print(pdfname);
    mycE_totsf_sec_dedxweighted_leak_2->Print(pdfname);
    tpdf->Clear();
    tpdf->cd();
    lat.DrawLatex(0.1,0.5,"Scenario: Total SF per section weigthed with dedx for each energy");
    lat.DrawLatex(0.1,0.4,"4. No fits for  #frac{E_{reco,totsfsecdedx}-E_{true}}{E_{true}} plots");
    tpdf->Update();
    tpdf->Print(pdfname);
    // mycE_overtrue_totsf_sec_dedxweighted_1->Print(pdfname);
    // mycE_overtrue_totsf_sec_dedxweighted_2->Print(pdfname);
    mycE_overtrue_totsf_sec_dedxweighted_leak_1->Print(pdfname);
    mycE_overtrue_totsf_sec_dedxweighted_leak_2->Print(pdfname);








    tpdf->Clear();
    tpdf->cd();
    lat.SetTextSize(0.06);
    lat.DrawLatex(0.1,0.5,"Scenario: Total SF per 5 sections for each energy");
    lat.DrawLatex(0.1,0.4,"1. Linearity plots");
    tpdf->Update();
    tpdf->Print(pdfname);
    //This is what I thought it was the linearity plots. Uncomment if needed. 
    // c273->Print(pdfname);
    // c274->Print(pdfname);
    // c->Print(pdfname);
    c270->Print(pdfname);
    tpdf->Clear();
    tpdf->cd();
    lat.DrawLatex(0.1,0.5,"Scenario: Total SF per 5 sections for each energy");
    lat.DrawLatex(0.1,0.4,"2. Resolution plots");
    tpdf->Update();
    tpdf->Print(pdfname);
    // resoFit_totsf_dedxsec->SetTitle("Relative energy resolution using the total sampling fraction per section for e- and all beam energies");
    // c->Update();
    // c->Print(pdfname);
    resoFit_totsf_dedxsec_leak->SetTitle("Relative energy resolution using the total sampling fraction per section adding leakage for e- and all beam energies");
    c278->Update();
    c278->Print(pdfname);
    tpdf->Clear();
    tpdf->cd();
    lat.DrawLatex(0.1,0.5,"Scenario: Total SF per 5 sections for each energy");
    lat.DrawLatex(0.1,0.4,"3. Fits for the previous plots");
    tpdf->Update();
    tpdf->Print(pdfname);
    // mycE_totsf_dedxsec_1->Print(pdfname);
    // mycE_totsf_dedxsec_2->Print(pdfname);
    mycE_totsf_dedxsec_leak_1->Print(pdfname);
    mycE_totsf_dedxsec_leak_2->Print(pdfname);
    tpdf->Clear();
    tpdf->cd();
    lat.DrawLatex(0.1,0.5,"Scenario: Total SF per 5 sections for each energy");
    lat.DrawLatex(0.1,0.4,"4. No fits for  #frac{E_{reco,totsfdedxsec}-E_{true}}{E_{true}} plots");
    tpdf->Update();
    tpdf->Print(pdfname);
    // mycE_overtrue_totsf_dedxsec_1->Print(pdfname);
    // mycE_overtrue_totsf_dedxsec_2->Print(pdfname);
    mycE_overtrue_totsf_dedxsec_leak_1->Print(pdfname);
    mycE_overtrue_totsf_dedxsec_leak_2->Print(pdfname);




    tpdf->Clear();
    tpdf->cd();
    lat.DrawLatex(0.1,0.5,"Scenario: dEdx");
    lat.DrawLatex(0.1,0.4,"1. Linearity plots");
    tpdf->Update();
    tpdf->Print(pdfname);
    //This is what I thought it was the linearity plots. Uncomment if needed. 
    // c24->Print(pdfname);
    c27->Print(pdfname);
    tpdf->Clear();
    tpdf->cd();
    lat.DrawLatex(0.1,0.5,"Scenario: dEdx");
    lat.DrawLatex(0.1,0.4,"2. Resolution plots");
    tpdf->Update();
    tpdf->Print(pdfname);
    c20->Print(pdfname);
    tpdf->Clear();
    tpdf->cd();
    lat.DrawLatex(0.1,0.5,"Scenario: dEdx");
    lat.DrawLatex(0.1,0.4,"5. Fits used for the previous pages");
    tpdf->Update();
    tpdf->Print(pdfname);
    mycE_val_1->Print(pdfname);
    mycE_val_2->Print(pdfname);
    tpdf->Clear();
    tpdf->cd();
    lat.DrawLatex(0.1,0.5,"Scenario: dEdx");
    lat.DrawLatex(0.1,0.4,"4. No fits for  #frac{E_{reco,dEdx}-E_{true}}{E_{true}} plots");
    tpdf->Update();
    tpdf->Print(pdfname);
    mycE_overtrue_val_1->Print(pdfname);
    mycE_overtrue_val_2->Print(pdfname);
    tpdf->Clear();
    tpdf->cd();
    lat.DrawLatex(0.1,0.5,"Combine resolution plots of: ");
    lat.DrawLatex(0.1,0.4,"1. Total SF ");
    lat.DrawLatex(0.1,0.3,"2. Total SF per section ");
    lat.DrawLatex(0.1,0.2,"3. dEdx ");
    tpdf->Update();
    tpdf->Print(pdfname);
    resoFit_totsf_leak->SetMarkerColor(1);resoFit_totsf_leak->SetLineColor(1);resoFit_totsf_leak->SetMarkerStyle(20);
    resoFit_totsf_sec_leak->SetMarkerColor(2);resoFit_totsf_sec_leak->SetLineColor(2);resoFit_totsf_sec_leak->SetMarkerStyle(21);
    resoFit_val->SetMarkerColor(3);resoFit_val->SetLineColor(3);resoFit_val->SetMarkerStyle(22);
    resoFit_totsf_leak->GetFunction("reso")->SetLineColor(1);
    resoFit_totsf_sec_leak->GetFunction("reso")->SetLineColor(2);
    resoFit_val->GetFunction("reso")->SetLineColor(3);
    resoFit_totsf_leak->SetTitle("Combining relative energy resolutions for e- and different calibration schemes");
    c237->Update();
    c237->Print(pdfname);

    tpdf->Clear();
    tpdf->cd();
    lat.DrawLatex(0.1,0.5,"Combine resolution plots of: ");
    lat.DrawLatex(0.1,0.4,"1. Total SF ");
    lat.DrawLatex(0.1,0.3,"2. Total SF per section weighted with dEdx");
    lat.DrawLatex(0.1,0.2,"3. dEdx ");
    tpdf->Update();
    tpdf->Print(pdfname);
    resoFit_totsf_leak->SetMarkerColor(1);resoFit_totsf_leak->SetLineColor(1);resoFit_totsf_leak->SetMarkerStyle(20);
    resoFit_totsf_sec_dedxweighted_leak->SetMarkerColor(2);resoFit_totsf_sec_dedxweighted_leak->SetLineColor(2);resoFit_totsf_sec_dedxweighted_leak->SetMarkerStyle(21);
    resoFit_val->SetMarkerColor(3);resoFit_val->SetLineColor(3);resoFit_val->SetMarkerStyle(22);
    resoFit_totsf_leak->GetFunction("reso")->SetLineColor(1);
    resoFit_totsf_sec_dedxweighted_leak->GetFunction("reso")->SetLineColor(2);
    resoFit_val->GetFunction("reso")->SetLineColor(3);
    resoFit_totsf_leak->SetTitle("Combining relative energy resolutions for e- and different calibration schemes");
    c298->Update();
    c298->Print(pdfname);

    tpdf->Clear();
    tpdf->cd();
    lat.DrawLatex(0.1,0.5,"Combine resolution plots of: ");
    lat.DrawLatex(0.1,0.4,"1. Total SF using dEdx weight");
    lat.DrawLatex(0.1,0.3,"2. Total SF per section weighted with dEdx");
    lat.DrawLatex(0.1,0.2,"3. dEdx ");
    tpdf->Update();
    tpdf->Print(pdfname);
    resoFit_totsf_dedxweighted_leak->SetMarkerColor(1);resoFit_totsf_dedxweighted_leak->SetLineColor(1);resoFit_totsf_dedxweighted_leak->SetMarkerStyle(20);
    resoFit_totsf_sec_dedxweighted_leak->SetMarkerColor(2);resoFit_totsf_sec_dedxweighted_leak->SetLineColor(2);resoFit_totsf_sec_dedxweighted_leak->SetMarkerStyle(21);
    resoFit_val->SetMarkerColor(3);resoFit_val->SetLineColor(3);resoFit_val->SetMarkerStyle(22);
    resoFit_totsf_dedxweighted_leak->GetFunction("reso")->SetLineColor(1);
    resoFit_totsf_sec_dedxweighted_leak->GetFunction("reso")->SetLineColor(2);
    resoFit_val->GetFunction("reso")->SetLineColor(3);
    resoFit_totsf_dedxweighted_leak->SetTitle("Combining relative energy resolutions for e- and different calibration schemes");
    c328->Update();
    c328->Print(pdfname);




    tpdf->Clear();
    tpdf->cd();
    lat.DrawLatex(0.1,0.5,"Combine resolution plots of: ");
    lat.DrawLatex(0.1,0.4,"1. Total SF ");
    lat.DrawLatex(0.1,0.3,"2. Total SF per 5 sections ");
    lat.DrawLatex(0.1,0.2,"3. dEdx ");
    tpdf->Update();
    tpdf->Print(pdfname);
    resoFit_totsf_leak->SetMarkerColor(1);resoFit_totsf_leak->SetLineColor(1);resoFit_totsf_leak->SetMarkerStyle(20);
    resoFit_totsf_dedxsec_leak->SetMarkerColor(2);resoFit_totsf_dedxsec_leak->SetLineColor(2);resoFit_totsf_dedxsec_leak->SetMarkerStyle(21);
    resoFit_val->SetMarkerColor(3);resoFit_val->SetLineColor(3);resoFit_val->SetMarkerStyle(22);
    resoFit_totsf_leak->GetFunction("reso")->SetLineColor(1);
    resoFit_totsf_dedxsec_leak->GetFunction("reso")->SetLineColor(2);
    resoFit_val->GetFunction("reso")->SetLineColor(3);
    resoFit_totsf_leak->SetTitle("Combining relative energy resolutions for e- and different calibration schemes");
    c282->Update();
    c282->Print(pdfname);

    tpdf->Clear();
    tpdf->cd();
    lat.SetTextSize(0.06);
    lat.DrawLatex(0.1,0.5,"Scenario: Total SF and shower max method");
    lat.DrawLatex(0.1,0.4,"1. Linearity plots");
    tpdf->Update();
    tpdf->Print(pdfname);
    //This is what I thought it was the linearity plots. Uncomment if needed. 
    // c211->Print(pdfname);
    // c212->Print(pdfname);
    // c209->Print(pdfname);
    c210->Print(pdfname);
    tpdf->Clear();
    tpdf->cd();
    lat.DrawLatex(0.1,0.5,"Scenario: Total SF and shower max method");
    lat.DrawLatex(0.1,0.4,"2. Resolution plots");
    tpdf->Update();
    tpdf->Print(pdfname);
    // resoFit_totsfnorm->SetTitle("Relative energy resolution using total SF with X_{0,max}/<X_{0,max}> for e- and all beam energies");
    // c213->Update();
    // c213->Print(pdfname);
    c214->Print(pdfname);
    tpdf->Clear();
    tpdf->cd();
    tpdf->Clear();
    tpdf->cd();
    lat.DrawLatex(0.1,0.5,"Scenario: Total SF and shower max method");
    lat.DrawLatex(0.1,0.4,"3. Fits for the previous plots");
    tpdf->Update();
    tpdf->Print(pdfname);
    mycE_totsfnorm_leak_1->Print(pdfname);
    mycE_totsfnorm_leak_2->Print(pdfname);
    tpdf->Clear();
    tpdf->cd();
    lat.DrawLatex(0.1,0.5,"Scenario: Total SF and shower max method");
    lat.DrawLatex(0.1,0.4,"4. No fits for  #frac{E_{reco,totsfshmax}-E_{true}}{E_{true}} plots");
    tpdf->Update();
    tpdf->Print(pdfname);
    mycE_overtrue_totsfnorm_leak_1->Print(pdfname);
    mycE_overtrue_totsfnorm_leak_2->Print(pdfname);
    
 

    tpdf->Clear();
    tpdf->cd();
    lat.SetTextSize(0.04);
    lat.DrawLatex(0.1,0.5,"Scenario: Total SF and shower max method using dedx weight");
    lat.DrawLatex(0.1,0.4,"1. Linearity plots");
    tpdf->Update();
    tpdf->Print(pdfname);
    //This is what I thought it was the linearity plots. Uncomment if needed. 
    // c->Print(pdfname);
    // c->Print(pdfname);
    // c337->Print(pdfname);
    c338->Print(pdfname);
    tpdf->Clear();
    tpdf->cd();
    lat.DrawLatex(0.1,0.5,"Scenario: Total SF and shower max method using dedx weight");
    lat.DrawLatex(0.1,0.4,"2. Resolution plots");
    tpdf->Update();
    tpdf->Print(pdfname);
    // resoFit_totsfnorm->SetTitle("Relative energy resolution using total SF with X_{0,max}/<X_{0,max}> for e- and all beam energies");
    // c->Update();
    // c->Print(pdfname);
    c342->Print(pdfname);
    tpdf->Clear();
    tpdf->cd();
    lat.DrawLatex(0.1,0.5,"Scenario: Total SF and shower max method using dedx weight");
    lat.DrawLatex(0.1,0.4,"3. Fits for the previous plots");
    tpdf->Update();
    tpdf->Print(pdfname);
    mycE_totsfnorm_dedxweighted_leak_1->Print(pdfname);
    mycE_totsfnorm_dedxweighted_leak_2->Print(pdfname);
    tpdf->Clear();
    tpdf->cd();
    lat.DrawLatex(0.1,0.5,"Scenario: Total SF and shower max method using dedx weight");
    lat.DrawLatex(0.1,0.4,"4. No fits for  #frac{E_{reco,totsfshmaxdedx}-E_{true}}{E_{true}} plots");
    tpdf->Update();
    tpdf->Print(pdfname);
    mycE_overtrue_totsfnorm_dedxweighted_1->Print(pdfname);
    mycE_overtrue_totsfnorm_dedxweighted_2->Print(pdfname);








    tpdf->Clear();
    tpdf->cd();
    lat.SetTextSize(0.04);
    lat.DrawLatex(0.1,0.5,"Scenario: Total SF per section and shower max method for each energy");
    lat.DrawLatex(0.1,0.4,"1. Linearity plots");
    tpdf->Update();
    tpdf->Print(pdfname);
    //This is what I thought it was the linearity plots. Uncomment if needed. 
    // c242->Print(pdfname);
    // c243->Print(pdfname);
    // c246->Print(pdfname);
    c247->Print(pdfname);
    tpdf->Clear();
    tpdf->cd();
    lat.SetTextSize(0.04);
    lat.DrawLatex(0.1,0.5,"Scenario: Total SF per section and shower max method for each energy");
    lat.DrawLatex(0.1,0.4,"2. Resolution plots");
    tpdf->Update();
    tpdf->Print(pdfname);
    // resoFit_totsfshmax_sec->SetTitle("Relative energy resolution using the total sampling fraction per section for e- and all beam energies");
    // c->Update();
    // c->Print(pdfname);
    // resoFit_totsf_sec_leak->SetTitle("Relative energy resolution using the total sampling fraction per section adding leakage for e- and all beam energies");
    // c250->Update();
    c251->Print(pdfname);
    tpdf->Clear();
    tpdf->cd();
    lat.DrawLatex(0.1,0.5,"Scenario: Total SF per section and shower max method for each energy");
    lat.DrawLatex(0.1,0.4,"3. Fits for the previous plots");
    tpdf->Update();
    tpdf->Print(pdfname);
    // mycE_totsfshmax_sec_1->Print(pdfname);
    // mycE_totsfshmax_sec_2->Print(pdfname);
    mycE_totsfshmax_sec_leak_1->Print(pdfname);
    mycE_totsfshmax_sec_leak_2->Print(pdfname);
    tpdf->Clear();
    tpdf->cd();
    lat.DrawLatex(0.1,0.5,"Scenario: Total SF per section and shower max method for each energy");
    lat.DrawLatex(0.1,0.4,"4. No fits for  #frac{E_{reco,totsfshmaxsec}-E_{true}}{E_{true}} plots");
    tpdf->Update();
    tpdf->Print(pdfname);
    // mycE_overtrue_totsfshmax_sec_1->Print(pdfname);
    // mycE_overtrue_totsfshmax_sec_2->Print(pdfname);
    mycE_overtrue_totsfshmax_sec_leak_1->Print(pdfname);
    mycE_overtrue_totsfshmax_sec_leak_2->Print(pdfname);

    //Total SF per dedx weighted section and shower max method
    tpdf->Clear();
    tpdf->cd();
    lat.SetTextSize(0.02);
    lat.DrawLatex(0.1,0.5,"Scenario: Total SF per dedx weighted section and shower max method for each energy");
    lat.DrawLatex(0.1,0.4,"1. Linearity plots");
    tpdf->Update();
    tpdf->Print(pdfname);
    //This is what I thought it was the linearity plots. Uncomment if needed. 
    // c309->Print(pdfname);
    // c310->Print(pdfname);
    // c307->Print(pdfname);
    c308->Print(pdfname);
    tpdf->Clear();
    tpdf->cd();
    lat.SetTextSize(0.02);
    lat.DrawLatex(0.1,0.5,"Scenario: Total SF per dedx weighted section and shower max method for each energy");
    lat.DrawLatex(0.1,0.4,"2. Resolution plots");
    tpdf->Update();
    tpdf->Print(pdfname);
    // resoFit_totsfshmax_sec_dedxweighted->SetTitle("Relative energy resolution using the total sampling fraction per section for e- and all beam energies");
    // c->Update();
    // c->Print(pdfname);
    // resoFit_totsf_sec_leak->SetTitle("Relative energy resolution using the total sampling fraction per section adding leakage for e- and all beam energies");
    // c250->Update();
    c312->Print(pdfname);
    tpdf->Clear();
    tpdf->cd();
    lat.DrawLatex(0.1,0.5,"Scenario: Total SF per dedx weighted section and shower max method for each energy");
    lat.DrawLatex(0.1,0.4,"3. Fits for the previous plots");
    tpdf->Update();
    tpdf->Print(pdfname);
    // mycE_totsfshmax_sec_dedxweighted_1->Print(pdfname);
    // mycE_totsfshmax_sec_dedxweighted_2->Print(pdfname);
    mycE_totsfshmax_sec_dedxweighted_leak_1->Print(pdfname);
    mycE_totsfshmax_sec_dedxweighted_leak_2->Print(pdfname);
    tpdf->Clear();
    tpdf->cd();
    lat.DrawLatex(0.1,0.5,"Scenario: Total SF per dedx weighted section and shower max method for each energy");
    lat.DrawLatex(0.1,0.4,"4. No fits for  #frac{E_{reco,totsfshmaxsecdedxweigthed}-E_{true}}{E_{true}} plots");
    tpdf->Update();
    tpdf->Print(pdfname);
    // mycE_overtrue_totsfshmax_sec_dedxweighted_1->Print(pdfname);
    // mycE_overtrue_totsfshmax_sec_dedxweighted_2->Print(pdfname);
    mycE_overtrue_totsfshmax_sec_dedxweighted_leak_1->Print(pdfname);
    mycE_overtrue_totsfshmax_sec_dedxweighted_leak_2->Print(pdfname);








    tpdf->Clear();
    tpdf->cd();
    lat.SetTextSize(0.02);
    lat.DrawLatex(0.1,0.5,"Scenario: Total SF per 5 sections and shower max method for each energy");
    lat.DrawLatex(0.1,0.4,"1. Linearity plots");
    tpdf->Update();
    tpdf->Print(pdfname);
    //This is what I thought it was the linearity plots. Uncomment if needed. 
    // c275->Print(pdfname);
    // c276->Print(pdfname);
    // c271->Print(pdfname);
    c272->Print(pdfname);
    tpdf->Clear();
    tpdf->cd();
    lat.SetTextSize(0.02);
    lat.DrawLatex(0.1,0.5,"Scenario: Total SF per 5 sections and shower max method for each energy");
    lat.DrawLatex(0.1,0.4,"2. Resolution plots");
    tpdf->Update();
    tpdf->Print(pdfname);
    // resoFit_totsfshmax_dedxsec->SetTitle("Relative energy resolution using the total sampling fraction per section for e- and all beam energies");
    // c279->Update();
    // c279->Print(pdfname);
    // resoFit_totsf_dedxsec_leak->SetTitle("Relative energy resolution using the total sampling fraction per section adding leakage for e- and all beam energies");
    // c280->Update();
    c280->Print(pdfname);
    tpdf->Clear();
    tpdf->cd();
    lat.DrawLatex(0.1,0.5,"Scenario: Total SF per 5 sections and shower max method for each energy");
    lat.DrawLatex(0.1,0.4,"3. Fits for the previous plots");
    tpdf->Update();
    tpdf->Print(pdfname);
    // mycE_totsfshmax_dedxsec_1->Print(pdfname);
    // mycE_totsfshmax_dedxsec_2->Print(pdfname);
    mycE_totsfshmax_dedxsec_leak_1->Print(pdfname);
    mycE_totsfshmax_dedxsec_leak_2->Print(pdfname);
    tpdf->Clear();
    tpdf->cd();
    lat.DrawLatex(0.1,0.5,"Scenario: Total SF per 5 sections and shower max method for each energy");
    lat.DrawLatex(0.1,0.4,"4. No fits for  #frac{E_{reco,totsfshmaxdedxsec}-E_{true}}{E_{true}} plots");
    tpdf->Update();
    tpdf->Print(pdfname);
    // mycE_overtrue_totsfshmax_dedxsec_1->Print(pdfname);
    // mycE_overtrue_totsfshmax_dedxsec_2->Print(pdfname);
    mycE_overtrue_totsfshmax_dedxsec_leak_1->Print(pdfname);
    mycE_overtrue_totsfshmax_dedxsec_leak_2->Print(pdfname);

    tpdf->Clear();
    tpdf->cd();
    lat.SetTextSize(0.04);
    lat.DrawLatex(0.1,0.5,"Combine resolution plots of: ");
    lat.DrawLatex(0.1,0.4,"1. Total SF and shower max method");
    lat.DrawLatex(0.1,0.3,"2. Total SF per section and shower max method");
    lat.DrawLatex(0.1,0.2,"3. dEdx");
    tpdf->Update();
    tpdf->Print(pdfname);
    resoFit_totsfnorm_leak->SetMarkerColor(1);resoFit_totsfnorm_leak->SetLineColor(1);resoFit_totsfnorm_leak->SetMarkerStyle(20);
    resoFit_totsfshmax_sec_leak->SetMarkerColor(2);resoFit_totsfshmax_sec_leak->SetLineColor(2);resoFit_totsfshmax_sec_leak->SetMarkerStyle(21);
    resoFit_val->SetMarkerColor(3);resoFit_val->SetLineColor(3);resoFit_val->SetMarkerStyle(22);
    resoFit_totsfnorm_leak->GetFunction("reso")->SetLineColor(1);
    resoFit_totsfshmax_sec_leak->GetFunction("reso")->SetLineColor(2);
    resoFit_val->GetFunction("reso")->SetLineColor(3);
    resoFit_totsfnorm_leak->SetTitle("Combining relative energy resolutions for e- and different calibration schemes");
    c252->Update();
    c252->Print(pdfname);
 
    tpdf->Clear();
    tpdf->cd();
    lat.DrawLatex(0.1,0.5,"Combine resolution plots of: ");
    lat.DrawLatex(0.1,0.4,"1. Total SF and shower max method using dedx weight");
    lat.DrawLatex(0.1,0.3,"2. Total SF per dEdx weighted section and shower max method");
    lat.DrawLatex(0.1,0.2,"3. dEdx");
    tpdf->Update();
    tpdf->Print(pdfname);
    resoFit_totsfnorm_dedxweighted_leak->SetMarkerColor(1);resoFit_totsfnorm_dedxweighted_leak->SetLineColor(1);resoFit_totsfnorm_dedxweighted_leak->SetMarkerStyle(20);
    resoFit_totsfshmax_sec_dedxweighted_leak->SetMarkerColor(2);resoFit_totsfshmax_sec_dedxweighted_leak->SetLineColor(2);resoFit_totsfshmax_sec_dedxweighted_leak->SetMarkerStyle(21);
    resoFit_val->SetMarkerColor(3);resoFit_val->SetLineColor(3);resoFit_val->SetMarkerStyle(22);
    resoFit_totsfnorm_dedxweighted_leak->GetFunction("reso")->SetLineColor(1);
    resoFit_totsfshmax_sec_dedxweighted_leak->GetFunction("reso")->SetLineColor(2);
    resoFit_val->GetFunction("reso")->SetLineColor(3);
    resoFit_totsfnorm_dedxweighted_leak->SetTitle("Combining relative energy resolutions for e- and different calibration schemes");
    c313->Update();
    c313->Print(pdfname);
 

    tpdf->Clear();
    tpdf->cd();
    lat.DrawLatex(0.1,0.5,"Combine resolution plots of: ");
    lat.DrawLatex(0.1,0.4,"1. Total SF and shower max method");
    lat.DrawLatex(0.1,0.3,"2. Total SF per 5 sections and shower max method");
    lat.DrawLatex(0.1,0.2,"3. dEdx");
    tpdf->Update();
    tpdf->Print(pdfname);
    resoFit_totsfnorm_leak->SetMarkerColor(1);resoFit_totsfnorm_leak->SetLineColor(1);resoFit_totsfnorm_leak->SetMarkerStyle(20);
    resoFit_totsfshmax_dedxsec_leak->SetMarkerColor(2);resoFit_totsfshmax_dedxsec_leak->SetLineColor(2);resoFit_totsfshmax_dedxsec_leak->SetMarkerStyle(21);
    resoFit_val->SetMarkerColor(3);resoFit_val->SetLineColor(3);resoFit_val->SetMarkerStyle(22);
    resoFit_totsfnorm_leak->GetFunction("reso")->SetLineColor(1);
    resoFit_totsfshmax_dedxsec_leak->GetFunction("reso")->SetLineColor(2);
    resoFit_val->GetFunction("reso")->SetLineColor(3);
    resoFit_totsfnorm_leak->SetTitle("Combining relative energy resolutions for e- and different calibration schemes");
    c283->Update();
    c283->Print(pdfname);
 









    tpdf->Clear();
    tpdf->cd();
    lat.SetTextSize(0.03);
    lat.DrawLatex(0.1,0.5,"Sampling fraction distributions for the various calibrations schemes:");
    lat.DrawLatex(0.1,0.4,"Total SF, Total SF vs energy");
    lat.DrawLatex(0.1,0.3,"Total SF with shower max method profile with/out sections based on dEdx ");
    lat.DrawLatex(0.1,0.2,"SF_{i} per layer with shower max method profile");
    tpdf->Update();
    tpdf->Print(pdfname);
    // tpdf->Clear();//c1->DrawClonePad(); 
    // tpdf->cd();
    c_totsfdistrib->Print(pdfname);
    c_totsfvsene->Print(pdfname);
    c_totsfvstprof->Print(pdfname);
    c_sf_sec1perlayervstprof->Print(pdfname);
    c_sf_sec2perlayervstprof->Print(pdfname);
    c_sf_sec3perlayervstprof->Print(pdfname);
    c_sf_dedxsec1perlayervstprof->Print(pdfname);
    c_sf_dedxsec2perlayervstprof->Print(pdfname);
    c_sf_dedxsec3perlayervstprof->Print(pdfname);
    c_sf_dedxsec4perlayervstprof->Print(pdfname);
    c_sf_dedxsec5perlayervstprof->Print(pdfname);
    c_sf_dedxsec6perlayervstprof->Print(pdfname);
    c_sfperlayervstprof->Print(pdfname);
    c_sfperlayervstprof_1->Print(pdfname);
    c_sfperlayervstprof_2->Print(pdfname+")");


   // tpdf->Clear();
    // tpdf->cd();
    // for(unsigned iL(startlayer); iL<(*ssvec).size()-startlayer; iL++){
    //   c_sfperlayervst_prof[(int) iL]->Print(pdfname);
    // }
    // tpdf->Clear();
    // tpdf->cd();
    // lat.DrawLatex(0.1,0.5,"Last page");
    // tpdf->Update();
    // tpdf->Print(pdfname+")");


    recoE_rawvsbeamenergy->Write();
    recoEvsbeamenergy->Write();
    recoE_leakvsbeamenergy->Write();
    recoE_leak_offsetvsbeamenergy->Write();
    if(upstream){
      recoE_corrvsbeamenergy->Write();
      recoE_corrovertrueEvsbeamenergy->Write();
    }
    recoEovertrueEvsbeamenergy->Write();
    //Reset histos
    // hist->Reset();

    results_com[k]->Close();
  
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

//===============================================================================
std::string DoubleToString (double number1)
{
  std::ostringstream oss1;
  oss1<< number1;
  return oss1.str();
}

//===============================================================================
void buildplot(TCanvas *c1, TH1F* hist, TLegend * leg, TString titleofplot, TString xtitle, TString ytitle, TString procNa, int numene){
  
  hist->SetTitle(titleofplot); 
  hist->GetXaxis()->SetTitle(xtitle); 
  hist->GetYaxis()->SetTitle(ytitle);
  // c1->SetLogy();
  if ( numene == 0){hist->SetLineColor(4);}
  if ( numene == 1){hist->SetLineColor(2);}
  if ( numene == 2){hist->SetLineColor(1);}
  if ( numene == 3){hist->SetLineColor(3);}
  if ( numene == 4){hist->SetLineColor(5);}
  if ( numene == 5){hist->SetLineColor(6);}
  if ( numene == 6){hist->SetLineColor(7);}
  if ( numene == 7){hist->SetLineColor(8);}
  if ( numene == 8){hist->SetLineColor(9);}
  if ( numene == 9){hist->SetLineColor(12);}
  if ( numene == 10){hist->SetLineColor(46);}
  if ( numene == 11){hist->SetLineColor(41);}
  if ( numene == 12){hist->SetLineColor(38);}
  if ( numene == 13){hist->SetLineColor(40);}
  if ( numene == 14){hist->SetLineColor(30);}
  if ( numene == 15){hist->SetLineColor(24);}
  if ( numene == 16){hist->SetLineColor(29);}
  if ( numene == 17){hist->SetLineColor(49);}

  leg->AddEntry(hist, procNa, "L");
  c1->cd();
  c1->Update(); 
  numene == 0 ? hist->Draw("HIST") : hist->Draw("HISTsame");
  numene == 0 ? leg->Draw() : leg->Draw("same");
  c1->Update(); 

}

//===============================================================================
FitResult extractMeanEnergy(TCanvas *mycE_1,TCanvas *mycE_2,  int color,
			    TH1F *histo, int rebin){
  
  FitResult res;
  // bool savecanvas = false;
  // TString savename = "";
  TString varname = histo->GetName();
  std::string varname1 = histo->GetName();
  // int p1 = varname1.find("E");
  // std::string number = varname1.substr(p1 + 1, varname1.length());
  std::size_t found = varname1.find_last_of('_');
  varname1 = varname1.substr(found+1);
  std::string st = boost::regex_replace(varname1,boost::regex("[^0-9]*([0-9]+).*"),std::string("\\1")); // \\1

  // std::cout << "TTTTTTTTTTTTTTT " << st << std::endl;
  double lmax = 0;
  // TCanvas *mycE_1 = new TCanvas("mycE_1_"+varname,"mycE_1",1500,1000);
  // TCanvas *mycE_2 = new TCanvas("mycE_2_"+varname,"mycE_2",1500,1000);
  double subcan1 = 0.;double subcan2 = 0.;
  switch (stoi(st)){ 
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
  if ( stof(st) < 90){mycE_1->cd(subcan1);}
  if ( stof(st) > 90){mycE_2->cd(subcan2);}


  // gStyle->SetOptStat("ksiourmen");  
  gStyle->SetOptStat();  
  histo->SetLineColor(color);
  histo->SetMarkerColor(color);
  // histo->SetLineColor(1);
  histo->Draw("");

  double meanE = histo->GetMean();
  double rmsE = histo->GetRMS();
  double meanEerr = histo->GetMeanError();
  double rmsEerr = histo->GetRMSError();
    
  std::cout << " --- " << varname << " = entries " << histo->GetEntries()
	    << " nbins " << histo->GetNbinsX()
	    << " " << histo->GetBinLowEdge(1)
	    << " " << histo->GetBinLowEdge(histo->GetNbinsX())
	    << " mean " << meanE
	    << " rms " << rmsE
	    << " MeanError " << meanEerr<< " RMSError " <<  rmsEerr
	    << " underflows " << histo->GetBinContent(0)
	    << " overflows " << histo->GetBinContent(histo->GetNbinsX()+1) 
	    << " check this  " << rmsE/meanE
	    << std::endl;
  res.mean = meanE; res.rms = rmsE; res.meanerr = meanEerr; res.rmserr = rmsEerr;

  if (histo->GetMaximum()>lmax) lmax = histo->GetMaximum();

  //fit
  // mycE_1->cd();
  // histo->SetMaximum(lmax);
  //==================
  //Rebin here
  histo->Rebin(rebin);
  //==================

  //fit
  histo->Fit("gaus","LR0","",
	     meanE-2*rmsE,
	     meanE+2*rmsE);
  // TF1 *fitResult = new TF1("fitResult","[0]*TMath::Gaus(x,[1],[2],0)",histo->GetXaxis()->GetXmin(),histo->GetXaxis()->GetXmax());
  TF1 *fitResult = histo->GetFunction("gaus");
  // TF1 *fitResult = new TF1("fitResult","[0]*TMath::Gaus(x,[1],[2],0)",histo->GetXaxis()->GetBinCenter(histo->GetMaximumBin()) - 2*rmsE,histo->GetXaxis()->GetBinCenter(histo->GetMaximumBin()) + 2*rmsE);
  // TF1 *fitResult = new TF1("fitResult","[0]*TMath::Gaus(x,[1],[2],0)",histo->GetXaxis()->GetXmin(),histo->GetXaxis()->GetXmax());
  // fitResult->SetParameters(histo->Integral(),
  // 			   histo->GetXaxis()->GetBinCenter(histo->GetMaximumBin()),
  // 			   histo->GetRMS());
  
  // histo->Fit("fitResult","LRQ+","same",
  // 	    fitResult->GetParameter(1)-2*fitResult->GetParameter(2),
  // 	    fitResult->GetParameter(1)+2*fitResult->GetParameter(2));

  fitResult->SetLineColor(2);

  double meanFitE = fitResult->GetParameter(1);
  double rmsFitE = fitResult->GetParameter(2);                                                                                  
  double meanFitEerr = fitResult->GetParError(1);
  double rmsFitEerr = fitResult->GetParError(2);    
  res.Fitmean = meanFitE; res.sigma = rmsFitE; res.Fitmeanerr = meanFitEerr; res.sigmaerr = rmsFitEerr;
    
  fitResult->Draw("same");
    
  std::cout <<  " Gauss fit: mean " << meanFitE << " RMS " << rmsFitE
	    << " Gauss fit: mean error " <<  meanFitEerr << " RMS error" << rmsFitEerr 
	    << std::endl;

    
  char buf[500];
  TLatex lat;
  //Limits of plots for viewing
  float range = 0.; double xmi = 0.; double xma = 0.; 
  //string for the reco case
  std::string recocase = "recoE_" + st;
  TString xtitle;
  //Putting again the first name here because we alter it in the way
  varname1 = histo->GetName();
  if( varname1.find("val") != std::string::npos && varname1.find("overtrue") == std::string::npos) {
    range = stof(st) < 50 ? 10 : 40;
    xmi = stof(st) - range;   
    xma = stof(st) + range;
    xtitle = "Erec";
  } else if ( varname1.find("raw") != std::string::npos && varname1.find("overtrue") == std::string::npos && varname1.find("leak") == std::string::npos){
    // range = meanE < 20000 ? 5000 : 10000;
    range = meanE < 20000 ? 5000 : 30000;
    xmi = meanE - range;
    xma = meanE + range;
    xtitle = "E_{uncalibrated}";
  } else if ( varname1.find(recocase) != std::string::npos && varname1.find("overtrue") == std::string::npos){
    range = stof(st) < 50 ? 10 : 40;
    xmi = stof(st) - range;   
    xma = stof(st) + range;
    xtitle = "Erec";
  } else if ( varname1.find("section") != std::string::npos && varname1.find("overtrue") == std::string::npos){
    range = stof(st) < 50 ? 10 : 40;
    xmi = stof(st) - range;   
    xma = stof(st) + range;
    xtitle = "Erec";
  } else if ( varname1.find("recoE") != std::string::npos && varname1.find("ch") != std::string::npos && varname1.find("overtrue") == std::string::npos){
    range = stof(st) < 50 ? 10 : 40;
    xmi = stof(st) - range;   
    xma = stof(st) + range;
    xtitle = "Erec";
  } else if ( varname1.find("recoE_leak") != std::string::npos && varname1.find("overtrue") == std::string::npos){
    range = stof(st) < 50 ? 10 : 40;
    xmi = stof(st) - range;   
    xma = stof(st) + range;
    xtitle = "Erec";
  } else if ( varname1.find("recoE_leak_offset") != std::string::npos && varname1.find("overtrue") == std::string::npos){
    range = stof(st) < 50 ? 10 : 40;
    xmi = stof(st) - range;   
    xma = stof(st) + range;
    xtitle = "Erec";
  } else if ( varname1.find("recoE_totsf") != std::string::npos && varname1.find("overtrue") == std::string::npos){
    range = stof(st) < 50 ? 10 : 40;
    xmi = stof(st) - range;   
    xma = stof(st) + range;
    xtitle = "Erec";
  } else if ( varname1.find("recoE_shmax") != std::string::npos && varname1.find("overtrue") == std::string::npos){
    range = stof(st) < 50 ? 10 : 40;
    xmi = stof(st) - range;   
    xma = stof(st) + range;
    xtitle = "Erec";
  } else if ( varname1.find("recoE_totsfnorm") != std::string::npos && varname1.find("overtrue") == std::string::npos){
    range = stof(st) < 50 ? 10 : 60;
    xmi = stof(st) - range;   
    xma = stof(st) + range;
    xtitle = "Erec";
  } else if ( varname1.find("recoE_sfnorm") != std::string::npos && varname1.find("overtrue") == std::string::npos){
    range = stof(st) < 50 ? 10 : 100;
    xmi = stof(st) - range;   
    xma = stof(st) + range;
    xtitle = "Erec";
  } else if ( varname1.find("recoE_sfnorm") != std::string::npos && varname1.find("overtrue") != std::string::npos){
    range = stof(st) > 15 ? 0.25 : 0.6;
    xmi = - range;   
    xma =  range;
    xtitle = "#frac{Erec-Etrue}{Etrue}";
  } else if ( varname1.find("recoE_totsfnorm") != std::string::npos && varname1.find("overtrue") == std::string::npos){
    range = stof(st) < 50 ? 30 : 200;
    xmi = stof(st) - range;   
    xma = stof(st) + range;
    xtitle = "Erec";
  } else if ( varname1.find("val") != std::string::npos && varname1.find("overtrue") != std::string::npos) {
    range = stof(st) > 15 ? 0.25 : 0.6;
    xmi = - range;   
    xma =  range;
    xtitle = "#frac{Erec-Etrue}{Etrue}";
  } else if ( varname1.find("recoEovertrue") != std::string::npos){
    range = stof(st) > 15 ? 0.25 : 0.6;
    xmi = - range;   
    xma =  range;
    xtitle = "#frac{Erec-Etrue}{Etrue}";
  } else if ( varname1.find("recoE") != std::string::npos && varname1.find("ch") != std::string::npos && varname1.find("overtrue") != std::string::npos){
    range = stof(st) > 15 ? 0.25 : 0.6;
    xmi = stof(st) - range;   
    xma = stof(st) + range;
    xtitle = "Erec";
  } else if ( varname1.find("raw") != std::string::npos && varname1.find("overtrue") == std::string::npos && varname1.find("leak") != std::string::npos){
    // range = meanE < 20000 ? 5000 : 10000;
    range = meanE < 20000 ? 5000 : 30000;
    if (meanE > 50000){range = 50000;}
    xmi = meanE - range;
    xma = meanE + range;
    xtitle = "E_{uncalibrated} with leakage";
  } else {
    xmi = histo->GetXaxis()->GetXmin();
    xma = histo->GetXaxis()->GetXmax();
    xtitle = "Erec";
  }
  if( (xmi < 0.) &&  varname1.find("overtrue") == std::string::npos){xmi = 0.;}

  // float range = stof(st) < 50 ? 10 : 40;
  // float range = meanE < 20000 ? 5000 : 10000;
  // double xmi = stof(st) - range;//-0.15;
  // double xmi = meanE - range;//-0.15;
  // if(xmi < 0.){xmi = 0.;}
  // double xma = stof(st)+range;//0.15;
  // double xma = meanE + range;//0.15;
  // double latx = histo->GetXaxis()->GetXmin()+(histo->GetXaxis()->GetXmax()-histo->GetXaxis()->GetXmin())/20.;
  double latx = xmi+(xma-xmi)/20.;
  double laty = histo->GetMaximum();
  sprintf(buf,"<Efit> = %3.3f +/- %3.3f",fitResult->GetParameter(1),fitResult->GetParError(1));
  lat.DrawLatex(latx,laty*0.9,buf);
  sprintf(buf,"RMSfit = %3.3f +/- %3.3f",fitResult->GetParameter(2),fitResult->GetParError(2));
  lat.DrawLatex(latx,laty*0.8,buf);
  sprintf(buf,"RMS/meanfit = %3.3f",fitResult->GetParameter(2)/fitResult->GetParameter(1));
  lat.DrawLatex(latx,laty*0.7,buf);
  
  sprintf(buf,"#chi^{2}/N = %3.3f/%d = %3.3f",fitResult->GetChisquare(),fitResult->GetNDF(),fitResult->GetChisquare()/fitResult->GetNDF());
  lat.DrawLatex(latx,laty*0.6,buf);
  if ( stof(st) < 90){mycE_1->Update();}  
  if ( stof(st) > 90){mycE_2->Update();}  

  // histo->GetXaxis()->SetTitle("#frac{Erec-Etrue}{Etrue}");
  histo->GetXaxis()->SetTitle(xtitle);
  histo->GetXaxis()->SetTitleOffset(1.1);
  histo->GetXaxis()->SetRangeUser(xmi,xma);
  // histo->GetYaxis()->SetRangeUser(0.,lmax);

  std::ostringstream saveName;
  saveName.str("");
  saveName << varname << "GeV_plusfit";
  if ( stof(st) < 90){mycE_1->Update();}  
  if ( stof(st) > 90){mycE_2->Update();}

  TString savename;
  stof(st) < 90 ? savename = mycE_1->GetName() : savename = mycE_2->GetName();
  stof(st) < 90 ? mycE_1->Print(savename+".png") : mycE_2->Print(savename+".png");
 

  // if ( stof(st) < 90){
  //   // mycE_1->Print((saveName.str()+".root").c_str());
  //   mycE_1->Print((saveName.str()+".png").c_str());
  // }
  // if ( stof(st) > 90){
  //   // mycE_2->Print((saveName.str()+".root").c_str());
  //   mycE_2->Print((saveName.str()+".png").c_str());
  // }

  return res;

};
//===============================================================================
//Fit for linearity
//Only assume that xaxis is GeV
void fitforlinearity(TCanvas *c1, TGraphErrors * gr, double fitenelow, double fitenemax, std::string yaxisunit){

  c1->cd();
  gr->Draw();
  TF1 *fitFunc=new TF1("linear","[0]+[1]*x",fitenelow,fitenemax);
  // fitFunc->SetLineColor(4);
  gr->Fit(fitFunc,"RMEQ");
  // gr->GetFunction("linear")->SetLineColor(4);
  char buf[500];
  TLatex lat;
  double latx = gr->GetHistogram()->GetXaxis()->GetXmin()+(gr->GetHistogram()->GetXaxis()->GetXmax()-gr->GetHistogram()->GetXaxis()->GetXmin())/20.;
  double laty = gr->GetHistogram()->GetMaximum();
  // lat.SetTextColor(4);
  sprintf(buf,"<E> #propto p0 + p1 #times E ");
  lat.DrawLatex(latx,laty*0.85,buf);
  sprintf(buf,"p0 = %3.3f #pm %3.3f %s",fitFunc->GetParameter(0),fitFunc->GetParError(0), yaxisunit.c_str());
  lat.DrawLatex(latx,laty*0.75,buf);
  sprintf(buf,"p1 = %3.3f #pm %3.3f %s/GeV",fitFunc->GetParameter(1),fitFunc->GetParError(1), yaxisunit.c_str());
  lat.DrawLatex(latx,laty*0.65,buf);
  sprintf(buf,"chi2/NDF = %3.3f/%d = %3.3f",fitFunc->GetChisquare(),fitFunc->GetNDF(),fitFunc->GetChisquare()/fitFunc->GetNDF());
  lat.DrawLatex(latx,laty*0.55,buf);
  c1->Update(); 
  TString savename = gr->GetName();
  c1->SaveAs(savename+".png");
}

//===============================================================================
//Fit for resolution
FitResultResolution fitforresolution(TCanvas *c1, TGraphErrors * res, int color, bool addNoiseTerm){

  FitResultResolution fitresol;
  char buf[500];
  TLatex lat;
  double latx = 0.;
  double laty = 0.;
  c1->cd();
  res->Draw();
  TF1 *fitFunc2;
  if (addNoiseTerm) {
    fitFunc2 = new TF1("reso","sqrt([0]/x+[1]+[2]/(x*x))",res->GetXaxis()->GetXmin(),res->GetXaxis()->GetXmax());
  } else {
    fitFunc2 = new TF1("reso","sqrt([0]/x+[1])",res->GetXaxis()->GetXmin(),res->GetXaxis()->GetXmax());
  }
  // fitFunc2->SetLineColor(4);
  res->Fit(fitFunc2,"RMEQ");
  res->GetFunction("reso")->SetLineColor(color);
  double sigmaStoch = sqrt(fitFunc2->GetParameter(0));
  double sigmaStochErr = fitFunc2->GetParError(0)/(2*sigmaStoch);
  double sigmaConst = sqrt(fitFunc2->GetParameter(1));
  double sigmaConstErr = fitFunc2->GetParError(1)/(2*sigmaConst);
  double sigmaNoise = 0;
  double sigmaNoiseErr = 0;
  if (addNoiseTerm) {
    sigmaNoise = sqrt(fitFunc2->GetParameter(2));
    sigmaNoiseErr = fitFunc2->GetParError(2)/(2*sigmaNoise);
  }
    
  // std::cout << "=================================" << std::endl;
  // std::cout << "Absolute error constant term " <<  fitFunc2->GetParError(1) << std::endl;

  latx = res->GetHistogram()->GetXaxis()->GetXmin()+(res->GetHistogram()->GetXaxis()->GetXmax()-res->GetHistogram()->GetXaxis()->GetXmin())/10.;
  laty = res->GetHistogram()->GetMaximum();
  // lat.SetTextColor(4);
  if (addNoiseTerm) {
    sprintf(buf,"#frac{#sigma}{E} #propto #frac{s}{#sqrt{E}} #oplus c #oplus #frac{n}{E}");
  } else {
    sprintf(buf,"#frac{#sigma}{E} #propto #frac{s}{#sqrt{E}} #oplus c");
  }
  lat.DrawLatex(latx,laty*0.85,buf);
  sprintf(buf,"s=%3.3f #pm %3.3f",sigmaStoch,sigmaStochErr);
  lat.DrawLatex(latx,laty*0.75,buf);
  sprintf(buf,"c=%3.3f #pm %3.3f",sigmaConst,sigmaConstErr);
  lat.DrawLatex(latx,laty*0.65,buf);
  sprintf(buf,"chi2/NDF = %3.3f/%d = %3.3f",fitFunc2->GetChisquare(),fitFunc2->GetNDF(),fitFunc2->GetChisquare()/fitFunc2->GetNDF());
  lat.DrawLatex(latx,laty*0.55,buf);
  if (addNoiseTerm) {
    sprintf(buf,"n=%3.2f #pm %3.2f",sigmaNoise,sigmaNoiseErr);
    lat.DrawLatex(latx,laty*0.45,buf);
  }

  c1->Update(); 

  fitresol.sigmaStoch = sigmaStoch;
  fitresol.sigmaStochErr = sigmaStochErr;
  fitresol.sigmaConst = sigmaConst;
  fitresol.sigmaConstErr = sigmaConstErr;
  fitresol.sigmaNoise = sigmaNoise;
  fitresol.sigmaNoiseErr = sigmaNoiseErr;

  return fitresol;

}


//===============================================================================
//Build the graph you want
void buildgraph(TCanvas *c1, TGraphErrors * res, TString titleofplot, TString xtitle, TString ytitle, int point, double energy, FitResult fitres, bool dogaussfit){

  double meanE = 0.; double rmsE = 0.; 
  // double meanEerr = 0.; double rmsEerr = 0.;
  if (dogaussfit){
    meanE = fabs(fitres.Fitmean); rmsE = fabs(fitres.sigma); 
    // meanEerr = fabs(fitres.Fitmeanerr); rmsEerr = fabs(fitres.sigmaerr); 
  } else {
    meanE = fabs(fitres.mean); rmsE = fabs(fitres.rms); 
    // meanEerr = fabs(fitres.meanerr); rmsEerr = fabs(fitres.rmserr); 
  }

  //For the y range and the overtrueE plots
  std::string varname = res->GetName();
  if( varname.find("overtrueE") != std::string::npos ){
    std::cout << "For the overtrueE plots taking meanE and meanEerr instead of rms" << std::endl;
    meanE = fabs(fitres.mean); rmsE = fabs(fitres.meanerr); 
  } 

  c1->cd();
  c1->Update(); 
  res->SetPoint( point , energy , meanE );
  res->SetPointError(point , 0. , rmsE );
  point==0 ? res->Draw("APL") : res->Draw("PSL");
  //For the y range we need to draw first the graph!!!
  if( varname.find("overtrueE") != std::string::npos && varname.find("overtrueE_sfnorm") == std::string::npos){
    res->GetYaxis()->SetRangeUser(-0.01,0.18);
  } else if ( varname.find("overtrueE") != std::string::npos && varname.find("overtrueE_sfnorm") != std::string::npos){ 
    res->GetYaxis()->SetRangeUser(-0.01,0.18);
  } else if ( varname.find("overtrueE") != std::string::npos && varname.find("overtrueE_totsfnorm_dedxweighted") != std::string::npos){ 
    res->GetYaxis()->SetRangeUser(-0.01,0.18);
  } else if ( varname.find("overtrueE") != std::string::npos && varname.find("overtrueE_totsf_dedxweighted") != std::string::npos){ 
    res->GetYaxis()->SetRangeUser(-0.01,0.18);
  } else {
    res->GetYaxis()->SetRangeUser(-0.01,520.);
  }
  res->GetXaxis()->SetRangeUser(0.,520.);
  res->SetTitle(titleofplot); 
  res->GetHistogram()->SetYTitle(ytitle); 
  res->GetHistogram()->SetXTitle(xtitle);
  //Line
  TLine *line = new TLine(0.,0.,520.,0.);
  line->SetLineColor(kRed);
  line->Draw();
  //For the Delta of <mean> 
  char buf[500];
  TLatex lat;
  lat.SetTextColor(kRed);
  double latx = 520./5.;
  double laty = 0.09;
  if(point==17 && varname.find("overtrueE") != std::string::npos){
    sprintf(buf,"#Delta<mean> = %3.3f ", res->GetMean(2) );//2 is for the y axis
    lat.DrawLatex(latx,laty*0.65,buf);
  }
  c1-> Update(); 

}

//===============================================================================
//Build the graph you want
void buildresolutiongraph(TCanvas *c1, TGraphErrors * res, TString titleofplot, int point, double energy, FitResult fitres, bool dogaussfit){

  double meanE = 0.; double rmsE = 0.; double meanEerr = 0.; double rmsEerr = 0.;
  if (dogaussfit){
    meanE = fabs(fitres.Fitmean); rmsE = fabs(fitres.sigma); meanEerr = fabs(fitres.Fitmeanerr); rmsEerr = fabs(fitres.sigmaerr); 
  } else {
    meanE = fabs(fitres.mean); rmsE = fabs(fitres.rms); meanEerr = fabs(fitres.meanerr); rmsEerr = fabs(fitres.rmserr); 
  }
  
  c1->cd();
  c1->Update(); 
  res->SetPoint( point , energy , rmsE/meanE );
  double err = rmsE/meanE*sqrt(pow(rmsEerr/rmsE,2)+pow(meanEerr/meanE,2));
  res->SetPointError(point , 0. , err );
  point==0 ? res->Draw("APL") : res->Draw("PSL");
  res->GetXaxis()->SetRangeUser(0.,520.);
  // res->GetHistogram()->SetMaximum(520.); 
  res->SetTitle(titleofplot); 
  res->GetHistogram()->SetYTitle("#sigma/E"); 
  res->GetHistogram()->SetXTitle("Energy (GeV)");
  c1-> Update(); 

}

//===============================================================================
//For the combined resolution in one plot
void combineresolution(TCanvas *c1, TGraphErrors * res1, int color1, FitResultResolution fitresol1, TGraphErrors * res2, int color2, FitResultResolution fitresol2,  TGraphErrors * res3, int color3, FitResultResolution fitresol3, TGraphErrors * res4, int color4, FitResultResolution fitresol4, TString particle, bool addNoiseTerm){
  
  c1->cd();
  res1->SetMarkerStyle(20); 
  // res1->SetMarkerSize(0.2); 
  res1->SetMarkerColor(color1);  
  res1->SetLineColor(color1);
  res1->Draw("AP");
  c1->Update();

  res2->SetMarkerStyle(21);
  // res2->SetMarkerSize(0.2); 
  res2->SetMarkerColor(color2);
  res2->SetLineColor(color2);
  res2->Draw("PS");
  c1->Update();

  res3->SetMarkerStyle(22);
  // res3->SetMarkerSize(0.2); 
  res3->SetMarkerColor(color3);
  res3->SetLineColor(color3);
  res3->Draw("PS");
  c1->Update();

  res4->SetMarkerStyle(23);
  // res4->SetMarkerSize(0.2); 
  res4->SetMarkerColor(color4);
  res4->SetLineColor(color4);
  res4->Draw("PS");
  c1->Update();

  TLegend *leg = new TLegend(0.25,0.55,0.89,0.89);  //coordinates are fractions
                                         //of pad dimensions

  double sigmaStoch = 0.; double sigmaStochErr = 0.; double sigmaConst = 0.; double sigmaConstErr = 0.;
  double sigmaNoise = 0; double sigmaNoiseErr = 0;

  sigmaStoch = fitresol1.sigmaStoch;
  sigmaStochErr = fitresol1.sigmaStochErr;
  sigmaConst = fitresol1.sigmaConst;
  sigmaConstErr = fitresol1.sigmaConstErr;
  sigmaNoise = fitresol1.sigmaNoise;
  sigmaNoiseErr = fitresol1.sigmaNoiseErr;

  char buf1[500];
  if (addNoiseTerm) {
    sprintf(buf1,"#frac{#sigma}{E} #propto #frac{%3.3f #pm %3.3f}{#sqrt{E}} #oplus %3.3f #pm %3.3f #oplus #frac{%3.2f #pm %3.2f}{E}",sigmaStoch,sigmaStochErr,sigmaConst,sigmaConstErr,sigmaNoise,sigmaNoiseErr);
  } else {
    sprintf(buf1,"#frac{#sigma}{E} #propto #frac{%3.3f #pm %3.3f}{#sqrt{E}} #oplus %3.3f #pm %3.3f",sigmaStoch,sigmaStochErr,sigmaConst,sigmaConstErr);
  }

  std::string res1str = "E_{reco}: ";
  res1str += buf1;

  sigmaStoch = fitresol2.sigmaStoch;
  sigmaStochErr = fitresol2.sigmaStochErr;
  sigmaConst = fitresol2.sigmaConst;
  sigmaConstErr = fitresol2.sigmaConstErr;
  sigmaNoise = fitresol2.sigmaNoise;
  sigmaNoiseErr = fitresol2.sigmaNoiseErr;

  char buf2[500];
  if (addNoiseTerm) {
    sprintf(buf2,"#frac{#sigma}{E} #propto #frac{%3.3f #pm %3.3f}{#sqrt{E}} #oplus %3.3f #pm %3.3f #oplus #frac{%3.2f #pm %3.2f}{E}",sigmaStoch,sigmaStochErr,sigmaConst,sigmaConstErr,sigmaNoise,sigmaNoiseErr);
  } else {
    sprintf(buf2,"#frac{#sigma}{E} #propto #frac{%3.3f #pm %3.3f}{#sqrt{E}} #oplus %3.3f #pm %3.3f",sigmaStoch,sigmaStochErr,sigmaConst,sigmaConstErr);
  }

  std::string res2str = "E_{reco,raw}: ";
  res2str += buf2;

  sigmaStoch = fitresol3.sigmaStoch;
  sigmaStochErr = fitresol3.sigmaStochErr;
  sigmaConst = fitresol3.sigmaConst;
  sigmaConstErr = fitresol3.sigmaConstErr;
  sigmaNoise = fitresol3.sigmaNoise;
  sigmaNoiseErr = fitresol3.sigmaNoiseErr;

  char buf3[500];
  if (addNoiseTerm) {
    sprintf(buf3,"#frac{#sigma}{E} #propto #frac{%3.3f #pm %3.3f}{#sqrt{E}} #oplus %3.3f #pm %3.3f #oplus #frac{%3.2f #pm %3.2f}{E}",sigmaStoch,sigmaStochErr,sigmaConst,sigmaConstErr,sigmaNoise,sigmaNoiseErr);
  } else {
    sprintf(buf3,"#frac{#sigma}{E} #propto #frac{%3.3f #pm %3.3f}{#sqrt{E}} #oplus %3.3f #pm %3.3f",sigmaStoch,sigmaStochErr,sigmaConst,sigmaConstErr);
  }

  std::string res3str = "E_{reco,material corr}: ";
  res3str += buf3;

  sigmaStoch = fitresol4.sigmaStoch;
  sigmaStochErr = fitresol4.sigmaStochErr;
  sigmaConst = fitresol4.sigmaConst;
  sigmaConstErr = fitresol4.sigmaConstErr;
  sigmaNoise = fitresol4.sigmaNoise;
  sigmaNoiseErr = fitresol4.sigmaNoiseErr;

  char buf4[500];
  if (addNoiseTerm) {
    sprintf(buf4,"#frac{#sigma}{E} #propto #frac{%3.3f #pm %3.3f}{#sqrt{E}} #oplus %3.3f #pm %3.3f #oplus #frac{%3.2f #pm %3.2f}{E}",sigmaStoch,sigmaStochErr,sigmaConst,sigmaConstErr,sigmaNoise,sigmaNoiseErr);
  } else {
    sprintf(buf4,"#frac{#sigma}{E} #propto #frac{%3.3f #pm %3.3f}{#sqrt{E}} #oplus %3.3f #pm %3.3f",sigmaStoch,sigmaStochErr,sigmaConst,sigmaConstErr);
  }

  std::string res4str = "E_{reco,dedx corr}: ";
  res4str += buf4;

  c1->Update(); 


  // leg->AddEntry(res1,"E_{reco}: " + buf1,"PL");  
  leg->AddEntry(res1,res1str.c_str(),"PL");  
  leg->AddEntry((TObject*)0, "", "");
  leg->AddEntry(res2,res2str.c_str(),"PL");  
  leg->AddEntry((TObject*)0, "", "");
  leg->AddEntry(res3,res3str.c_str(),"PL");  
  leg->AddEntry((TObject*)0, "", "");
  leg->AddEntry(res4,res4str.c_str(),"PL");  
  leg->SetFillColor(0);
  // leg->SetHeader(" Channels ");
                                            
  leg->Draw("PS");
  c1->Update();

  TString titleofplot, titleofplot1;
  titleofplot1 = "Relative energy resolution for "; 
  titleofplot = titleofplot1 + particle;
  
  res1-> SetTitle( titleofplot  );

  res1->GetHistogram()->SetXTitle("Energy (GeV)");
  res1->GetHistogram()->SetYTitle("#sigma/E");
  res1->GetXaxis()->SetTitleOffset(1.2);
  
  // res1->GetYaxis()->SetRangeUser(0.,4000.);
  // res1->GetXaxis()->SetRangeUser(0.,5000.);
  c1->Update();
}

//===============================================================================
//For the combined resolution in one plot : No upstream material in front
void combineresolution_noupstream(TCanvas *c1, TGraphErrors * res1, std::string res1str, int color1, FitResultResolution fitresol1, TGraphErrors * res2,  std::string res2str, int color2, FitResultResolution fitresol2,  TGraphErrors * res3, std::string res3str, int color3, FitResultResolution fitresol3, TString particle, bool addNoiseTerm){
  
  c1->cd();
  res1->SetMarkerStyle(20); 
  // res1->SetMarkerSize(0.2); 
  res1->SetMarkerColor(color1);  
  res1->SetLineColor(color1);
  res1->Draw("AP");
  c1->Update();

  res2->SetMarkerStyle(21);
  // res2->SetMarkerSize(0.2); 
  res2->SetMarkerColor(color2);
  res2->SetLineColor(color2);
  res2->Draw("PS");
  c1->Update();

  res3->SetMarkerStyle(22);
  // res3->SetMarkerSize(0.2); 
  res3->SetMarkerColor(color3);
  res3->SetLineColor(color3);
  res3->Draw("PS");
  c1->Update();

  TLegend *leg = new TLegend(0.25,0.55,0.89,0.89);  //coordinates are fractions
                                         //of pad dimensions

  double sigmaStoch = 0.; double sigmaStochErr = 0.; double sigmaConst = 0.; double sigmaConstErr = 0.;
  double sigmaNoise = 0; double sigmaNoiseErr = 0;

  sigmaStoch = fitresol1.sigmaStoch;
  sigmaStochErr = fitresol1.sigmaStochErr;
  sigmaConst = fitresol1.sigmaConst;
  sigmaConstErr = fitresol1.sigmaConstErr;
  sigmaNoise = fitresol1.sigmaNoise;
  sigmaNoiseErr = fitresol1.sigmaNoiseErr;

  char buf1[500];
  if (addNoiseTerm) {
    sprintf(buf1,"#frac{#sigma}{E} #propto #frac{%3.3f #pm %3.3f}{#sqrt{E}} #oplus %3.3f #pm %3.3f #oplus #frac{%3.2f #pm %3.2f}{E}",sigmaStoch,sigmaStochErr,sigmaConst,sigmaConstErr,sigmaNoise,sigmaNoiseErr);
  } else {
    sprintf(buf1,"#frac{#sigma}{E} #propto #frac{%3.3f #pm %3.3f}{#sqrt{E}} #oplus %3.3f #pm %3.3f",sigmaStoch,sigmaStochErr,sigmaConst,sigmaConstErr);
  }

  // std::string res1str = "E_{sim}: ";
  res1str += buf1;

  sigmaStoch = fitresol2.sigmaStoch;
  sigmaStochErr = fitresol2.sigmaStochErr;
  sigmaConst = fitresol2.sigmaConst;
  sigmaConstErr = fitresol2.sigmaConstErr;
  sigmaNoise = fitresol2.sigmaNoise;
  sigmaNoiseErr = fitresol2.sigmaNoiseErr;

  char buf2[500];
  if (addNoiseTerm) {
    sprintf(buf2,"#frac{#sigma}{E} #propto #frac{%3.3f #pm %3.3f}{#sqrt{E}} #oplus %3.3f #pm %3.3f #oplus #frac{%3.2f #pm %3.2f}{E}",sigmaStoch,sigmaStochErr,sigmaConst,sigmaConstErr,sigmaNoise,sigmaNoiseErr);
  } else {
    sprintf(buf2,"#frac{#sigma}{E} #propto #frac{%3.3f #pm %3.3f}{#sqrt{E}} #oplus %3.3f #pm %3.3f",sigmaStoch,sigmaStochErr,sigmaConst,sigmaConstErr);
  }

  // std::string res2str = "E_{reco}: ";
  res2str += buf2;

  sigmaStoch = fitresol3.sigmaStoch;
  sigmaStochErr = fitresol3.sigmaStochErr;
  sigmaConst = fitresol3.sigmaConst;
  sigmaConstErr = fitresol3.sigmaConstErr;
  sigmaNoise = fitresol3.sigmaNoise;
  sigmaNoiseErr = fitresol3.sigmaNoiseErr;

  char buf3[500];
  if (addNoiseTerm) {
    sprintf(buf3,"#frac{#sigma}{E} #propto #frac{%3.3f #pm %3.3f}{#sqrt{E}} #oplus %3.3f #pm %3.3f #oplus #frac{%3.2f #pm %3.2f}{E}",sigmaStoch,sigmaStochErr,sigmaConst,sigmaConstErr,sigmaNoise,sigmaNoiseErr);
  } else {
    sprintf(buf3,"#frac{#sigma}{E} #propto #frac{%3.3f #pm %3.3f}{#sqrt{E}} #oplus %3.3f #pm %3.3f",sigmaStoch,sigmaStochErr,sigmaConst,sigmaConstErr);
  }

  // std::string res3str = "E_{reco,dEdx corr}: ";
  res3str += buf3;

  c1->Update(); 


  // leg->AddEntry(res1,"E_{reco}: " + buf1,"PL");  
  leg->AddEntry(res1,res1str.c_str(),"PL");  
  leg->AddEntry((TObject*)0, "", "");
  leg->AddEntry(res2,res2str.c_str(),"PL");  
  leg->AddEntry((TObject*)0, "", "");
  leg->AddEntry(res3,res3str.c_str(),"PL");  
  leg->SetFillColor(0);
  // leg->SetHeader(" Channels ");
                                            
  leg->Draw("PS");
  c1->Update();

  TString titleofplot, titleofplot1;
  titleofplot1 = "Relative energy resolution for "; 
  titleofplot = titleofplot1 + particle;
  
  res1-> SetTitle( titleofplot  );

  res1->GetHistogram()->SetXTitle("Energy (GeV)");
  res1->GetHistogram()->SetYTitle("#sigma/E");
  res1->GetXaxis()->SetTitleOffset(1.2);
  res1->GetYaxis()->SetRangeUser(0.,2.);

  res1->GetYaxis()->SetRangeUser(0.,0.15);
  // res1->GetXaxis()->SetRangeUser(0.,5000.);
  c1->Update();


}
//===============================================================================
//For the combined resolution in one plot : No upstream material in front
//4 plots for 4 sections
void combineresolution_noupstream_sections(TCanvas *c1, TGraphErrors * res1, std::string res1str, int color1, FitResultResolution fitresol1, TGraphErrors * res2,  std::string res2str, int color2, FitResultResolution fitresol2,  TGraphErrors * res3, std::string res3str, int color3, FitResultResolution fitresol3, TGraphErrors * res4, std::string res4str, int color4, FitResultResolution fitresol4, TString particle, bool addNoiseTerm){
  
  c1->cd();
  res1->SetMarkerStyle(20); 
  // res1->SetMarkerSize(0.2); 
  res1->SetMarkerColor(color1);  
  res1->SetLineColor(color1);
  res1->Draw("AP");
  c1->Update();

  res2->SetMarkerStyle(21);
  // res2->SetMarkerSize(0.2); 
  res2->SetMarkerColor(color2);
  res2->SetLineColor(color2);
  res2->Draw("PS");
  c1->Update();

  res3->SetMarkerStyle(22);
  // res3->SetMarkerSize(0.2); 
  res3->SetMarkerColor(color3);
  res3->SetLineColor(color3);
  res3->Draw("PS");
  c1->Update();

  res4->SetMarkerStyle(29);
  // res4->SetMarkerSize(0.2); 
  res4->SetMarkerColor(color4);
  res4->SetLineColor(color4);
  res4->Draw("PS");
  c1->Update();

  TLegend *leg = new TLegend(0.2,0.5,0.89,0.89);  //coordinates are fractions
                                         //of pad dimensions

  double sigmaStoch = 0.; double sigmaStochErr = 0.; double sigmaConst = 0.; double sigmaConstErr = 0.;
  double sigmaNoise = 0; double sigmaNoiseErr = 0;

  sigmaStoch = fitresol1.sigmaStoch;
  sigmaStochErr = fitresol1.sigmaStochErr;
  sigmaConst = fitresol1.sigmaConst;
  sigmaConstErr = fitresol1.sigmaConstErr;
  sigmaNoise = fitresol1.sigmaNoise;
  sigmaNoiseErr = fitresol1.sigmaNoiseErr;

  char buf1[500];
  if (addNoiseTerm) {
    sprintf(buf1,"#frac{#sigma}{E} #propto #frac{%3.3f #pm %3.3f}{#sqrt{E}} #oplus %3.3f #pm %3.3f #oplus #frac{%3.2f #pm %3.2f}{E}",sigmaStoch,sigmaStochErr,sigmaConst,sigmaConstErr,sigmaNoise,sigmaNoiseErr);
  } else {
    sprintf(buf1,"#frac{#sigma}{E} #propto #frac{%3.3f #pm %3.3f}{#sqrt{E}} #oplus %3.3f #pm %3.3f",sigmaStoch,sigmaStochErr,sigmaConst,sigmaConstErr);
  }

  // std::string res1str = "E_{sim}: ";
  res1str += buf1;

  sigmaStoch = fitresol2.sigmaStoch;
  sigmaStochErr = fitresol2.sigmaStochErr;
  sigmaConst = fitresol2.sigmaConst;
  sigmaConstErr = fitresol2.sigmaConstErr;
  sigmaNoise = fitresol2.sigmaNoise;
  sigmaNoiseErr = fitresol2.sigmaNoiseErr;

  char buf2[500];
  if (addNoiseTerm) {
    sprintf(buf2,"#frac{#sigma}{E} #propto #frac{%3.3f #pm %3.3f}{#sqrt{E}} #oplus %3.3f #pm %3.3f #oplus #frac{%3.2f #pm %3.2f}{E}",sigmaStoch,sigmaStochErr,sigmaConst,sigmaConstErr,sigmaNoise,sigmaNoiseErr);
  } else {
    sprintf(buf2,"#frac{#sigma}{E} #propto #frac{%3.3f #pm %3.3f}{#sqrt{E}} #oplus %3.3f #pm %3.3f",sigmaStoch,sigmaStochErr,sigmaConst,sigmaConstErr);
  }

  // std::string res2str = "E_{reco}: ";
  res2str += buf2;

  sigmaStoch = fitresol3.sigmaStoch;
  sigmaStochErr = fitresol3.sigmaStochErr;
  sigmaConst = fitresol3.sigmaConst;
  sigmaConstErr = fitresol3.sigmaConstErr;
  sigmaNoise = fitresol3.sigmaNoise;
  sigmaNoiseErr = fitresol3.sigmaNoiseErr;

  char buf3[500];
  if (addNoiseTerm) {
    sprintf(buf3,"#frac{#sigma}{E} #propto #frac{%3.3f #pm %3.3f}{#sqrt{E}} #oplus %3.3f #pm %3.3f #oplus #frac{%3.2f #pm %3.2f}{E}",sigmaStoch,sigmaStochErr,sigmaConst,sigmaConstErr,sigmaNoise,sigmaNoiseErr);
  } else {
    sprintf(buf3,"#frac{#sigma}{E} #propto #frac{%3.3f #pm %3.3f}{#sqrt{E}} #oplus %3.3f #pm %3.3f",sigmaStoch,sigmaStochErr,sigmaConst,sigmaConstErr);
  }

  // std::string res3str = "E_{reco,dEdx corr}: ";
  res3str += buf3;

  sigmaStoch = fitresol4.sigmaStoch;
  sigmaStochErr = fitresol4.sigmaStochErr;
  sigmaConst = fitresol4.sigmaConst;
  sigmaConstErr = fitresol4.sigmaConstErr;
  sigmaNoise = fitresol4.sigmaNoise;
  sigmaNoiseErr = fitresol4.sigmaNoiseErr;

  char buf4[500];
  if (addNoiseTerm) {
    sprintf(buf4,"#frac{#sigma}{E} #propto #frac{%3.3f #pm %3.3f}{#sqrt{E}} #oplus %3.3f #pm %3.3f #oplus #frac{%3.2f #pm %3.2f}{E}",sigmaStoch,sigmaStochErr,sigmaConst,sigmaConstErr,sigmaNoise,sigmaNoiseErr);
  } else {
    sprintf(buf4,"#frac{#sigma}{E} #propto #frac{%3.3f #pm %3.3f}{#sqrt{E}} #oplus %3.3f #pm %3.3f",sigmaStoch,sigmaStochErr,sigmaConst,sigmaConstErr);
  }

  // std::string res4str = "E_{reco,dEdx corr}: ";
  res4str += buf4;

  c1->Update(); 


  // leg->AddEntry(res1,"E_{reco}: " + buf1,"PL");  
  leg->AddEntry(res1,res1str.c_str(),"PL");  
  leg->AddEntry((TObject*)0, "", "");
  leg->AddEntry(res2,res2str.c_str(),"PL");  
  leg->AddEntry((TObject*)0, "", "");
  leg->AddEntry(res3,res3str.c_str(),"PL");  
  leg->AddEntry((TObject*)0, "", "");
  leg->AddEntry(res4,res4str.c_str(),"PL");  
  leg->SetFillColor(0);
  // leg->SetHeader(" Channels ");
                                            
  leg->Draw("PS");
  c1->Update();

  TString titleofplot, titleofplot1;
  titleofplot1 = "Relative energy resolution for "; 
  titleofplot = titleofplot1 + particle;
  
  res1-> SetTitle( titleofplot  );

  res1->GetHistogram()->SetXTitle("Energy (GeV)");
  res1->GetHistogram()->SetYTitle("#sigma/E");
  res1->GetXaxis()->SetTitleOffset(1.2);
  res1->GetYaxis()->SetRangeUser(0.,2.);

  res1->GetYaxis()->SetRangeUser(0.,0.15);
  // res1->GetXaxis()->SetRangeUser(0.,5000.);
  c1->Update();


}


