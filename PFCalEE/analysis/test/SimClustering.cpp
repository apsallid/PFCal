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

//-----------------------------------------------------------------------------------
//Declaration here definition after main
std::string IntToString(int number1);
std::string DoubleToString (double number2);
void buildplot(TCanvas *c1, TH1F* hist, TLegend * leg, TString titleofplot, TString xtitle, TString ytitle, TString procNa, int numene);

//======================================================================================
int main(){//main
  
  //Particle type and energy
  TString particle = "e-"; 
  TString particle_id = "11"; 
  const int numberofenergies = 7; 
  std::vector<int> particleenergies;
  particleenergies.clear();
  particleenergies.push_back(10);
  particleenergies.push_back(20);
  particleenergies.push_back(30);
  particleenergies.push_back(35);
  particleenergies.push_back(50);
  particleenergies.push_back(100);
  particleenergies.push_back(300);
 //======================================================================================
  //Read all files and trees
  TString filename[numberofenergies];
  for (int l=0; l<numberofenergies; l++){
    TString fp = "/tmp/apsallid/Single_" + particle_id;	
    // TString fp = "Single_" + particle_id;	
    TString tp = "_" + IntToString(particleenergies[l]);
    TString sp = ".root";
    filename[l] = fp + tp + sp;
  }

  //The tree in the files that we want to get
  TChain *mytree_simhits[numberofenergies];
  TChain *mytree_clustertrack[numberofenergies];
  TChain *mytree_rechits[numberofenergies];
  
  for (int l=0; l<numberofenergies; l++){
    mytree_simhits[l] = new TChain("FakeClusterCaloFace/mytree_simhits");
    mytree_clustertrack[l] = new TChain("FakeClusterCaloFace/mytree_clustertrack");
    mytree_rechits[l] = new TChain("FakeClusterCaloFace/mytree_rechits");
    mytree_simhits[l]->Add(filename[l]);
    mytree_clustertrack[l]->Add(filename[l]);
    mytree_rechits[l]->Add(filename[l]);
  }

  //Set the branches here
  //Variables in the tree : mytree_simhits and dir: FakeClusterCaloFace directory
  float ene_simhit_matched; //, eta_simhit_matched, phi_simhit_matched; 
  //float ene_simhit_unmatched, eta_simhit_unmatched, phi_simhit_unmatched; 

  //Variables in the tree : mytree_clustertrack and dir: FakeClusterCaloFace directory
  float ene_track, eta_track, phi_track, pt_track, ene_cluster, eta_cluster, phi_cluster, pt_cluster, cluster_x, cluster_y, cluster_z, dr_track_cluster, de_track_cluster;
  float total_ene_rechit_matched, total_ene_rechit_unmatched;
  float total_ene_simhit_matched, total_ene_simhit_unmatched;
  int totalnumofrechits, totalnumofrechits_matched, totalnumofrechits_unmatched;
  int totalnumofsimhits, totalnumofsimhits_matched, totalnumofsimhits_unmatched;
  float hitefficiency; ///totalnumofrechits_matched/totalnumofrechits_matched+unmatched
  float simhitefficiency; // match/matched+unmatched

  //Variables in the tree : mytree_rechits and dir: FakeClusterCaloFace directory
  float ene_rechit_matched, eta_rechit_matched, phi_rechit_matched; 
  float ene_rechit_unmatched, eta_rechit_unmatched, phi_rechit_unmatched; 

  //Set the branches here
  for (int l=0; l<numberofenergies; l++){
    mytree_simhits[l]->SetBranchAddress("ene_simhit_matched",&ene_simhit_matched);

    mytree_clustertrack[l]->SetBranchAddress("ene_track",&ene_track);
    mytree_clustertrack[l]->SetBranchAddress("eta_track",&eta_track);
    mytree_clustertrack[l]->SetBranchAddress("phi_track",&phi_track);
    mytree_clustertrack[l]->SetBranchAddress("pt_track",&pt_track);
    mytree_clustertrack[l]->SetBranchAddress("ene_cluster",&ene_cluster);
    mytree_clustertrack[l]->SetBranchAddress("eta_cluster",&eta_cluster);
    mytree_clustertrack[l]->SetBranchAddress("phi_cluster",&phi_cluster);
    mytree_clustertrack[l]->SetBranchAddress("pt_cluster",&pt_cluster);
    mytree_clustertrack[l]->SetBranchAddress("cluster_x",&cluster_x);
    mytree_clustertrack[l]->SetBranchAddress("cluster_y",&cluster_y);
    mytree_clustertrack[l]->SetBranchAddress("cluster_z",&cluster_z);
    mytree_clustertrack[l]->SetBranchAddress("dr_track_cluster",&dr_track_cluster);
    mytree_clustertrack[l]->SetBranchAddress("de_track_cluster",&de_track_cluster);
    mytree_clustertrack[l]->SetBranchAddress("total_ene_rechit_matched",&total_ene_rechit_matched);
    mytree_clustertrack[l]->SetBranchAddress("total_ene_rechit_unmatched",&total_ene_rechit_unmatched);
    mytree_clustertrack[l]->SetBranchAddress("total_ene_simhit_matched",&total_ene_simhit_matched);
    mytree_clustertrack[l]->SetBranchAddress("total_ene_simhit_unmatched",&total_ene_simhit_unmatched);
    mytree_clustertrack[l]->SetBranchAddress("totalnumofrechits",&totalnumofrechits);
    mytree_clustertrack[l]->SetBranchAddress("totalnumofrechits_matched",&totalnumofrechits_matched);
    mytree_clustertrack[l]->SetBranchAddress("totalnumofrechits_unmatched",&totalnumofrechits_unmatched);
    mytree_clustertrack[l]->SetBranchAddress("totalnumofsimhits",&totalnumofsimhits);
    mytree_clustertrack[l]->SetBranchAddress("totalnumofsimhits_matched",&totalnumofsimhits_matched);
    mytree_clustertrack[l]->SetBranchAddress("totalnumofsimhits_unmatched",&totalnumofsimhits_unmatched);
    mytree_clustertrack[l]->SetBranchAddress("hitefficiency",&hitefficiency);
    mytree_clustertrack[l]->SetBranchAddress("simhitefficiency",&simhitefficiency);

    mytree_rechits[l]->SetBranchAddress("ene_rechit_matched",&ene_rechit_matched);  
    mytree_rechits[l]->SetBranchAddress("eta_rechit_matched",&eta_rechit_matched);  
    mytree_rechits[l]->SetBranchAddress("phi_rechit_matched",&phi_rechit_matched);  
    mytree_rechits[l]->SetBranchAddress("ene_rechit_unmatched",&ene_rechit_unmatched);  
    mytree_rechits[l]->SetBranchAddress("eta_rechit_unmatched",&eta_rechit_unmatched);  
    mytree_rechits[l]->SetBranchAddress("phi_rechit_unmatched",&phi_rechit_unmatched);  
  }
  
  //======================================================================================
  // Histos
  std::vector<TH1F*> h_ene_simhit_matched; h_ene_simhit_matched.clear();

  std::vector<TH1F*> h_ene_track; h_ene_track.clear();
  std::vector<TH1F*> h_eta_track; h_eta_track.clear();
  std::vector<TH1F*> h_phi_track; h_phi_track.clear();
  std::vector<TH1F*> h_pt_track; h_pt_track.clear();
  std::vector<TH1F*> h_ene_cluster; h_ene_cluster.clear();
  std::vector<TH1F*> h_eta_cluster; h_eta_cluster.clear();
  std::vector<TH1F*> h_phi_cluster; h_phi_cluster.clear();
  std::vector<TH1F*> h_pt_cluster; h_pt_cluster.clear();
  std::vector<TH1F*> h_cluster_x; h_cluster_x.clear();
  std::vector<TH1F*> h_cluster_y; h_cluster_y.clear();
  std::vector<TH1F*> h_cluster_z; h_cluster_z.clear();
  std::vector<TH1F*> h_dr_track_cluster; h_dr_track_cluster.clear();
  std::vector<TH1F*> h_de_track_cluster; h_de_track_cluster.clear();
  std::vector<TH1F*> h_total_ene_rechit_matched; h_total_ene_rechit_matched.clear();
  std::vector<TH1F*> h_total_ene_rechit_unmatched; h_total_ene_rechit_unmatched.clear();
  std::vector<TH1F*> h_total_ene_simhit_matched; h_total_ene_simhit_matched.clear();
  std::vector<TH1F*> h_total_ene_simhit_unmatched; h_total_ene_simhit_unmatched.clear();
  std::vector<TH1F*> h_totalnumofrechits; h_totalnumofrechits.clear();
  std::vector<TH1F*> h_totalnumofrechits_matched; h_totalnumofrechits_matched.clear();
  std::vector<TH1F*> h_totalnumofrechits_unmatched; h_totalnumofrechits_unmatched.clear();
  std::vector<TH1F*> h_totalnumofsimhits; h_totalnumofsimhits.clear();
  std::vector<TH1F*> h_totalnumofsimhits_matched; h_totalnumofsimhits_matched.clear();
  std::vector<TH1F*> h_totalnumofsimhits_unmatched; h_totalnumofsimhits_unmatched.clear();
  std::vector<TH1F*> h_hitefficiency; h_hitefficiency.clear();
  std::vector<TH1F*> h_simhitefficiency; h_simhitefficiency.clear();

  std::vector<TH1F*> h_ene_rechit_matched; h_ene_rechit_matched.clear();  
  std::vector<TH1F*> h_eta_rechit_matched; h_eta_rechit_matched.clear();  
  std::vector<TH1F*> h_phi_rechit_matched; h_phi_rechit_matched.clear();  
  std::vector<TH1F*> h_ene_rechit_unmatched; h_ene_rechit_unmatched.clear();  
  std::vector<TH1F*> h_eta_rechit_unmatched; h_eta_rechit_unmatched.clear();  
  std::vector<TH1F*> h_phi_rechit_unmatched; h_phi_rechit_unmatched.clear();  

  std::vector<TH1F*> h_dptoverpt; h_dptoverpt.clear();  
  std::vector<TH1F*> h_deta_track_cluster; h_deta_track_cluster.clear();  
  std::vector<TH1F*> h_dphi_track_cluster; h_dphi_track_cluster.clear();  

  for (Int_t k=0; k<numberofenergies; k++){
    h_ene_simhit_matched.push_back(new TH1F(("h_ene_simhit_matched_" + IntToString(particleenergies[k])).c_str(),";;Events",100,0.,0.00035)); 
    h_ene_track.push_back(new TH1F(("h_ene_track_" + IntToString(particleenergies[k])).c_str(),";;Events",70,0.,700.)); 
    h_eta_track.push_back(new TH1F(("h_eta_track_" + IntToString(particleenergies[k])).c_str(),";;Events",120,-6.,6.)); 
    h_phi_track.push_back(new TH1F(("h_phi_track_" + IntToString(particleenergies[k])).c_str(),";;Events",120,-6.,6.)); 
    h_pt_track.push_back(new TH1F(("h_pt_track_" + IntToString(particleenergies[k])).c_str(),";;Events",500,0.,500.)); 
    h_ene_cluster.push_back(new TH1F(("h_ene_cluster_" + IntToString(particleenergies[k])).c_str(),";;Events",100,0.,0.5)); 
    h_eta_cluster.push_back(new TH1F(("h_eta_cluster_" + IntToString(particleenergies[k])).c_str(),";;Events",120,-6.,6.)); 
    h_phi_cluster.push_back(new TH1F(("h_phi_cluster_" + IntToString(particleenergies[k])).c_str(),";;Events",120,-6.,6.)); 
    h_pt_cluster.push_back(new TH1F(("h_pt_cluster_" + IntToString(particleenergies[k])).c_str(),";;Events",1000,0.,0.05)); 
    h_cluster_x.push_back(new TH1F(("h_cluster_x_" + IntToString(particleenergies[k])).c_str(),";;Events",400,-200.,200.)); 
    h_cluster_y.push_back(new TH1F(("h_cluster_y_" + IntToString(particleenergies[k])).c_str(),";;Events",400,-200,200.)); 
    h_cluster_z.push_back(new TH1F(("h_cluster_z_" + IntToString(particleenergies[k])).c_str(),";;Events",600,-600,600.)); 
    h_dr_track_cluster.push_back(new TH1F(("h_dr_track_cluster_" + IntToString(particleenergies[k])).c_str(),";dR(Track,Cluster);Events",1000,0.,10.)); 
    h_de_track_cluster.push_back(new TH1F(("h_de_track_cluster_" + IntToString(particleenergies[k])).c_str(),";;Events",70,0.,700.)); 
    h_total_ene_rechit_matched.push_back(new TH1F(("h_total_ene_rechit_matched_" + IntToString(particleenergies[k])).c_str(),";;Events",100,0.,0.5)); 
    h_total_ene_rechit_unmatched.push_back(new TH1F(("h_total_ene_rechit_unmatched_" + IntToString(particleenergies[k])).c_str(),";;Events",100,0.,0.5)); 
    h_total_ene_simhit_matched.push_back(new TH1F(("h_total_ene_simhit_matched_" + IntToString(particleenergies[k])).c_str(),";;Events",100,0.,0.01)); 
    h_total_ene_simhit_unmatched.push_back(new TH1F(("h_total_ene_simhit_unmatched_" + IntToString(particleenergies[k])).c_str(),";;Events",100,0.,0.01)); 
    h_totalnumofrechits.push_back(new TH1F(("h_totalnumofrechits_" + IntToString(particleenergies[k])).c_str(),";;Events",100,35000,45000)); 
    h_totalnumofrechits_matched.push_back(new TH1F(("h_totalnumofrechits_matched_" + IntToString(particleenergies[k])).c_str(),";;Events",400,0.,4000.)); 
    h_totalnumofrechits_unmatched.push_back(new TH1F(("h_totalnumofrechits_unmatched_" + IntToString(particleenergies[k])).c_str(),";;Events",500,0.,10000.)); 
    h_totalnumofsimhits.push_back(new TH1F(("h_totalnumofsimhits_" + IntToString(particleenergies[k])).c_str(),";;Events",800,0.,80000.)); 
    h_totalnumofsimhits_matched.push_back(new TH1F(("h_totalnumofsimhits_matched_" + IntToString(particleenergies[k])).c_str(),";;Events",250.,0.,25000.)); 
    h_totalnumofsimhits_unmatched.push_back(new TH1F(("h_totalnumofsimhits_unmatched_" + IntToString(particleenergies[k])).c_str(),";;Events",400,0.,150000.)); 
    h_hitefficiency.push_back(new TH1F(("h_hitefficiency_" + IntToString(particleenergies[k])).c_str(),";;Events",60,0.,0.6)); 
    h_simhitefficiency.push_back(new TH1F(("h_simhitefficiency_" + IntToString(particleenergies[k])).c_str(),";;Events",150,0.,1.5)); 

    h_ene_rechit_matched.push_back(new TH1F(("h_ene_rechit_matched_" + IntToString(particleenergies[k])).c_str(),";;Events",100,0.,0.5)); 
    h_eta_rechit_matched.push_back(new TH1F(("h_eta_rechit_matched_" + IntToString(particleenergies[k])).c_str(),";;Events",64,-3.2,3.2)); 
    h_phi_rechit_matched.push_back(new TH1F(("h_phi_rechit_matched_" + IntToString(particleenergies[k])).c_str(),";;Events",64,-3.2,3.2)); 
    h_ene_rechit_unmatched.push_back(new TH1F(("h_ene_rechit_unmatched_" + IntToString(particleenergies[k])).c_str(),";;Events",100,0.,0.5));   
    h_eta_rechit_unmatched.push_back(new TH1F(("h_eta_rechit_unmatched_" + IntToString(particleenergies[k])).c_str(),";;Events",64,-3.2,3.2));   
    h_phi_rechit_unmatched.push_back(new TH1F(("h_phi_rechit_unmatched_" + IntToString(particleenergies[k])).c_str(),";;Events",64,-3.2,3.2));   


    h_dptoverpt.push_back(new TH1F(("h_dptoverpt_" + IntToString(particleenergies[k])).c_str(),";;Events",100,-3.,3.));   
    h_deta_track_cluster.push_back(new TH1F(("h_deta_track_cluster_" + IntToString(particleenergies[k])).c_str(),";;Events",300,-0.2,0.2));   
    h_dphi_track_cluster.push_back(new TH1F(("h_dphi_track_cluster_" + IntToString(particleenergies[k])).c_str(),";;Events",300,-0.2,0.2));   

  }
  
  //======================================================================================
  //The file that we will store the results of this analysis
  TString res = "SimClustering_output.root";
  TFile* results = new TFile(res,"recreate");

  std::cout << "============================================================"<< std::endl;
  std::cout << "Results file " << res << std::endl;


  //======================================================================================
  //Loop on energies
  for (int l=0; l<numberofenergies; l++){
    //======================================================================================
    // Loop on entries
    for (Int_t ievt=0; ievt<mytree_clustertrack[l]->GetEntries(); ievt++) {
      mytree_clustertrack[l]->GetEntry(ievt);
      /* if (ievt%(10000)==0) std::cout << "Entry " << ievt << std::endl; */
      // if (ievt==100){break;}
      // std::cout << "Entry " << ievt << std::endl;

      h_ene_track[l]->Fill(ene_track); 
      h_eta_track[l]->Fill(eta_track); 
      h_phi_track[l]->Fill(phi_track); 
      h_pt_track[l]->Fill(pt_track); 
      h_ene_cluster[l]->Fill(ene_cluster); 
      h_eta_cluster[l]->Fill(eta_cluster); 
      h_phi_cluster[l]->Fill(phi_cluster); 
      h_pt_cluster[l]->Fill(pt_cluster); 
      h_cluster_x[l]->Fill(cluster_x); 
      h_cluster_y[l]->Fill(cluster_y); 
      h_cluster_z[l]->Fill(cluster_z); 
      h_dr_track_cluster[l]->Fill(dr_track_cluster);                  
      h_de_track_cluster[l]->Fill(de_track_cluster);                  
      h_total_ene_rechit_matched[l]->Fill(total_ene_rechit_matched);          
      h_total_ene_rechit_unmatched[l]->Fill(total_ene_rechit_unmatched);        
      h_total_ene_simhit_matched[l]->Fill(total_ene_simhit_matched);          
      h_total_ene_simhit_unmatched[l]->Fill(total_ene_simhit_unmatched);        
      h_totalnumofrechits[l]->Fill(totalnumofrechits);                 
      h_totalnumofrechits_matched[l]->Fill(totalnumofrechits_matched);         
      h_totalnumofrechits_unmatched[l]->Fill(totalnumofrechits_unmatched);       
      h_totalnumofsimhits[l]->Fill(totalnumofsimhits);                 
      h_totalnumofsimhits_matched[l]->Fill(totalnumofsimhits_matched);         
      h_totalnumofsimhits_unmatched[l]->Fill(totalnumofsimhits_unmatched);       
      h_hitefficiency[l]->Fill(hitefficiency);   
      h_simhitefficiency[l]->Fill(  (float) totalnumofsimhits_matched/ (float) totalnumofsimhits );   
      
      //deta = eta_track - eta_cluster
      h_deta_track_cluster[l]->Fill(eta_track - eta_cluster);         
      //dphi = phi_track - phi_cluster
      h_dphi_track_cluster[l]->Fill(phi_track - phi_cluster); 

    }// loop on entries

    //======================================================================================
    //Write histos 
    h_dr_track_cluster[l]->Write();
    h_ene_track[l]->Write(); 
    h_eta_track[l]->Write(); 
    h_phi_track[l]->Write(); 
    h_pt_track[l]->Write(); 
    h_ene_cluster[l]->Write(); 
    h_eta_cluster[l]->Write(); 
    h_phi_cluster[l]->Write(); 
    h_pt_cluster[l]->Write(); 
    h_cluster_x[l]->Write(); 
    h_cluster_y[l]->Write(); 
    h_cluster_z[l]->Write(); 
    h_dr_track_cluster[l]->Write();                  
    h_de_track_cluster[l]->Write();                  
    h_total_ene_rechit_matched[l]->Write();          
    h_total_ene_rechit_unmatched[l]->Write();        
    h_total_ene_simhit_matched[l]->Write();          
    h_total_ene_simhit_unmatched[l]->Write();        
    h_totalnumofrechits[l]->Write();                 
    h_totalnumofrechits_matched[l]->Write();         
    h_totalnumofrechits_unmatched[l]->Write();       
    h_totalnumofsimhits[l]->Write();                 
    h_totalnumofsimhits_matched[l]->Write();         
    h_totalnumofsimhits_unmatched[l]->Write();       
    h_hitefficiency[l]->Write();                     
    h_simhitefficiency[l]->Write();                
    h_deta_track_cluster[l]->Write();     
    h_dphi_track_cluster[l]->Write();     

    //======================================================================================
    //We should here clear the histograms because we want them empty for the next file. 
    //Reset histos
    h_dr_track_cluster[l]->Reset();
    h_ene_track[l]->Reset(); 
    h_eta_track[l]->Reset(); 
    h_phi_track[l]->Reset(); 
    h_pt_track[l]->Reset(); 
    h_ene_cluster[l]->Reset(); 
    h_eta_cluster[l]->Reset(); 
    h_phi_cluster[l]->Reset(); 
    h_pt_cluster[l]->Reset(); 
    h_cluster_x[l]->Reset(); 
    h_cluster_y[l]->Reset(); 
    h_cluster_z[l]->Reset(); 
    h_dr_track_cluster[l]->Reset();                  
    h_de_track_cluster[l]->Reset();                  
    h_total_ene_rechit_matched[l]->Reset();          
    h_total_ene_rechit_unmatched[l]->Reset();        
    h_total_ene_simhit_matched[l]->Reset();          
    h_total_ene_simhit_unmatched[l]->Reset();        
    h_totalnumofrechits[l]->Reset();                 
    h_totalnumofrechits_matched[l]->Reset();         
    h_totalnumofrechits_unmatched[l]->Reset();       
    h_totalnumofsimhits[l]->Reset();                 
    h_totalnumofsimhits_matched[l]->Reset();         
    h_totalnumofsimhits_unmatched[l]->Reset();       
    h_hitefficiency[l]->Reset();                     
    h_simhitefficiency[l]->Reset(); 
    h_deta_track_cluster[l]->Reset(); 
    h_dphi_track_cluster[l]->Reset(); 


  }//loop on energies

  results->Close();

  //======================================================================================
  //Make plots
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

  //For the legend
  TLegend* leg1 = new TLegend(0.5, 0.7, 0.8, 0.9); leg1->SetHeader("P_{T}"); leg1->SetFillColor(17);
  TLegend* leg2 = new TLegend(0.5, 0.7, 0.8, 0.9); leg2->SetHeader("P_{T}"); leg2->SetFillColor(17);
  TLegend* leg3 = new TLegend(0.5, 0.7, 0.8, 0.9); leg3->SetHeader("P_{T}"); leg3->SetFillColor(17);
  TLegend* leg4 = new TLegend(0.5, 0.7, 0.8, 0.9); leg4->SetHeader("P_{T}"); leg4->SetFillColor(17);
  TLegend* leg5 = new TLegend(0.5, 0.7, 0.8, 0.9); leg5->SetHeader("P_{T}"); leg5->SetFillColor(17);
  TLegend* leg6 = new TLegend(0.5, 0.7, 0.8, 0.9); leg6->SetHeader("P_{T}"); leg6->SetFillColor(17);
  TLegend* leg7 = new TLegend(0.5, 0.7, 0.8, 0.9); leg7->SetHeader("P_{T}"); leg7->SetFillColor(17);
  TLegend* leg8 = new TLegend(0.5, 0.7, 0.8, 0.9); leg8->SetHeader("P_{T}"); leg8->SetFillColor(17);
  TLegend* leg9 = new TLegend(0.5, 0.7, 0.8, 0.9); leg9->SetHeader("P_{T}"); leg9->SetFillColor(17);
  TLegend* leg10 = new TLegend(0.5, 0.7, 0.8, 0.9); leg10->SetHeader("P_{T}"); leg10->SetFillColor(17);
  TLegend* leg11 = new TLegend(0.5, 0.7, 0.8, 0.9); leg11->SetHeader("P_{T}"); leg11->SetFillColor(17);
  TLegend* leg12 = new TLegend(0.5, 0.7, 0.8, 0.9); leg12->SetHeader("P_{T}"); leg12->SetFillColor(17);
  TLegend* leg13 = new TLegend(0.5, 0.7, 0.8, 0.9); leg13->SetHeader("P_{T}"); leg13->SetFillColor(17);
  TLegend* leg14 = new TLegend(0.5, 0.7, 0.8, 0.9); leg14->SetHeader("P_{T}"); leg14->SetFillColor(17);
  TLegend* leg15 = new TLegend(0.5, 0.7, 0.8, 0.9); leg15->SetHeader("P_{T}"); leg15->SetFillColor(17);
  TLegend* leg16 = new TLegend(0.5, 0.7, 0.8, 0.9); leg16->SetHeader("P_{T}"); leg16->SetFillColor(17);
  TLegend* leg17 = new TLegend(0.5, 0.7, 0.8, 0.9); leg17->SetHeader("P_{T}"); leg17->SetFillColor(17);
  TLegend* leg18 = new TLegend(0.5, 0.7, 0.8, 0.9); leg18->SetHeader("P_{T}"); leg18->SetFillColor(17);
  TLegend* leg19 = new TLegend(0.5, 0.7, 0.8, 0.9); leg19->SetHeader("P_{T}"); leg19->SetFillColor(17);
  TLegend* leg20 = new TLegend(0.5, 0.7, 0.8, 0.9); leg20->SetHeader("P_{T}"); leg20->SetFillColor(17);
  TLegend* leg21 = new TLegend(0.5, 0.7, 0.8, 0.9); leg21->SetHeader("P_{T}"); leg21->SetFillColor(17);
  TLegend* leg22 = new TLegend(0.5, 0.7, 0.8, 0.9); leg22->SetHeader("P_{T}"); leg22->SetFillColor(17);
  TLegend* leg23 = new TLegend(0.5, 0.7, 0.8, 0.9); leg23->SetHeader("P_{T}"); leg23->SetFillColor(17);
  TLegend* leg24 = new TLegend(0.5, 0.7, 0.8, 0.9); leg24->SetHeader("P_{T}"); leg24->SetFillColor(17);
  TLegend* leg25 = new TLegend(0.5, 0.7, 0.8, 0.9); leg25->SetHeader("P_{T}"); leg25->SetFillColor(17);
  TLegend* leg26 = new TLegend(0.5, 0.7, 0.8, 0.9); leg26->SetHeader("P_{T}"); leg26->SetFillColor(17);
  TLegend* leg27 = new TLegend(0.5, 0.7, 0.8, 0.9); leg27->SetHeader("P_{T}"); leg27->SetFillColor(17);

  TString procNa[numberofenergies];
  procNa[0] = "10 GeV";
  procNa[1] = "20 GeV";
  procNa[2] = "30 GeV";
  procNa[3] = "35 GeV";
  procNa[4] = "50 GeV";
  procNa[5] = "100 GeV";
  procNa[6] = "300 GeV";

  //The file that we store the histos
  results= TFile::Open(res,"read");
  std::cout << "Opened results file " << res << std::endl;
  
  TH1F *hist1, *hist2, *hist3, *hist4, *hist5, *hist6, *hist7, *hist8, *hist9, *hist10, *hist11, *hist12, *hist13, *hist14, *hist15, *hist16, *hist17, *hist18, *hist19, *hist20, *hist21, *hist22, *hist23, *hist24, *hist25, *hist26, *hist27;

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  TString titleofplot1,titleofplot2,titleofplot; 
  //======================================================================================
  //Loop on energies
  for (int l=0; l<numberofenergies; l++){
    //-------------------------------------------------------------------------------------------------
    titleofplot = "DR(track,cluster) for " + particle; 
    hist1 = (TH1F*)results->Get( ("h_dr_track_cluster_" + IntToString(particleenergies[l])).c_str() );
    hist1->GetYaxis()->SetRangeUser(0.0001,100000.);
    buildplot(c1, hist1, leg1, titleofplot, "DR(track,cluster)", "Events/0.01", procNa[l], l);
    std::cout << "Working on plot: " << hist1->GetName() << std::endl;
    //-------------------------------------------------------------------------------------------------
    titleofplot = "Energy of the track for " + particle; 
    hist2 = (TH1F*)results->Get( ("h_ene_track_" + IntToString(particleenergies[l])).c_str() );
    buildplot(c2, hist2, leg2, titleofplot, "Track Energy", "Events", procNa[l], l);
    std::cout << "Working on plot: " << hist2->GetName() << std::endl;
    //-------------------------------------------------------------------------------------------------
    titleofplot = "Eta of the track for " + particle; 
    hist3 = (TH1F*)results->Get( ("h_eta_track_" + IntToString(particleenergies[l])).c_str() );
    hist3->GetYaxis()->SetRangeUser(0.0001,5000.);
    buildplot(c3, hist3, leg3, titleofplot, "#eta", "Events", procNa[l], l);
    std::cout << "Working on plot: " << hist3->GetName() << std::endl;
    //-------------------------------------------------------------------------------------------------
    titleofplot = "Phi of the track for " + particle; 
    hist4 = (TH1F*)results->Get( ("h_phi_track_" + IntToString(particleenergies[l])).c_str() );
    hist4->GetYaxis()->SetRangeUser(0.0001,5000.);
    buildplot(c4, hist4, leg4, titleofplot, "#phi", "Events", procNa[l], l);
    std::cout << "Working on plot: " << hist4->GetName() << std::endl;
    //-------------------------------------------------------------------------------------------------
    titleofplot = "P_{T} of the track for " + particle; 
    hist5 = (TH1F*)results->Get( ("h_pt_track_" + IntToString(particleenergies[l])).c_str() );
    buildplot(c5, hist5, leg5, titleofplot, "", "Events", procNa[l], l);
    std::cout << "Working on plot: " << hist5->GetName() << std::endl;
    //-------------------------------------------------------------------------------------------------
    titleofplot = "Energy of the cluster for " + particle; 
    hist6 = (TH1F*)results->Get( ("h_ene_cluster_" + IntToString(particleenergies[l])).c_str() );
    buildplot(c6, hist6, leg6, titleofplot, "", "Events", procNa[l], l);
    std::cout << "Working on plot: " << hist6->GetName() << std::endl;
    //-------------------------------------------------------------------------------------------------
    titleofplot = "Eta of the cluster for " + particle; 
    hist7 = (TH1F*)results->Get( ("h_eta_cluster_" + IntToString(particleenergies[l])).c_str() );
    hist7->GetYaxis()->SetRangeUser(0.0001,5000.);
    buildplot(c7, hist7, leg7, titleofplot, "#eta", "Events", procNa[l], l);
    std::cout << "Working on plot: " << hist7->GetName() << std::endl;
    //-------------------------------------------------------------------------------------------------
    titleofplot = "Phi of the cluster for " + particle; 
    hist8 = (TH1F*)results->Get( ("h_phi_cluster_" + IntToString(particleenergies[l])).c_str() );
    hist8->GetYaxis()->SetRangeUser(0.0001,5000.);
    buildplot(c8, hist8, leg8, titleofplot, "#phi", "Events", procNa[l], l);
    std::cout << "Working on plot: " << hist8->GetName() << std::endl;
    //-------------------------------------------------------------------------------------------------
    titleofplot = "P_{T} of the cluster for " + particle; 
    hist9 = (TH1F*)results->Get( ("h_pt_cluster_" + IntToString(particleenergies[l])).c_str() );
    buildplot(c9, hist9, leg9, titleofplot, "", "Events", procNa[l], l);
    std::cout << "Working on plot: " << hist9->GetName() << std::endl;
    //-------------------------------------------------------------------------------------------------
    titleofplot = "Cluster x for " + particle; 
    hist10 = (TH1F*)results->Get( ("h_cluster_x_" + IntToString(particleenergies[l])).c_str() );
    buildplot(c10, hist10, leg10, titleofplot, "", "Events", procNa[l], l);
    std::cout << "Working on plot: " << hist10->GetName() << std::endl;
    //-------------------------------------------------------------------------------------------------
    titleofplot = "Cluster y for " + particle; 
    hist11 = (TH1F*)results->Get( ("h_cluster_y_" + IntToString(particleenergies[l])).c_str() );
    buildplot(c11, hist11, leg11, titleofplot, "", "Events", procNa[l], l);
    std::cout << "Working on plot: " << hist11->GetName() << std::endl;
    //-------------------------------------------------------------------------------------------------
    titleofplot = "Cluster z for " + particle; 
    hist12 = (TH1F*)results->Get( ("h_cluster_z_" + IntToString(particleenergies[l])).c_str() );
    buildplot(c12, hist12, leg12, titleofplot, "", "Events", procNa[l], l);
    std::cout << "Working on plot: " << hist12->GetName() << std::endl;
    //-------------------------------------------------------------------------------------------------
    titleofplot = "DE(track,cluster) for " + particle; 
    hist13 = (TH1F*)results->Get( ("h_de_track_cluster_" + IntToString(particleenergies[l])).c_str() );
    buildplot(c13, hist13, leg13, titleofplot, "DE(track,cluster)", "Events", procNa[l], l);
    std::cout << "Working on plot: " << hist13->GetName() << std::endl;
    //-------------------------------------------------------------------------------------------------
    titleofplot = "Total energy of matched rechits for " + particle; 
    hist14 = (TH1F*)results->Get( ("h_total_ene_rechit_matched_" + IntToString(particleenergies[l])).c_str() );
    buildplot(c14, hist14, leg14, titleofplot, "E_{tot,rechit,matched}", "Events", procNa[l], l);
    std::cout << "Working on plot: " << hist14->GetName() << std::endl;
    //-------------------------------------------------------------------------------------------------
    titleofplot = "Total energy of unmatched rechits for " + particle; 
    hist15 = (TH1F*)results->Get( ("h_total_ene_rechit_unmatched_" + IntToString(particleenergies[l])).c_str() );
    buildplot(c15, hist15, leg15, titleofplot, "E_{tot,rechit,unmatched}", "Events", procNa[l], l);
    std::cout << "Working on plot: " << hist15->GetName() << std::endl;
    //-------------------------------------------------------------------------------------------------
    titleofplot = "Total energy of matched simhits for " + particle; 
    hist16 = (TH1F*)results->Get( ("h_total_ene_simhit_matched_" + IntToString(particleenergies[l])).c_str() );
    buildplot(c16, hist16, leg16, titleofplot, "E_{tot,simhits,matched}", "Events", procNa[l], l);
    std::cout << "Working on plot: " << hist16->GetName() << std::endl;
    //-------------------------------------------------------------------------------------------------
    titleofplot = "Total energy of unmatched simhits for " + particle; 
    hist17 = (TH1F*)results->Get( ("h_total_ene_simhit_unmatched_" + IntToString(particleenergies[l])).c_str() );
    buildplot(c17, hist17, leg17, titleofplot, "E_{tot,simhits,unmatched}", "Events", procNa[l], l);
    std::cout << "Working on plot: " << hist17->GetName() << std::endl;
    //-------------------------------------------------------------------------------------------------
    titleofplot = "Total number of rechits for " + particle; 
    hist18 = (TH1F*)results->Get( ("h_totalnumofrechits_" + IntToString(particleenergies[l])).c_str() );
    hist18->GetYaxis()->SetRangeUser(0.0001,7000.);
    buildplot(c18, hist18, leg18, titleofplot, "N_{tot,rechits}", "Events", procNa[l], l);
    std::cout << "Working on plot: " << hist18->GetName() << std::endl;
    //-------------------------------------------------------------------------------------------------
    titleofplot = "Total number of matched rechits for " + particle; 
    hist19 = (TH1F*)results->Get( ("h_totalnumofrechits_matched_" + IntToString(particleenergies[l])).c_str() );
    hist19->GetYaxis()->SetRangeUser(0.0001,20000.);
    buildplot(c19, hist19, leg19, titleofplot, "N_{tot,matched rechits}", "Events", procNa[l], l);
    std::cout << "Working on plot: " << hist19->GetName() << std::endl;
    //-------------------------------------------------------------------------------------------------
    titleofplot = "Total number of unmatched rechits for " + particle; 
    hist20 = (TH1F*)results->Get( ("h_totalnumofrechits_unmatched_" + IntToString(particleenergies[l])).c_str() );
    buildplot(c20, hist20, leg20, titleofplot, "N_{tot,unmatched rechits}", "Events", procNa[l], l);
    std::cout << "Working on plot: " << hist20->GetName() << std::endl;
    //-------------------------------------------------------------------------------------------------
    titleofplot = "Total number of simhits for " + particle; 
    hist21 = (TH1F*)results->Get( ("h_totalnumofsimhits_" + IntToString(particleenergies[l])).c_str() );
    hist21->GetYaxis()->SetRangeUser(0.0001,10000.);
    buildplot(c21, hist21, leg21, titleofplot, "N_{tot,simhits}", "Events", procNa[l], l);
    std::cout << "Working on plot: " << hist21->GetName() << std::endl;
    //-------------------------------------------------------------------------------------------------
    titleofplot = "Total number of matched simhits for " + particle; 
    hist22 = (TH1F*)results->Get( ("h_totalnumofsimhits_matched_" + IntToString(particleenergies[l])).c_str() );
    hist22->GetYaxis()->SetRangeUser(0.0001,20000.);
    buildplot(c22, hist22, leg22, titleofplot, "N_{tot,matched simhits}", "Events", procNa[l], l);
    std::cout << "Working on plot: " << hist22->GetName() << std::endl;
    //-------------------------------------------------------------------------------------------------
    titleofplot = "Total number of unmatched simhits for " + particle; 
    hist23 = (TH1F*)results->Get( ("h_totalnumofsimhits_unmatched_" + IntToString(particleenergies[l])).c_str() );
    buildplot(c23, hist23, leg23, titleofplot, "N_{tot,unmatched simhits}", "Events", procNa[l], l);
    std::cout << "Working on plot: " << hist23->GetName() << std::endl;
    //-------------------------------------------------------------------------------------------------
    titleofplot = "Hit efficiency for " + particle; 
    hist24 = (TH1F*)results->Get( ("h_hitefficiency_" + IntToString(particleenergies[l])).c_str() );
    hist24->GetYaxis()->SetRangeUser(0.0001,80000.);
    buildplot(c24, hist24, leg24, titleofplot, "", "Events", procNa[l], l);
    std::cout << "Working on plot: " << hist24->GetName() << std::endl;
    //-------------------------------------------------------------------------------------------------
    titleofplot = "Simhit efficiency for " + particle; 
    hist25 = (TH1F*)results->Get( ("h_simhitefficiency_" + IntToString(particleenergies[l])).c_str() );
    hist25->GetYaxis()->SetRangeUser(0.0001,20000.);
    buildplot(c25, hist25, leg25, titleofplot, "", "Events", procNa[l], l);
    std::cout << "Working on plot: " << hist25->GetName() << std::endl;
    //-------------------------------------------------------------------------------------------------
    titleofplot = "#Delta#eta(track,cluster) for " + particle; 
    hist26 = (TH1F*)results->Get( ("h_deta_track_cluster_" + IntToString(particleenergies[l])).c_str() );
    hist26->GetYaxis()->SetRangeUser(0.0001,50000.);
    buildplot(c26, hist26, leg26, titleofplot, "#Delta#eta(track,cluster)", "Events", procNa[l], l);
    std::cout << "Working on plot: " << hist26->GetName() << std::endl;
    // //-------------------------------------------------------------------------------------------------
    titleofplot = "#Delta#phi(track,cluster) for " + particle; 
    hist27 = (TH1F*)results->Get( ("h_dphi_track_cluster_" + IntToString(particleenergies[l])).c_str() );
    hist27->GetYaxis()->SetRangeUser(0.0001,50000.);
    buildplot(c27, hist27, leg27, titleofplot, "#Delta#phi(track,cluster)", "Events", procNa[l], l);
    std::cout << "Working on plot: " << hist27->GetName() << std::endl;
    
    // titleofplot = " for " + particle; 
    // hist = (TH1F*)results->Get( ("_" + IntToString(particleenergies[l])).c_str() );
    // buildplot(c, hist, leg, titleofplot, "", "Events", procNa[l], l);
    // std::cout << "Working on plot: " << hist->GetName() << std::endl;

  }

  //======================================================================================
  //The files that we will store the results of this analysis for the combined plot
  TString res_com = "SimClustering_combined_plots.root";
  TFile* results_com = new TFile(res_com,"recreate");
     
  //======================================================================================
  //Write canvas with combined plot 
  c1->Write();
  c2->Write();
  c3->Write();
  c4->Write();
  c5->Write();
  c6->Write();
  c7->Write();
  c8->Write();
  c9->Write();
  c10->Write();
  c11->Write();
  c12->Write();
  c13->Write();
  c14->Write();
  c15->Write();
  c16->Write();
  c17->Write();
  c18->Write();
  c19->Write();
  c20->Write();
  c21->Write();
  c22->Write();
  c23->Write();
  c24->Write();
  c25->Write();
  c26->Write();
  c27->Write();

  results_com->Close();

  
}



//-----------------------------------------------------------------------------------
//Definitions
std::string IntToString (int number1)
{
  std::ostringstream oss1;
  oss1<< number1;
  return oss1.str();
}

std::string DoubleToString (double number2)
{
  std::ostringstream oss2;
  oss2<< number2;
  return oss2.str();
}

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

  leg->AddEntry(hist, procNa, "L");
  //hist->GetYaxis()->SetRangeUser(0.,900.);
  //hist->GetYaxis()->SetRangeUser(0.,1500.);
  //hist->GetXaxis()->SetRangeUser(0.,520.);
  c1->cd();
  c1->Update(); 
  // if (l==5){mg->Draw("ap same");}
  numene == 0 ? hist->Draw("HIST") : hist->Draw("HISTsame");
  numene == 0 ? leg->Draw() : leg->Draw("same");
  c1->Update(); 

}







  
  

