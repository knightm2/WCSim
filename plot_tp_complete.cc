#include <algorithm>
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <TObject.h>
#include <TCanvas.h>
#include <TGraph2D.h>
#include <TGraph.h>
#include <TF1.h>
#include <TF2.h>
#include <TStyle.h>

void plot_tp_complete(){
// makes 2D plots for theta-phi effective area simulations, with and without reflectors, for a single mPMT configuration
  gROOT->SetBatch(kTRUE);

  TString filepath = "/neut/datasrv2a/tknight/WCSim/";
  TString mPMT_filepath = filepath + "z_NP1_jobs_Dec13_COC/"; 
  TString single_pmt_filepath = filepath + "z_NP1_jobs_Dec13_COC/"; //if simply gating on single PMTs, use same filepath as mPMT_filepath
  TString default_reflector_filepath = filepath + "z_NP1_jobs_Nov30_COC/"; //filepath to files with default reflectors. Only important if using custom reflectors
  TString default_reflector_single_pmt_filepath = filepath + "z_NP1_jobs_Nov30_COC/"; //filepath to files with default reflectors. only important if using custom reflectors
  TString output_path = filepath + "plots/";

  //// parameters to set before running code
  TString mpmttype          = "NP1";  // NP1, NP2
  int nevents               = 250; // depending on value in x_mPMT_tp.mac
  bool single_pmt           = false; //Set to true to gate on a single pmt of the mpmt
  int pmt                   = 19;      // the pmt of the mPMT that is of interest, if using single pmt
  bool custom_reflectors    = false; //If using reflectors that are not the ones currently used in design (As of December 2017)
  TString reflector_angle   = "23deg"; //If using custom reflectors
  TString reflector_height  = "13p5mm"; //If using custom reflectors

  static const int thetamin = -90;     // minimum theta value simulated
  static const int thetamax = 90;     // maximum theta value simulated
  static const int thetaint = 10;     // interval of theta values
  static const int phimin   = -90;     // minimum phi value simulated
  static const int phimax   = 90;     // maximum phi value simulated
  static const int phiint   = 10;     // interval of phi values
//  double radius             = 0.0001;    // radius of source (cm)
  
// The direction of the photons from the source are sweeped through in steps of theta and phi (from the reference of the source)
  static const int s_thetamin = -60;
  static const int s_thetamax = 60;
  static const int s_thetaint = 4;
  static const int s_phimin   = -60;
  static const int s_phimax   = 60;
  static const int s_phiint   = 4;
//  static const int s_theta_steps = 5;  //
  
//  static const int s_phi_steps   = 5;  //
//  int step_angle           = 5; //in degrees

  static const int ntheta = (thetamax-thetamin)/thetaint + 1;  // number of theta values simulated
  static const int nphi = (phimax-phimin)/phiint + 1;   // number of phi values simulated
  static const int n = ntheta*nphi; // number of .root files for refl/norefl
  static const int ns_theta = (s_thetamax-s_thetamin)/s_thetaint + 1;
  static const int ns_phi = (s_phimax-s_phimin)/s_phiint + 1;
  static const int ns = ns_theta*ns_phi;
  double n_pi = 3.14159265359;
  double theta_r, phi_r;
 

//  Light coverage is a rectangle spanned by 1/2 the distance between the current point and all other neighbouring points
// NOTE: although this is only true for angles close to theta = 90,phi = anything, the number of photons generated corrects for it
  double AngCovPerStep = (s_thetaint*n_pi/180.0)*(s_phiint*n_pi/180.0); //Using small angle approx, taking source as 1m from mPMT 
  
  cout << "Generating angular coverage plots for " << mpmttype << "..." << endl;
  
  //// define arrays for theta and phi
  int theta[ntheta];
  TString thetas[ntheta];
  cout << "theta values (degrees):" << endl;
  for (int i=0;i<ntheta;i++){
    theta[i] = thetamin + i*thetaint;
    thetas[i] = Form("%d",theta[i]);
    cout << thetas[i] << ", ";
  }
  cout << endl;
  
  int phi[nphi];
  TString phis[nphi];
  cout << "phi values (degrees):" << endl;
  for (int i=0;i<nphi;i++){
    phi[i] = phimin + i*phiint;
    phis[i] = Form("%d",phi[i]);
    cout << phis[i] << ", ";
  }
  cout << endl;
  
  int s_theta[ns_theta];
  TString s_thetas[ns_theta];
  cout << "source theta values (degrees):" << endl;
  for (int i=0;i<ns_theta;i++){
    s_theta[i] = s_thetamin + i*s_thetaint;
    s_thetas[i] = Form("%d",s_theta[i]);
    cout << s_thetas[i] << ", ";
  }
  cout << endl;
  
  int s_phi[ns_phi];
  TString s_phis[ns_phi];
  cout << "phi values (degrees):" << endl;
  for (int i=0;i<ns_phi;i++){
    s_phi[i] = s_phimin + i*s_phiint;
    s_phis[i] = Form("%d",s_phi[i]);
    cout << s_phis[i] << ", ";
  }
  cout << endl;

  //// define arrays for effective area
//  int nhits_norefl[ntheta][nphi];  // NHits (no reflectors)
//  int nhits_refl[ntheta][nphi];    // NHits (with reflectors)
//  double ea_norefl[ntheta][nphi];  // effective area (no reflectors)
//  double ea_refl[ntheta][nphi];    // effective area (with reflectors)
  double percinc[ntheta][nphi];    // percent increase of reflectors
  double perc_cov_refl[ntheta][nphi];
  double perc_cov_norefl[ntheta][nphi];
  double any_direction_refl[ntheta][nphi];
  double any_direction_norefl[ntheta][nphi];

  int NHits;
  int mPMT_pmt;
  int mPMT;
  int Hits_tp;
  int nphotons;

  TH2D *AngCovReflHist = new TH2D("AngCovReflHist","AngCovRefHist",ntheta,thetamin,thetamax,nphi,phimin,phimax);
  TH2D *AngCovNoReflHist = new TH2D("AngCovNoReflHist","AngCovNoRefHist",ntheta,thetamin,thetamax,nphi,phimin,phimax);
  TH2D *PercIncHist = new TH2D("PercIncHist","PercIncHist",ntheta,thetamin,thetamax,nphi,phimin,phimax);
  TH2D *DirsHitReflHist = new TH2D("DirsHitReflHist","DirsHitReflHist",ntheta,thetamin,thetamax,nphi,phimin,phimax);
  TH2D *DirsHitNoReflHist = new TH2D("DirsHitNoReflHist","DirsHitNoReflHist",ntheta,thetamin,thetamax,nphi,phimin,phimax);
 
  //// open files, get nhits, and calculate EA (with refl)

  //Calculate angular coverage with reflector at each point
//  TString filepath = "/neut/datasrv2a/tknight/WCSim/";
  int ncount = -1;
  for (int itheta=0;itheta<ntheta;itheta++){ for (int iphi=0;iphi<nphi;iphi++){
    theta_r = (double)theta[itheta]*n_pi/180.0;
    phi_r = (double)phi[iphi]*n_pi/180.0;
//    nhits_refl[itheta][iphi] = 0;
//    ea_refl[itheta][iphi] = 0.0;
    perc_cov_refl[itheta][iphi] = 0.0;
    any_direction_refl[itheta][iphi] = 0.0;
    for (int istheta=0;istheta<ns_theta;istheta++){ for (int isphi=0;isphi<ns_phi;isphi++){
//    for (int iy=0;iy<=1;iy++){ for (int iz=0;iz<1;iz++){

//      AngCovPerStep = (s_phiint*n_pi/180.0)*fabs(sin(((double)s_theta[istheta]-(double)phi[iphi]-s_thetaint)*n_pi/180.0) - sin(((double)s_theta[istheta]-(double)phi[iphi]+s_thetaint)*n_pi/180.0));
      ncount = ncount+1;
      cout << "getting NHits from reflector file " << ncount+1 << " of " << n*ns << "..." << endl;
      if(/*abs(s_theta[istheta]+phi[iphi]) > 90 ||*/ abs(s_theta[istheta]-phi[iphi]) > 90){continue;}

      // define filename
      TString filename;
      if(custom_reflectors){
        if(single_pmt){filename = single_pmt_filepath+"refl_"+thetas[itheta]+"/x_mPMT_"+mpmttype+"_refl_"+reflector_angle+"_"+reflector_height+"_"+thetas[itheta]+"_"+phis[iphi]+"_"+s_thetas[istheta]+"_"+s_phis[isphi]+"_flat.root";}  
        else{filename = mPMT_filepath+"refl_"+thetas[itheta]+"/x_mPMT_"+mpmttype+"_refl_"+reflector_angle+"_"+reflector_height+"_"+thetas[itheta]+"_"+phis[iphi]+"_"+s_thetas[istheta]+"_"+s_phis[isphi]+"_flat.root";}
      }
      else{
        if(single_pmt){filename = single_pmt_filepath+"refl_"+thetas[itheta]+"/x_mPMT_"+mpmttype+"_refl_"+thetas[itheta]+"_"+phis[iphi]+"_"+s_thetas[istheta]+"_"+s_phis[isphi]+"_flat.root";}   
        else{filename = mPMT_filepath+"refl_"+thetas[itheta]+"/x_mPMT_"+mpmttype+"_refl_"+thetas[itheta]+"_"+phis[iphi]+"_"+s_thetas[istheta]+"_"+s_phis[isphi]+"_flat.root";}
      }
     // get NHits from file
      TFile *f = new TFile(filename);
      if(f->IsZombie()){
   //     nhits_refl[itheta][iphi] = 0;
   //     ea_refl[itheta][iphi] = 0;
        delete f;
        continue;
      }
      TTree *t = (TTree*)f->Get("CherenkovHits");
      
//      int NHits;
//      int mPMT_pmt;
//      int mPMT;
      Hits_tp = 0;
  
      t->SetBranchAddress("NHits",&NHits);    
      t->SetBranchAddress("mPMT_pmt",&mPMT_pmt);
//      nhits_refl[itheta][iphi] = 0;
      nphotons = t->GetEntries();
//      cout << t->GetEntries() << endl;
      for (int j=0; j<t->GetEntries(); j++){
      //  t->SetBranchAddress("NHits",&NHits);
        t->GetEntry(j);
        if(NHits==1){
          if(single_pmt){
      //      t->SetBranchAddress("mPMT_pmt",&mPMT_pmt);
      //      t->GetEntry(j);
            if(mPMT_pmt != pmt){continue;}
          }
          Hits_tp += 1;        
//          nhits_refl[itheta][iphi] += NHits;
//          t->SetBranchAddress("NHits",&NHits);    
//          }
        }
      }
    // convert NHits to effective area
      //ea_refl[itheta][iphi] = (double)nhits_refl[itheta][iphi]/(double)nevents*n_pi*radius*radius;
//      if(abs(theta_r)>abs(phi_r)){if(iz == 0){AngCovPerStep = (step_distance/100.0)*(step_distance/100.0); }}
//      else if(abs(theta_r)<abs(phi_r)){if(iy == 0){AngCovPerStep = (step_distance/100.0)*(step_distance/100.0);}}
//      else if(iy == 0 && iz == 0){AngCovPerStep = (step_distance/100.0)*(step_distance/100.0);} //Central hit gets full coverage
      perc_cov_refl[itheta][iphi] += (double)Hits_tp*AngCovPerStep/((double)nevents);

//      AngCovPerStep = (step_distance/100.0)*(step_distance/100.0)*cos(theta_r)*cos(phi_r); //Back to normal
      
  //double photon_count_factor = (cos());
      if(Hits_tp >= nphotons/2){any_direction_refl[itheta][iphi] += 1.0;}
      delete t;
      f->Close();
      delete f;
    }}
    AngCovReflHist->SetBinContent(itheta+1,iphi+1,perc_cov_refl[itheta][iphi]);
    DirsHitReflHist->SetBinContent(itheta+1,iphi+1,any_direction_refl[itheta][iphi]);
    cout << "Angular coverage = " << perc_cov_refl[itheta][iphi] << "    Directions Hit = " << any_direction_refl[itheta][iphi] << endl;

  }}


  //// open files, get nhits, and calculate EA + percent increase (no refl)

  ncount = -1;
  for (int itheta=0;itheta<ntheta;itheta++){ for (int iphi=0;iphi<nphi;iphi++){
    theta_r = (double)theta[itheta]*n_pi/180.0;
    phi_r = (double)phi[iphi]*n_pi/180.0;
//    nhits_norefl[itheta][iphi] = 0;
//    ea_norefl[itheta][iphi] = 0.0;
    perc_cov_norefl[itheta][iphi] = 0.0;
    any_direction_norefl[itheta][iphi] = 0.0;
  
  for (int istheta=0;istheta<ns_theta;istheta++){ for (int isphi=0;isphi<ns_phi;isphi++){
//    for (int iy=0;iy<=1;iy++){ for (int iz=0;iz<1;iz++){
       
//      AngCovPerStep = (s_phiint*n_pi/180.0)*fabs(sin(((double)s_theta[istheta]-(double)phi[iphi]-s_thetaint)*n_pi/180.0) - sin(((double)s_theta[istheta]-(double)phi[iphi]+s_thetaint)*n_pi/180.0));
      ncount = ncount+1;
      cout << "getting NHits from no-reflector file " << ncount+1 << " of " << n*ns << "..." << endl;
      if(/*abs(s_theta[istheta]+phi[iphi]) > 90 ||*/ abs(s_theta[istheta]-phi[iphi]) > 90){continue;}

      // define filename
      TString filename;

      if(custom_reflectors){
        if(single_pmt){filename = default_reflector_single_pmt_filepath+"refl_"+thetas[itheta]+"/x_mPMT_"+mpmttype+"_refl_"+thetas[itheta]+"_"+phis[iphi]+"_"+s_thetas[istheta]+"_"+s_phis[isphi]+"_flat.root";}   
        else{filename = default_reflector_filepath+"refl_"+thetas[itheta]+"/x_mPMT_"+mpmttype+"_refl_"+thetas[itheta]+"_"+phis[iphi]+"_"+s_thetas[istheta]+"_"+s_phis[isphi]+"_flat.root";}
      }
      else{
        if(single_pmt){filename = single_pmt_filepath+"norefl_"+thetas[itheta]+"/x_mPMT_"+mpmttype+"_norefl_"+thetas[itheta]+"_"+phis[iphi]+"_"+s_thetas[istheta]+"_"+s_phis[isphi]+"_flat.root";}   
        else{filename = mPMT_filepath+"norefl_"+thetas[itheta]+"/x_mPMT_"+mpmttype+"_norefl_"+thetas[itheta]+"_"+phis[iphi]+"_"+s_thetas[istheta]+"_"+s_phis[isphi]+"_flat.root";}
      }     
// get NHits from file
      TFile *f = new TFile(filename);
      if(f->IsZombie()){
   //     nhits_refl[itheta][iphi] = 0;
   //     ea_refl[itheta][iphi] = 0;
   //
        delete f;
        continue;
      }

      TTree *t = (TTree*)f->Get("CherenkovHits");
      
//      int NHits;
//      int mPMT_pmt;
//      int mPMT;
      Hits_tp = 0;

      t->SetBranchAddress("NHits",&NHits);    
      t->SetBranchAddress("mPMT_pmt",&mPMT_pmt);
//      nhits_norefl[itheta][iphi] = 0;
      nphotons = t->GetEntries(); 
      for (int j=0; j<t->GetEntries(); j++){
      //  t->SetBranchAddress("NHits",&NHits);
        t->GetEntry(j);
        if(NHits==1){
          if(single_pmt){
      //      t->SetBranchAddress("mPMT_pmt",&mPMT_pmt);
      //      t->GetEntry(j);
            if(mPMT_pmt != pmt){continue;}
          }
          Hits_tp += 1;        
        }
      }

    // convert NHits to effective area
      //ea_refl[itheta][iphi] = (double)nhits_refl[itheta][iphi]/(double)nevents*n_pi*radius*radius;
//      if(abs(theta_r)>abs(phi_r)){if(iz == 0){AngCovPerStep = (step_distance/100.0)*(step_distance/100.0); }}
//      else if(abs(theta_r)<abs(phi_r)){if(iy == 0){AngCovPerStep = (step_distance/100.0)*(step_distance/100.0);}}
//      else if(iy == 0 && iz == 0){AngCovPerStep = (step_distance/100.0)*(step_distance/100.0);} //Central hit gets full coverage
      perc_cov_norefl[itheta][iphi] += (double)Hits_tp*AngCovPerStep/((double)nevents);
//      AngCovPerStep = (step_distance/100.0)*(step_distance/100.0)*cos(theta_r)*cos(phi_r); //Back to normal
      if(Hits_tp >= nphotons/2){any_direction_norefl[itheta][iphi] += 1.0;}
      delete t;
      f->Close();
      delete f; 
    }}
    AngCovNoReflHist->SetBinContent(itheta+1,iphi+1,perc_cov_norefl[itheta][iphi]);
    DirsHitNoReflHist->SetBinContent(itheta+1,iphi+1,any_direction_norefl[itheta][iphi]);
    cout << "Angular Coverage = " << perc_cov_norefl[itheta][iphi] << "    Directions Hit = " << any_direction_norefl[itheta][iphi] << endl;

    // calculate percent increase
//    percinc[itheta][iphi] = (ea_refl[itheta][iphi]-ea_norefl[itheta][iphi])/ea_norefl[itheta][iphi]*100.;
//    if(perc_cov_norefl[itheta][iphi] <= 0.0){perc_cov_norefl = 0.000001;}
    percinc[itheta][iphi] = (perc_cov_refl[itheta][iphi] - perc_cov_norefl[itheta][iphi])/perc_cov_norefl[itheta][iphi]*100.;   
    
    PercIncHist->SetBinContent(itheta+1,iphi+1,percinc[itheta][iphi]);
 
 }}

  //// make 2D graphs
  double thetag[n];
  double phig[n];
  double percincg[n];
  double perc_cov_reflg[n];
  double perc_cov_noreflg[n];
  double any_direction_reflg[n];
  double any_direction_noreflg[n];
  ncount = -1;
  for (int itheta=0;itheta<ntheta;itheta++){ for (int iphi=0;iphi<nphi;iphi++){
    ncount = ncount + 1;
    thetag[ncount]     = theta[itheta];
    phig[ncount]       = phi[iphi];
//    ea_reflg[ncount]   = ea_refl[itheta][iphi];
//    ea_noreflg[ncount] = ea_norefl[itheta][iphi];
    percincg[ncount]   = percinc[itheta][iphi];
    perc_cov_reflg[ncount] = perc_cov_refl[itheta][iphi];
    perc_cov_noreflg[ncount] = perc_cov_norefl[itheta][iphi];
    any_direction_reflg[ncount] = any_direction_refl[itheta][iphi];
    any_direction_noreflg[ncount] = any_direction_norefl[itheta][iphi];
  }}
  
  
  TGraph2D *gper_cov_refl   = new TGraph2D(n,thetag,phig,perc_cov_reflg);
  TGraph2D *gper_cov_norefl = new TGraph2D(n,thetag,phig,perc_cov_noreflg);
  TGraph2D *gpercinc   = new TGraph2D(n,thetag,phig,percincg);
  TGraph2D *gany_direction_refl = new TGraph2D(n,thetag,phig,any_direction_reflg);
  TGraph2D *gany_direction_norefl = new TGraph2D(n,thetag,phig,any_direction_noreflg);
  
  //// draw angular coverage plots
  // with reflector
  TCanvas *c1 = new TCanvas("c1",mpmttype+" Effective Area (refl)",660,600);
  gStyle->SetPalette(1);
  c1->SetRightMargin(0.2);
  gper_cov_refl->GetZaxis()->SetTitleOffset(1.5);
  gper_cov_refl->SetMaximum(*max_element(perc_cov_reflg,perc_cov_reflg+n));
  gper_cov_refl->SetMinimum(*min_element(perc_cov_reflg,perc_cov_reflg+n));

  if(custom_reflectors){
    gper_cov_refl->SetTitle("Angular coverage seen (custom refl);#theta (degrees);#phi (degrees);Solid angle");
  }
  else{
    gper_cov_refl->SetTitle("Angular coverage seen (refl);#theta (degrees);#phi (degrees);Solid angle");
  }
  gper_cov_refl->Draw("COLZ");


  // without reflector
  TCanvas *c2 = new TCanvas("c2",mpmttype+" Effective Area (norefl)",660,600);
  gStyle->SetPalette(1);
  c2->SetRightMargin(0.2);
//  gea_norefl->SetTitle(mpmttype+" Effective Area (no refl)");
  gper_cov_norefl->GetZaxis()->SetTitleOffset(1.5);
  gper_cov_norefl->SetMaximum(*max_element(perc_cov_reflg,perc_cov_reflg+n));
  gper_cov_norefl->SetMinimum(*min_element(perc_cov_reflg,perc_cov_reflg+n));

  if(custom_reflectors){
    gper_cov_norefl->SetTitle("Angular coverage seen (default refl);#theta (degrees);#phi (degrees);Solid angle");
  }
  else{
    gper_cov_norefl->SetTitle("Angular coverage seen (no refl);#theta (degrees);#phi (degrees);Solid angle");
  }
  gper_cov_norefl->Draw("COLZ");

  // with/without reflector percent increase
  TCanvas *c3 = new TCanvas("c3",mpmttype+" Effective Area Percent Increase",660,600);
  gStyle->SetPalette(1);
  c3->SetRightMargin(0.2);
  gpercinc->SetTitle("Angular Coverage Percent Increase;#theta (degrees);#phi (degrees);Angular Coverage Increase (%)");
  
//  double maximum = max_element(percinc,percinc+n);

//  gpercinc->SetMaximum(*max_element(percinc,percinc+n));

  if(gpercinc->GetMaximum()>30.0){
    gpercinc->SetMaximum(30.0);
  }
  if(gpercinc->GetMinimum()<-30.0){
    gpercinc->SetMinimum(-30.0);
  }
//  if(*max_element(percinc,percinc+n)<30.0){gpercinc->SetMaximum(*max_element(percinc,percinc+n));}
//  else {gpercinc->SetMaximum(30.0);}
//  if(*min_element(percinc,percinc+n)>-30.0){
//  gpercinc->SetMinimum(-18.0);
//  gpercinc->SetMaximum(18.0);
  //}
  //else {gpercinc->SetMinimum(-30.0);}
  gpercinc->GetZaxis()->SetTitleOffset(1.5);
//  gpercinc->SetMaximum(min(*max_element(percinc,percinc+n),30.));
  gpercinc->Draw("COLZ");

  //Directions hit with reflector
  TCanvas *c4 = new TCanvas("c4",mpmttype+" Directions hit by light at (#theta,#phi)",660,600);
  gStyle->SetPalette(1);
  c4->SetRightMargin(0.2);
  gany_direction_refl->SetTitle(mpmttype+" Directions hit by light at (#theta,#phi) (refl);#theta (degrees);#phi (degrees);Directions Hit");
  gany_direction_refl->GetZaxis()->SetTitleOffset(1.5);
  gany_direction_refl->Draw("COLZ");

  //Directions hit without reflector
  TCanvas *c5 = new TCanvas("c5",mpmttype+" Directions hit by light at (#theta,#phi)",660,600);
  gStyle->SetPalette(1);
  c5->SetRightMargin(0.2);
  gany_direction_norefl->SetTitle(mpmttype+" Directions hit by light at (#theta,#phi) (no refl);#theta (degrees); #phi (degrees);Directions Hit");
  gany_direction_norefl->GetZaxis()->SetTitleOffset(1.5);
  gany_direction_norefl->Draw("COLZ");

  //// write graphs to file
  TString pmt_str;
  pmt_str.Form("%d",pmt);

  if(custom_reflectors){
    if(single_pmt){
      c1->Print(output_path+mpmttype+"_tp_plots_complete_pmt"+pmt_str+"_"+reflector_angle+"_"+reflector_height+".pdf(","pdf");
      c2->Print(output_path+mpmttype+"_tp_plots_complete_pmt"+pmt_str+"_"+reflector_angle+"_"+reflector_height+".pdf","pdf");
      c3->Print(output_path+mpmttype+"_tp_plots_complete_pmt"+pmt_str+"_"+reflector_angle+"_"+reflector_height+".pdf","pdf");
      c4->Print(output_path+mpmttype+"_tp_plots_complete_pmt"+pmt_str+"_"+reflector_angle+"_"+reflector_height+".pdf","pdf");
      c5->Print(output_path+mpmttype+"_tp_plots_complete_pmt"+pmt_str+"_"+reflector_angle+"_"+reflector_height+".pdf)","pdf");
    }
    else{
      c1->Print(output_path+mpmttype+"_tp_plots_complete_"+reflector_angle+"_"+reflector_height+".pdf(","pdf");
      c2->Print(output_path+mpmttype+"_tp_plots_complete_"+reflector_angle+"_"+reflector_height+".pdf","pdf");
      c3->Print(output_path+mpmttype+"_tp_plots_complete_"+reflector_angle+"_"+reflector_height+".pdf","pdf");
      c4->Print(output_path+mpmttype+"_tp_plots_complete_"+reflector_angle+"_"+reflector_height+".pdf","pdf");
      c5->Print(output_path+mpmttype+"_tp_plots_complete_"+reflector_angle+"_"+reflector_height+".pdf)","pdf");
    }
  }
  else{
    if(single_pmt){
      c1->Print(output_path+mpmttype+"_tp_plots_complete_pmt"+pmt_str+".pdf(","pdf");
      c2->Print(output_path+mpmttype+"_tp_plots_complete_pmt"+pmt_str+".pdf","pdf");
      c3->Print(output_path+mpmttype+"_tp_plots_complete_pmt"+pmt_str+".pdf","pdf");
      c4->Print(output_path+mpmttype+"_tp_plots_complete_pmt"+pmt_str+".pdf","pdf");
      c5->Print(output_path+mpmttype+"_tp_plots_complete_pmt"+pmt_str+".pdf)","pdf");
    }
    else{
      c1->Print(output_path+mpmttype+"_tp_plots_complete.pdf(","pdf");
      c2->Print(output_path+mpmttype+"_tp_plots_complete.pdf","pdf");
      c3->Print(output_path+mpmttype+"_tp_plots_complete.pdf","pdf");
      c4->Print(output_path+mpmttype+"_tp_plots_complete.pdf","pdf");
      c5->Print(output_path+mpmttype+"_tp_plots_complete.pdf)","pdf");
    }
  }

  TFile *HistogramFile;

  if(custom_reflectors){
    if(single_pmt){TFile * HistogramFile = new TFile(output_path+"CompleteGraphsAndHistograms_pmt"+pmt_str+"_"+reflector_angle+"_"+reflector_height+".root","RECREATE");}
    else{TFile * HistogramFile = new TFile(output_path+"CompleteGraphsAndHistograms_"+reflector_angle+"_"+reflector_height+".root","RECREATE");}  
  }
  else{
    if(single_pmt){TFile * HistogramFile = new TFile(output_path+"CompleteGraphsAndHistograms_pmt"+pmt_str+".root","RECREATE");}
    else{TFile * HistogramFile = new TFile(output_path+"CompleteGraphsAndHistograms.root","RECREATE");}
  }
  if(HistogramFile->IsOpen()){cout << "Histogram ROOT File opened correctly." << endl;}
  else{cout << "Error: Histogram ROOT File did not open correctly." << endl;}

  //Writing graphs and histograms to root file
  gper_cov_refl->Write("AngularCoverageReflector");
  gper_cov_norefl->Write("AngularCoverageNoReflector");
  gpercinc->Write("AngularCoverageIncrease");
  gany_direction_refl->Write("NumberOfDirectionsReflector");
  gany_direction_norefl->Write("NumberOfDirectionsNoReflector");

  //fitting Histograms with 2D gaussians and writing to root file
  TF2 *AngCovReflHist_fit = new TF2("AngCovReflHist_fit","[0]*TMath::Gaus(x,[1],[2])*TMath::Gaus(y,[3],[4])",-90,90,-90,90);
  AngCovReflHist_fit->SetParameters(*max_element(perc_cov_reflg,perc_cov_reflg+n),0,20,0,20);
  AngCovReflHist_fit->SetContour(2);
  AngCovReflHist_fit->SetLineColor(1);
//  AngCovReflHist->Fit(AngCovReflHist_fit);
//  AngCovReflHist_fit->SetContourLevel(1,AngCovReflHist_fit->GetParameter(0)*0.5);
  AngCovReflHist->Write("AngCovReflHist");

  AngCovReflHist->Fit(AngCovReflHist_fit);
  AngCovReflHist_fit->SetContourLevel(1,AngCovReflHist_fit->GetParameter(0)*0.5);
  AngCovReflHist->Write("AngCovReflHist_Fit");
  AngCovReflHist_fit->Write("AngCovReflHist_FitFunction");

  TF2 *AngCovNoReflHist_fit = new TF2("AngCovReflNoHist_fit","[0]*TMath::Gaus(x,[1],[2])*TMath::Gaus(y,[3],[4])",-90,90,-90,90);
  AngCovNoReflHist_fit->SetParameters(*max_element(perc_cov_noreflg,perc_cov_noreflg+n),0,20,0,20);
  AngCovNoReflHist_fit->SetContour(2);
  AngCovNoReflHist_fit->SetLineColor(1);
//  AngCovNoReflHist->Fit(AngCovReflHist_fit);
//  AngCovNoReflHist_fit->SetContourLevel(1,AngCovReflNoHist_fit->GetParameter(0)*0.5);
  AngCovNoReflHist->Write("AngCovNoReflHist");

  AngCovNoReflHist->Fit(AngCovNoReflHist_fit);
  AngCovNoReflHist_fit->SetContourLevel(1,AngCovNoReflHist_fit->GetParameter(0)*0.5);
  AngCovNoReflHist->Write("AngCovNoReflHist_Fit");
  AngCovNoReflHist_fit->Write("AngCovNoRefl_FitFunction");

  TF2 *PercIncHist_fit = new TF2("PercIncHist_fit","[0]*TMath::Gaus(x,[1],[2])*TMath::Gaus(y,[3],[4])",-90,90,-90,90);
  PercIncHist_fit->SetParameters(*max_element(percincg,percincg+n),0,20,0,20);
  PercIncHist_fit->SetContour(2);
  PercIncHist_fit->SetLineColor(1);
  PercIncHist->Write("PercIncHist");
 
  PercIncHist->Fit(AngCovReflHist_fit);
  PercIncHist_fit->SetContourLevel(1,PercIncHist_fit->GetParameter(0)*0.5);
  PercIncHist->Write("PercIncHist_Fit");
  PercIncHist_fit->Write("PercIncHist_FitFunction");
//  AngCovReflHist->Write("AngCovReflHist");
//  AngCovNoReflHist->Write("AngCovNoReflHist");
//  PercIncHist->Write("PercIncHist");

  DirsHitReflHist->Write("DirsHitReflHist");
 
  DirsHitNoReflHist->Write("DirsHitNoReflHist");

  HistogramFile->Close();  

}












