#define GridGen_C

#include "NueAna/MultiBinAna/GridGen.h"
#include <cstdlib>

GridGen::GridGen()
{
  SetOutputFile();

  SetNDeltaSteps();
  SetNSinSq2Th13Steps();
  SetDeltaRange();
  SetSinSq2Th13Range();
  SetNormalHierarchy();

  SetTheta12();
  SetTheta23();
  SetAbsValDeltaMSq23();
  SetDeltaMSq12();

  SetTheta13();

  SetNExperiments();

  FreezeTheta23(false);

  return;
}
GridGen::~GridGen()
{
}
void GridGen::SetDeltaRange(double l,double h)
{
  DeltaLow = l;
  DeltaHigh = h;

  return;
}
void GridGen::SetSinSq2Th13Range(double l,double h)
{
  if(l<0 || l>1)
  {
    cout<<"Unphysical value of SinSq2Th13.  Setting low value to 0."<<endl;
    l=0;
  }
  if(h<0 || h>1)
  {
    cout<<"Unphysical value of SinSq2Th13.  Setting high value to 1."<<endl;
    h=1;
  }

  SinSq2Th13Low = l;
  SinSq2Th13High = h;

  return;
}
void GridGen::SetTheta12(double val, double errup, double errdn)
{
  Theta12 = val;
  dTheta12_up = errup;
  dTheta12_dn = errdn;
  return;
}
void GridGen::SetTheta23(double val, double errup, double errdn)
{
  Theta23 = val;
  dTheta23_up = errup;
  dTheta23_dn = errdn;
  return;
}
void GridGen::SetAbsValDeltaMSq23(double val, double errup, double errdn)
{
  DeltaMSq23 = TMath::Abs(val);
  dDeltaMSq23_up = errup;
  dDeltaMSq23_dn = errdn;
  return;
}
void GridGen::SetDeltaMSq12(double val, double errup, double errdn)
{
  DeltaMSq12 = val;
  dDeltaMSq12_up = errup;
  dDeltaMSq12_dn = errdn;
  return;
}
void GridGen::SetTheta13(double val, double errup, double errdn)
{
  Theta13 = val;
  dTheta13_up = errup;
  dTheta13_dn = errdn;
  return;
}
void GridGen::AddExtrap(Extrapolate2D* E)
{
  Extrap.push_back(E);
  return;
}
Double_t GridGen::AsymGaus(Double_t sm, Double_t sp)
{
  // sm = sigma-minus
  // sp = sigma-plus
  Double_t u = gRandom->Uniform();
  Double_t g = gRandom->Gaus(0,1);
  if (u>sp/(sp+sm)) {
    return -TMath::Abs(g*sm);
  }
  else {
    return TMath::Abs(g*sp);
  }
}
double GridGen::DrawTheta23(double dtheta)
{
  if (FrozenTheta23) {
    return Theta23;
  }
  else {
    Double_t sstt = 1-TMath::Abs(gRandom->Gaus(0,dtheta));
    Double_t stt;
    Double_t tt;
    if (gRandom->Uniform()>0.5) {
      stt  = TMath::Sqrt(sstt);
      tt   = TMath::Pi()/2+(gRandom->Uniform()>0.5?-1:1)*(TMath::Pi()/2-TMath::ASin(stt));
    }
    else {
      stt  = -TMath::Sqrt(sstt);
      tt   = -TMath::Pi()/2+(gRandom->Uniform()>0.5?-1:1)*(-TMath::Pi()/2-TMath::ASin(stt));
    }
    double theta23new = tt/2;

    return theta23new;
  }
}
void GridGen::Run()
{
  if(Extrap.size()==0)
  {
    cout<<"No Extrapolate2D input.  Quitting..."<<endl;
    return;
  }

  unsigned int ie;
  for(ie=0;ie<Extrap.size();ie++)
  {
    Extrap[ie]->GetPrediction();//this initializes everything (note: if creating Extrapolate2D object from a file, the number of bins, etc. doesn't get set until GetPrediction() is called)
  }

  for(ie=0;ie<Extrap.size();ie++)
  {
    Extrap[ie]->SetOscPar(OscPar::kTh12,Theta12);
    Extrap[ie]->SetOscPar(OscPar::kTh23,Theta23);
    Extrap[ie]->SetOscPar(OscPar::kDeltaM12,DeltaMSq12);
    Extrap[ie]->SetOscPar(OscPar::kDeltaM23,DeltaMSq23);
    if(!NormalHier) Extrap[ie]->InvertMassHierarchy();
    Extrap[ie]->OscillatePrediction();
  }

  int i;
  unsigned int j;
  int nbins=Extrap[0]->Pred_TotalBkgd_VsBinNumber->GetNbinsX();
  double delta,ssq2th13;
  double sig, bkgd;
  double nc,numucc,bnue,tau,nue;

  vector<TTree*> ftree;
  vector< vector<TTree*> > ftree2;
  TTree *ttmp;
  for(i=0;i<nbins;i++)
  {
    ttmp = new TTree(Form("Bin_%i",i),Form("Bin_%i",i));
    ttmp->Branch("Delta",&delta,"Delta/D");
    ttmp->Branch("Th13Axis",&ssq2th13,"Th13Axis/D");
    ttmp->Branch("Signal",&sig,"Signal/D");
    ttmp->Branch("Background",&bkgd,"Background/D");

    ftree.push_back(ttmp);

    ttmp->Reset();

    ftree2.push_back( vector<TTree*>() );

    for(j=0;j<Extrap.size();j++)
    {
      ttmp = new TTree(Form("Bin_%i_Run_%i",i,j),Form("Bin_%i_Run_%i",i,j));
      ttmp->Branch("Delta",&delta,"Delta/D");
      ttmp->Branch("Th13Axis",&ssq2th13,"Th13Axis/D");
      ttmp->Branch("Signal",&nue,"Signal/D");
      ttmp->Branch("NC",&nc,"NC/D");
      ttmp->Branch("NuMuCC",&numucc,"NuMuCC/D");
      ttmp->Branch("BNueCC",&bnue,"BNueCC/D");
      ttmp->Branch("NuTauCC",&tau,"NuTauCC/D");

      ftree2[i].push_back(ttmp);

      ttmp->Reset();
    }
  }

  double delta_increment = 0;
  if(nDeltaSteps>0) delta_increment = (DeltaHigh - DeltaLow)/(nDeltaSteps);
  double ssq2th13_increment = 0;
  if(nSinSq2Th13Steps>0) ssq2th13_increment = (SinSq2Th13High - SinSq2Th13Low)/(nSinSq2Th13Steps);

  int id,is,l;
  int ir,ip;
  int nPID = Extrap[0]->GetNPID();

  l=0;
  for(id=0;id<nDeltaSteps+1;id++)
  {
    delta = (id*delta_increment + DeltaLow)*TMath::Pi();

    for(is=0;is<nSinSq2Th13Steps+1;is++)
    {
      ssq2th13 = (is*ssq2th13_increment + SinSq2Th13Low);
      for(ie=0;ie<Extrap.size();ie++)
      {
        Extrap[ie]->SetDeltaCP(delta);
        Extrap[ie]->SetSinSq2Th13(ssq2th13);
        Extrap[ie]->OscillatePrediction();
      }

      for(i=0;i<nbins;i++)
      {
	sig=0;
	bkgd=0;
	ir = int(i/nPID);
        ip = i%nPID;
	for(ie=0;ie<Extrap.size();ie++)
	{
	  sig += Extrap[ie]->Pred_Signal_VsBinNumber->GetBinContent(i+1);
	  bkgd += Extrap[ie]->Pred_TotalBkgd_VsBinNumber->GetBinContent(i+1);

          nc = Extrap[ie]->Pred[Background::kNC]->GetBinContent(ip+1,ir+1);
          numucc = Extrap[ie]->Pred[Background::kNuMuCC]->GetBinContent(ip+1,ir+1);
          bnue = Extrap[ie]->Pred[Background::kBNueCC]->GetBinContent(ip+1,ir+1);
          tau = Extrap[ie]->Pred[Background::kNuTauCC]->GetBinContent(ip+1,ir+1);
          nue = Extrap[ie]->Pred[Background::kNueCC]->GetBinContent(ip+1,ir+1);

          ftree2[i][ie]->Fill();
	}
	ftree[i]->Fill();
      }

      if(l%100==0) cout<<100.*l/((nDeltaSteps+1)*(nSinSq2Th13Steps+1))<<"% complete"<<endl;
      l++;
    }
  }

  double nPOTNear,nPOTFar;
  double dm23;

  TTree *paramtree = new TTree("paramtree","paramtree");
  paramtree->Branch("nearPOT",&nPOTNear,"nearPOT/D");
  paramtree->Branch("farPOT",&nPOTFar,"farPOT/D");
  paramtree->Branch("Theta12",&Theta12,"Theta12/D");
  paramtree->Branch("Theta23",&Theta23,"Theta23/D");
  paramtree->Branch("DeltaMSq23",&dm23,"DeltaMSq23/D");
  paramtree->Branch("DeltaMSq12",&DeltaMSq12,"DeltaMSq12/D");

  dm23 = DeltaMSq23;
  if(!NormalHier) dm23 = -1.*DeltaMSq23;

  for(ie=0;ie<Extrap.size();ie++)
  {
    nPOTNear = Extrap[ie]->GetNearPOT();
    nPOTFar = Extrap[ie]->GetFarPOT();
    paramtree->Fill();
  }

  TFile *fout = new TFile(gSystem->ExpandPathName(outFileName.c_str()),"RECREATE");
  for(i=0;i<nbins;i++)
  {
    ftree[i]->Write();
    for(j=0;j<Extrap.size();j++)
    {
      ftree2[i][j]->Write();
    }
  }
  paramtree->Write();
  fout->Close();

  return;
}
void GridGen::RunWithOscParErrs(string s)
{
  if(Extrap.size()==0)
  {
    cout<<"No Extrapolate2D input.  Quitting..."<<endl;
    return;
  }

  unsigned int ie;
  for(ie=0;ie<Extrap.size();ie++)
  {
//     Extrap[ie]->SetPrintResult();
    Extrap[ie]->GetPrediction();
  }

  if(Extrap[0]->Pred_TotalBkgd_VsBinNumber->GetNbinsX()!=1)
  {
    cout<<"RunWithOscParErrs() only works for 1 bin right now.  Quitting..."<<endl;
    return;
  }

  for(ie=0;ie<Extrap.size();ie++)
  {
    Extrap[ie]->SetOscPar(OscPar::kTh12,Theta12);
    Extrap[ie]->SetOscPar(OscPar::kTh23,Theta23);
    Extrap[ie]->SetOscPar(OscPar::kDeltaM12,DeltaMSq12);
    Extrap[ie]->SetOscPar(OscPar::kDeltaM23,DeltaMSq23);
    if(!NormalHier) Extrap[ie]->InvertMassHierarchy();
    Extrap[ie]->OscillatePrediction();
  }

  double delta,t13axis;
  double sig, bkgd;
  double oscparerr;

  TTree *ftree = new TTree("Bin_0","Bin_0");
  ftree->Branch("Delta",&delta,"Delta/D");
  ftree->Branch("Th13Axis",&t13axis,"Th13Axis/D");
  ftree->Branch("Signal",&sig,"Signal/D");
  ftree->Branch("Background",&bkgd,"Background/D");
  ftree->Branch("DNExp_DOscPars",&oscparerr,"DNExp_DOscPars/D");

  double delta_increment = 0;
  if(nDeltaSteps>0) delta_increment = (DeltaHigh - DeltaLow)/(nDeltaSteps);
  double ssq2th13_increment = 0;
  if(nSinSq2Th13Steps>0) ssq2th13_increment = (SinSq2Th13High - SinSq2Th13Low)/(nSinSq2Th13Steps);

  int id,is,l,u;
  double theta23,theta12,dm21,dm32,ssq2th13;
  double nexp,nobs;
  double delnexp;
  TH1D *delnexphist = new TH1D("delnexphist","",400,-1,1);

  gRandom->SetSeed(0);

  TFile *f = new TFile(gSystem->ExpandPathName(s.c_str()),"RECREATE");//save the osc par err distributions to this file

  l=0;
  for(id=0;id<nDeltaSteps+1;id++)
  {
    delta = (id*delta_increment + DeltaLow)*TMath::Pi();

    for(is=0;is<nSinSq2Th13Steps+1;is++)
    {
      t13axis = (is*ssq2th13_increment + SinSq2Th13Low);
      ssq2th13 = t13axis/(2.*TMath::Sin(Theta23)*TMath::Sin(Theta23));

      //get nominal prediction
      sig=0;
      bkgd=0;
      for(ie=0;ie<Extrap.size();ie++)
      {
        Extrap[ie]->SetOscPar(OscPar::kTh12,Theta12);
        Extrap[ie]->SetOscPar(OscPar::kTh23,Theta23);
        Extrap[ie]->SetOscPar(OscPar::kDeltaM12,DeltaMSq12);
        Extrap[ie]->SetOscPar(OscPar::kDeltaM23,DeltaMSq23);//note that DeltaMSq23 is always positive
        if(!NormalHier) Extrap[ie]->InvertMassHierarchy();
        Extrap[ie]->SetDeltaCP(delta);
        Extrap[ie]->SetSinSq2Th13(ssq2th13);
        Extrap[ie]->OscillatePrediction();
        bkgd += (Extrap[ie]->Pred_TotalBkgd_VsBinNumber->GetBinContent(1));
        sig += (Extrap[ie]->Pred_Signal_VsBinNumber->GetBinContent(1));
      }
      nexp = sig+bkgd;

      //do pseudo experiments
      delnexphist->Reset();
      delnexphist->SetName(Form("DeltaNexp_%i_%i",id,is));
      for(u=0;u<NumExpts;u++)
      {
//         theta12 = gRandom->Gaus(Theta12,dTheta12);
//         dm21 = gRandom->Gaus(DeltaMSq12,dDeltaMSq12);
//         dm32 = gRandom->Gaus(DeltaMSq23,dDeltaMSq23);
//         theta23 = gRandom->Gaus(Theta23,dTheta23);
        theta12 = Theta12 + AsymGaus(dTheta12_dn,dTheta12_up);
        dm21 = DeltaMSq12 + AsymGaus(dDeltaMSq12_dn,dDeltaMSq12_up);
        dm32 = DeltaMSq23 + AsymGaus(dDeltaMSq23_dn,dDeltaMSq23_up);
        //theta23 = DrawTheta23(dTheta23_up);
        theta23 = Theta23 + AsymGaus(dTheta23_dn,dTheta23_up);

        ssq2th13 = t13axis/(2.*TMath::Sin(theta23)*TMath::Sin(theta23));

        nobs=0;
        for(ie=0;ie<Extrap.size();ie++)
        {
          Extrap[ie]->SetOscPar(OscPar::kTh12,theta12);
          Extrap[ie]->SetOscPar(OscPar::kDeltaM12,dm21);
          Extrap[ie]->SetOscPar(OscPar::kDeltaM23,dm32);
          if(!NormalHier) Extrap[ie]->InvertMassHierarchy();
          Extrap[ie]->SetOscPar(OscPar::kTh23,theta23);
          Extrap[ie]->SetSinSq2Th13(ssq2th13);
          Extrap[ie]->OscillatePrediction();
          nobs += (Extrap[ie]->Pred_TotalBkgd_VsBinNumber->GetBinContent(1));
          nobs += (Extrap[ie]->Pred_Signal_VsBinNumber->GetBinContent(1));
        }
        delnexp = (nobs-nexp)/nexp;
        delnexphist->Fill(delnexp);
      }
      oscparerr = delnexphist->GetRMS();
      delnexphist->Write();
      f->Close();

      gROOT->cd("/");
      ftree->Fill();

      f = new TFile(gSystem->ExpandPathName(s.c_str()),"UPDATE");

      if(l%100==0) cout<<100.*l/((nDeltaSteps+1)*(nSinSq2Th13Steps+1))<<"% complete"<<endl;
      l++;
    }
  }

  f->Close();

  double nPOTNear,nPOTFar;

  TTree *paramtree = new TTree("paramtree","paramtree");
  paramtree->Branch("nearPOT",&nPOTNear,"nearPOT/D");
  paramtree->Branch("farPOT",&nPOTFar,"farPOT/D");
  paramtree->Branch("Theta12",&Theta12,"Theta12/D");
  paramtree->Branch("Theta23",&Theta23,"Theta23/D");
  paramtree->Branch("DeltaMSq23",&dm32,"DeltaMSq23/D");
  paramtree->Branch("DeltaMSq12",&DeltaMSq12,"DeltaMSq12/D");

  dm32 = DeltaMSq23;
  if(!NormalHier) dm32 = -1.*DeltaMSq23;

  for(ie=0;ie<Extrap.size();ie++)
  {
    nPOTNear = Extrap[ie]->GetNearPOT();
    nPOTFar = Extrap[ie]->GetFarPOT();
    paramtree->Fill();
  }

  TFile *fout = new TFile(gSystem->ExpandPathName(outFileName.c_str()),"RECREATE");
  ftree->Write();
  paramtree->Write();
  fout->Close();

  return;
}
void GridGen::RunMultiBinOscParErrs(string s)
{
  if(Extrap.size()==0)
  {
    cout<<"No Extrapolate2D input.  Quitting..."<<endl;
    return;
  }

  unsigned int ie;
  for(ie=0;ie<Extrap.size();ie++)
  {
    Extrap[ie]->GetPrediction();
  }

  for(ie=0;ie<Extrap.size();ie++)
  {
    Extrap[ie]->SetOscPar(OscPar::kTh12,Theta12);
    Extrap[ie]->SetOscPar(OscPar::kTh23,Theta23);
    Extrap[ie]->SetOscPar(OscPar::kDeltaM12,DeltaMSq12);
    Extrap[ie]->SetOscPar(OscPar::kDeltaM23,DeltaMSq23);
    if(!NormalHier) Extrap[ie]->InvertMassHierarchy();
    Extrap[ie]->OscillatePrediction();
  }

  int i,k;
  unsigned int j;
  int nbins=Extrap[0]->Pred_TotalBkgd_VsBinNumber->GetNbinsX();
  int nPID = Extrap[0]->GetNPID();
  double delta,t13axis;
  vector<double> sig, bkgd;
  vector< vector<double> > nc,numucc,bnue,tau,nue;
  vector<double> oscparerr;
  vector<double> oscparerr_offdiag;
  int noff;

  for(i=0;i<nbins;i++)
  {
    sig.push_back(0);
    bkgd.push_back(0);
    oscparerr.push_back(0);
    nc.push_back( vector<double>() );
    numucc.push_back( vector<double>() );
    bnue.push_back( vector<double>() );
    tau.push_back( vector<double>() );
    nue.push_back( vector<double>() );
    for(k=0;k<nbins;k++)
    {
      if(k>i)
      {
        oscparerr_offdiag.push_back(0);
      }
    }
    for(j=0;j<Extrap.size();j++)
    {
      nc[i].push_back(0);
      numucc[i].push_back(0);
      bnue[i].push_back(0);
      tau[i].push_back(0);
      nue[i].push_back(0);
    }
  }

  vector<TTree*> ftree;
  vector< vector<TTree*> > ftree2;
  TTree *ttmp;
  noff=0;
  for(i=0;i<nbins;i++)
  {
    ttmp = new TTree(Form("Bin_%i",i),Form("Bin_%i",i));
    ttmp->Branch("Delta",&delta,"Delta/D");
    ttmp->Branch("Th13Axis",&t13axis,"Th13Axis/D");
    ttmp->Branch("Signal",&sig[i],"Signal/D");
    ttmp->Branch("Background",&bkgd[i],"Background/D");
    ttmp->Branch("DNExp_DOscPars",&oscparerr[i],"DNExp_DOscPars/D");
    for(k=0;k<nbins;k++)
    {
      if(k>i)
      {
        ttmp->Branch(Form("Bin_%i_Bin_%i",i,k),&oscparerr_offdiag[noff],Form("Bin_%i_Bin_%i/D",i,k));
	noff++;
      }
    }
    ftree.push_back(ttmp);

    ttmp->Reset();

    ftree2.push_back( vector<TTree*>() );

    for(j=0;j<Extrap.size();j++)
    {
      ttmp = new TTree(Form("Bin_%i_Run_%i",i,j),Form("Bin_%i_Run_%i",i,j));
      ttmp->Branch("Delta",&delta,"Delta/D");
      ttmp->Branch("Th13Axis",&t13axis,"Th13Axis/D");
      ttmp->Branch("Signal",&nue[i][j],"Signal/D");
      ttmp->Branch("NC",&nc[i][j],"NC/D");
      ttmp->Branch("NuMuCC",&numucc[i][j],"NuMuCC/D");
      ttmp->Branch("BNueCC",&bnue[i][j],"BNueCC/D");
      ttmp->Branch("NuTauCC",&tau[i][j],"NuTauCC/D");

      ftree2[i].push_back(ttmp);

      ttmp->Reset();
    }
  }

  double delta_increment = 0;
  if(nDeltaSteps>0) delta_increment = (DeltaHigh - DeltaLow)/(nDeltaSteps);
  double ssq2th13_increment = 0;
  if(nSinSq2Th13Steps>0) ssq2th13_increment = (SinSq2Th13High - SinSq2Th13Low)/(nSinSq2Th13Steps);

  int id,is,l,u;
  int ip,ir;
  double theta23,theta12,dm21,dm32,ssq2th13;
  vector<double> nexp,nobs;
  vector<TH1D*> delnexphist;
  vector<TH1D*> delnexphist_offdiag;
  for(i=0;i<nbins;i++)
  {
    delnexphist.push_back(new TH1D(Form("delnexphist_%i",i),"",400,-1,1));
    nexp.push_back(0);
    nobs.push_back(0);
    for(k=0;k<nbins;k++)
    {
      if(k>i)
      {
	delnexphist_offdiag.push_back(new TH1D(Form("delnexphist_%i_%i",i,k),"",400,-1,1));
      }
    }
  }

  gRandom->SetSeed(0);

  TFile *f = new TFile(gSystem->ExpandPathName(s.c_str()),"RECREATE");//save the osc par err distributions to this file

  l=0;
  for(id=0;id<nDeltaSteps+1;id++)
  {
    delta = (id*delta_increment + DeltaLow)*TMath::Pi();

    for(is=0;is<nSinSq2Th13Steps+1;is++)
    {
      t13axis = (is*ssq2th13_increment + SinSq2Th13Low);
      ssq2th13 = t13axis/(2.*TMath::Sin(Theta23)*TMath::Sin(Theta23));

      //get nominal prediction
      for(ie=0;ie<Extrap.size();ie++)
      {
        Extrap[ie]->SetOscPar(OscPar::kTh12,Theta12);
        Extrap[ie]->SetOscPar(OscPar::kTh23,Theta23);
        Extrap[ie]->SetOscPar(OscPar::kDeltaM12,DeltaMSq12);
        Extrap[ie]->SetOscPar(OscPar::kDeltaM23,DeltaMSq23);//note that DeltaMSq23 is always positive
        if(!NormalHier) Extrap[ie]->InvertMassHierarchy();
        Extrap[ie]->SetDeltaCP(delta);
        Extrap[ie]->SetSinSq2Th13(ssq2th13);
        Extrap[ie]->OscillatePrediction();
      }

      for(i=0;i<nbins;i++)
      {
        sig[i]=0;
        bkgd[i]=0;
	ir = int(i/nPID);
        ip = i%nPID;
        for(ie=0;ie<Extrap.size();ie++)
        {
          bkgd[i] += (Extrap[ie]->Pred_TotalBkgd_VsBinNumber->GetBinContent(i+1));
          sig[i] += (Extrap[ie]->Pred_Signal_VsBinNumber->GetBinContent(i+1));

	  nc[i][ie] = Extrap[ie]->Pred[Background::kNC]->GetBinContent(ip+1,ir+1);
          numucc[i][ie] = Extrap[ie]->Pred[Background::kNuMuCC]->GetBinContent(ip+1,ir+1);
          bnue[i][ie] = Extrap[ie]->Pred[Background::kBNueCC]->GetBinContent(ip+1,ir+1);
          tau[i][ie] = Extrap[ie]->Pred[Background::kNuTauCC]->GetBinContent(ip+1,ir+1);
          nue[i][ie] = Extrap[ie]->Pred[Background::kNueCC]->GetBinContent(ip+1,ir+1);
        }
        nexp[i] = sig[i]+bkgd[i];
      }

      //do pseudo experiments
      noff=0;
      for(i=0;i<nbins;i++)
      {
        delnexphist[i]->Reset();
        delnexphist[i]->SetName(Form("DeltaNexp_%i_%i_Diag_%i",id,is,i));
	for(k=0;k<nbins;k++)
	{
	  if(k>i)
	  {
	    delnexphist_offdiag[noff]->Reset();
	    delnexphist_offdiag[noff]->SetName(Form("DeltaNexp_%i_%i_OffDiag_%i_%i",id,is,i,k));
	    noff++;
	  }
	}
      }
      for(u=0;u<NumExpts;u++)
      {
        theta12 = Theta12 + AsymGaus(dTheta12_dn,dTheta12_up);
        dm21 = DeltaMSq12 + AsymGaus(dDeltaMSq12_dn,dDeltaMSq12_up);
        dm32 = DeltaMSq23 + AsymGaus(dDeltaMSq23_dn,dDeltaMSq23_up);
        //theta23 = DrawTheta23(dTheta23_up);
        theta23 = Theta23 + AsymGaus(dTheta23_dn,dTheta23_up);

        ssq2th13 = t13axis/(2.*TMath::Sin(theta23)*TMath::Sin(theta23));

	for(ie=0;ie<Extrap.size();ie++)
        {
          Extrap[ie]->SetOscPar(OscPar::kTh12,theta12);
          Extrap[ie]->SetOscPar(OscPar::kDeltaM12,dm21);
          Extrap[ie]->SetOscPar(OscPar::kDeltaM23,dm32);
          if(!NormalHier) Extrap[ie]->InvertMassHierarchy();
          Extrap[ie]->SetOscPar(OscPar::kTh23,theta23);
          Extrap[ie]->SetSinSq2Th13(ssq2th13);
          Extrap[ie]->OscillatePrediction();
        }

	noff=0;
	for(i=0;i<nbins;i++)
	{
	  nobs[i]=0;
	  for(ie=0;ie<Extrap.size();ie++)
          {
            nobs[i] += (Extrap[ie]->Pred_TotalBkgd_VsBinNumber->GetBinContent(i+1));
            nobs[i] += (Extrap[ie]->Pred_Signal_VsBinNumber->GetBinContent(i+1));
          }
          delnexphist[i]->Fill((nobs[i]-nexp[i])/(nexp[i]));
	}
	for(i=0;i<nbins;i++)
	{
	  for(k=0;k<nbins;k++)
	  {
	    if(k>i)
	    {
	      delnexphist_offdiag[noff]->Fill((nobs[i]-nexp[i])*(nobs[k]-nexp[k])/(nexp[i]*nexp[k]));
	      noff++;
	    }
	  }
	}
      }

      noff=0;
      for(i=0;i<nbins;i++)
      {
	oscparerr[i] = delnexphist[i]->GetRMS();
	delnexphist[i]->Write();
	for(k=0;k<nbins;k++)
	{
	  if(k>i)
	  {
	    oscparerr_offdiag[noff] = delnexphist_offdiag[noff]->GetRMS();
            if(delnexphist_offdiag[noff]->GetMean()<0) oscparerr_offdiag[noff] = -1.*oscparerr_offdiag[noff];
	    delnexphist_offdiag[noff]->Write();
	    noff++;
	  }
	}
      }
      f->Close();

      gROOT->cd("/");

      for(i=0;i<nbins;i++)
      {
	for(ie=0;ie<Extrap.size();ie++)
	{
	  ftree2[i][ie]->Fill();
	}
        ftree[i]->Fill();
      }

      f = new TFile(gSystem->ExpandPathName(s.c_str()),"UPDATE");

      if(l%100==0) cout<<100.*l/((nDeltaSteps+1)*(nSinSq2Th13Steps+1))<<"% complete"<<endl;
      l++;
    }
  }

  f->Close();

  double nPOTNear,nPOTFar;

  TTree *paramtree = new TTree("paramtree","paramtree");
  paramtree->Branch("nearPOT",&nPOTNear,"nearPOT/D");
  paramtree->Branch("farPOT",&nPOTFar,"farPOT/D");
  paramtree->Branch("Theta12",&Theta12,"Theta12/D");
  paramtree->Branch("Theta23",&Theta23,"Theta23/D");
  paramtree->Branch("DeltaMSq23",&dm32,"DeltaMSq23/D");
  paramtree->Branch("DeltaMSq12",&DeltaMSq12,"DeltaMSq12/D");

  dm32 = DeltaMSq23;
  if(!NormalHier) dm32 = -1.*DeltaMSq23;

  for(ie=0;ie<Extrap.size();ie++)
  {
    nPOTNear = Extrap[ie]->GetNearPOT();
    nPOTFar = Extrap[ie]->GetFarPOT();
    paramtree->Fill();
  }

  TFile *fout = new TFile(gSystem->ExpandPathName(outFileName.c_str()),"RECREATE");
  for(i=0;i<nbins;i++)
  {
    ftree[i]->Write();
    for(ie=0;ie<Extrap.size();ie++)
    {
      ftree2[i][ie]->Write();
    }
  }
  paramtree->Write();
  fout->Close();

  return;
}
void GridGen::RunMultiBin_VaryTheta13(string s)
{
  if(Extrap.size()==0)
  {
    cout<<"No Extrapolate2D input.  Quitting..."<<endl;
    return;
  }

  unsigned int ie;
  for(ie=0;ie<Extrap.size();ie++)
  {
    Extrap[ie]->GetPrediction();
  }

  for(ie=0;ie<Extrap.size();ie++)
  {
    Extrap[ie]->SetOscPar(OscPar::kTh13,Theta13);
    Extrap[ie]->SetOscPar(OscPar::kTh12,Theta12);
    Extrap[ie]->SetOscPar(OscPar::kTh23,Theta23);
    Extrap[ie]->SetOscPar(OscPar::kDeltaM12,DeltaMSq12);
    Extrap[ie]->SetOscPar(OscPar::kDeltaM23,DeltaMSq23);
    if(!NormalHier) Extrap[ie]->InvertMassHierarchy();
    Extrap[ie]->OscillatePrediction();
  }

  int i,k;
  unsigned int j;
  int nbins=Extrap[0]->Pred_TotalBkgd_VsBinNumber->GetNbinsX();
  int nPID = Extrap[0]->GetNPID();
  double delta,t13axis;
  vector<double> sig, bkgd;
  vector< vector<double> > nc,numucc,bnue,tau,nue;
  vector<double> oscparerr;
  vector<double> oscparerr_offdiag;
  int noff;

  t13axis = TMath::Sin(2*Theta13)*TMath::Sin(2*Theta13);

  for(i=0;i<nbins;i++)
  {
    sig.push_back(0);
    bkgd.push_back(0);
    oscparerr.push_back(0);
    nc.push_back( vector<double>() );
    numucc.push_back( vector<double>() );
    bnue.push_back( vector<double>() );
    tau.push_back( vector<double>() );
    nue.push_back( vector<double>() );
    for(k=0;k<nbins;k++)
    {
      if(k>i)
      {
        oscparerr_offdiag.push_back(0);
      }
    }
    for(j=0;j<Extrap.size();j++)
    {
      nc[i].push_back(0);
      numucc[i].push_back(0);
      bnue[i].push_back(0);
      tau[i].push_back(0);
      nue[i].push_back(0);
    }
  }

  vector<TTree*> ftree;
  vector< vector<TTree*> > ftree2;
  TTree *ttmp;
  noff=0;
  for(i=0;i<nbins;i++)
  {
    ttmp = new TTree(Form("Bin_%i",i),Form("Bin_%i",i));
    ttmp->Branch("Delta",&delta,"Delta/D");
    ttmp->Branch("Th13Axis",&t13axis,"Th13Axis/D");
    ttmp->Branch("Signal",&sig[i],"Signal/D");
    ttmp->Branch("Background",&bkgd[i],"Background/D");
    ttmp->Branch("DNExp_DOscPars",&oscparerr[i],"DNExp_DOscPars/D");
    for(k=0;k<nbins;k++)
    {
      if(k>i)
      {
        ttmp->Branch(Form("Bin_%i_Bin_%i",i,k),&oscparerr_offdiag[noff],Form("Bin_%i_Bin_%i/D",i,k));
        noff++;
      }
    }
    ftree.push_back(ttmp);

    ttmp->Reset();

    ftree2.push_back( vector<TTree*>() );

    for(j=0;j<Extrap.size();j++)
    {
      ttmp = new TTree(Form("Bin_%i_Run_%i",i,j),Form("Bin_%i_Run_%i",i,j));
      ttmp->Branch("Delta",&delta,"Delta/D");
      ttmp->Branch("Th13Axis",&t13axis,"Th13Axis/D");
      ttmp->Branch("Signal",&nue[i][j],"Signal/D");
      ttmp->Branch("NC",&nc[i][j],"NC/D");
      ttmp->Branch("NuMuCC",&numucc[i][j],"NuMuCC/D");
      ttmp->Branch("BNueCC",&bnue[i][j],"BNueCC/D");
      ttmp->Branch("NuTauCC",&tau[i][j],"NuTauCC/D");

      ftree2[i].push_back(ttmp);

      ttmp->Reset();
    }
  }

  double delta_increment = 0;
  if(nDeltaSteps>0) delta_increment = (DeltaHigh - DeltaLow)/(nDeltaSteps);

  int id,l,u;
  int ip,ir;
  double theta23,theta12,dm21,dm32,theta13;
  vector<double> nexp,nobs;
  vector<TH1D*> delnexphist;
  vector<TH1D*> delnexphist_offdiag;
  for(i=0;i<nbins;i++)
  {
    delnexphist.push_back(new TH1D(Form("delnexphist_%i",i),"",400,-1,1));
    nexp.push_back(0);
    nobs.push_back(0);
    for(k=0;k<nbins;k++)
    {
      if(k>i)
      {
        delnexphist_offdiag.push_back(new TH1D(Form("delnexphist_%i_%i",i,k),"",400,-1,1));
      }
    }
  }

  gRandom->SetSeed(0);

  TFile *f = new TFile(gSystem->ExpandPathName(s.c_str()),"RECREATE");//save the osc par err distributions to this file

  l=0;
  for(id=0;id<nDeltaSteps+1;id++)
  {
    delta = (id*delta_increment + DeltaLow)*TMath::Pi();

    //get nominal prediction
    for(ie=0;ie<Extrap.size();ie++)
    {
      Extrap[ie]->SetOscPar(OscPar::kTh13,Theta13);
      Extrap[ie]->SetOscPar(OscPar::kTh12,Theta12);
      Extrap[ie]->SetOscPar(OscPar::kTh23,Theta23);
      Extrap[ie]->SetOscPar(OscPar::kDeltaM12,DeltaMSq12);
      Extrap[ie]->SetOscPar(OscPar::kDeltaM23,DeltaMSq23);//note that DeltaMSq23 is always positive
      if(!NormalHier) Extrap[ie]->InvertMassHierarchy();
      Extrap[ie]->SetDeltaCP(delta);
      Extrap[ie]->OscillatePrediction();
    }

    for(i=0;i<nbins;i++)
    {
      sig[i]=0;
      bkgd[i]=0;
      ir = int(i/nPID);
      ip = i%nPID;
      for(ie=0;ie<Extrap.size();ie++)
      {
        bkgd[i] += (Extrap[ie]->Pred_TotalBkgd_VsBinNumber->GetBinContent(i+1));
        sig[i] += (Extrap[ie]->Pred_Signal_VsBinNumber->GetBinContent(i+1));

        nc[i][ie] = Extrap[ie]->Pred[Background::kNC]->GetBinContent(ip+1,ir+1);
        numucc[i][ie] = Extrap[ie]->Pred[Background::kNuMuCC]->GetBinContent(ip+1,ir+1);
        bnue[i][ie] = Extrap[ie]->Pred[Background::kBNueCC]->GetBinContent(ip+1,ir+1);
        tau[i][ie] = Extrap[ie]->Pred[Background::kNuTauCC]->GetBinContent(ip+1,ir+1);
        nue[i][ie] = Extrap[ie]->Pred[Background::kNueCC]->GetBinContent(ip+1,ir+1);
      }
      nexp[i] = sig[i]+bkgd[i];
    }

    //do pseudo experiments
    noff=0;
    for(i=0;i<nbins;i++)
    {
      delnexphist[i]->Reset();
      delnexphist[i]->SetName(Form("DeltaNexp_%i_Diag_%i",id,i));
      for(k=0;k<nbins;k++)
      {
        if(k>i)
        {
          delnexphist_offdiag[noff]->Reset();
          delnexphist_offdiag[noff]->SetName(Form("DeltaNexp_%i_OffDiag_%i_%i",id,i,k));
          noff++;
        }
      }
    }
    for(u=0;u<NumExpts;u++)
    {
      theta13 = Theta13 + AsymGaus(dTheta13_dn,dTheta13_up);
      theta12 = Theta12 + AsymGaus(dTheta12_dn,dTheta12_up);
      dm21 = DeltaMSq12 + AsymGaus(dDeltaMSq12_dn,dDeltaMSq12_up);
      dm32 = DeltaMSq23 + AsymGaus(dDeltaMSq23_dn,dDeltaMSq23_up);
      theta23 = Theta23 + AsymGaus(dTheta23_dn,dTheta23_up);

      for(ie=0;ie<Extrap.size();ie++)
      {
        Extrap[ie]->SetOscPar(OscPar::kTh13,theta13);
        Extrap[ie]->SetOscPar(OscPar::kTh12,theta12);
        Extrap[ie]->SetOscPar(OscPar::kDeltaM12,dm21);
        Extrap[ie]->SetOscPar(OscPar::kDeltaM23,dm32);
        if(!NormalHier) Extrap[ie]->InvertMassHierarchy();
        Extrap[ie]->SetOscPar(OscPar::kTh23,theta23);
        Extrap[ie]->OscillatePrediction();
      }

      noff=0;
      for(i=0;i<nbins;i++)
      {
        nobs[i]=0;
        for(ie=0;ie<Extrap.size();ie++)
        {
          nobs[i] += (Extrap[ie]->Pred_TotalBkgd_VsBinNumber->GetBinContent(i+1));
          nobs[i] += (Extrap[ie]->Pred_Signal_VsBinNumber->GetBinContent(i+1));
        }
        delnexphist[i]->Fill((nobs[i]-nexp[i])/(nexp[i]));
      }
      for(i=0;i<nbins;i++)
      {
        for(k=0;k<nbins;k++)
        {
          if(k>i)
          {
            delnexphist_offdiag[noff]->Fill((nobs[i]-nexp[i])*(nobs[k]-nexp[k])/(nexp[i]*nexp[k]));
            noff++;
          }
        }
      }
    }

    noff=0;
    for(i=0;i<nbins;i++)
    {
      oscparerr[i] = delnexphist[i]->GetRMS();
      delnexphist[i]->Write();
      for(k=0;k<nbins;k++)
      {
        if(k>i)
        {
          oscparerr_offdiag[noff] = delnexphist_offdiag[noff]->GetRMS();
          if(delnexphist_offdiag[noff]->GetMean()<0) oscparerr_offdiag[noff] = -1.*oscparerr_offdiag[noff];
          delnexphist_offdiag[noff]->Write();
          noff++;
        }
      }
    }
    f->Close();

    gROOT->cd("/");

    for(i=0;i<nbins;i++)
    {
      for(ie=0;ie<Extrap.size();ie++)
      {
        ftree2[i][ie]->Fill();
      }
      ftree[i]->Fill();
    }

    f = new TFile(gSystem->ExpandPathName(s.c_str()),"UPDATE");

    if(l%100==0) cout<<100.*l/(nDeltaSteps+1)<<"% complete"<<endl;
    l++;
  }

  f->Close();

  double nPOTNear,nPOTFar;

  TTree *paramtree = new TTree("paramtree","paramtree");
  paramtree->Branch("nearPOT",&nPOTNear,"nearPOT/D");
  paramtree->Branch("farPOT",&nPOTFar,"farPOT/D");
  paramtree->Branch("Theta13",&Theta13,"Theta13/D");
  paramtree->Branch("Theta12",&Theta12,"Theta12/D");
  paramtree->Branch("Theta23",&Theta23,"Theta23/D");
  paramtree->Branch("DeltaMSq23",&dm32,"DeltaMSq23/D");
  paramtree->Branch("DeltaMSq12",&DeltaMSq12,"DeltaMSq12/D");

  dm32 = DeltaMSq23;
  if(!NormalHier) dm32 = -1.*DeltaMSq23;

  for(ie=0;ie<Extrap.size();ie++)
  {
    nPOTNear = Extrap[ie]->GetNearPOT();
    nPOTFar = Extrap[ie]->GetFarPOT();
    paramtree->Fill();
  }

  TFile *fout = new TFile(gSystem->ExpandPathName(outFileName.c_str()),"RECREATE");
  for(i=0;i<nbins;i++)
  {
    ftree[i]->Write();
    for(ie=0;ie<Extrap.size();ie++)
    {
      ftree2[i][ie]->Write();
    }
  }
  paramtree->Write();
  fout->Close();

  return;
}




/*
	Dung's Testing Code Area
*/


void GridGen::SetNStepSinSqTh14SterileFit(int n) {
  NStepSinSqTh14 = n;
  return;
}

void GridGen::SetNStepSinSqTh24SterileFit(int n) {
  NStepSinSqTh24 = n;
  return;
}

void GridGen::SetSinSqTh14SterileFit(double val) {
  SinSqTh14GridValue = val;
  return;
}

void GridGen::SetSinSqTh24SterileFit(double val) {
  SinSqTh24GridValue = val;
  return;
}

void GridGen::RunMultiBinOscParErrsSterileFit(string s, double SetDM41) {
  if (Extrap.size() == 0) {
    cout << "No Extrapolate2D input.  Quitting..." << endl;
    return;
  }

  for (unsigned int ie = 0; ie < Extrap.size(); ie++) {
    Extrap[ie]->GetPrediction();
  }

  for (unsigned int ie = 0; ie < Extrap.size(); ie++) {
    Extrap[ie]->SetOscPar(OscPar::kTh12, Theta12);
    Extrap[ie]->SetOscPar(OscPar::kTh13, Theta13);
    Extrap[ie]->SetOscPar(OscPar::kTh23, Theta23);
    Extrap[ie]->SetOscPar(OscPar::kTh14, 0);
    Extrap[ie]->SetOscPar(OscPar::kTh24, 0);
    Extrap[ie]->SetOscPar(OscPar::kTh34, 0);
    Extrap[ie]->SetOscPar(OscPar::kDeltaM12, DeltaMSq12);
    Extrap[ie]->SetOscPar(OscPar::kDeltaM23, DeltaMSq23);
    Extrap[ie]->SetOscPar(OscPar::kDm41, 0);
    Extrap[ie]->SetOscPar(OscPar::kDelta, 0);
    Extrap[ie]->SetOscPar(OscPar::kDelta14, 0);
    Extrap[ie]->SetOscPar(OscPar::kDelta24, 0);
    if (!NormalHier) {
      Extrap[ie]->InvertMassHierarchy();
    }
    Extrap[ie]->OscillatePrediction();
  }

  int nbins = Extrap[0]->Pred_TotalBkgd_VsBinNumber->GetNbinsX();
  int nPID = Extrap[0]->GetNPID();
  vector<double> sig, bkgd;
  vector< vector<double> > nc, numucc, bnue, tau, nue;
  vector<double> oscparerr;
  vector<double> oscparerr_offdiag;
  int noff;

  for (unsigned int i = 0; i < nbins; i++) {
    sig.push_back(0);
    bkgd.push_back(0);
    oscparerr.push_back(0);
    nc.push_back(vector<double>() );
    numucc.push_back(vector<double>() );
    bnue.push_back(vector<double>() );
    tau.push_back(vector<double>() );
    nue.push_back(vector<double>() );
    for (unsigned int k = 0; k < nbins; k++) {
      if (k > i) {
        oscparerr_offdiag.push_back(0);
      }
    }
    for (unsigned int j = 0; j < Extrap.size(); j++) {
      nc[i].push_back(0);
      numucc[i].push_back(0);
      bnue[i].push_back(0);
      tau[i].push_back(0);
      nue[i].push_back(0);
    }
  }

  double dmsq41;
  double ssq2th14;
  double ssq2th24;
  double delta24;
  double delta14;
  double delta13;
  double th34;
  vector<TTree *> ftree;
  vector< vector<TTree *> > ftree2;
  TTree * ttmp;
  noff = 0;
  for (unsigned int i = 0; i < nbins; i++) {
    ttmp = new TTree(Form("Bin_%i", i), Form("Bin_%i", i));
    ttmp->Branch("Dmsq41", &dmsq41, "Dmsq41/D");
    ttmp->Branch("SinSq2Th14", &ssq2th14, "SinSq2Th14/D");
    ttmp->Branch("SinSq2Th24", &ssq2th24, "SinSq2Th24/D");
    ttmp->Branch("Delta13", &delta13, "Delta13/D");
    ttmp->Branch("Delta14", &delta14, "Delta14/D");
    ttmp->Branch("Delta24", &delta24, "Delta24/D");
    ttmp->Branch("Theta34", &th34, "Theta34/D");
    ttmp->Branch("Signal", &sig[i], "Signal/D");
    ttmp->Branch("Background", &bkgd[i], "Background/D");
    ttmp->Branch("DNExp_DOscPars", &oscparerr[i], "DNExp_DOscPars/D");
    for (unsigned int k = 0; k < nbins; k++) {
      if (k > i) {
        ttmp->Branch(Form("Bin_%i_Bin_%i", i, k), &oscparerr_offdiag[noff], Form("Bin_%i_Bin_%i/D", i, k));
        noff++;
      }
    }
    ftree.push_back(ttmp);
    ttmp->Reset();
    ftree2.push_back(vector<TTree *>() );
    for (unsigned int j = 0; j < Extrap.size(); j++) {
      ttmp = new TTree(Form("Bin_%i_Run_%i", i, j), Form("Bin_%i_Run_%i", i, j));
      ttmp->Branch("Dmsq41", &dmsq41, "Dmsq41/D");
      ttmp->Branch("SinSq2Th14", &ssq2th14, "SinSq2Th14/D");
      ttmp->Branch("SinSq2Th24", &ssq2th24, "SinSq2Th24/D");
      ttmp->Branch("Delta13", &delta13, "Delta13/D");
      ttmp->Branch("Delta14", &delta14, "Delta14/D");
      ttmp->Branch("Delta24", &delta24, "Delta14/D");
      ttmp->Branch("Theta34", &th34, "Theta34/D");
      ttmp->Branch("Signal", &nue[i][j], "Signal/D");
      ttmp->Branch("NC", &nc[i][j], "NC/D");
      ttmp->Branch("NuMuCC", &numucc[i][j], "NuMuCC/D");
      ttmp->Branch("BNueCC", &bnue[i][j], "BNueCC/D");
      ttmp->Branch("NuTauCC", &tau[i][j], "NuTauCC/D");
      ftree2[i].push_back(ttmp);
      ttmp->Reset();
    }
  }

  vector<double> nexp, nobs;
  vector<TH1D *> delnexphist;
  vector<TH1D *> delnexphist_offdiag;
  for (unsigned int i = 0; i < nbins; i++) {
    delnexphist.push_back(new TH1D(Form("delnexphist_%i", i), "", 400, -1, 1));
    nexp.push_back(0);
    nobs.push_back(0);
    for (unsigned int k = 0; k < nbins; k++) {
      if (k > i) {
        delnexphist_offdiag.push_back(new TH1D(Form("delnexphist_%i_%i", i, k), "", 400, -1, 1));
      }
    }
  }

  gRandom->SetSeed(0);

  TFile * f = new TFile(gSystem->ExpandPathName(s.c_str()), "RECREATE");

  dmsq41 = SetDM41;
  ssq2th14 = SinSqTh14GridValue;
  ssq2th24 = SinSqTh24GridValue;
  for (unsigned int ie = 0; ie < Extrap.size(); ie++) {
    Extrap[ie]->SetOscPar(OscPar::kDm41, dmsq41);
    Extrap[ie]->SetOscPar(OscPar::kTh14, TMath::ASin(TMath::Sqrt(ssq2th14)));
    Extrap[ie]->SetOscPar(OscPar::kTh24, TMath::ASin(TMath::Sqrt(ssq2th24)));
  }
  for (int iTh34 = 0; iTh34 < 3; iTh34++) {
    th34 = (double)iTh34 * (TMath::Pi() / 4);
    for (unsigned int ie = 0; ie < Extrap.size(); ie++) {
      Extrap[ie]->SetOscPar(OscPar::kTh34, th34);
    }
    for (int iDelCP = 0; iDelCP < 3; iDelCP++) {
      delta13 = (double)iDelCP * (2 * TMath::Pi() / 3);
      for (unsigned int ie = 0; ie < Extrap.size(); ie++) {
        Extrap[ie]->SetOscPar(OscPar::kDelta, delta13);
      }
      for (int iDelEff = 0; iDelEff < 8; iDelEff++) {
        delta14 = (double)(((double)(rand() % 100)) / 100) * 2 * TMath::Pi();
        delta24 = (double)iDelEff * (TMath::Pi() / 4) + delta14;
        for (unsigned int ie = 0; ie < Extrap.size(); ie++) { // Get nominal prediction
          Extrap[ie]->SetOscPar(OscPar::kDelta14, delta14);
          Extrap[ie]->SetOscPar(OscPar::kDelta24, delta24);
          Extrap[ie]->SetOscPar(OscPar::kTh12, Theta12);
          Extrap[ie]->SetOscPar(OscPar::kTh13, Theta13);
          Extrap[ie]->SetOscPar(OscPar::kTh23, Theta23);
          Extrap[ie]->SetOscPar(OscPar::kDeltaM12, DeltaMSq12);
          Extrap[ie]->SetOscPar(OscPar::kDeltaM23, DeltaMSq23);
          //if (!NormalHier) { Extrap[ie]->InvertMassHierarchy(); }
          Extrap[ie]->OscillatePrediction();
        }

        for (unsigned int i = 0; i < nbins; i++) {
          sig[i] = 0;
          bkgd[i] = 0;
          int ir = int(i / nPID);
          int ip = i % nPID;
          for (unsigned int ie = 0; ie < Extrap.size(); ie++) {
            bkgd[i] += (Extrap[ie]->Pred_TotalBkgd_VsBinNumber->GetBinContent(i + 1));
            sig[i] += (Extrap[ie]->Pred_Signal_VsBinNumber->GetBinContent(i + 1));
            nc[i][ie] = Extrap[ie]->Pred[Background::kNC]->GetBinContent(ip + 1, ir + 1);
            numucc[i][ie] = Extrap[ie]->Pred[Background::kNuMuCC]->GetBinContent(ip + 1, ir + 1);
            bnue[i][ie] = Extrap[ie]->Pred[Background::kBNueCC]->GetBinContent(ip + 1, ir + 1);
            tau[i][ie] = Extrap[ie]->Pred[Background::kNuTauCC]->GetBinContent(ip + 1, ir + 1);
            nue[i][ie] = Extrap[ie]->Pred[Background::kNueCC]->GetBinContent(ip + 1, ir + 1);
          }
          nexp[i] = sig[i] + bkgd[i];
        }

        // do pseudo experiments
        noff = 0;
        for (unsigned int i = 0; i < nbins; i++) {
          delnexphist[i]->Reset();
          delnexphist[i]->SetName(Form("DeltaNexp_%i_%i_%i_Diag_%i", iTh34, iDelCP, iDelEff, i));
          for (unsigned int k = 0; k < nbins; k++) {
            if (k > i) {
              delnexphist_offdiag[noff]->Reset();
              delnexphist_offdiag[noff]->SetName(Form("DeltaNexp_%i_%i_%i_OffDiag_%i_%i", iTh34, iDelCP, iDelEff, i, k));
              noff++;
            }
          }
        }

        for (unsigned int u = 0; u < NumExpts; u++) {
          double theta12 = Theta12 + AsymGaus(dTheta12_dn, dTheta12_up);
          double theta13 = Theta13 + AsymGaus(dTheta13_dn, dTheta13_up);
          double theta23 = Theta23 + AsymGaus(dTheta23_dn, dTheta23_up);
          double dm21 = DeltaMSq12 + AsymGaus(dDeltaMSq12_dn, dDeltaMSq12_up);
          double dm32 = DeltaMSq23 + AsymGaus(dDeltaMSq23_dn, dDeltaMSq23_up);

          for (unsigned int ie = 0; ie < Extrap.size(); ie++) {
            Extrap[ie]->SetOscPar(OscPar::kTh12, theta12);
            Extrap[ie]->SetOscPar(OscPar::kTh13, theta13);
            Extrap[ie]->SetOscPar(OscPar::kTh23, theta23);
            Extrap[ie]->SetOscPar(OscPar::kDeltaM12, dm21);
            Extrap[ie]->SetOscPar(OscPar::kDeltaM23, dm32);
            //if (!NormalHier) { Extrap[ie]->InvertMassHierarchy(); }
            Extrap[ie]->OscillatePrediction();
          }

          noff = 0;
          for (unsigned int i = 0; i < nbins; i++) {
            nobs[i] = 0;
            for (unsigned int ie = 0; ie < Extrap.size(); ie++) {
              nobs[i] += (Extrap[ie]->Pred_TotalBkgd_VsBinNumber->GetBinContent(i + 1));
              nobs[i] += (Extrap[ie]->Pred_Signal_VsBinNumber->GetBinContent(i + 1));
            }
            delnexphist[i]->Fill((nobs[i] - nexp[i]) / (nexp[i]));
          }
          for (unsigned int i = 0; i < nbins; i++) {
            for (unsigned int k = 0; k < nbins; k++) {
              if (k > i) {
                delnexphist_offdiag[noff]->Fill((nobs[i] - nexp[i]) * (nobs[k] - nexp[k]) / (nexp[i] * nexp[k]));
                noff++;
              }
            }
          }
        }

        noff = 0;
        for (unsigned int i = 0; i < nbins; i++) {
          oscparerr[i] = delnexphist[i]->GetRMS();
          delnexphist[i]->Write();
          for (unsigned int k = 0; k < nbins; k++) {
            if (k > i) {
              oscparerr_offdiag[noff] = delnexphist_offdiag[noff]->GetRMS();
              if (delnexphist_offdiag[noff]->GetMean() < 0) {
                oscparerr_offdiag[noff] = -1. * oscparerr_offdiag[noff];
              }
              delnexphist_offdiag[noff]->Write();
              noff++;
            }
          }
        }
        f->Close();

        gROOT->cd("/");

        for (unsigned int i = 0; i < nbins; i++) {
          for (unsigned int ie = 0; ie < Extrap.size(); ie++) {
            ftree2[i][ie]->Fill();
          }
          ftree[i]->Fill();
        }

        f = new TFile(gSystem->ExpandPathName(s.c_str()), "UPDATE");
      }
    }
  }

  f->Close();

  double nPOTNear, nPOTFar;

  TTree * paramtree = new TTree("paramtree", "paramtree");
  paramtree->Branch("nearPOT", &nPOTNear, "nearPOT/D");
  paramtree->Branch("farPOT", &nPOTFar, "farPOT/D");
  paramtree->Branch("Theta12", &Theta12, "Theta12/D");
  paramtree->Branch("Theta13", &Theta13, "Theta13/D");
  paramtree->Branch("Theta23", &Theta23, "Theta23/D");
  paramtree->Branch("DeltaMSq41", &dmsq41, "DeltaMSq41/D");
  paramtree->Branch("DeltaMSq23", &DeltaMSq23, "DeltaMSq23/D");
  paramtree->Branch("DeltaMSq12", &DeltaMSq12, "DeltaMSq12/D");
  paramtree->Branch("SinSq2Th14", &ssq2th14, "SinSq2Th14/D");
  paramtree->Branch("SinSq2Th24", &ssq2th24, "SinSq2Th24/D");

  if (!NormalHier) {
    DeltaMSq23 = -1. * DeltaMSq23;
  }

  for (unsigned int ie = 0; ie < Extrap.size(); ie++) {
    nPOTNear = Extrap[ie]->GetNearPOT();
    nPOTFar = Extrap[ie]->GetFarPOT();
    paramtree->Fill();
  }

  TFile * fout = new TFile(gSystem->ExpandPathName(outFileName.c_str()), "RECREATE");
  for (unsigned int i = 0; i < nbins; i++) {
    ftree[i]->Write();
    for (unsigned int ie = 0; ie < Extrap.size(); ie++) {
      ftree2[i][ie]->Write();
    }
  }
  paramtree->Write();
  fout->Close();

  return;
}
