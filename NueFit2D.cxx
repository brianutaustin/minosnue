#define NueFit2D_C

#include "NueAna/MultiBinAna/NueFit2D.h"
#include <fstream>

// ////////////RBT Start edits///////////////
static void dlnlFunction(int &npar, double * /*gin*/, double &result, double * par, int /*iflag*/) {

  std::vector<double> f;
  NueFit2D * nf2d = (NueFit2D *) gMinuit->GetObjectFit();

  for (Int_t i = 0; i < npar; i++) {
    f.push_back(par[i]);
  }
  result = nf2d->StdLikeComparison(f);

}

static void globalGridFunction(int &, double * /*gin*/, double &likelihoodValue, double * gridPar, int /*iflag*/) {

  std::vector<double> f;
  NueFit2D * nf2d = (NueFit2D *) gMinuit->GetObjectFit();

  for (Int_t i = 0; i < 6; i++) {
    f.push_back(gridPar[i]);
  }
  likelihoodValue = nf2d->StdLikeComparisonForGlobalMinSterileFit(f);

}

double NueFit2D::StdLikeComparisonForGlobalMinSterileFit(vector<double> npar) {
  if (NObs == 0) {
    cout << "NObs not set.  Quitting..." << endl;
    return 0;
  }
  if (FitMethod == 1 && (FracErr_Bkgd == 0 || FracErr_Sig == 0)) {
    cout << "FracErr_Bkgd and FracErr_Sig need to be set for ScaledChi2.  Quitting..." << endl;
    return 0;
  }
  if (Extrap.size() == 0) {
    cout << "No Extrapolate2D input.  Quitting..." << endl;
    return 0;
  }

  if (FitMethod == 3) {
    DefineStdDlnLMinuit();
  }
  if (FitMethod == 4) {
    DefineBinDlnLMinuit();
  }
  if (FitMethod == 3 || FitMethod == 4) {
    if (ErrCalc != 0) {
      ErrCalc->SetUseGrid(false);
    }
  }

  Bkgd = (TH1D *) NObs->Clone("Bkgd");
  Bkgd->Reset();
  Sig = (TH1D *) NObs->Clone("Sig");
  Sig->Reset();
  TH1D * NExp = (TH1D *) NObs->Clone("NExp");
  NExp->Reset();

  for (unsigned int ie = 0; ie < Extrap.size(); ie++) {
    Extrap[ie]->GetPrediction();
    Extrap[ie]->SetOscPar(OscPar::kTh14, TMath::ASin(TMath::Sqrt(npar.at(0))));
    Extrap[ie]->SetOscPar(OscPar::kTh24, TMath::ASin(TMath::Sqrt(npar.at(1))));
    Extrap[ie]->SetOscPar(OscPar::kTh34, npar.at(2));
    Extrap[ie]->SetOscPar(OscPar::kDelta, npar.at(3));
    Extrap[ie]->SetOscPar(OscPar::kDelta14, npar.at(4));
    Extrap[ie]->SetOscPar(OscPar::kDelta24, npar.at(5) + npar.at(4));
    Extrap[ie]->SetOscPar(OscPar::kTh12, grid_n_th12);
    Extrap[ie]->SetOscPar(OscPar::kTh13, grid_n_th13);
    Extrap[ie]->SetOscPar(OscPar::kTh23, grid_n_th23);
    Extrap[ie]->SetOscPar(OscPar::kDeltaM12, grid_n_dm2_21);
    Extrap[ie]->SetOscPar(OscPar::kDeltaM23, grid_n_dm2_32);
    Extrap[ie]->SetOscPar(OscPar::kDm41, grid_dmsq41);
    Extrap[ie]->OscillatePrediction();
  }

  if (ErrorMatrix == 0) {
    ErrorMatrix = new TH2D("ErrorMatrix", "", nBins, -0.5, nBins - 0.5, nBins, -0.5, nBins - 0.5);
  }
  if (InvErrorMatrix == 0) {
    InvErrorMatrix = new TH2D("InvErrorMatrix", "", nBins, -0.5, nBins - 0.5, nBins, -0.5, nBins - 0.5);
  }

  Bkgd->Reset();
  Sig->Reset();
  NExp->Reset();
  for (unsigned int ie = 0; ie < Extrap.size(); ie++) {
    Bkgd->Add(Extrap[ie]->Pred_TotalBkgd_VsBinNumber);
    Sig->Add(Extrap[ie]->Pred_Signal_VsBinNumber);
  }
  NExp->Add(Bkgd);
  NExp->Add(Sig);

  double c2;
  vector<double> emptyPar;
  switch (FitMethod) {
    case 0: {
      c2 = PoissonChi2(NExp);
      break;
    }
    case 1: {
      c2 = ScaledChi2(Bkgd, Sig);
      break;
    }
    case 2: {
      c2 = StandardChi2(NExp);
      break;
    }
    case 3: {
      c2 = StandardLikelihood();
      break;
    }
    case 4: {
      c2 = BinLikelihood();
      break;
    }
    default: {
      cout << "Error in GetMinLikelihoodSterileFit(): Unknown 'FitMethod'. BinLikelihood() is used." << endl;
      c2 = BinLikeComparison(emptyPar);
      break;
    }
  }

  return c2;
}

double NueFit2D::DoGlobalMinSearchSterileFit(vector<double> SinStart) {
  TMinuit *ssgMinuit = new TMinuit(6);
  gMinuit = ssgMinuit;
  ssgMinuit->SetObjectFit(this);
  ssgMinuit->SetFCN(globalGridFunction);
  ssgMinuit->SetPrintLevel(1);

  Double_t * startValue = new Double_t[6];
  Double_t * stepValue = new Double_t[6];
  Double_t * Bmin = new Double_t[6];
  Double_t * Bmax = new Double_t[6];
  Int_t iErflg = 0;
  double arglist[1];
  arglist[0] = 1.e-2;

  // Set the parameters
  for (Int_t i = 0; i < 2; i++) { // SinSqTh14, SinSqTh24
    startValue[i] = SinStart.at(i);
    stepValue[i] = 0.01;
    Bmin[i] = 0.0;
    Bmax[i] = 1.0;
    ssgMinuit->mnparm(i, Form("Sin%i", i), startValue[i], stepValue[i], Bmin[i], Bmax[i], iErflg);
  }
  // Theta34
  startValue[2] = 0.0;
  stepValue[2] = TMath::Pi() / 4;
  Bmin[2] = 0.0;
  Bmax[2] = TMath::Pi() / 2;
  ssgMinuit->mnparm(2, Form("Sin%i", 2), startValue[2], stepValue[2], Bmin[2], Bmax[2], iErflg);

  // Delta13
  startValue[3] = 0.0;
  stepValue[3] = 2 * TMath::Pi() / 3;
  Bmin[3] = 0.0;
  Bmax[3] = 4 * TMath::Pi() / 3;
  ssgMinuit->mnparm(3, Form("Sin%i", 3), startValue[3], stepValue[3], Bmin[3], Bmax[3], iErflg);

  // DeltaEff
  startValue[4] = 0.0;
  stepValue[4] = 0.01;
  Bmin[4] = 0.0;
  Bmax[4] = 0.0;
  ssgMinuit->mnparm(4, Form("Sin%i", 4), startValue[4], stepValue[4], Bmin[4], Bmax[4], iErflg);
  startValue[5] = 0.0;
  stepValue[5] = TMath::Pi() / 4;
  Bmin[5] = 0.0;
  Bmax[5] = 7 * TMath::Pi() / 4;
  ssgMinuit->mnparm(5, Form("Sin%i", 5), startValue[5], stepValue[5], Bmin[5], Bmax[5], iErflg);

  ssgMinuit->mnexcm("SET EPS", arglist, 1, iErflg);

  // Resets function value and errors to UNDEFINED:
  ssgMinuit->mnrset(1);

  // 1 std dev for dlnl:
  ssgMinuit->SetErrorDef(1.0);

  // Max iterations:
  ssgMinuit->SetMaxIterations(500);
  // Go minimize!
  cout << "Start minimizing..." << endl;
  ssgMinuit->Migrad();
  cout << "Done minimizing!" << endl;
  // Get the minimum for the function
  double minpoint = ssgMinuit->fAmin;

  return minpoint;
}

static void binFunction(int &npar, double * /*gin*/, double &result, double * par, int /*iflag*/) {

  std::vector<double> f;
  NueFit2D * nf2d = (NueFit2D *) gMinuit->GetObjectFit();

  for (Int_t i = 0; i < npar; i++) {
    f.push_back(par[i]);
  }
  result = nf2d->BinLikeComparison(f);

}
// //////////End edits///////////////////////

double Poisson(double mean, int n, bool numOnly) {
  double numerator = 1;

  if (mean > 100 || n > 100) {
    // use roots protected function -
    if (!numOnly) {
      return TMath::Poisson(n, mean);
    }

    // else     pull a trick to prevent explosions
    double logN = n * TMath::Log(mean) - mean;
    numerator = TMath::Exp(logN);
  } else {
    numerator = TMath::Power(mean, n) * TMath::Exp(-mean);
  }
  double denom = 1.0;
  if (!numOnly) {
    denom = TMath::Factorial(n);
  }
  return numerator / denom;
}

double Gaussian(double x, double mean, double sig) {
  double front = 1 / TMath::Sqrt(2 * 3.1415926 * sig);
  double exp = (x - mean) * (x - mean) / (2 * sig);

  return front * TMath::Exp(-exp);
}

NueFit2D::NueFit2D() {

  NObs = 0;
  ErrCalc = 0;
  ErrorMatrix = 0;
  InvErrorMatrix = 0;
  ExternalErrorMatrix = 0;

  SetOutputFile();
  SetNExperiments();

  SetGridNorm();
  GridScale_Normal = 1.;
  GridScale_Inverted = 1.;

  nBins = 0;

  IncludeOscParErrs = false;

  FormatChi2Hists();

  SetNDeltaSteps();
  SetNSinSq2Th13Steps();
  SetDeltaRange();
  SetSinSq2Th13Range();
  SetNepsetauSteps();
  SetepsetauRange();

  SetNSinSq2Th14Steps();
  SetSinSq2Th14Range();

  SetNSinSqTh14Steps();
  SetNSinSqTh24Steps();
  SetSinSqTh14Range();
  SetSinSqTh24Range();

  SetFitMethod();

  return;
}
NueFit2D::~NueFit2D() {
}
void NueFit2D::SetNObs(TH1D * n) {
  NObs = (TH1D *) n->Clone("NObs");
  nBins = NObs->GetNbinsX();
  return;
}
void NueFit2D::SetFracBkgdError(TH1D * n) {
  FracErr_Bkgd = (TH1D *) n->Clone("FracErr_Bkgd");
  nBins = FracErr_Bkgd->GetNbinsX();
  return;
}
void NueFit2D::SetFracSigError(TH1D * n) {
  FracErr_Sig = (TH1D *) n->Clone("FracErr_Sig");
  nBins = FracErr_Sig->GetNbinsX();
  return;
}

void NueFit2D::SetSystErrorMatrix(TH2D * n) {
  ExternalErrorMatrix = (TH2D *) n->Clone("ExternalErrorMatrix");
  nBins = ExternalErrorMatrix->GetNbinsX();

  return;
}
void NueFit2D::SetOutputFile(string s) {
  outFileName = s;
  return;
}
void NueFit2D::SetGridFiles(string snorm, string sinvt) {
  GridFileName_Normal = snorm;
  GridFileName_Inverted = sinvt;
  return;
}
void NueFit2D::AddExtrap(Extrapolate2D * E) {
  Extrap.push_back(E);
  return;
}

void NueFit2D::ReadGridFiles() {
  double fp;
  unsigned int i;
  TTree * temp = 0;
  unsigned int j;

  for (i = 0; i < nBins; i++) {
    grid_bin_oscparerr.push_back(0);
  }

  TFile * fnorm;
  GridTree_Normal.clear();
  GridTree_2_Normal.clear();
  nPts_Normal = 0;
  if (!gSystem->AccessPathName(gSystem->ExpandPathName(GridFileName_Normal.c_str()))) { // if file exists
    fnorm = new TFile(gSystem->ExpandPathName(GridFileName_Normal.c_str()), "READ");
    for (i = 0; i < nBins; i++) {
      temp = (TTree *) fnorm->Get(Form("Bin_%i", i));
      temp->SetName(Form("Bin_%i_Normal", i));
      temp->SetBranchAddress("Background", &grid_background);
      temp->SetBranchAddress("Signal", &grid_signal);
      temp->SetBranchAddress("Delta", &grid_delta);
      temp->SetBranchAddress("Th13Axis", &grid_sinsq2th13);
      if (IncludeOscParErrs) {
        temp->SetBranchAddress("DNExp_DOscPars", &grid_bin_oscparerr[i]);
      }
      GridTree_Normal.push_back(temp);

      GridTree_2_Normal.push_back(vector<TTree *>() );

      for (j = 0; j < Extrap.size(); j++) {
        temp = (TTree *) fnorm->Get(Form("Bin_%i_Run_%i", i, j));
        temp->SetName(Form("Bin_%i_Run_%i_Normal", i, j));
        temp->SetBranchAddress("NC", &grid_nc);
        temp->SetBranchAddress("NuMuCC", &grid_numucc);
        temp->SetBranchAddress("BNueCC", &grid_bnuecc);
        temp->SetBranchAddress("NuTauCC", &grid_nutaucc);
        temp->SetBranchAddress("Signal", &grid_nue);
        temp->SetBranchAddress("Delta", &grid_delta);
        temp->SetBranchAddress("Th13Axis", &grid_sinsq2th13);
        GridTree_2_Normal[i].push_back(temp);
      }
    }
    nPts_Normal = GridTree_Normal[0]->GetEntries();

    paramtree_Normal = (TTree *) fnorm->Get("paramtree");
    paramtree_Normal->SetName("paramtree_Normal");
    paramtree_Normal->SetBranchAddress("farPOT", &fp);
    paramtree_Normal->SetBranchAddress("Theta12", &grid_n_th12);
    paramtree_Normal->SetBranchAddress("Theta23", &grid_n_th23);
    paramtree_Normal->SetBranchAddress("DeltaMSq23", &grid_n_dm2_32);
    paramtree_Normal->SetBranchAddress("DeltaMSq12", &grid_n_dm2_21);
    paramtree_Normal->SetBranchAddress("Theta13", &grid_n_th13);
    paramtree_Normal->GetEntry(0);

    if (GridNorm > 0) {
      GridScale_Normal = GridNorm / fp;
    }
  } else {
    cout << "Grid file (normal hierarchy) doesn't exist." << endl;
    return;
  }

  TFile * finvt;
  GridTree_Inverted.clear();
  nPts_Inverted = 0;
  if (!gSystem->AccessPathName(gSystem->ExpandPathName(GridFileName_Inverted.c_str()))) { // if file exists
    finvt = new TFile(gSystem->ExpandPathName(GridFileName_Inverted.c_str()), "READ");
    for (i = 0; i < nBins; i++) {
      temp = (TTree *) finvt->Get(Form("Bin_%i", i));
      temp->SetName(Form("Bin_%i_Inverted", i));
      temp->SetBranchAddress("Background", &grid_background);
      temp->SetBranchAddress("Signal", &grid_signal);
      temp->SetBranchAddress("Delta", &grid_delta);
      temp->SetBranchAddress("Th13Axis", &grid_sinsq2th13);
      if (IncludeOscParErrs) {
        temp->SetBranchAddress("DNExp_DOscPars", &grid_bin_oscparerr[i]);
      }
      GridTree_Inverted.push_back(temp);

      GridTree_2_Inverted.push_back(vector<TTree *>() );

      for (j = 0; j < Extrap.size(); j++) {
        temp = (TTree *) finvt->Get(Form("Bin_%i_Run_%i", i, j));
        temp->SetName(Form("Bin_%i_Run_%i_Inverted", i, j));
        temp->SetBranchAddress("NC", &grid_nc);
        temp->SetBranchAddress("NuMuCC", &grid_numucc);
        temp->SetBranchAddress("BNueCC", &grid_bnuecc);
        temp->SetBranchAddress("NuTauCC", &grid_nutaucc);
        temp->SetBranchAddress("Signal", &grid_nue);
        temp->SetBranchAddress("Delta", &grid_delta);
        temp->SetBranchAddress("Th13Axis", &grid_sinsq2th13);
        GridTree_2_Inverted[i].push_back(temp);
      }
    }
    nPts_Inverted = GridTree_Inverted[0]->GetEntries();

    paramtree_Inverted = (TTree *) finvt->Get("paramtree");
    paramtree_Inverted->SetName("paramtree_Inverted");
    paramtree_Inverted->SetBranchAddress("farPOT", &fp);
    paramtree_Inverted->SetBranchAddress("Theta12", &grid_i_th12);
    paramtree_Inverted->SetBranchAddress("Theta23", &grid_i_th23);
    paramtree_Inverted->SetBranchAddress("DeltaMSq23", &grid_i_dm2_32);
    paramtree_Inverted->SetBranchAddress("DeltaMSq12", &grid_i_dm2_21);
    paramtree_Inverted->SetBranchAddress("Theta13", &grid_i_th13);
    paramtree_Inverted->GetEntry(0);

    if (GridNorm > 0) {
      GridScale_Inverted = GridNorm / fp;
    }
  } else {
    cout << "Grid file (inverted hierarchy) doesn't exist." << endl;
    return;
  }

  return;
}


void NueFit2D::DefineBinDlnLMinuit() {
  // Function sets the maximum size of Minuit.

  // Clear things first:

  int npar = 0;

  if (nBins != 0) {
    npar = nBins;
  } else if (NObs != 0) {
    npar = NObs->GetNbinsX();
  } else if (Extrap.size() != 0) {
    npar = Extrap[0]->Pred_TotalBkgd_VsBinNumber->GetNbinsX();
  } else {
    cout << "ERROR: Add extrapolation or NObs before initializing Minuit size"
         << endl;
    return;
  }

  // Make new minuit:
  minuit = new TMinuit(npar);
  minuit->SetPrintLevel(-1);

  double arglist[1];
  int ierflg = 0;

  arglist[0] = 1.e-3;
  minuit->mnexcm("SET EPS", arglist, 1, ierflg);


}

void NueFit2D::DefineStdDlnLMinuit() {
  // Function sets the maximum size of Minuit.

  // Clear things first:
  // minuit->mncler();

  int npar = 0;

  // Number of systematics inputted:
  if (FracErr_Bkgd_List.size() != 0) {
    npar = FracErr_Bkgd_List.size();
  }

  int nb = 0;
  // Number of bins (for HOOHE):
  if (nBins != 0) {
    nb = nBins;
  } else if (NObs != 0) {
    nb = NObs->GetNbinsX();
  } else if (Extrap.size() != 0) {
    nb = Extrap[0]->Pred_TotalBkgd_VsBinNumber->GetNbinsX();
  } else {
    cout << "ERROR: Add extrapolation or NObs before initializing Minuit size"
         << endl;
    return;
  }

  npar += nb;

  // make new minuit
  minuit = new TMinuit(npar + nb);
  minuit->SetPrintLevel(-1);

}

void NueFit2D::AddError(ErrorCalc * Err) {
  // References the ErrorCalc object
  ErrCalc = Err;
  return;
}

void NueFit2D::CalculateDlnLMatrix(TH2D * SysMatrix = 0, TH2D * HOOHEMatrix = 0) {
  // Calculate the covariance matrix for the "bin by bin" method

  // Make new inverted matrix if not done already
  if (InvErrorMatrix == 0) {
    if (SysMatrix) {
      InvErrorMatrix = (TH2D *) SysMatrix->Clone("InvErrorMatrix");
    } else if (HOOHEMatrix) {
      InvErrorMatrix = (TH2D *) HOOHEMatrix->Clone("InvErrorMatrix");
    } else {
      cout << "Didn't give me any matrices to invert!" << endl;
      return;
    }
  }
  // Reset everything to zero:
  InvErrorMatrix->Reset();
  int totbins = nBins;
  int i = 0;
  int j = 0;
  int k = 0;

  // Add together error elements into a single array
  double * Varray = new double[totbins * totbins];
  for (i = 0; i < totbins; i++) {
    for (j = 0; j < totbins; j++) {
      k = i * totbins + j;
      Varray[k] = 0;
      if (SysMatrix != 0) {
        Varray[k] += SysMatrix->GetBinContent(i + 1, j + 1);
      }

      if (HOOHEMatrix != 0) {
        Varray[k] += HOOHEMatrix->GetBinContent(i + 1, j + 1);
      }
    }
  }

  // Hand elements to a TMatrix
  TMatrixD * Vmatrix = new TMatrixD(totbins, totbins, Varray);

  // Insert it (determ found for debugging purposes)!
  double determ;
  TMatrixD Vinvmatrix = Vmatrix->Invert(&determ);
  // cout << "DETERM=" << determ << endl;

  // Extract the array of the inverted matrix:
  Double_t * Vinvarray = Vinvmatrix.GetMatrixArray();

  // Turn it into the inverted covariance matrix

  // make Vinv out of array
  for (i = 0; i < totbins; i++) {
    for (j = 0; j < totbins; j++) {
      InvErrorMatrix->SetBinContent(i + 1, j + 1, Vinvarray[i * totbins + j]);
    }
  }

  delete [] Varray;

  return;

}

void NueFit2D::SetFracError(vector<TH1D *> bkglist, vector<TH1D *> siglist) {
  // Set the vectors of fractional errors for the "Standard" likelihood

  if (bkglist.size() == 0) {
    cout << "Forgot to add a Background systematic list!" << endl;
    cout << "->Nothing will be set." << endl;
    return;
  }

  if (siglist.size() == 0) {
    cout << "Forgot to add a Signal systematic list!" << endl;
    cout << "->Nothing will be set." << endl;
    return;
  }

  if (bkglist.size() != siglist.size()) {
    cout << "Signal and Background error lists are different sizes." << endl;
    cout << "->Nothing will be set." << endl;
    return;
  }

  // Reset everything:
  FracErr_Bkgd_List.clear();
  FracErr_Sig_List.clear();

  for (unsigned int i = 0; i < bkglist.size(); i++) {

    FracErr_Bkgd_List.push_back((TH1D *) (bkglist[i]->Clone(Form("FracErr_Bkgd_%i", i))));
    FracErr_Sig_List.push_back((TH1D *) (siglist[i]->Clone(Form("FracErr_Sig_%i", i))));

  }

}

double NueFit2D::DoStdMinParam() {

  // Size of systematics.  The +nBins will be for HOOHE.
  int sys_size = FracErr_Bkgd_List.size() + nBins;

  if (sys_size > minuit->fMaxpar) {
    cout << "WARNING: WRONG MINUIT SIZE" << endl;
  }

  Double_t * vstrt = new Double_t[sys_size];
  Double_t * stp = new Double_t[sys_size];
  Double_t * bmin = new Double_t[sys_size];
  Double_t * bmax = new Double_t[sys_size];
  Int_t ierflg = 0;

  // Pass on the object and set the static function:
  gMinuit = minuit;
  minuit->SetObjectFit(this);
  minuit->SetFCN(dlnlFunction);

  // Set the parameters
  for (Int_t i = 0; i < sys_size; i++) {
    vstrt[i] = 0.0;
    stp[i] = 1.0;
    bmin[i] = 0.0;
    bmax[i] = 0.0;

    // HOO disabled for now:
    // if (InvErrorMatrix==0&&(unsigned int)i>=FracErr_Bkgd_List.size()) vstrt[i] = 1.0;
    if ((unsigned int) i >= FracErr_Bkgd_List.size()) {
      vstrt[i] = 1.0;
    }

    minuit->mnparm(i, Form("f%i", i), vstrt[i], stp[i], bmin[i], bmax[i], ierflg);

    // HOO disabled for now:
    // if (InvErrorMatrix==0&&(unsigned int)i>=FracErr_Bkgd_List.size()) minuit->FixParameter(i);
    if ((unsigned int) i >= FracErr_Bkgd_List.size()) {
      minuit->FixParameter(i);
    }
  }

  // Resets function value and errors to UNDEFINED:
  minuit->mnrset(1);

  // 2*dlnl style error definition
  minuit->SetErrorDef(1.0);

  // Max iterations:
  minuit->SetMaxIterations(500);

  // Go minimize!
  minuit->Migrad();

  // // Uncomment this for examining nuisance parameters.
  //  for (Int_t i=0;i<sys_size;i++){
  //    double opt_f;
  //    double err_f;
  //
  //    minuit->GetParameter(i,opt_f,err_f);
  //
  //    cout << opt_f << " ";
  //  }
  //  cout << endl;

  double minpoint = minuit->fAmin;

  return minpoint;

}

// RBT: more exception handling?
double NueFit2D::StandardLikelihood() {

  // HOO Matrix does not converge.  Don't use me for now.

  // Define the HOO matrix
  // if (ErrCalc!=0){

  // Calculate the decomposition systematic covariance matrix
  // ErrCalc->CalculateHOOError();

  // Invert HOO matrix and put it in inverror
  // CalculateDlnLMatrix(0,ErrCalc->CovMatrix_Decomp);

  // }

  // Call the likelihood equation either with or without systematics.

  double dlnl = 1.0e10;

  // empty vector
  vector<double> npar;

  unsigned int sys_size = FracErr_Bkgd_List.size();

  // If we are using systematics, do a minimization
  if (sys_size > 0) {

    dlnl = DoStdMinParam();

    // Check if the fit converged.
    Double_t fmin, fedm, errdef;
    Int_t npari, nparx, istat;
    minuit->mnstat(fmin, fedm, errdef, npari, nparx, istat);
    if (istat != 3 && fedm > 0.01) {
      cout << "Fit failed with status " << istat << endl;
      cout << "---> FMIN=" << fmin << " FEDM=" << fedm
           << " ERRDEF=" << errdef << " NPARI=" << npari
           << " NPARX=" << nparx << endl;
    }

  } else {
    // This will bypass the minimization and systematics.
    dlnl = StdLikeComparison(npar);
  }
  return dlnl;

}

double NueFit2D::StdLikeComparison(vector<double> npar) {

  // Actual likelihood comparison happens here.  This gets minimized.

  float dlnl = 0.0;
  double temp = 0.0;

  double sig = 0.0;
  double bkg = 0.0;
  double nobs = 0.0;
  double nexp = 0.0;
  double sigma_bkg = 0.0;
  double sigma_sig = 0.0;

  // unsigned int size_syst = npar.size();

  // Check there is one nuisance parameter per systematic:
  // (The "+nBins" is for the HOO background)
  // if (size_syst>0&&size_syst!=FracErr_Bkgd_List.size()+nBins){
  //  cout << "Error: different number of systematics and nuisance parameters!"
  //  << endl;
  //  return 1.0e10;
  // }


  for (unsigned int i = 0; i < nBins; i++) {
    // Signal and background predictions:
    sig = Sig->GetBinContent(i + 1);
    bkg = Bkgd->GetBinContent(i + 1);

    // Loop through systematics and add contributions
    for (unsigned int j = 0; j < FracErr_Bkgd_List.size(); j++) {

      // Get the nuisance parameter:
      double f = npar.at(j);

      // Total shift on background:
      sigma_bkg = FracErr_Bkgd_List[j]->GetBinContent(i + 1);
      bkg += f * sigma_bkg * bkg;

      // Total shift on signal:
      sigma_sig = FracErr_Sig_List[j]->GetBinContent(i + 1);
      sig += f * sigma_sig * sig;

    }

    // RBT: this is disabled for now:
    // if (InvErrorMatrix!=0){
    //  //RBT: Add a shift due to HOO error in this bin.
    //  unsigned int par_hoo = FracErr_Bkgd_List.size() + i;
    //  double g = npar.at(par_hoo);
    //  bkg += g;
    // }

    // "Data" distribution:
    nobs = NObs->GetBinContent(i + 1);

    // Expected distribution (with appropriate nuisance shifts)
    nexp = sig + bkg;
    temp = 0;

    // Likelihood comparison:
    if (nobs > 0 && nexp > 0) {
      // Regular (both distributions positive):
      temp = nexp - nobs + nobs * TMath::Log(nobs) - nobs * TMath::Log(nexp);
    } else if (nobs == 0 && nexp > 0) {
      // No data was seen in this bin:
      temp = nexp;
    } else if (nexp == 0 && nobs == 0) {
      // Nothing was seen, nothing expected:
      temp = 0;
    } else {
      // Something weird happened:
      return 1.0e10;
    }

    // Turn it into a chi2:
    dlnl  += 2.0 * temp;

  }

  // Add penalty term for nuisance parameters:
  for (unsigned int j = 0; j < FracErr_Bkgd_List.size(); j++) {
    dlnl += npar.at(j) * npar.at(j);
  }

  // HOO gets a covariance term as well:

  // if (InvErrorMatrix!=0){
  //  for (unsigned int i=FracErr_Bkgd_List.size(); i<FracErr_Bkgd_List.size()+nBins;i++){
  //    for (unsigned int j=FracErr_Bkgd_List.size(); j<FracErr_Bkgd_List.size()+nBins;j++){
  //	//This needs to be defined.  Will be zero for now!
  //	dlnl += npar.at(i)*InvErrorMatrix->GetBinContent(i+1,j+1)*npar.at(j);
  //    }
  //  }
  // }

  // Return the log likelihood:
  return dlnl;

}

// RBT: Should be fine here, but check the initialization, etc.
double NueFit2D::DoBinMinParam() {

  int sys_size = nBins;
  if (sys_size > minuit->fMaxpar) {
    cout << "WARNING: WRONG MINUIT SIZE" << endl;
  }
  Double_t * vstrt = new Double_t[sys_size];
  Double_t * stp = new Double_t[sys_size];
  Double_t * bmin = new Double_t[sys_size];
  Double_t * bmax = new Double_t[sys_size];
  Int_t ierflg = 0;
  // Pass on the object and set the static function:
  gMinuit = minuit;
  minuit->SetObjectFit(this);
  minuit->SetFCN(binFunction);

  // Set the parameters
  for (Int_t i = 0; i < sys_size; i++) {
    vstrt[i] = 0.0;
    stp[i] = 1.0;
    bmin[i] = 0.0;
    bmax[i] = 0.0;
    minuit->mnparm(i, Form("f%i", i), vstrt[i], stp[i], bmin[i], bmax[i], ierflg);
  }

  // Resets function value and errors to UNDEFINED:
  minuit->mnrset(1);

  // 1 std dev for dlnl:
  minuit->SetErrorDef(1.0);

  // Max iterations:
  minuit->SetMaxIterations(500);
  // Go minimize!
  minuit->SetPrintLevel(-1);
  cout << "Migrading..." << endl;
  minuit->Migrad();
  cout << "Done Migrading..." << endl;
  // Get the minimum for the function
  double minpoint = minuit->fAmin;

  return minpoint;

}

double NueFit2D::BinLikelihood() {

  // Define the covariance matrix
  if (ErrCalc != 0) {

    // Calculate the systematic error covariance matrix:
    ErrCalc->CalculateSystErrorMatrix();

    // Calculate the decomposition systematic covariance matrix
    ErrCalc->CalculateHOOError();

    // Combine HOOHE and Syst into matrix, and invert
    CalculateDlnLMatrix(ErrCalc->CovMatrix, ErrCalc->CovMatrix_Decomp);
  }


  double dlnl = 1.0e10;
  vector<double> npar;
  if (ErrCalc == 0) {
    // Just use the empty parameter array if no error matrix set:
    dlnl = BinLikeComparison(npar);
  } else {
    // Otherwise, do a Minuit minimization
    dlnl = DoBinMinParam();

    Double_t fmin, fedm, errdef;
    Int_t npari, nparx, istat;
    minuit->mnstat(fmin, fedm, errdef, npari, nparx, istat);
    if (istat != 3 && fedm > 0.01) {
      cout << "Fit failed with status " << istat << endl;
      cout << "---> FMIN=" << fmin << " FEDM=" << fedm
           << " ERRDEF=" << errdef << " NPARI=" << npari
           << " NPARX=" << nparx << endl;
    }

  }

  return dlnl;
}

double NueFit2D::BinLikeComparison(vector<double> npar) {

  float dlnl = 0.0;
  double temp = 0.0;

  double sig = 0.0;
  double bkg = 0.0;
  double nobs = 0.0;
  double nexp = 0.0;

  // Loop over fit bins:
  for (unsigned int i = 0; i < nBins; i++) {

    // Get Signal and background for prediction:
    sig = Sig->GetBinContent(i + 1);
    bkg = Bkgd->GetBinContent(i + 1);

    // Prediction:
    nexp = sig + bkg;

    // Add nuisance parameter shift if using matrix:
    if (npar.size() > 0 && ErrCalc != 0 && InvErrorMatrix != 0) {
      double f = npar.at(i);
      nexp += f;
    }

    // Observed distribution:
    nobs = NObs->GetBinContent(i + 1);

    temp = 0;

    // Likelihood comparison:
    if (nobs > 0 && nexp > 0) {
      // Regular (both distributions positive):
      temp = nexp - nobs + nobs * TMath::Log(nobs) - nobs * TMath::Log(nexp);
    } else if (nobs == 0 && nexp > 0) {
      // No data was seen in this bin:
      temp = nexp;
    } else if (nexp == 0 && nobs == 0) {
      // Nothing was seen, nothing expected:
      temp = 0;
    } else {
      // Something weird happened:
      return 1.0e10;
    }

    dlnl  += 2.0 * temp;

  }
  if (npar.size() > 0 && InvErrorMatrix != 0) {
    for (unsigned int i = 0; i < npar.size(); i++) {
      for (unsigned int j = 0; j < npar.size(); j++) {
        // Covariance terms
        dlnl += npar.at(i) * InvErrorMatrix->GetBinContent(i + 1, j + 1) * npar.at(j);
      }
    }
  }
  /*
     for (unsigned int j=0; j< npar.size();j++){
     cout << "NPAR"<<j<<" = "<< npar.at(j)<< endl;
     }
   */

  return dlnl;

}

void NueFit2D::FormatChi2Hists(int nx, double xlow, double xhigh, int ny, double ylow, double yhigh) {
  nSinSq2Th13_Chi2Hist = nx;
  SinSq2Th13Low_Chi2Hist = xlow;
  SinSq2Th13High_Chi2Hist = xhigh;
  nDelta_Chi2Hist = ny;
  DeltaLow_Chi2Hist = ylow;
  DeltaHigh_Chi2Hist = yhigh;

  return;
}
void NueFit2D::SetupChi2Hists() {
  Chi2_Normal = new TH2D("Chi2_Normal", "", nSinSq2Th13_Chi2Hist, SinSq2Th13Low_Chi2Hist, SinSq2Th13High_Chi2Hist, nDelta_Chi2Hist, DeltaLow_Chi2Hist, DeltaHigh_Chi2Hist);
  Chi2_Inverted = new TH2D("Chi2_Inverted", "", nSinSq2Th13_Chi2Hist, SinSq2Th13Low_Chi2Hist, SinSq2Th13High_Chi2Hist, nDelta_Chi2Hist, DeltaLow_Chi2Hist, DeltaHigh_Chi2Hist);

  return;
}
void NueFit2D::SetDeltaRange(double l, double h) {
  DeltaLow = l;
  DeltaHigh = h;

  return;
}
void NueFit2D::SetSinSq2Th13Range(double l, double h) {
  if (l < 0 || l > 1) {
    cout << "Unphysical value of SinSq2Th13.  Setting low value to 0." << endl;
    l = 0;
  }
  if (h < 0 || h > 1) {
    cout << "Unphysical value of SinSq2Th13.  Setting high value to 1." << endl;
    h = 1;
  }

  SinSq2Th13Low = l;
  SinSq2Th13High = h;

  return;
}
void NueFit2D::SetSinSq2Th14Range(double l, double h) {
  if (l < 0 || l > 1) {
    cout << "Unphysical value of SinSq2Th14.  Setting low value to 0." << endl;
    l = 0;
  }
  if (h < 0 || h > 1) {
    cout << "Unphysical value of SinSq2Th14.  Setting high value to 1." << endl;
    h = 1;
  }

  SinSq2Th14Low = l;
  SinSq2Th14High = h;

  return;
}

void NueFit2D::SetSinSqTh14Range(double l, double h) {
  if (l < 0 || l > 1) {
    cout << "Unphysical value of SinSqTh14.  Setting low value to 0." << endl;
    l = 0;
  }
  if (h < 0 || h > 1) {
    cout << "Unphysical value of SinSqTh14.  Setting high value to 1." << endl;
    h = 1;
  }

  SinSqTh14Low = l;
  SinSqTh14High = h;

  return;
}

void NueFit2D::SetSinSqTh24Range(double l, double h) {
  if (l < 0 || l > 1) {
    cout << "Unphysical value of SinSqTh24.  Setting low value to 0." << endl;
    l = 0;
  }
  if (h < 0 || h > 1) {
    cout << "Unphysical value of SinSqTh24.  Setting high value to 1." << endl;
    h = 1;
  }

  SinSqTh24Low = l;
  SinSqTh24High = h;

  return;
}

void NueFit2D::SetepsetauRange(double l, double h) {
  /*
     if(l<0 || l>1)
     {
     cout<<"Unphysical value of Eps_etau.  Setting low value to 0."<<endl;
     l=0;
     }
     if(h<0 || h>1)
     {
     cout<<"Unphysical value of Eps_etau.  Setting high value to 1."<<endl;
     h=1;
     }
   */
  epsetauLow = l;
  epsetauHigh = h;

  return;
}

void NueFit2D::RunStandardChi2Sensitivity() {
  if (NObs == 0) {
    cout << "NObs not set.  Quitting..." << endl;
    return;
  }
  if (Extrap.size() == 0) {
    cout << "No Extrapolate2D input.  Quitting..." << endl;
    return;
  }

  Bkgd = (TH1D *) NObs->Clone("Bkgd");
  Bkgd->Reset();
  Sig = (TH1D *) NObs->Clone("Sig");
  Sig->Reset();
  TH1D * NExp = (TH1D *) NObs->Clone("NExp");
  NExp->Reset();

  unsigned int ie;
  for (ie = 0; ie < Extrap.size(); ie++) {
    Extrap[ie]->GetPrediction();
  }

  if (ErrorMatrix == 0) {
    ErrorMatrix = new TH2D("ErrorMatrix", "", nBins, -0.5, nBins - 0.5, nBins, -0.5, nBins - 0.5);
  }
  if (InvErrorMatrix == 0) {
    InvErrorMatrix = new TH2D("InvErrorMatrix", "", nBins, -0.5, nBins - 0.5, nBins, -0.5, nBins - 0.5);
  }

  Int_t i, j;

  Double_t ss2th13 = 0;
  Double_t delta = 0;
  Double_t chi2 = 0, chi2prev = 0;
  Double_t contourlvl = 2.71;
  Double_t prev = 0;
  Double_t best = 0;
  Double_t Deltaincrement = 0;
  if (nDeltaSteps > 0) {
    Deltaincrement = (DeltaHigh - DeltaLow) / (nDeltaSteps);
  }
  Double_t Th13increment = 0;
  if (nSinSq2Th13Steps > 0) {
    Th13increment = (SinSq2Th13High - SinSq2Th13Low) / (nSinSq2Th13Steps);
  }

  vector<Double_t> limit90_norm;
  vector<Double_t> limit90_invt;
  vector<double> deltapts_norm;
  vector<double> deltapts_invt;

  for (j = 0; j < nDeltaSteps + 1; j++) {
    delta = j * Deltaincrement * TMath::Pi() + DeltaLow;
    for (ie = 0; ie < Extrap.size(); ie++) {
      Extrap[ie]->SetDeltaCP(delta);
    }

    contourlvl = 2.71;

    best = 0;
    chi2prev = 0;
    prev = 0;
    for (i = 0; i < nSinSq2Th13Steps + 1; i++) {
      ss2th13 = i * Th13increment + SinSq2Th13Low;
      for (ie = 0; ie < Extrap.size(); ie++) {
        Extrap[ie]->SetSinSq2Th13(ss2th13);
        Extrap[ie]->OscillatePrediction();
      }

      Bkgd->Reset();
      Sig->Reset();
      NExp->Reset();

      for (ie = 0; ie < Extrap.size(); ie++) {
        Bkgd->Add(Extrap[ie]->Pred_TotalBkgd_VsBinNumber);
        Sig->Add(Extrap[ie]->Pred_Signal_VsBinNumber);
      }
      NExp->Add(Bkgd);
      NExp->Add(Sig);

      chi2 = StandardChi2(NExp);

      if (chi2 > contourlvl && chi2prev < contourlvl) {
        best = ss2th13 + ((ss2th13 - prev) / (chi2 - chi2prev)) * (contourlvl - chi2);
        limit90_norm.push_back(best);
        deltapts_norm.push_back(delta / TMath::Pi());
        cout << chi2prev << ", " << prev << ", " << chi2 << ", " << ss2th13 << ", " << best << endl;
        break;
      }
      chi2prev = chi2;
      prev = ss2th13;
    }
  }

  for (ie = 0; ie < Extrap.size(); ie++) {
    Extrap[ie]->InvertMassHierarchy();
  }

  for (j = 0; j < nDeltaSteps + 1; j++) {
    delta = j * Deltaincrement * TMath::Pi() + DeltaLow;
    for (ie = 0; ie < Extrap.size(); ie++) {
      Extrap[ie]->SetDeltaCP(delta);
    }

    contourlvl = 2.71;

    best = 0;
    chi2prev = 0;
    prev = 0;
    for (i = 0; i < nSinSq2Th13Steps + 1; i++) {
      ss2th13 = i * Th13increment + SinSq2Th13Low;
      for (ie = 0; ie < Extrap.size(); ie++) {
        Extrap[ie]->SetSinSq2Th13(ss2th13);
        Extrap[ie]->OscillatePrediction();
      }

      Bkgd->Reset();
      Sig->Reset();
      NExp->Reset();

      for (ie = 0; ie < Extrap.size(); ie++) {
        Bkgd->Add(Extrap[ie]->Pred_TotalBkgd_VsBinNumber);
        Sig->Add(Extrap[ie]->Pred_Signal_VsBinNumber);
      }
      NExp->Add(Bkgd);
      NExp->Add(Sig);

      chi2 = StandardChi2(NExp);

      if (chi2 > contourlvl && chi2prev < contourlvl) {
        best = ss2th13 + ((ss2th13 - prev) / (chi2 - chi2prev)) * (contourlvl - chi2);
        limit90_invt.push_back(best);
        deltapts_invt.push_back(delta / TMath::Pi());
        cout << chi2prev << ", " << prev << ", " << chi2 << ", " << ss2th13 << ", " << best << endl;
        break;
      }
      chi2prev = chi2;
      prev = ss2th13;
    }
  }

  if (limit90_norm.size() == 0) {
    cout << chi2 << endl;
    cout << "Didn't find limit.  Quitting..." << endl;
    return;
  }

  double * s_n = new double[nDeltaSteps + 1];
  double * s_i = new double[nDeltaSteps + 1];
  double * de_n = new double[nDeltaSteps + 1];
  double * de_i = new double[nDeltaSteps + 1];
  Int_t n_n, n_i;
  n_n = limit90_norm.size();
  n_i = limit90_invt.size();
  for (i = 0; i < n_n; i++) {
    de_n[i] = deltapts_norm.at(i);
    s_n[i] = limit90_norm.at(i);
  }
  for (i = 0; i < n_i; i++) {
    de_i[i] = deltapts_invt.at(i);
    s_i[i] = limit90_invt.at(i);
  }
  TGraph * gn = new TGraph(n_n, s_n, de_n);
  gn->SetMarkerStyle(20);
  gn->SetTitle("");
  gn->GetYaxis()->SetTitle("#delta");
  gn->GetXaxis()->SetTitle("sin^{2}2#theta_{13}");
  gn->GetXaxis()->SetLimits(0, 0.6);
  gn->SetLineWidth(2);
  gn->SetLineColor(kBlue);
  gn->SetMarkerColor(kBlue);
  gn->SetMaximum(2);
  gn->SetName("Limit_Normal");
  TGraph * gi = new TGraph(n_i, s_i, de_i);
  gi->SetMarkerStyle(20);
  gi->SetTitle("");
  gi->GetYaxis()->SetTitle("#delta");
  gi->GetXaxis()->SetTitle("sin^{2}2#theta_{13}");
  gi->GetXaxis()->SetLimits(0, 0.6);
  gi->SetLineWidth(2);
  gi->SetLineColor(kRed);
  gi->SetMarkerColor(kRed);
  gi->SetMaximum(2);
  gi->SetName("Limit_Inverted");

  cout << "90% confidence level limit = " << limit90_norm.at(0) << ", " << limit90_invt.at(0) << endl;

  TFile * fout = new TFile(gSystem->ExpandPathName(outFileName.c_str()), "RECREATE");
  gn->Write();
  gi->Write();
  fout->Write();
  fout->Close();

  delete [] s_n;
  delete [] s_i;
  delete [] de_n;
  delete [] de_i;

  return;
}
void NueFit2D::RunScaledChi2Sensitivity() {
  if (NObs == 0) {
    cout << "NObs not set.  Quitting..." << endl;
    return;
  }

  if (FracErr_Bkgd == 0) {
    cout << "FracErr_Bkgd not set. Quitting..." << endl;
    return;
  }

  if (Extrap.size() == 0) {
    cout << "No Extrapolate2D input.  Quitting..." << endl;
    return;
  }

  Bkgd = (TH1D *) NObs->Clone("Bkgd");
  Bkgd->Reset();
  Sig = (TH1D *) NObs->Clone("Sig");
  Sig->Reset();
  FracErr_Sig = (TH1D *) NObs->Clone("FracErr_Sig");
  FracErr_Sig->Reset();
  // for sensitivity, signal err is ignored

  unsigned int ie;
  for (ie = 0; ie < Extrap.size(); ie++) {
    Extrap[ie]->GetPrediction();
  }

  Int_t i, j;

  Double_t ss2th13 = 0;
  Double_t delta = 0;
  Double_t chi2 = 0, chi2prev = 0;
  Double_t contourlvl = 2.71;
  Double_t prev = 0;
  Double_t best = 0;
  Double_t Deltaincrement = 0;
  if (nDeltaSteps > 0) {
    Deltaincrement = (DeltaHigh - DeltaLow) / (nDeltaSteps);
  }
  Double_t Th13increment = 0;
  if (nSinSq2Th13Steps > 0) {
    Th13increment = (SinSq2Th13High - SinSq2Th13Low) / (nSinSq2Th13Steps);
  }

  vector<Double_t> limit90_norm;
  vector<Double_t> limit90_invt;
  vector<double> deltapts_norm;
  vector<double> deltapts_invt;

  for (j = 0; j < nDeltaSteps + 1; j++) {
    delta = j * Deltaincrement * TMath::Pi() + DeltaLow;
    for (ie = 0; ie < Extrap.size(); ie++) {
      Extrap[ie]->SetDeltaCP(delta);
    }

    contourlvl = 2.71;

    best = 0;
    chi2prev = 0;
    prev = 0;
    for (i = 0; i < nSinSq2Th13Steps + 1; i++) {
      ss2th13 = i * Th13increment + SinSq2Th13Low;
      for (ie = 0; ie < Extrap.size(); ie++) {
        Extrap[ie]->SetSinSq2Th13(ss2th13);
        Extrap[ie]->OscillatePrediction();
      }

      Bkgd->Reset();
      Sig->Reset();

      for (ie = 0; ie < Extrap.size(); ie++) {
        Bkgd->Add(Extrap[ie]->Pred_TotalBkgd_VsBinNumber);
        Sig->Add(Extrap[ie]->Pred_Signal_VsBinNumber);
      }

      chi2 = ScaledChi2(Bkgd, Sig);

      if (chi2 > contourlvl && chi2prev < contourlvl) {
        best = ss2th13 + ((ss2th13 - prev) / (chi2 - chi2prev)) * (contourlvl - chi2);
        limit90_norm.push_back(best);
        deltapts_norm.push_back(delta / TMath::Pi());
        cout << chi2prev << ", " << prev << ", " << chi2 << ", " << ss2th13 << ", " << best << endl;
        break;
      }
      chi2prev = chi2;
      prev = ss2th13;
    }
  }

  for (ie = 0; ie < Extrap.size(); ie++) {
    Extrap[ie]->InvertMassHierarchy();
  }

  for (j = 0; j < nDeltaSteps + 1; j++) {
    delta = j * Deltaincrement * TMath::Pi() + DeltaLow;
    for (ie = 0; ie < Extrap.size(); ie++) {
      Extrap[ie]->SetDeltaCP(delta);
    }

    contourlvl = 2.71;

    best = 0;
    chi2prev = 0;
    prev = 0;
    for (i = 0; i < nSinSq2Th13Steps + 1; i++) {
      ss2th13 = i * Th13increment + SinSq2Th13Low;
      for (ie = 0; ie < Extrap.size(); ie++) {
        Extrap[ie]->SetSinSq2Th13(ss2th13);
        Extrap[ie]->OscillatePrediction();
      }

      Bkgd->Reset();
      Sig->Reset();

      for (ie = 0; ie < Extrap.size(); ie++) {
        Bkgd->Add(Extrap[ie]->Pred_TotalBkgd_VsBinNumber);
        Sig->Add(Extrap[ie]->Pred_Signal_VsBinNumber);
      }

      chi2 = ScaledChi2(Bkgd, Sig);

      if (chi2 > contourlvl && chi2prev < contourlvl) {
        best = ss2th13 + ((ss2th13 - prev) / (chi2 - chi2prev)) * (contourlvl - chi2);
        limit90_invt.push_back(best);
        deltapts_invt.push_back(delta / TMath::Pi());
        cout << chi2prev << ", " << prev << ", " << chi2 << ", " << ss2th13 << ", " << best << endl;
        break;
      }
      chi2prev = chi2;
      prev = ss2th13;
    }
  }

  if (limit90_norm.size() == 0) {
    cout << chi2 << endl;
    cout << "Didn't find limit.  Quitting..." << endl;
    return;
  }

  double * s_n = new double[nDeltaSteps + 1];
  double * s_i = new double[nDeltaSteps + 1];
  double * de_n = new double[nDeltaSteps + 1];
  double * de_i = new double[nDeltaSteps + 1];
  Int_t n_n, n_i;
  n_n = limit90_norm.size();
  n_i = limit90_invt.size();
  for (i = 0; i < n_n; i++) {
    de_n[i] = deltapts_norm.at(i);
    s_n[i] = limit90_norm.at(i);
  }
  for (i = 0; i < n_i; i++) {
    de_i[i] = deltapts_invt.at(i);
    s_i[i] = limit90_invt.at(i);
  }
  TGraph * gn = new TGraph(n_n, s_n, de_n);
  gn->SetMarkerStyle(20);
  gn->SetTitle("");
  gn->GetYaxis()->SetTitle("#delta");
  gn->GetXaxis()->SetTitle("sin^{2}2#theta_{13}");
  gn->GetXaxis()->SetLimits(0, 0.6);
  gn->SetLineWidth(2);
  gn->SetLineColor(kBlue);
  gn->SetMarkerColor(kBlue);
  gn->SetMaximum(2);
  gn->SetName("Limit_Normal");
  TGraph * gi = new TGraph(n_i, s_i, de_i);
  gi->SetMarkerStyle(20);
  gi->SetTitle("");
  gi->GetYaxis()->SetTitle("#delta");
  gi->GetXaxis()->SetTitle("sin^{2}2#theta_{13}");
  gi->GetXaxis()->SetLimits(0, 0.6);
  gi->SetLineWidth(2);
  gi->SetLineColor(kRed);
  gi->SetMarkerColor(kRed);
  gi->SetMaximum(2);
  gi->SetName("Limit_Inverted");

  cout << "90% confidence level limit = " << limit90_norm.at(0) << ", " << limit90_invt.at(0) << endl;

  TFile * fout = new TFile(gSystem->ExpandPathName(outFileName.c_str()), "RECREATE");
  gn->Write();
  gi->Write();
  fout->Write();
  fout->Close();

  delete [] s_n;
  delete [] s_i;
  delete [] de_n;
  delete [] de_i;

  return;
}
void NueFit2D::RunFCAnalytical() {
  if (NObs == 0) {
    cout << "NObs not set.  Quitting..." << endl;
    return;
  }
  if (FracErr_Bkgd == 0) {
    cout << "FracErr_Bkgd not set.  Quitting..." << endl;
    return;
  }
  if (FracErr_Sig == 0) {
    cout << "FracErr_Sig not set.  Quitting..." << endl;
    return;
  }

  if (nBins > 1) {
    cout << "Analytical Feldman-Cousins Method should only be used for a single bin fit.  Quitting..." << endl;
    return;
  }

  ReadGridFiles();
  SetupChi2Hists();

  if (nPts_Normal == 0 || nPts_Inverted == 0) {
    return;
  }

  int i;
  double omega;

  for (i = 0; i < nPts_Normal; i++) {
    if (i % 100 == 0) {
      cout << 100. * i / nPts_Normal << "% complete for normal hierarchy" << endl;
    }
    GridTree_Normal[0]->GetEntry(i);
    omega = EvaluateOmega(grid_signal * GridScale_Normal, grid_background * GridScale_Normal);
    Chi2_Normal->Fill(grid_sinsq2th13, grid_delta / TMath::Pi(), omega);
  }

  for (i = 0; i < nPts_Inverted; i++) {
    if (i % 100 == 0) {
      cout << 100. * i / nPts_Inverted << "% complete for inverted hierarchy" << endl;
    }
    GridTree_Inverted[0]->GetEntry(i);
    omega = EvaluateOmega(grid_signal * GridScale_Inverted, grid_background * GridScale_Inverted);
    Chi2_Inverted->Fill(grid_sinsq2th13, grid_delta / TMath::Pi(), omega);
  }

  TFile * f = new TFile(gSystem->ExpandPathName(outFileName.c_str()), "RECREATE");
  Chi2_Normal->Write();
  Chi2_Inverted->Write();
  f->Close();

  return;
}
void NueFit2D::RunFCTraditional() {
  if (NObs == 0) {
    cout << "NObs not set.  Quitting..." << endl;
    return;
  }
  if (Extrap.size() == 0) {
    cout << "No Extrapolate2D input.  Quitting..." << endl;
    return;
  }

  ReadGridFiles();
  SetupChi2Hists();

  if (nPts_Normal == 0 || nPts_Inverted == 0) {
    return;
  }

  TH1D * nexp_bkgd = new TH1D("nexp_bkgd", "", nBins, -0.5, nBins - 0.5);
  TH1D * nexp_signal = new TH1D("nexp_signal", "", nBins, -0.5, nBins - 0.5);
  TH1D * nexp = new TH1D("nexp", "", nBins, -0.5, nBins - 0.5);

  bool calcchi2min; // if nobs<nexpmin for any bin, then chi2min!=0 and needs to be calculated
  TH1D * nexpmin = 0;

  int i;
  unsigned int j;

  if (gSystem->AccessPathName(gSystem->ExpandPathName(PseudoExpFile.c_str()))) {
    cout << "Pseudo-experiment file doesn't exist." << endl;
    return;
  }

  TFile * f = new TFile(gSystem->ExpandPathName(PseudoExpFile.c_str()), "READ");

  double chi2data;
  int chi2databin;
  TH1D * chi2hist = (TH1D *) f->Get("Chi2Hist_Normal_0");
  double chi2binwidth = chi2hist->GetBinWidth(1);
  double chi2start = chi2hist->GetXaxis()->GetBinLowEdge(1);
  double frac;
  delete chi2hist;

  // normal hierarchy
  for (j = 0; j < Extrap.size(); j++) {
    Extrap[j]->GetPrediction();
    Extrap[j]->SetSinSq2Th13(0);
    Extrap[j]->OscillatePrediction();
    if (j == 0) {
      nexpmin = (TH1D *) Extrap[j]->Pred_TotalBkgd_VsBinNumber->Clone("nexpmin");
      nexpmin->Add(Extrap[j]->Pred_Signal_VsBinNumber);
    } else {
      nexpmin->Add(Extrap[j]->Pred_TotalBkgd_VsBinNumber);
      nexpmin->Add(Extrap[j]->Pred_Signal_VsBinNumber);
    }
  }

  for (i = 0; i < nPts_Normal; i++) {
    if (i % 100 == 0) {
      cout << 100. * i / nPts_Normal << "% complete for normal hierarchy" << endl;
    }

    nexp_bkgd->Reset();
    nexp_signal->Reset();
    nexp->Reset();
    for (j = 0; j < nBins; j++) {
      GridTree_Normal[j]->GetEntry(i);
      nexp_bkgd->SetBinContent(j + 1, grid_background * GridScale_Normal);
      nexp_signal->SetBinContent(j + 1, grid_signal * GridScale_Normal);
    }
    nexp->Add(nexp_bkgd, nexp_signal, 1, 1);
    chi2data = PoissonChi2(nexp);
    calcchi2min = false;
    for (j = 0; j < nBins; j++) {
      if (NObs->GetBinContent(j + 1) < nexpmin->GetBinContent(j + 1)) {
        calcchi2min = true;
      }
    }
    if (calcchi2min) {
      chi2data = chi2data - PoissonChi2(nexpmin);
    }
    chi2databin = int((chi2data - chi2start) / chi2binwidth) + 1;

    TH1D * chi2hist = (TH1D *) f->Get(Form("Chi2Hist_Normal_%i", i));
    frac = chi2hist->Integral(1, chi2databin - 1) / chi2hist->Integral();
    if (chi2databin == 1) {
      frac = 0;
    }

    Chi2_Normal->Fill(grid_sinsq2th13, grid_delta / TMath::Pi(), frac);

    delete chi2hist;
  }

  // inverted hierarchy
  for (j = 0; j < Extrap.size(); j++) {
    Extrap[j]->SetSinSq2Th13(0);
    Extrap[j]->InvertMassHierarchy();
    Extrap[j]->OscillatePrediction();
    if (j == 0) {
      nexpmin = (TH1D *) Extrap[j]->Pred_TotalBkgd_VsBinNumber->Clone("nexpmin");
      nexpmin->Add(Extrap[j]->Pred_Signal_VsBinNumber);
    } else {
      nexpmin->Add(Extrap[j]->Pred_TotalBkgd_VsBinNumber);
      nexpmin->Add(Extrap[j]->Pred_Signal_VsBinNumber);
    }
  }

  for (i = 0; i < nPts_Inverted; i++) {
    if (i % 100 == 0) {
      cout << 100. * i / nPts_Inverted << "% complete for inverted hierarchy" << endl;
    }

    nexp_bkgd->Reset();
    nexp_signal->Reset();
    nexp->Reset();
    for (j = 0; j < nBins; j++) {
      GridTree_Inverted[j]->GetEntry(i);
      nexp_bkgd->SetBinContent(j + 1, grid_background * GridScale_Inverted);
      nexp_signal->SetBinContent(j + 1, grid_signal * GridScale_Inverted);
    }
    nexp->Add(nexp_bkgd, nexp_signal, 1, 1);
    chi2data = PoissonChi2(nexp);
    calcchi2min = false;
    for (j = 0; j < nBins; j++) {
      if (NObs->GetBinContent(j + 1) < nexpmin->GetBinContent(j + 1)) {
        calcchi2min = true;
      }
    }
    if (calcchi2min) {
      chi2data = chi2data - PoissonChi2(nexpmin);
    }
    chi2databin = int((chi2data - chi2start) / chi2binwidth) + 1;

    TH1D * chi2hist = (TH1D *) f->Get(Form("Chi2Hist_Inverted_%i", i));
    frac = chi2hist->Integral(1, chi2databin - 1) / chi2hist->Integral();
    if (chi2databin == 1) {
      frac = 0;
    }

    Chi2_Inverted->Fill(grid_sinsq2th13, grid_delta / TMath::Pi(), frac);

    delete chi2hist;
  }

  f->Close();

  TFile * fout = new TFile(gSystem->ExpandPathName(outFileName.c_str()), "RECREATE");
  Chi2_Normal->Write();
  Chi2_Inverted->Write();
  fout->Close();

  return;
}
void NueFit2D::RunMultiBinFC() {
  if (NObs == 0) {
    cout << "NObs not set.  Quitting..." << endl;
    return;
  }
  if (Extrap.size() == 0) {
    cout << "No Extrapolate2D input.  Quitting..." << endl;
    return;
  }
  for (unsigned int ie = 0; ie < Extrap.size(); ie++) {
    Extrap[ie]->GetPrediction();
  }
  if (FitMethod == 1 && (FracErr_Bkgd == 0 || FracErr_Sig == 0)) {
    cout << "FracErr_Bkgd and FracErr_Sig need to be set for ScaledChi2.  Quitting..." << endl;
    return;
  }
  if (ErrCalc == 0) {
    cout << "No ErrorCalc object set!  Quitting..." << endl;
    return;
  }
  if (FitMethod == 3) {
    DefineStdDlnLMinuit();
  }
  if (FitMethod == 4) {
    DefineBinDlnLMinuit();
  }

  nBins = NObs->GetNbinsX();
  if (ErrorMatrix == 0) {
    ErrorMatrix = new TH2D("ErrorMatrix", "", nBins, -0.5, nBins - 0.5, nBins, -0.5, nBins - 0.5);
  }
  if (InvErrorMatrix == 0) {
    InvErrorMatrix = new TH2D("InvErrorMatrix", "", nBins, -0.5, nBins - 0.5, nBins, -0.5, nBins - 0.5);
  }

  ReadGridFiles();
  SetupChi2Hists();

  if (nPts_Normal == 0 || nPts_Inverted == 0) {
    return;
  }

  TH1D * nexp_bkgd = new TH1D("nexp_bkgd", "", nBins, -0.5, nBins - 0.5);
  TH1D * nexp_signal = new TH1D("nexp_signal", "", nBins, -0.5, nBins - 0.5);
  TH1D * nexp = new TH1D("nexp", "", nBins, -0.5, nBins - 0.5);

  int i;
  unsigned int j, k;

  if (gSystem->AccessPathName(gSystem->ExpandPathName(PseudoExpFile.c_str()))) {
    cout << "Pseudo-experiment file doesn't exist." << endl;
    return;
  }

  TFile * f = new TFile(gSystem->ExpandPathName(PseudoExpFile.c_str()), "READ");

  double chi2data, chi2min;
  int chi2databin;
  TH1D * chi2hist = (TH1D *) f->Get("Chi2Hist_Normal_0");
  double chi2binwidth = chi2hist->GetBinWidth(1);
  double chi2start = chi2hist->GetXaxis()->GetBinLowEdge(1);
  double frac;
  delete chi2hist;

  vector< vector<double> > nc, numucc, bnuecc, nutaucc, sig;
  for (j = 0; j < nBins; j++) {
    nc.push_back(vector<double>() );
    numucc.push_back(vector<double>() );
    bnuecc.push_back(vector<double>() );
    nutaucc.push_back(vector<double>() );
    sig.push_back(vector<double>() );
    for (k = 0; k < Extrap.size(); k++) {
      nc[j].push_back(0);
      numucc[j].push_back(0);
      bnuecc[j].push_back(0);
      nutaucc[j].push_back(0);
      sig[j].push_back(0);
    }
  }

  Bkgd = (TH1D *) NObs->Clone("Bkgd");
  Bkgd->Reset();
  Sig = (TH1D *) NObs->Clone("Sig");
  Sig->Reset();

  // normal hierarchy

  for (i = 0; i < nPts_Normal; i++) {
    if (i % 100 == 0) {
      cout << 100. * i / nPts_Normal << "% complete for normal hierarchy" << endl;
    }

    nexp_bkgd->Reset();
    nexp_signal->Reset();
    nexp->Reset();

    for (j = 0; j < nBins; j++) {
      GridTree_Normal[j]->GetEntry(i);
      nexp_bkgd->SetBinContent(j + 1, grid_background * GridScale_Normal);
      nexp_signal->SetBinContent(j + 1, grid_signal * GridScale_Normal);

      for (k = 0; k < Extrap.size(); k++) {
        GridTree_2_Normal[j][k]->GetEntry(i);
        nc[j][k] = grid_nc * GridScale_Normal;
        numucc[j][k] = grid_numucc * GridScale_Normal;
        bnuecc[j][k] = grid_bnuecc * GridScale_Normal;
        nutaucc[j][k] = grid_nutaucc * GridScale_Normal;
        sig[j][k] = grid_nue * GridScale_Normal;
      }
    }
    nexp->Add(nexp_bkgd, nexp_signal, 1, 1);
    ErrCalc->SetGridPred(nBins, nc, numucc, bnuecc, nutaucc, sig);

    chi2min = GetMinLikelihood(grid_delta, true);

    ErrCalc->SetUseGrid(true); // use grid predictions set up above
    chi2data = 1e10;
    if (FitMethod == 0) {
      chi2data = PoissonChi2(nexp) - chi2min;
    } else if (FitMethod == 1) {
      chi2data = ScaledChi2(Bkgd, Sig) - chi2min;
    } else if (FitMethod == 2) {
      chi2data = StandardChi2(nexp) - chi2min;
    } else if (FitMethod == 3) {
      // Likelihood: "Standard" (N syst, N nuisance)
      // Calculate the likelihood (x2 for chi)
      Bkgd->Reset();
      Bkgd->Add(nexp_bkgd);
      Sig->Reset();
      Sig->Add(nexp_signal);
      chi2data = StandardLikelihood() - chi2min;
    } else if (FitMethod == 4) {
      // Likelihood: Bin by Bin Calculation of Systematics
      // Calculate the likelihood (x2 for chi)
      Bkgd->Reset();
      Bkgd->Add(nexp_bkgd);
      Sig->Reset();
      Sig->Add(nexp_signal);
      chi2data = BinLikelihood() - chi2min;
    } else {
      cout << "Error in RunMultiBinFC(): Unknown 'FitMethod'." << endl;
    }
    chi2databin = int((chi2data - chi2start) / chi2binwidth) + 1;

    TH1D * chi2hist = (TH1D *) f->Get(Form("Chi2Hist_Normal_%i", i));
    if (chi2hist->Integral() < 1) {
      cout << "Warning, chi2hist is empty." << endl;
      frac = 0;
    } else {
      frac = chi2hist->Integral(1, chi2databin - 1) / chi2hist->Integral();
    }

    if (chi2databin == 1) {
      frac = 0;
    }

    Chi2_Normal->Fill(grid_sinsq2th13, grid_delta / TMath::Pi(), frac);

    delete chi2hist;
  }

  // inverted hierarchy

  for (i = 0; i < nPts_Inverted; i++) {
    if (i % 100 == 0) {
      cout << 100. * i / nPts_Inverted << "% complete for inverted hierarchy" << endl;
    }

    nexp_bkgd->Reset();
    nexp_signal->Reset();
    nexp->Reset();

    for (j = 0; j < nBins; j++) {
      GridTree_Inverted[j]->GetEntry(i);
      nexp_bkgd->SetBinContent(j + 1, grid_background * GridScale_Inverted);
      nexp_signal->SetBinContent(j + 1, grid_signal * GridScale_Inverted);

      for (k = 0; k < Extrap.size(); k++) {
        GridTree_2_Inverted[j][k]->GetEntry(i);
        nc[j][k] = grid_nc * GridScale_Inverted;
        numucc[j][k] = grid_numucc * GridScale_Inverted;
        bnuecc[j][k] = grid_bnuecc * GridScale_Inverted;
        nutaucc[j][k] = grid_nutaucc * GridScale_Inverted;
        sig[j][k] = grid_nue * GridScale_Inverted;
      }
    }
    nexp->Add(nexp_bkgd, nexp_signal, 1, 1);
    ErrCalc->SetGridPred(nBins, nc, numucc, bnuecc, nutaucc, sig);

    chi2min = GetMinLikelihood(grid_delta, false);

    ErrCalc->SetUseGrid(true); // use grid predictions set up above
    chi2data = 1e10;
    if (FitMethod == 0) {
      chi2data = PoissonChi2(nexp) - chi2min;
    } else if (FitMethod == 1) {
      chi2data = ScaledChi2(Bkgd, Sig) - chi2min;
    } else if (FitMethod == 2) {
      chi2data = StandardChi2(nexp) - chi2min;
    } else if (FitMethod == 3) {
      // Likelihood: "Standard" (N syst, N nuisance)
      // Calculate the likelihood (x2 for chi)
      Bkgd->Reset();
      Bkgd->Add(nexp_bkgd);
      Sig->Reset();
      Sig->Add(nexp_signal);
      chi2data = StandardLikelihood() - chi2min;
    } else if (FitMethod == 4) {
      // Likelihood: Bin by Bin Calculation of Systematics
      // Calculate the likelihood (x2 for chi)
      Bkgd->Reset();
      Bkgd->Add(nexp_bkgd);
      Sig->Reset();
      Sig->Add(nexp_signal);
      chi2data = BinLikelihood() - chi2min;
    } else {
      cout << "Error in RunMultiBinFC(): Unknown 'FitMethod'." << endl;
    }
    chi2databin = int((chi2data - chi2start) / chi2binwidth) + 1;

    TH1D * chi2hist = (TH1D *) f->Get(Form("Chi2Hist_Inverted_%i", i));
    if (chi2hist->Integral() < 1) {
      cout << "Warning, chi2hist is empty." << endl;
      frac = 0;
    } else {
      frac = chi2hist->Integral(1, chi2databin - 1) / chi2hist->Integral();
    }

    if (chi2databin == 1) {
      frac = 0;
    }

    Chi2_Inverted->Fill(grid_sinsq2th13, grid_delta / TMath::Pi(), frac);

    delete chi2hist;
  }

  f->Close();

  TFile * fout = new TFile(gSystem->ExpandPathName(outFileName.c_str()), "RECREATE");
  Chi2_Normal->Write();
  Chi2_Inverted->Write();
  fout->Close();

  return;
}
void NueFit2D::RunDataGrid(string filenorm, string fileinvt) {
  if (NObs == 0) {
    cout << "NObs not set.  Quitting..." << endl;
    return;
  }
  if (Extrap.size() == 0) {
    cout << "No Extrapolate2D input.  Quitting..." << endl;
    return;
  }
  for (unsigned int ie = 0; ie < Extrap.size(); ie++) {
    Extrap[ie]->GetPrediction();
  }
  if (FitMethod == 1 && (FracErr_Bkgd == 0 || FracErr_Sig == 0)) {
    cout << "FracErr_Bkgd and FracErr_Sig need to be set for ScaledChi2.  Quitting..." << endl;
    return;
  }
  if (ErrCalc == 0) {
    cout << "No ErrorCalc object set!  Quitting..." << endl;
    return;
  }
  if (FitMethod == 3) {
    DefineStdDlnLMinuit();
  }
  if (FitMethod == 4) {
    DefineBinDlnLMinuit();
  }

  nBins = NObs->GetNbinsX();
  if (ErrorMatrix == 0) {
    ErrorMatrix = new TH2D("ErrorMatrix", "", nBins, -0.5, nBins - 0.5, nBins, -0.5, nBins - 0.5);
  }
  if (InvErrorMatrix == 0) {
    InvErrorMatrix = new TH2D("InvErrorMatrix", "", nBins, -0.5, nBins - 0.5, nBins, -0.5, nBins - 0.5);
  }

  ReadGridFiles();

  if (nPts_Normal == 0 || nPts_Inverted == 0) {
    return;
  }

  TH1D * nexp_bkgd = new TH1D("nexp_bkgd", "", nBins, -0.5, nBins - 0.5);
  TH1D * nexp_signal = new TH1D("nexp_signal", "", nBins, -0.5, nBins - 0.5);
  TH1D * nexp = new TH1D("nexp", "", nBins, -0.5, nBins - 0.5);

  int i;
  unsigned int j, k;

  double chi2data, chi2min;

  vector< vector<double> > nc, numucc, bnuecc, nutaucc, sig;
  for (j = 0; j < nBins; j++) {
    nc.push_back(vector<double>() );
    numucc.push_back(vector<double>() );
    bnuecc.push_back(vector<double>() );
    nutaucc.push_back(vector<double>() );
    sig.push_back(vector<double>() );
    for (k = 0; k < Extrap.size(); k++) {
      nc[j].push_back(0);
      numucc[j].push_back(0);
      bnuecc[j].push_back(0);
      nutaucc[j].push_back(0);
      sig[j].push_back(0);
    }
  }

  Bkgd = (TH1D *) NObs->Clone("Bkgd");
  Bkgd->Reset();
  Sig = (TH1D *) NObs->Clone("Sig");
  Sig->Reset();

  // normal hierarchy

  ofstream myfile;
  myfile.open(gSystem->ExpandPathName(filenorm.c_str()));

  for (i = 0; i < nPts_Normal; i++) {
    if (i % 100 == 0) {
      cout << 100. * i / nPts_Normal << "% complete for normal hierarchy" << endl;
    }

    nexp_bkgd->Reset();
    nexp_signal->Reset();
    nexp->Reset();

    for (j = 0; j < nBins; j++) {
      GridTree_Normal[j]->GetEntry(i);
      nexp_bkgd->SetBinContent(j + 1, grid_background * GridScale_Normal);
      nexp_signal->SetBinContent(j + 1, grid_signal * GridScale_Normal);

      for (k = 0; k < Extrap.size(); k++) {
        GridTree_2_Normal[j][k]->GetEntry(i);
        nc[j][k] = grid_nc * GridScale_Normal;
        numucc[j][k] = grid_numucc * GridScale_Normal;
        bnuecc[j][k] = grid_bnuecc * GridScale_Normal;
        nutaucc[j][k] = grid_nutaucc * GridScale_Normal;
        sig[j][k] = grid_nue * GridScale_Normal;
      }
    }
    nexp->Add(nexp_bkgd, nexp_signal, 1, 1);
    ErrCalc->SetGridPred(nBins, nc, numucc, bnuecc, nutaucc, sig);

    chi2min = GetMinLikelihood(grid_delta, true);

    ErrCalc->SetUseGrid(true); // use grid predictions set up above
    chi2data = 1e10;
    if (FitMethod == 0) {
      chi2data = PoissonChi2(nexp) - chi2min;
    } else if (FitMethod == 1) {
      chi2data = ScaledChi2(Bkgd, Sig) - chi2min;
    } else if (FitMethod == 2) {
      chi2data = StandardChi2(nexp) - chi2min;
    } else if (FitMethod == 3) {
      // Likelihood: "Standard" (N syst, N nuisance)
      // Calculate the likelihood (x2 for chi)
      Bkgd->Reset();
      Bkgd->Add(nexp_bkgd);
      Sig->Reset();
      Sig->Add(nexp_signal);
      chi2data = StandardLikelihood() - chi2min;
    } else if (FitMethod == 4) {
      // Likelihood: Bin by Bin Calculation of Systematics
      // Calculate the likelihood (x2 for chi)
      Bkgd->Reset();
      Bkgd->Add(nexp_bkgd);
      Sig->Reset();
      Sig->Add(nexp_signal);
      chi2data = BinLikelihood() - chi2min;
    } else {
      cout << "Error in RunMultiBinFC(): Unknown 'FitMethod'." << endl;
    }

    myfile << grid_sinsq2th13 << " " << grid_delta << " " << nexp_signal->Integral() << " ";
    for (j = 0; j < nBins; j++) {
      myfile << nexp_signal->GetBinContent(j + 1) << " ";
    }
    myfile << chi2data << endl;
  }

  myfile.close();

  // inverted hierarchy

  myfile.open(gSystem->ExpandPathName(fileinvt.c_str()));

  for (i = 0; i < nPts_Inverted; i++) {
    if (i % 100 == 0) {
      cout << 100. * i / nPts_Inverted << "% complete for inverted hierarchy" << endl;
    }

    nexp_bkgd->Reset();
    nexp_signal->Reset();
    nexp->Reset();

    for (j = 0; j < nBins; j++) {
      GridTree_Inverted[j]->GetEntry(i);
      nexp_bkgd->SetBinContent(j + 1, grid_background * GridScale_Inverted);
      nexp_signal->SetBinContent(j + 1, grid_signal * GridScale_Inverted);

      for (k = 0; k < Extrap.size(); k++) {
        GridTree_2_Inverted[j][k]->GetEntry(i);
        nc[j][k] = grid_nc * GridScale_Inverted;
        numucc[j][k] = grid_numucc * GridScale_Inverted;
        bnuecc[j][k] = grid_bnuecc * GridScale_Inverted;
        nutaucc[j][k] = grid_nutaucc * GridScale_Inverted;
        sig[j][k] = grid_nue * GridScale_Inverted;
      }
    }
    nexp->Add(nexp_bkgd, nexp_signal, 1, 1);
    ErrCalc->SetGridPred(nBins, nc, numucc, bnuecc, nutaucc, sig);

    chi2min = GetMinLikelihood(grid_delta, false);

    ErrCalc->SetUseGrid(true); // use grid predictions set up above
    chi2data = 1e10;
    if (FitMethod == 0) {
      chi2data = PoissonChi2(nexp) - chi2min;
    } else if (FitMethod == 1) {
      chi2data = ScaledChi2(Bkgd, Sig) - chi2min;
    } else if (FitMethod == 2) {
      chi2data = StandardChi2(nexp) - chi2min;
    } else if (FitMethod == 3) {
      // Likelihood: "Standard" (N syst, N nuisance)
      // Calculate the likelihood (x2 for chi)
      Bkgd->Reset();
      Bkgd->Add(nexp_bkgd);
      Sig->Reset();
      Sig->Add(nexp_signal);
      chi2data = StandardLikelihood() - chi2min;
    } else if (FitMethod == 4) {
      // Likelihood: Bin by Bin Calculation of Systematics
      // Calculate the likelihood (x2 for chi)
      Bkgd->Reset();
      Bkgd->Add(nexp_bkgd);
      Sig->Reset();
      Sig->Add(nexp_signal);
      chi2data = BinLikelihood() - chi2min;
    } else {
      cout << "Error in RunMultiBinFC(): Unknown 'FitMethod'." << endl;
    }

    myfile << grid_sinsq2th13 << " " << grid_delta << " " << nexp_signal->Integral() << " ";
    for (j = 0; j < nBins; j++) {
      myfile << nexp_signal->GetBinContent(j + 1) << " ";
    }
    myfile << chi2data << endl;
  }

  myfile.close();

  return;
}
void NueFit2D::RunPseudoExperiments() {
  if (FracErr_Bkgd == 0) {
    cout << "FracErr_Bkgd not set.  Quitting..." << endl;
    return;
  }
  if (FracErr_Sig == 0) {
    cout << "FracErr_Sig not set.  Quitting..." << endl;
    return;
  }
  if (Extrap.size() == 0) {
    cout << "No Extrapolate2D input.  Quitting..." << endl;
    return;
  }

  gRandom->SetSeed(0);

  ReadGridFiles();

  if (nPts_Normal == 0 || nPts_Inverted == 0) {
    return;
  }

  int i, k;
  unsigned int j;
  TH1D * nexp_bkgd = new TH1D("nexp_bkgd", "", nBins, -0.5, nBins - 0.5);
  TH1D * nexp_signal = new TH1D("nexp_signal", "", nBins, -0.5, nBins - 0.5);
  TH1D * dnexp_oscpar = new TH1D("dnexp_oscpar", "", nBins, -0.5, nBins - 0.5);
  TH1D * nexp = new TH1D("nexp", "", nBins, -0.5, nBins - 0.5);
  double delchi2;
  bool calcchi2min; // if nobs<nexpmin for any bin, then chi2min!=0 and needs to be calculated
  TH1D * nexpmin = 0;
  TH1D * chi2hist = new TH1D("chi2hist", "", 110000, -10, 100);

  // normal hierarchy
  for (j = 0; j < Extrap.size(); j++) {
    Extrap[j]->GetPrediction();
    Extrap[j]->SetSinSq2Th13(0);
    Extrap[j]->OscillatePrediction();
    if (j == 0) {
      nexpmin = (TH1D *) Extrap[j]->Pred_TotalBkgd_VsBinNumber->Clone("nexpmin");
      nexpmin->Add(Extrap[j]->Pred_Signal_VsBinNumber);
    } else {
      nexpmin->Add(Extrap[j]->Pred_TotalBkgd_VsBinNumber);
      nexpmin->Add(Extrap[j]->Pred_Signal_VsBinNumber);
    }
  }

  TFile * f = new TFile(gSystem->ExpandPathName(outFileName.c_str()), "RECREATE");

  for (i = 0; i < nPts_Normal; i++) {
    if (i % 100 == 0) {
      cout << 100. * i / nPts_Normal << "% complete for normal hierarchy" << endl;
    }

    nexp_bkgd->Reset();
    nexp_signal->Reset();
    nexp->Reset();
    dnexp_oscpar->Reset();
    for (j = 0; j < nBins; j++) {
      GridTree_Normal[j]->GetEntry(i);
      nexp_bkgd->SetBinContent(j + 1, grid_background * GridScale_Normal);
      nexp_signal->SetBinContent(j + 1, grid_signal * GridScale_Normal);
      dnexp_oscpar->SetBinContent(j + 1, grid_oscparerr);
    }
    nexp->Add(nexp_bkgd, nexp_signal, 1, 1);

    chi2hist->Reset();
    chi2hist->SetName(Form("Chi2Hist_Normal_%i", i));

    for (k = 0; k < NumExpts; k++) {
      if (IncludeOscParErrs) {
        GenerateOneExperiment(nexp_bkgd, nexp_signal, dnexp_oscpar);
      } else {    GenerateOneExperiment(nexp_bkgd, nexp_signal); }
      delchi2 = PoissonChi2(nexp);
      calcchi2min = false;
      for (j = 0; j < nBins; j++) {
        if (NObs->GetBinContent(j + 1) < nexpmin->GetBinContent(j + 1)) {
          calcchi2min = true;
        }
      }
      if (calcchi2min) {
        delchi2 = delchi2 - PoissonChi2(nexpmin);
      }
      chi2hist->Fill(delchi2);
    }
    chi2hist->Write();
    f->Close();
    f = new TFile(gSystem->ExpandPathName(outFileName.c_str()), "UPDATE");
  }

  for (j = 0; j < nBins; j++) {
    GridTree_Normal[j]->Write();
  }

  f->Close();

  // inverted hierarchy
  for (j = 0; j < Extrap.size(); j++) {
    Extrap[j]->SetSinSq2Th13(0);
    Extrap[j]->InvertMassHierarchy();
    Extrap[j]->OscillatePrediction();
    if (j == 0) {
      nexpmin = (TH1D *) Extrap[j]->Pred_TotalBkgd_VsBinNumber->Clone("nexpmin");
      nexpmin->Add(Extrap[j]->Pred_Signal_VsBinNumber);
    } else {
      nexpmin->Add(Extrap[j]->Pred_TotalBkgd_VsBinNumber);
      nexpmin->Add(Extrap[j]->Pred_Signal_VsBinNumber);
    }
  }

  f = new TFile(gSystem->ExpandPathName(outFileName.c_str()), "UPDATE");

  for (i = 0; i < nPts_Inverted; i++) {
    if (i % 100 == 0) {
      cout << 100. * i / nPts_Inverted << "% complete for inverted hierarchy" << endl;
    }

    nexp_bkgd->Reset();
    nexp_signal->Reset();
    nexp->Reset();
    dnexp_oscpar->Reset();
    for (j = 0; j < nBins; j++) {
      GridTree_Inverted[j]->GetEntry(i);
      nexp_bkgd->SetBinContent(j + 1, grid_background * GridScale_Inverted);
      nexp_signal->SetBinContent(j + 1, grid_signal * GridScale_Inverted);
      dnexp_oscpar->SetBinContent(j + 1, grid_oscparerr);
    }
    nexp->Add(nexp_bkgd, nexp_signal, 1, 1);

    chi2hist->Reset();
    chi2hist->SetName(Form("Chi2Hist_Inverted_%i", i));

    for (k = 0; k < NumExpts; k++) {
      if (IncludeOscParErrs) {
        GenerateOneExperiment(nexp_bkgd, nexp_signal, dnexp_oscpar);
      } else {    GenerateOneExperiment(nexp_bkgd, nexp_signal); }
      delchi2 = PoissonChi2(nexp);
      calcchi2min = false;
      for (j = 0; j < nBins; j++) {
        if (NObs->GetBinContent(j + 1) < nexpmin->GetBinContent(j + 1)) {
          calcchi2min = true;
        }
      }
      if (calcchi2min) {
        delchi2 = delchi2 - PoissonChi2(nexpmin);
      }
      chi2hist->Fill(delchi2);
    }
    chi2hist->Write();
    f->Close();
    f = new TFile(gSystem->ExpandPathName(outFileName.c_str()), "UPDATE");
  }

  for (j = 0; j < nBins; j++) {
    GridTree_Inverted[j]->Write();
  }

  f->Close();

  return;
}

void NueFit2D::RunMultiBinPseudoExpts(bool Print) {
  if (FitMethod == 1 && (FracErr_Bkgd == 0 || FracErr_Sig == 0)) {
    cout << "FracErr_Bkgd and FracErr_Sig need to be set for ScaledChi2.  Quitting..." << endl;
    return;
  }
  if (Extrap.size() == 0) {
    cout << "No Extrapolate2D input.  Quitting..." << endl;
    return;
  }
  for (unsigned int ie = 0; ie < Extrap.size(); ie++) {
    Extrap[ie]->GetPrediction();
  }
  if (ErrCalc == 0) {
    cout << "Need to set ErrorCalc object!  Quitting..." << endl;
    return;
  }
  if (FitMethod == 3) {
    DefineStdDlnLMinuit();
  }
  if (FitMethod == 4) {
    DefineBinDlnLMinuit();
  }

  nBins = Extrap[0]->Pred_TotalBkgd_VsBinNumber->GetNbinsX();
  if (ErrorMatrix == 0) {
    ErrorMatrix = new TH2D("ErrorMatrix", "", nBins, -0.5, nBins - 0.5, nBins, -0.5, nBins - 0.5);
  }
  if (InvErrorMatrix == 0) {
    InvErrorMatrix = new TH2D("InvErrorMatrix", "", nBins, -0.5, nBins - 0.5, nBins, -0.5, nBins - 0.5);
  }

  TH2D * Error4Expts = new TH2D("Error4Expts", "", nBins, -0.5, nBins - 0.5, nBins, -0.5, nBins - 0.5);

  ReadGridFiles();

  if (nPts_Normal == 0 || nPts_Inverted == 0) {
    return;
  }

  gRandom->SetSeed(0);

  int i, u;
  unsigned int j, k;
  TH1D * nexp_bkgd = new TH1D("nexp_bkgd", "", nBins, -0.5, nBins - 0.5);
  TH1D * nexp_signal = new TH1D("nexp_signal", "", nBins, -0.5, nBins - 0.5);
  TH1D * nexp = new TH1D("nexp", "", nBins, -0.5, nBins - 0.5);
  double delchi2, chi2min;
  TH1D * chi2hist = new TH1D("chi2hist", "", 110000, -10, 100);
  double ele;
  int noff;

  vector< vector<double> > nc, numucc, bnuecc, nutaucc, sig;
  for (j = 0; j < nBins; j++) {
    nc.push_back(vector<double>() );
    numucc.push_back(vector<double>() );
    bnuecc.push_back(vector<double>() );
    nutaucc.push_back(vector<double>() );
    sig.push_back(vector<double>() );
    for (k = 0; k < Extrap.size(); k++) {
      nc[j].push_back(0);
      numucc[j].push_back(0);
      bnuecc[j].push_back(0);
      nutaucc[j].push_back(0);
      sig[j].push_back(0);
    }
  }

  Bkgd = (TH1D *) nexp->Clone("Bkgd");
  Bkgd->Reset();
  Sig = (TH1D *) nexp->Clone("Sig");
  Sig->Reset();

  ofstream myfile;
  string file, ofile;

  // normal hierarchy

  TFile * f = new TFile(gSystem->ExpandPathName(outFileName.c_str()), "RECREATE");

  if (Print) {
    ofile = gSystem->ExpandPathName(outFileName.c_str());
    file = ofile.substr(0, ofile.length() - 5) + "_Normal.dat";
    myfile.open(gSystem->ExpandPathName(file.c_str()));
  }

  for (i = 0; i < nPts_Normal; i++) {
    cout << "point " << (i + 1) << "/" << nPts_Normal << " (normal hierarchy)" << endl;

    nexp_bkgd->Reset();
    nexp_signal->Reset();
    nexp->Reset();

    for (j = 0; j < nBins; j++) {
      GridTree_Normal[j]->GetEntry(i);
      nexp_bkgd->SetBinContent(j + 1, grid_background * GridScale_Normal);
      nexp_signal->SetBinContent(j + 1, grid_signal * GridScale_Normal);

      for (k = 0; k < Extrap.size(); k++) {
        GridTree_2_Normal[j][k]->GetEntry(i);
        nc[j][k] = grid_nc * GridScale_Normal;
        numucc[j][k] = grid_numucc * GridScale_Normal;
        bnuecc[j][k] = grid_bnuecc * GridScale_Normal;
        nutaucc[j][k] = grid_nutaucc * GridScale_Normal;
        sig[j][k] = grid_nue * GridScale_Normal;
      }
    }
    nexp->Add(nexp_bkgd, nexp_signal, 1, 1);
    ErrCalc->SetGridPred(nBins, nc, numucc, bnuecc, nutaucc, sig);

    Error4Expts->Reset();
    ErrCalc->SetUseGrid(true);
    ErrCalc->CalculateSystErrorMatrix();
    Error4Expts->Add(ErrCalc->CovMatrix);
    ErrCalc->CalculateHOOError();
    Error4Expts->Add(ErrCalc->CovMatrix_Decomp);
    if (IncludeOscParErrs) {
      noff = 0;
      for (j = 0; j < nBins; j++) {
        ele = Error4Expts->GetBinContent(j + 1, j + 1);
        ele += (grid_bin_oscparerr[j] * grid_bin_oscparerr[j] * nexp->GetBinContent(j + 1) * nexp->GetBinContent(j + 1));
        Error4Expts->SetBinContent(j + 1, j + 1, ele);

        for (k = 0; k < nBins; k++) {
          if (k > j) {
            ele = Error4Expts->GetBinContent(j + 1, k + 1);
            ele += (grid_bin_oscparerr[j] * grid_bin_oscparerr[k] * nexp->GetBinContent(j + 1) * nexp->GetBinContent(k + 1));
            Error4Expts->SetBinContent(j + 1, k + 1, ele);

            ele = Error4Expts->GetBinContent(k + 1, j + 1);
            ele += (grid_bin_oscparerr[j] * grid_bin_oscparerr[k] * nexp->GetBinContent(j + 1) * nexp->GetBinContent(k + 1));
            Error4Expts->SetBinContent(k + 1, j + 1, ele);

            noff++;
          }
        }
      }
    }

    chi2hist->Reset();
    chi2hist->SetName(Form("Chi2Hist_Normal_%i", i));

    for (u = 0; u < NumExpts; u++) {
      cout << "expt " << (u + 1) << "/" << NumExpts << endl;

      GenerateOneCorrelatedExp(nexp, Error4Expts);
      if (Print) {
        myfile << grid_sinsq2th13 << " " << grid_delta << " ";
        for (j = 0; j < nBins; j++) {
          myfile << NObs->GetBinContent(j + 1) << " ";
        }
        myfile << endl;
      }

      chi2min = GetMinLikelihood(grid_delta, true);

      ErrCalc->SetUseGrid(true); // will use the grid predictions set above
      delchi2 = 1e10;
      if (FitMethod == 0) {
        delchi2 = PoissonChi2(nexp) - chi2min;
      } else if (FitMethod == 1) {
        delchi2 = ScaledChi2(nexp_bkgd, nexp_signal) - chi2min;
      } else if (FitMethod == 2) {
        delchi2 = StandardChi2(nexp) - chi2min;
      } else if (FitMethod == 3) {
        // Likelihood: "Standard" (N syst, N nuisance)
        // Calculate the likelihood (x2 for chi)
        Bkgd->Reset();
        Bkgd->Add(nexp_bkgd);
        Sig->Reset();
        Sig->Add(nexp_signal);
        delchi2 = StandardLikelihood() - chi2min;
      } else if (FitMethod == 4) {
        // Likelihood: Bin by Bin Calculation of Systematics
        // Calculate the likelihood (x2 for chi)
        Bkgd->Reset();
        Bkgd->Add(nexp_bkgd);
        Sig->Reset();
        Sig->Add(nexp_signal);
        delchi2 = BinLikelihood() - chi2min;
      } else {
        cout << "Error in RunMultiBinPseudoExpts(): Unknown 'FitMethod'." << endl;
      }
      chi2hist->Fill(delchi2);

    }
    f->cd();
    chi2hist->Write();
    f->Close();

    f = new TFile(gSystem->ExpandPathName(outFileName.c_str()), "UPDATE");
  }

  if (Print) {
    myfile.close();
  }

  for (j = 0; j < nBins; j++) {
    GridTree_Normal[j]->Write();

    for (k = 0; k < Extrap.size(); k++) {
      GridTree_2_Normal[j][k]->Write();
    }
  }

  f->Close();

  // inverted hierarchy

  f = new TFile(gSystem->ExpandPathName(outFileName.c_str()), "UPDATE");

  if (Print) {
    ofile = gSystem->ExpandPathName(outFileName.c_str());
    file = ofile.substr(0, ofile.length() - 5) + "_Inverted.dat";
    myfile.open(gSystem->ExpandPathName(file.c_str()));
  }

  for (i = 0; i < nPts_Inverted; i++) {
    cout << "point " << (i + 1) << "/" << nPts_Inverted << " (inverted hierarchy)" << endl;

    nexp_bkgd->Reset();
    nexp_signal->Reset();
    nexp->Reset();

    for (j = 0; j < nBins; j++) {
      GridTree_Inverted[j]->GetEntry(i);
      nexp_bkgd->SetBinContent(j + 1, grid_background * GridScale_Inverted);
      nexp_signal->SetBinContent(j + 1, grid_signal * GridScale_Inverted);

      for (k = 0; k < Extrap.size(); k++) {
        GridTree_2_Inverted[j][k]->GetEntry(i);
        nc[j][k] = grid_nc * GridScale_Inverted;
        numucc[j][k] = grid_numucc * GridScale_Inverted;
        bnuecc[j][k] = grid_bnuecc * GridScale_Inverted;
        nutaucc[j][k] = grid_nutaucc * GridScale_Inverted;
        sig[j][k] = grid_nue * GridScale_Inverted;
      }
    }
    nexp->Add(nexp_bkgd, nexp_signal, 1, 1);
    ErrCalc->SetGridPred(nBins, nc, numucc, bnuecc, nutaucc, sig);

    Error4Expts->Reset();
    ErrCalc->SetUseGrid(true);
    ErrCalc->CalculateSystErrorMatrix();
    Error4Expts->Add(ErrCalc->CovMatrix);
    ErrCalc->CalculateHOOError();
    Error4Expts->Add(ErrCalc->CovMatrix_Decomp);
    if (IncludeOscParErrs) {
      noff = 0;
      for (j = 0; j < nBins; j++) {
        ele = Error4Expts->GetBinContent(j + 1, j + 1);
        ele += (grid_bin_oscparerr[j] * grid_bin_oscparerr[j] * nexp->GetBinContent(j + 1) * nexp->GetBinContent(j + 1));
        Error4Expts->SetBinContent(j + 1, j + 1, ele);

        for (k = 0; k < nBins; k++) {
          if (k > j) {
            ele = Error4Expts->GetBinContent(j + 1, k + 1);
            ele += (grid_bin_oscparerr[j] * grid_bin_oscparerr[k] * nexp->GetBinContent(j + 1) * nexp->GetBinContent(k + 1));
            Error4Expts->SetBinContent(j + 1, k + 1, ele);

            ele = Error4Expts->GetBinContent(k + 1, j + 1);
            ele += (grid_bin_oscparerr[j] * grid_bin_oscparerr[k] * nexp->GetBinContent(j + 1) * nexp->GetBinContent(k + 1));
            Error4Expts->SetBinContent(k + 1, j + 1, ele);

            noff++;
          }
        }
      }
    }

    chi2hist->Reset();
    chi2hist->SetName(Form("Chi2Hist_Inverted_%i", i));

    for (u = 0; u < NumExpts; u++) {
      cout << "expt " << (u + 1) << "/" << NumExpts << endl;

      GenerateOneCorrelatedExp(nexp, Error4Expts);
      if (Print) {
        myfile << grid_sinsq2th13 << " " << grid_delta << " ";
        for (j = 0; j < nBins; j++) {
          myfile << NObs->GetBinContent(j + 1) << " ";
        }
        myfile << endl;
      }

      chi2min = GetMinLikelihood(grid_delta, false);

      ErrCalc->SetUseGrid(true); // will use the grid predictions set above
      delchi2 = 1e10;
      if (FitMethod == 0) {
        delchi2 = PoissonChi2(nexp) - chi2min;
      } else if (FitMethod == 1) {
        delchi2 = ScaledChi2(nexp_bkgd, nexp_signal) - chi2min;
      } else if (FitMethod == 2) {
        delchi2 = StandardChi2(nexp) - chi2min;
      } else if (FitMethod == 3) {
        // Likelihood: "Standard" (N syst, N nuisance)
        // Calculate the likelihood (x2 for chi)
        Bkgd->Reset();
        Bkgd->Add(nexp_bkgd);
        Sig->Reset();
        Sig->Add(nexp_signal);
        delchi2 = StandardLikelihood() - chi2min;
      } else if (FitMethod == 4) {
        // Likelihood: Bin by Bin Calculation of Systematics
        // Calculate the likelihood (x2 for chi)
        Bkgd->Reset();
        Bkgd->Add(nexp_bkgd);
        Sig->Reset();
        Sig->Add(nexp_signal);
        delchi2 = BinLikelihood() - chi2min;
      } else {
        cout << "Error in RunMultiBinPseudoExpts(): Unknown 'FitMethod'." << endl;
      }
      chi2hist->Fill(delchi2);
    }
    f->cd();
    chi2hist->Write();
    f->Close();
    f = new TFile(gSystem->ExpandPathName(outFileName.c_str()), "UPDATE");
  }

  if (Print) {
    myfile.close();
  }

  for (j = 0; j < nBins; j++) {
    GridTree_Inverted[j]->Write();
    for (k = 0; k < Extrap.size(); k++) {
      GridTree_2_Inverted[j][k]->Write();
    }
  }

  f->Close();

  return;
}
void NueFit2D::GenerateOneExperiment(TH1D * nexp_bkgd, TH1D * nexp_signal, TH1D * dnexp_oscpar) {
  NObs = (TH1D *) nexp_bkgd->Clone("NObs");
  NObs->Reset();
  nBins = NObs->GetNbinsX();

  int obs;
  double b, s, eb, es, bgauss, sgauss;
  double tot, etot;
  double temp;

  for (unsigned int i = 0; i < nBins; i++) {
    b = nexp_bkgd->GetBinContent(i + 1);
    s = nexp_signal->GetBinContent(i + 1);
    eb = FracErr_Bkgd->GetBinContent(i + 1) * nexp_bkgd->GetBinContent(i + 1);
    es = FracErr_Sig->GetBinContent(i + 1) * nexp_signal->GetBinContent(i + 1);

    bgauss = gRandom->Gaus(b, eb);
    sgauss = gRandom->Gaus(s, es);
    tot = bgauss + sgauss;

    if (IncludeOscParErrs) {
      etot = dnexp_oscpar->GetBinContent(i + 1) * tot;
      temp = gRandom->Gaus(tot, etot);
      tot = temp;
    }

    obs = gRandom->Poisson(tot);

    NObs->SetBinContent(i + 1, obs);
  }

  return;
}
void NueFit2D::GenerateOneCorrelatedExp(TH1D * nexp, TH2D * err) {
  NObs = (TH1D *) nexp->Clone("NObs");
  NObs->Reset();
  nBins = NObs->GetNbinsX();

  unsigned int i, j, k;

  // error matrix
  double * Varray = new double[nBins * nBins];
  for (i = 0; i < nBins; i++) {
    for (j = 0; j < nBins; j++) {
      k = i * nBins + j;
      Varray[k] = err->GetBinContent(i + 1, j + 1);
    }
  }
  const TMatrixD M(nBins, nBins, Varray);
//   cout<<"M:"<<endl;
//   for(i=0;i<nBins;i++)
//   {
//     for(j=0;j<nBins;j++)
//     {
//       cout<<M[i][j]<<" ";
//     }
//     cout<<endl;
//   }

  // get eigenvectors of M
  TMatrixDEigen V(M);
  // construct matrix of eigenvectors
  TMatrixD AT = V.GetEigenVectors();
//   cout<<"AT:"<<endl;
//   for(i=0;i<nBins;i++)
//   {
//     for(j=0;j<nBins;j++)
//     {
//       cout<<AT[i][j]<<" ";
//     }
//     cout<<endl;
//   }

  // transpose to get transformation matrix A
  TMatrixD A = AT;
  A.Transpose(A);
//   cout<<"A:"<<endl;
//   for(i=0;i<nBins;i++)
//   {
//     for(j=0;j<nBins;j++)
//     {
//       cout<<A[i][j]<<" ";
//     }
//     cout<<endl;
//   }

  // diagonalize M
  TMatrixD D(nBins, nBins);
  D = A * M * AT;
//   cout<<"D:"<<endl;
//   for(i=0;i<nBins;i++)
//   {
//     for(j=0;j<nBins;j++)
//     {
//       cout<<D[i][j]<<" ";
//     }
//     cout<<endl;
//   }

  // transform prediction to diagonal basis using A
  vector<double> nnew;
  double t;
  for (i = 0; i < nBins; i++) {
    t = 0;
    for (j = 0; j < nBins; j++) {
      t += A[i][j] * nexp->GetBinContent(j + 1);
    }
    nnew.push_back(t);
  }

  // randomize prediction in diagonal basis
  vector<double> nrand;
  for (i = 0; i < nBins; i++) {
//     cout<<"nexp_ = "<<nexp_->GetBinContent(i+1)<<" "<<sqrt(D[i][i])<<endl;
    if (D[i][i] < 0) {
      cout << "Warning in NueFit2D::GenerateOneCorrelatedExp(): Negative element in diagonalized systematic error matrix ... setting to 0." << endl;
      D[i][i] = 0;
    }
    nrand.push_back(gRandom->Gaus(nnew[i], sqrt(D[i][i])));
  }

  // transform randomized prediction back to original basis
  for (i = 0; i < nBins; i++) {
    t = 0;
    for (j = 0; j < nBins; j++) {
      t += AT[i][j] * nrand[j];
    }
    t = gRandom->Poisson(t);
    NObs->SetBinContent(i + 1, t);
//     cout<<t<<endl;
  }

  delete [] Varray;

  return;
}
double NueFit2D::ScaledChi2(TH1D * nexp_bkgd, TH1D * nexp_signal) {
  double chi2 = 0;
  double nb, ns, nobs, eb, es, nexp, errscale;

  for (unsigned int i = 0; i < nBins; i++) {
    nb = nexp_bkgd->GetBinContent(i + 1);
    ns = nexp_signal->GetBinContent(i + 1);
    nexp = nb + ns;
    nobs = NObs->GetBinContent(i + 1);
    eb = FracErr_Bkgd->GetBinContent(i + 1);
    es = FracErr_Sig->GetBinContent(i + 1);
    errscale = nexp / (eb * eb * nb * nb + es * es * ns * ns + nexp);
    if (nobs > 0) {
      chi2 += (2 * (nexp - nobs + nobs * TMath::Log(nobs / nexp)) * errscale);
    } else if (nobs == 0) {
      chi2 += (2 * nexp * errscale);
    }
    // chi2 is undefined if nexp is 0
  }

  return chi2;
}
double NueFit2D::PoissonChi2(TH1D * nexp) {
  double chi2 = 0;
  double no, ne;

  for (unsigned int i = 0; i < nBins; i++) {
    ne = nexp->GetBinContent(i + 1);
    no = NObs->GetBinContent(i + 1);
    if (no > 0) {
      chi2 += (2 * (ne - no + no * TMath::Log(no / ne)));
    } else if (no == 0) {
      chi2 += (2 * ne);
    }
    // chi2 is undefined if ne is 0
  }

  return chi2;
}
double NueFit2D::StandardChi2(TH1D * nexp) {
  double chi2 = 0;

  CalculateErrorMatrixInv(nexp);
  for (unsigned int i = 0; i < nBins; i++) {
    for (unsigned int j = 0; j < nBins; j++) {
      chi2 += ((NObs->GetBinContent(i + 1) - nexp->GetBinContent(i + 1)) * (NObs->GetBinContent(j + 1) - nexp->GetBinContent(j + 1)) * InvErrorMatrix->GetBinContent(i + 1, j + 1));
    }
  }

  return chi2;
}
void NueFit2D::CalculateErrorMatrixInv(TH1D * nexp) {
  ErrorMatrix->Reset();
  InvErrorMatrix->Reset();

  if (ErrCalc != 0) { // if ErrorCalc object has been added, use it to calculate covariance matrix
    ErrCalc->CalculateSystErrorMatrix();
    ErrorMatrix->Add(ErrCalc->CovMatrix);
    ErrCalc->CalculateHOOError();
    ErrorMatrix->Add(ErrCalc->CovMatrix_Decomp);
  } else if (ExternalErrorMatrix != 0) { // if setting a systematic matrix externally, use it
    ErrorMatrix->Add(ExternalErrorMatrix);
  }
  // otherwise, it'll be statistics only

  unsigned int i, j, k;
  double ele;
  for (j = 0; j < nBins; j++) {
    ele = ErrorMatrix->GetBinContent(j + 1, j + 1);
    ele += nexp->GetBinContent(j + 1);
    ErrorMatrix->SetBinContent(j + 1, j + 1, ele);
  }

  // take CovMatrix and make array appropiate for TMatrix
  double * Varray = new double[nBins * nBins];
  for (i = 0; i < nBins; i++) {
    for (j = 0; j < nBins; j++) {
      k = i * nBins + j;
      Varray[k] = ErrorMatrix->GetBinContent(i + 1, j + 1);
    }
  }

  // make TMatrix
  TMatrixD * Vmatrix = new TMatrixD(nBins, nBins, Varray);

  // get determinant
  //   Double_t det = Vmatrix->Determinant();
  //   cout<<"Covariance Matrix det = "<<det<<endl;

  // invert
  TMatrixD Vinvmatrix = Vmatrix->Invert();

  // make array out of inverse
  Double_t * Vinvarray = Vinvmatrix.GetMatrixArray();

  // make Vinv out of array
  for (i = 0; i < nBins; i++) {
    for (j = 0; j < nBins; j++) {
      InvErrorMatrix->SetBinContent(i + 1, j + 1, Vinvarray[i * nBins + j]);
    }
  }

  delete [] Varray;

  return;
}
double NueFit2D::EvaluateOmega(double signal, double background) {
  // this will only be called for a single bin fit

  int n0 = (int) NObs->GetBinContent(1, 1);
  double b0 = background;
  double s0 = signal;
  double k0 = 1.0;

  double errbkgd = FracErr_Bkgd->GetBinContent(1, 1);
  double errsig = FracErr_Sig->GetBinContent(1, 1);

  double errK = 1.0;

  if (s0 > 0) {
    errK = (errsig) * (errsig);
  }
  if (errK < 1e-5) {
    errK = 1.0;  std::cout << "can't do perfect signal" << std::endl;
  }
  double errBg = errbkgd * b0 * errbkgd * b0;

  double rank0 = CalculateRank(n0, s0, b0, errBg, k0, errK);

  double results[3];
  FindBestFit(n0, s0, b0, errBg, k0, errK, results);

  double betaHat = results[0];
  double kHat = results[1];

  errBg = (betaHat * errbkgd) * (betaHat * errbkgd);
  errK = (errsig) * (errsig) * kHat * kHat;

//   std::cout<<n0<<"  "<<s0<<"  "<<b0<<"  "<<k0<<"  "<<betaHat<<" "<<kHat<<"  "<<errBg<<" "<<errK<<std::endl;
//   std::cout<<s0*kHat + betaHat<<"  "<<n0<<std::endl;

  double omega = 0;
  double omegaBar = Poisson(s0 * kHat + betaHat, n0);

  int nBase = int(s0 * kHat + betaHat);

  int nHigh = nBase;
  int nLow  = nBase - 1;

  double delta = 1.0;
  bool filled = false;

  const int NUMB = 30;
  const int NUMK = 30;

  static double bVal[NUMB];
  static double kVal[NUMK];
  static double errBVar[NUMB];
  static double errKVar[NUMK];

  static bool first = true;
  static double scale = 1.0;

  double bStart = 0.4 * b0;   double bStop = 2.2 * b0;
  double kStart = 0.4;   double kStop = 1.6;

  if (first) {
    for (int i = 0; i < NUMB; i++) {
      bVal[i] = bStart + i * (bStop - bStart) / float(NUMB);
      errBVar[i] = errbkgd * errbkgd * bVal[i] * bVal[i];
    }

    for (int i = 0; i < NUMK; i++) {
      kVal[i] = kStart + i * (kStop - kStart) / float(NUMK);
      errKVar[i] = (errsig) * (errsig) * kVal[i] * kVal[i];
    }
    first = false;
    scale = (bStop - bStart) * (kStop - kStart) / (NUMK * NUMB);
  }
  double bkProb[NUMB][NUMK];

  /*   so lets be clear about this, the values for these gaussians are the same
     but its the rank that changes for any given value of n (the values themselve differ with mu)
     so the first time through i calculate the contribution for each point in the space
     then when looping through just look it up in a giant array             */

  // better yet - I can save cycles by building the arrays separately on the first pass

  double bGauss[NUMB];
  double kGauss[NUMK];

  for (int i = 0; i < NUMB; i++) {
    bGauss[i] = Gaussian(betaHat, bVal[i], errBg);
  }
  for (int i = 0; i < NUMB; i++) {
    kGauss[i] = Gaussian(kHat, kVal[i], errK);
  }

  bool risingH = false, risingL = false;
  bool doneH = false;  // some calc savers
  bool doneL = false;  // some calc savers

  double ThreshHold = 1e-5;
  double pfThresh = ThreshHold * 1e-2;

  double slip = 0;   // Error on the numerical integral

  while (delta > ThreshHold || (1 - (omega + omegaBar) > 2 * ThreshHold) ) {
    if (nHigh == n0) {
      nHigh++;
    }
    if (nLow == n0) {
      nLow--;
    }

    double dOmegaH = 0.0, dOmegaBarH = 0.0;
    double dOmegaL = 0.0, dOmegaBarL = 0.0;

    double PrefactorH = Poisson(s0 * kHat + betaHat, nHigh);
    double PrefactorL = Poisson(s0 * kHat + betaHat, nLow);

    if (PrefactorH > pfThresh) {
      risingH = true;
    }
    if (PrefactorL > pfThresh) {
      risingL = true;
    }
    if (PrefactorH < pfThresh && risingH) {
      doneH = true;
    }
    if (PrefactorL < pfThresh && risingL) {
      doneL = true;
    }
    if (doneH && doneL) {
      // just a sanity check that the code doesn't think its done too early, this will cause an inf loop
      if (1 - (omega + omegaBar) < 1.5 * slip) {
        break;
      }
      std::cout << "Thats unexpected:  " << 1 - (omega + omegaBar) << "  " << slip << std::endl;
    }

    if (PrefactorH > pfThresh || PrefactorL > pfThresh) { // No need to loop if contribution too small
      for (int i = 0; i < NUMB; i++) { //  0 to infinity
        double b = bVal[i];
        double eb = errBVar[i];
        for (int j = 0; j < NUMK; j++) { //  -inf to inf
          double k = kVal[j];
          if (!filled) {
            bkProb[i][j] = bGauss[i] * kGauss[j];
          }

          double val = bkProb[i][j];
          double eK = errKVar[j];
          if (val < ThreshHold * 1e-4) {
            continue;                         // no point in using these contributions

          }
          double rank = 1.0;

          if (PrefactorH > pfThresh) {
            rank = CalculateRank(nHigh, s0, b, eb, k, eK);
            if (rank > rank0) {
              dOmegaH += val * scale;
            } else {                   dOmegaBarH += val * scale; }
          }

          if (PrefactorL > pfThresh) {
            if (nLow >= 0) {
              rank = CalculateRank(nLow, s0, b, eb, k, eK);
              if (rank > rank0) {
                dOmegaL += val * scale;
              } else {                     dOmegaBarL += val * scale; }
            }
          }
        }
      }
      if (!filled) {
        filled = true;
      }
      if (PrefactorH > pfThresh && 1 - (dOmegaH + dOmegaBarH) > slip) {
        slip = 1 - (dOmegaH + dOmegaBarH);
      }
      if (PrefactorL > pfThresh && 1 - (dOmegaL + dOmegaBarL) > slip) {
        slip = 1 - (dOmegaL + dOmegaBarL);
      }
    }
    delta = TMath::Max(PrefactorH * dOmegaH + PrefactorL * dOmegaL, PrefactorH * dOmegaBarH + PrefactorL * dOmegaBarL);
    omega += PrefactorH * dOmegaH + PrefactorL * dOmegaL;
    omegaBar += PrefactorH * dOmegaBarH + PrefactorL * dOmegaBarL;
    nHigh++; nLow--;
  }

  if (slip > ThreshHold) {
    std::cout << "Integral error past threshold:  " << omega << "  " << omegaBar << "  " << omega + omegaBar << "  " << slip << std::endl;
    std::cout << n0 << "  " << s0 << "  " << b0 << "  " << k0 << "  " << betaHat << " " << kHat << "  " << errBg << " " << errK << std::endl;
    std::cout << s0 * kHat + betaHat << "  " << n0 << std::endl;
  }

  return omega;
}
double NueFit2D::CalculateRank(int n, double s, double b, double errBg, double k, double errK) {
  double results[3];

  // Calculate the numerator
  FindBestFit(n, s, b, errBg, k, errK, results);

//   double numerator = results[2];

  double sBest, kBest, bBest;

  // Now calculate the denominator
  if (n >= b) { // then we have n = sk+b, Beta = b; k = k0// don't have to do anything
    sBest = (n - b) / k; kBest = k; bBest = b;
  } else { // if(n < b) {
    bBest = 0.5 * (b - errBg + TMath::Sqrt( (b - errBg) * (b - errBg) + 4 * n * errBg));
    // then we have no signal at best fit point
    sBest = 0;  kBest = k;
  }

  double betaHat = results[0];
  double kHat = results[1];
  double logN = n * TMath::Log(s * kHat + betaHat) - (s * kHat + betaHat) - (betaHat - b) * (betaHat - b) / (2 * errBg)
    - (kHat - k) * (kHat - k) / (2 * errK);
  double logD = n * TMath::Log(sBest * kBest + bBest) - (sBest * kBest + bBest) - (bBest - b) * (bBest - b) / (2 * errBg)
    - (kBest - k) * (kBest - k) / (2 * errK);

  double rank = 1.0;
  if (logD != 0) {
    rank = TMath::Exp(logN - logD);
  } else { std::cout << "odd" << std::endl; }

  return rank;
}
bool NueFit2D::FindBestFit(int n,  double s, double b, double errBg, double k, double errK, double * res) {
//  At the moment this function is solvable in ~closed form (we ignore the dependance of errBg, errK on khat, betaHat
//  this lets us avoid having a numerical minimization, but the second half of this function is the code that allows for that

  double sol_A = 1 + errK / errBg * s * s;
  double sol_B = errBg + s * k - 2 * errK * s * s * b / errBg - b + s * s * errK;
  double sol_C = -n * errBg + s * k * errBg - s * s * errK * b - s * k * b + s * s * b * b * errK / errBg;

  double rad = sol_B * sol_B - 4 * sol_A * sol_C;

  if (rad < 0) {
    std::cout << "negative radical - not sure about this...." << std::endl;
  }

  double betaHat = (-sol_B + TMath::Sqrt(rad)) / (2 * sol_A);

  // solving for kHat given betaHat
  double kHat = k - s * errK / errBg * (b - betaHat);

  res[0] = betaHat;
  res[1] = kHat;

/*
   fFixedFitPar[0] = n;
   fFixedFitPar[1] = s;
   fFixedFitPar[2] = b;
   fFixedFitPar[3] = errBg;
   fFixedFitPar[4] = k;
   fFixedFitPar[5] = errK;
 */
//   res[2] = fResult = Poisson(s*kHat+betaHat, n, true)*Gaussian(betaHat, b, errBg)*Gaussian(kHat, k, errK);
/*
   static bool first = true;
   static TMinuit *min = 0;

   if(first){
   min = new TMinuit(2);
   gMinuit = min;
   min->SetFCN(WrapperFunction);
   min->SetObjectFit(this);
   min->DefineParameter(0,"bfit",fObserved,1, 0,60);
   min->DefineParameter(1,"kfit",1,0.99, 0,22);

   const double ERRDEF=1.;
   min->SetErrorDef(ERRDEF);
   min->SetPrintLevel(-1);
   first = false;
   std::cout<<"Built"<<std::endl;
   }
   Double_t arglist[10];
   Int_t ierflg = 0;

   arglist[0] = 10000;       //max calls
   arglist[1] = 0.01;         //tolerance

   min->mnexcm("SIMPLEX",arglist,2,ierflg);

   double errs[2];
   for(int i=0;i<2;i++){
   min->GetParameter(i,res[i],errs[i]);
   //     std::cout<<res[i]<<"  "<<errs[i]<<std::endl;
   }
 */
//   res[2] = fResult;

  return true;
}
double NueFit2D::GetSensitivityAt(double delta, bool normalhier) {

  if (NObs == 0) {
    cout << "NObs not set.  Quitting..." << endl;
    return 0;
  }
  if (FitMethod == 1 && (FracErr_Bkgd == 0 || FracErr_Sig == 0)) {
    cout << "FracErr_Bkgd and FracErr_Sig need to be set for ScaledChi2.  Quitting..." << endl;
    return 0;
  }
  if (Extrap.size() == 0) {
    cout << "No Extrapolate2D input.  Quitting..." << endl;
    return 0;
  }

  if (FitMethod == 3) {
    DefineStdDlnLMinuit();
  }
  if (FitMethod == 4) {
    DefineBinDlnLMinuit();
  }

  Bkgd = (TH1D *) NObs->Clone("Bkgd");
  Bkgd->Reset();
  Sig = (TH1D *) NObs->Clone("Sig");
  Sig->Reset();
  TH1D * NExp = (TH1D *) NObs->Clone("NExp");
  NExp->Reset();

  unsigned int ie;
  for (ie = 0; ie < Extrap.size(); ie++) {
    Extrap[ie]->GetPrediction();
  }

  Int_t i;

  if (ErrorMatrix == 0) {
    ErrorMatrix = new TH2D("ErrorMatrix", "", nBins, -0.5, nBins - 0.5, nBins, -0.5, nBins - 0.5);
  }
  if (InvErrorMatrix == 0) {
    InvErrorMatrix = new TH2D("InvErrorMatrix", "", nBins, -0.5, nBins - 0.5, nBins, -0.5, nBins - 0.5);
  }

  Double_t ss2th13 = 0;
  Double_t delchi2 = 0, delchi2prev = 0;
  Double_t contourlvl = 2.71;
  Double_t prev = 0;
  Double_t best = 0;
  Double_t Th13increment = 0;
  if (nSinSq2Th13Steps > 0) {
    Th13increment = (SinSq2Th13High - SinSq2Th13Low) / (nSinSq2Th13Steps);
  }
  double sens = -1;
  double chi2[3], val[3];
  TGraph * g3;
  TF1 * fit;
  double minchi2;
  double limit = 0;

  for (ie = 0; ie < Extrap.size(); ie++) {
    Extrap[ie]->SetDeltaCP(delta);
    if (!normalhier) {
      Extrap[ie]->InvertMassHierarchy();
    }
  }

  contourlvl = 2.71;

  best = -1.;
  minchi2 = 100;
  for (i = 0; i < 3; i++) {
    chi2[i] = -1.;
    val[i] = -1.;
  }
  for (i = 0; i < nSinSq2Th13Steps + 1; i++) {
    chi2[0] = chi2[1];
    chi2[1] = chi2[2];
    val[0] = val[1];
    val[1] = val[2];
    ss2th13 = i * Th13increment + SinSq2Th13Low;
    for (ie = 0; ie < Extrap.size(); ie++) {
      Extrap[ie]->SetSinSq2Th13(ss2th13);
      Extrap[ie]->OscillatePrediction();
    }
    Bkgd->Reset();
    Sig->Reset();
    NExp->Reset();
    for (ie = 0; ie < Extrap.size(); ie++) {
      Bkgd->Add(Extrap[ie]->Pred_TotalBkgd_VsBinNumber);
      Sig->Add(Extrap[ie]->Pred_Signal_VsBinNumber);
    }
    NExp->Add(Bkgd);
    NExp->Add(Sig);

    chi2[2] = 1e10;
    if (FitMethod == 0) {
      chi2[2] = PoissonChi2(NExp);
    } else if (FitMethod == 1) {
      chi2[2] = ScaledChi2(Bkgd, Sig);
    } else if (FitMethod == 2) {
      chi2[2] = StandardChi2(NExp);
    } else if (FitMethod == 3) {
      // Likelihood: "Standard" (N syst, N nuisance)
      // Calculate the likelihood (x2 for chi)
      chi2[2] = StandardLikelihood();
    } else if (FitMethod == 4) {
      // Likelihood: Bin by Bin Calculation of Systematics
      // Calculate the likelihood (x2 for chi)
      chi2[2] = BinLikelihood();
    } else {
      cout << "Error in GetSensitivityAt(): Unknown 'FitMethod'." << endl;
    }

    val[2] = ss2th13;

    if (i < 2) {
      continue;
    }

    if (i < 3 && chi2[2] > chi2[1] && chi2[1] > chi2[0]) { // first three points are increasing, first point is minimum.
      best = val[0];
      minchi2 = chi2[0];
      cout << "minimum at 1st point: " << minchi2 << " at " << best << endl;
      break;
    }

    if (chi2[2] > chi2[1] && chi2[0] > chi2[1]) { // found minimum
      g3 = new TGraph(3, val, chi2);
      fit = new TF1("pol2", "pol2");
      g3->Fit(fit, "Q"); // fit to second order polynominal
      if (fit->GetParameter(2) > 0) { // if the x^2 term is nonzero
        best = -fit->GetParameter(1) / (2 * fit->GetParameter(2)); // the location of the minimum is -p1/(2*p2)
        minchi2 = fit->GetParameter(0) + fit->GetParameter(1) * best + fit->GetParameter(2) * best * best;
        cout << "minimum with fit: " << minchi2 << " at " << best << endl;
      } else { // if the x^2 term is zero, then just use the minimum you got by scanning
        best = val[1];
        minchi2 = chi2[1];
        cout << "minimum with scan: " << minchi2 << " at " << best << endl;
      }
      break;
    }
  }

  limit = -1;
  delchi2prev = 1000;
  prev = 0;
  for (i = 0; i < nSinSq2Th13Steps + 1; i++) {
    ss2th13 = i * Th13increment + SinSq2Th13Low;
    for (ie = 0; ie < Extrap.size(); ie++) {
      Extrap[ie]->SetSinSq2Th13(ss2th13);
      Extrap[ie]->OscillatePrediction();
    }

    Bkgd->Reset();
    Sig->Reset();
    NExp->Reset();

    for (ie = 0; ie < Extrap.size(); ie++) {
      Bkgd->Add(Extrap[ie]->Pred_TotalBkgd_VsBinNumber);
      Sig->Add(Extrap[ie]->Pred_Signal_VsBinNumber);
    }
    NExp->Add(Bkgd);
    NExp->Add(Sig);

    delchi2 = 1e10;
    if (FitMethod == 0) {
      delchi2 = PoissonChi2(NExp) - minchi2;
    } else if (FitMethod == 1) {
      delchi2 = ScaledChi2(Bkgd, Sig) - minchi2;
    } else if (FitMethod == 2) {
      delchi2 = StandardChi2(NExp) - minchi2;
    } else if (FitMethod == 3) {
      // Likelihood: "Standard" (N syst, N nuisance)
      // Calculate the likelihood (x2 for chi)
      delchi2 = StandardLikelihood() - minchi2;
    } else if (FitMethod == 4) {
      // Likelihood: Bin by Bin Calculation of Systematics
      // Calculate the likelihood (x2 for chi)
      delchi2 = BinLikelihood() - minchi2;
    } else {
      cout << "Error in GetSensitivityAt(): Unknown 'FitMethod'." << endl;
    }

    // RBT: How to run the separate options

//     if (opt==0){
//       //Chi2 standard
//       chi2 = StandardChi2(NExp);
//     } else if (opt==1){
//       //Likelihood: Bin by Bin Calculation of Systematics
//       //Calculate the likelihood (x2 for chi)
//       chi2 = BinLikelihood();
//     } else if (opt==2){
//       //Likelihood: "Standard" (N syst, N nuisance)
//       //Calculate the likelihood (x2 for chi)
//       chi2 = StandardLikelihood();
//     }

    if (delchi2 > contourlvl && delchi2prev < contourlvl && TMath::Abs(delchi2prev - delchi2) < 0.1) {
      limit = ss2th13 + ((ss2th13 - prev) / (delchi2 - delchi2prev)) * (contourlvl - delchi2);
      sens = limit;
      break;
    }
    delchi2prev = delchi2;
    prev = ss2th13;
  }

  return sens;
}
void NueFit2D::RunDeltaChi2Contour(int cl) {
  if (NObs == 0) {
    cout << "NObs not set.  Quitting..." << endl;
    return;
  }
  if (FitMethod == 1 && (FracErr_Bkgd == 0 || FracErr_Sig == 0)) {
    cout << "FracErr_Bkgd and FracErr_Sig need to be set for ScaledChi2.  Quitting..." << endl;
    return;
  }
  if (Extrap.size() == 0) {
    cout << "No Extrapolate2D input.  Quitting..." << endl;
    return;
  }

  if (FitMethod == 3) {
    DefineStdDlnLMinuit();
  }
  if (FitMethod == 4) {
    DefineBinDlnLMinuit();
  }

  Bkgd = (TH1D *) NObs->Clone("Bkgd");
  Bkgd->Reset();
  Sig = (TH1D *) NObs->Clone("Sig");
  Sig->Reset();
  NExp = (TH1D *) NObs->Clone("NExp");
  NExp->Reset();

  unsigned int ie;
  for (ie = 0; ie < Extrap.size(); ie++) {
    Extrap[ie]->GetPrediction();
  }

  if (ErrorMatrix == 0) {
    ErrorMatrix = new TH2D("ErrorMatrix", "", nBins, -0.5, nBins - 0.5, nBins, -0.5, nBins - 0.5);
  }
  if (InvErrorMatrix == 0) {
    InvErrorMatrix = new TH2D("InvErrorMatrix", "", nBins, -0.5, nBins - 0.5, nBins, -0.5, nBins - 0.5);
  }

  Int_t i, j;

  Double_t ss2th13 = 0;
  Double_t delta = 0;
  Double_t Deltaincrement = 0;
  if (nDeltaSteps > 0) {
    Deltaincrement = (DeltaHigh - DeltaLow) / (nDeltaSteps);
  }
  Double_t Th13increment = 0;
  if (nSinSq2Th13Steps > 0) {
    Th13increment = (SinSq2Th13High - SinSq2Th13Low) / (nSinSq2Th13Steps);
  }

  double chi2[3], val[3];
  TGraph * g3;
  TF1 * fit;
  double best;
  double minchi2;

  double limit;
  double delchi2;
  double sprev, delchi2prev;
  double contourlvl = 0;
  if (cl == 0) { // 90% CL
    contourlvl = 2.71;
  } else if (cl == 1) { // 68.3% cL
    contourlvl = 1.0;
  } else {
    cout << "Error in RunDeltaChi2Contour(): Input value should be 0 or 1 for 90% or 68.3%.  Quitting..." << endl;
    return;
  }

  cout << "Seeking ";
  if (cl == 0) {
    cout << "90% ";
  } else { cout << "68% "; }
  cout << " CL upper limit" << endl;

  vector<double> deltapts;
  vector<double> bestfit_norm;
  vector<double> bestfit_invt;
  vector<double> limit_norm;
  vector<double> limit_invt;
  vector<double> lowlimit_norm;
  vector<double> lowlimit_invt;


  for (j = 0; j < nDeltaSteps + 1; j++) {
    delta = j * Deltaincrement * TMath::Pi() + DeltaLow;
    for (ie = 0; ie < Extrap.size(); ie++) {
      Extrap[ie]->SetDeltaCP(delta);

    }
    deltapts.push_back(delta / TMath::Pi());

    best = -1.;
    minchi2 = 100;
    for (i = 0; i < 3; i++) {
      chi2[i] = -1.;
      val[i] = -1.;
    }
    for (i = 0; i < nSinSq2Th13Steps + 1; i++) {
      chi2[0] = chi2[1];
      chi2[1] = chi2[2];
      val[0] = val[1];
      val[1] = val[2];
      ss2th13 = i * Th13increment + SinSq2Th13Low;
      for (ie = 0; ie < Extrap.size(); ie++) {
        Extrap[ie]->SetSinSq2Th13(ss2th13);
        Extrap[ie]->OscillatePrediction();
      }
      Bkgd->Reset();
      Sig->Reset();
      NExp->Reset();
      for (ie = 0; ie < Extrap.size(); ie++) {
        Bkgd->Add(Extrap[ie]->Pred_TotalBkgd_VsBinNumber);
        Sig->Add(Extrap[ie]->Pred_Signal_VsBinNumber);
      }
      NExp->Add(Bkgd);
      NExp->Add(Sig);

      chi2[2] = 1e10;
      if (FitMethod == 0) {
        chi2[2] = PoissonChi2(NExp);
      } else if (FitMethod == 1) {
        chi2[2] = ScaledChi2(Bkgd, Sig);
      } else if (FitMethod == 2) {
        chi2[2] = StandardChi2(NExp);
      } else if (FitMethod == 3) {
        // Likelihood: "Standard" (N syst, N nuisance)
        // Calculate the likelihood (x2 for chi)
        chi2[2] = StandardLikelihood();
      } else if (FitMethod == 4) {
        // Likelihood: Bin by Bin Calculation of Systematics
        // Calculate the likelihood (x2 for chi)
        chi2[2] = BinLikelihood();
      } else {
        cout << "Error in RunDeltaChi2Contour(): Unknown 'FitMethod'." << endl;
      }

      val[2] = ss2th13;

      if (i < 2) {
        continue;
      }

      if (i < 3 && chi2[2] > chi2[1] && chi2[1] > chi2[0]) { // first three points are increasing, first point is minimum.
        best = val[0];
        minchi2 = chi2[0];
        cout << "minimum at 1st point: " << minchi2 << " at " << best << endl;
        break;
      }

      if (chi2[2] > chi2[1] && chi2[0] > chi2[1]) { // found minimum
        g3 = new TGraph(3, val, chi2);
        fit = new TF1("pol2", "pol2");
        g3->Fit(fit, "Q"); // fit to second order polynominal
        if (fit->GetParameter(2) > 0) { // if the x^2 term is nonzero
          best = -fit->GetParameter(1) / (2 * fit->GetParameter(2)); // the location of the minimum is -p1/(2*p2)
          minchi2 = fit->GetParameter(0) + fit->GetParameter(1) * best + fit->GetParameter(2) * best * best;
          cout << "minimum with fit: " << minchi2 << " at " << best << endl;
        } else { // if the x^2 term is zero, then just use the minimum you got by scanning
          best = val[1];
          minchi2 = chi2[1];
          cout << "minimum with scan: " << minchi2 << " at " << best << endl;
        }
        break;
      }
    }
    bestfit_norm.push_back(best);

    limit = -1.;
    delchi2prev = 1000;
    sprev = 0;
    for (i = 0; i < (SinSq2Th13High - best) / Th13increment; i++) {

      ss2th13 = i * Th13increment + best; // SinSq2Th13Low;
      for (ie = 0; ie < Extrap.size(); ie++) {
        Extrap[ie]->SetSinSq2Th13(ss2th13);
        Extrap[ie]->OscillatePrediction();
      }
      Bkgd->Reset();
      Sig->Reset();
      NExp->Reset();
      for (ie = 0; ie < Extrap.size(); ie++) {
        Bkgd->Add(Extrap[ie]->Pred_TotalBkgd_VsBinNumber);
        Sig->Add(Extrap[ie]->Pred_Signal_VsBinNumber);
      }
      NExp->Add(Bkgd);
      NExp->Add(Sig);

      delchi2 = 1e10;
      if (FitMethod == 0) {
        delchi2 = PoissonChi2(NExp) - minchi2;
      } else if (FitMethod == 1) {
        delchi2 = ScaledChi2(Bkgd, Sig) - minchi2;
      } else if (FitMethod == 2) {
        delchi2 = StandardChi2(NExp) - minchi2;
      } else if (FitMethod == 3) {
        // Likelihood: "Standard" (N syst, N nuisance)
        // Calculate the likelihood (x2 for chi)
        delchi2 = StandardLikelihood() - minchi2;
      } else if (FitMethod == 4) {
        // Likelihood: Bin by Bin Calculation of Systematics
        // Calculate the likelihood (x2 for chi)
        delchi2 = BinLikelihood() - minchi2;
      } else {
        cout << "Error in RunDeltaChi2Contour(): Unknown 'FitMethod'." << endl;
      }

      if (i == 1) {
        continue;
      }

      if (delchi2 > contourlvl && delchi2prev < contourlvl) {
        limit = ss2th13 + ((ss2th13 - sprev) / (delchi2 - delchi2prev)) * (contourlvl - delchi2);
        cout << delchi2prev << ", " << sprev << ", " << delchi2 << ", " << ss2th13 << ", " << limit << endl;
        break;
      }
      delchi2prev = delchi2;
      sprev = ss2th13;
    }
    limit_norm.push_back(limit);


    limit = 0.;
    delchi2prev = 1000;
    sprev = 0;
    for (i = 0; i < (best - SinSq2Th13Low) / Th13increment; i++) {

      ss2th13 = -i * Th13increment + best;
      for (ie = 0; ie < Extrap.size(); ie++) {
        Extrap[ie]->SetSinSq2Th13(ss2th13);
        Extrap[ie]->OscillatePrediction();
      }
      Bkgd->Reset();
      Sig->Reset();
      NExp->Reset();
      for (ie = 0; ie < Extrap.size(); ie++) {
        Bkgd->Add(Extrap[ie]->Pred_TotalBkgd_VsBinNumber);
        Sig->Add(Extrap[ie]->Pred_Signal_VsBinNumber);
      }
      NExp->Add(Bkgd);
      NExp->Add(Sig);

      delchi2 = 1e10;
      if (FitMethod == 0) {
        delchi2 = PoissonChi2(NExp) - minchi2;
      } else if (FitMethod == 1) {
        delchi2 = ScaledChi2(Bkgd, Sig) - minchi2;
      } else if (FitMethod == 2) {
        delchi2 = StandardChi2(NExp) - minchi2;
      } else if (FitMethod == 3) {
        // Likelihood: "Standard" (N syst, N nuisance)
        // Calculate the likelihood (x2 for chi)
        delchi2 = StandardLikelihood() - minchi2;
      } else if (FitMethod == 4) {
        // Likelihood: Bin by Bin Calculation of Systematics
        // Calculate the likelihood (x2 for chi)
        delchi2 = BinLikelihood() - minchi2;
      } else {
        cout << "Error in RunDeltaChi2Contour(): Unknown 'FitMethod'." << endl;
      }

      if (i == 1) {
        continue;
      }

      if (delchi2 > contourlvl && delchi2prev < contourlvl) {
        limit = ss2th13 + ((ss2th13 - sprev) / (delchi2 - delchi2prev)) * (contourlvl - delchi2);
        cout << delchi2prev << ", " << sprev << ", " << delchi2 << ", " << ss2th13 << ", " << limit << endl;
        break;
      }
      delchi2prev = delchi2;
      sprev = ss2th13;
    }
    lowlimit_norm.push_back(limit);

  }

  for (ie = 0; ie < Extrap.size(); ie++) {
    Extrap[ie]->InvertMassHierarchy();
  }

  for (j = 0; j < nDeltaSteps + 1; j++) {
    delta = j * Deltaincrement * TMath::Pi() + DeltaLow;
    for (ie = 0; ie < Extrap.size(); ie++) {
      Extrap[ie]->SetDeltaCP(delta);
    }

    best = -1.;
    minchi2 = 100;
    for (i = 0; i < 3; i++) {
      chi2[i] = -1.;
      val[i] = -1.;
    }
    for (i = 0; i < nSinSq2Th13Steps + 1; i++) {
      chi2[0] = chi2[1];
      chi2[1] = chi2[2];
      val[0] = val[1];
      val[1] = val[2];
      ss2th13 = i * Th13increment + SinSq2Th13Low;
      for (ie = 0; ie < Extrap.size(); ie++) {
        Extrap[ie]->SetSinSq2Th13(ss2th13);
        Extrap[ie]->OscillatePrediction();
      }
      Bkgd->Reset();
      Sig->Reset();
      NExp->Reset();
      for (ie = 0; ie < Extrap.size(); ie++) {
        Bkgd->Add(Extrap[ie]->Pred_TotalBkgd_VsBinNumber);
        Sig->Add(Extrap[ie]->Pred_Signal_VsBinNumber);
      }
      NExp->Add(Bkgd);
      NExp->Add(Sig);

      chi2[2] = 1e10;
      if (FitMethod == 0) {
        chi2[2] = PoissonChi2(NExp);
      } else if (FitMethod == 1) {
        chi2[2] = ScaledChi2(Bkgd, Sig);
      } else if (FitMethod == 2) {
        chi2[2] = StandardChi2(NExp);
      } else if (FitMethod == 3) {
        // Likelihood: "Standard" (N syst, N nuisance)
        // Calculate the likelihood (x2 for chi)
        chi2[2] = StandardLikelihood();
      } else if (FitMethod == 4) {
        // Likelihood: Bin by Bin Calculation of Systematics
        // Calculate the likelihood (x2 for chi)
        chi2[2] = BinLikelihood();
      } else {
        cout << "Error in RunDeltaChi2Contour(): Unknown 'FitMethod'." << endl;
      }

      val[2] = ss2th13;

      if (i < 2) {
        continue;
      }

      if (i < 3 && chi2[2] > chi2[1] && chi2[1] > chi2[0]) { // first three points are increasing, first point is minimum.
        best = val[0];
        minchi2 = chi2[0];
        cout << "minimum at 1st point: " << minchi2 << " at " << best << endl;
        break;
      }

      if (chi2[2] > chi2[1] && chi2[0] > chi2[1]) { // found minimum
        g3 = new TGraph(3, val, chi2);
        fit = new TF1("pol2", "pol2");
        g3->Fit(fit, "Q"); // fit to second order polynominal
        if (fit->GetParameter(2) > 0) { // if the x^2 term is nonzero
          best = -fit->GetParameter(1) / (2 * fit->GetParameter(2)); // the location of the minimum is -p1/(2*p2)
          minchi2 = fit->GetParameter(0) + fit->GetParameter(1) * best + fit->GetParameter(2) * best * best;
          cout << "minimum with fit: " << minchi2 << " at " << best << endl;
        } else { // if the x^2 term is zero, then just use the minimum you got by scanning
          best = val[1];
          minchi2 = chi2[1];
          cout << "minimum with scan: " << minchi2 << " at " << best << endl;
        }
        break;
      }
    }
    bestfit_invt.push_back(best);

    limit = -1.;
    delchi2prev = 1000;
    sprev = 0;
    for (i = 0; i < (SinSq2Th13High - best) / Th13increment; i++) {
      ss2th13 = i * Th13increment + best; // SinSq2Th13Low;
      for (ie = 0; ie < Extrap.size(); ie++) {
        Extrap[ie]->SetSinSq2Th13(ss2th13);
        Extrap[ie]->OscillatePrediction();
      }
      Bkgd->Reset();
      Sig->Reset();
      NExp->Reset();
      for (ie = 0; ie < Extrap.size(); ie++) {
        Bkgd->Add(Extrap[ie]->Pred_TotalBkgd_VsBinNumber);
        Sig->Add(Extrap[ie]->Pred_Signal_VsBinNumber);
      }
      NExp->Add(Bkgd);
      NExp->Add(Sig);

      delchi2 = 1e10;
      if (FitMethod == 0) {
        delchi2 = PoissonChi2(NExp) - minchi2;
      } else if (FitMethod == 1) {
        delchi2 = ScaledChi2(Bkgd, Sig) - minchi2;
      } else if (FitMethod == 2) {
        delchi2 = StandardChi2(NExp) - minchi2;
      } else if (FitMethod == 3) {
        // Likelihood: "Standard" (N syst, N nuisance)
        // Calculate the likelihood (x2 for chi)
        delchi2 = StandardLikelihood() - minchi2;
      } else if (FitMethod == 4) {
        // Likelihood: Bin by Bin Calculation of Systematics
        // Calculate the likelihood (x2 for chi)
        delchi2 = BinLikelihood() - minchi2;
      } else {
        cout << "Error in RunDeltaChi2Contour(): Unknown 'FitMethod'." << endl;
      }

      if (i == 1) {
        continue;
      }

      if (delchi2 > contourlvl && delchi2prev < contourlvl) {
        limit = ss2th13 + ((ss2th13 - sprev) / (delchi2 - delchi2prev)) * (contourlvl - delchi2);
        cout << delchi2prev << ", " << sprev << ", " << delchi2 << ", " << ss2th13 << ", " << limit << endl;
        break;
      }
      delchi2prev = delchi2;
      sprev = ss2th13;
    }
    limit_invt.push_back(limit);

    limit = 0.;
    delchi2prev = 1000;
    sprev = 0;
    for (i = 0; i < (best - SinSq2Th13Low) / Th13increment; i++) {
      ss2th13 = -i * Th13increment + best; // SinSq2Th13Low;
      for (ie = 0; ie < Extrap.size(); ie++) {
        Extrap[ie]->SetSinSq2Th13(ss2th13);
        Extrap[ie]->OscillatePrediction();
      }
      Bkgd->Reset();
      Sig->Reset();
      NExp->Reset();
      for (ie = 0; ie < Extrap.size(); ie++) {
        Bkgd->Add(Extrap[ie]->Pred_TotalBkgd_VsBinNumber);
        Sig->Add(Extrap[ie]->Pred_Signal_VsBinNumber);
      }
      NExp->Add(Bkgd);
      NExp->Add(Sig);

      delchi2 = 1e10;
      if (FitMethod == 0) {
        delchi2 = PoissonChi2(NExp) - minchi2;
      } else if (FitMethod == 1) {
        delchi2 = ScaledChi2(Bkgd, Sig) - minchi2;
      } else if (FitMethod == 2) {
        delchi2 = StandardChi2(NExp) - minchi2;
      } else if (FitMethod == 3) {
        // Likelihood: "Standard" (N syst, N nuisance)
        // Calculate the likelihood (x2 for chi)
        delchi2 = StandardLikelihood() - minchi2;
      } else if (FitMethod == 4) {
        // Likelihood: Bin by Bin Calculation of Systematics
        // Calculate the likelihood (x2 for chi)
        delchi2 = BinLikelihood() - minchi2;
      } else {
        cout << "Error in RunDeltaChi2Contour(): Unknown 'FitMethod'." << endl;
      }

      if (i == 1) {
        continue;
      }

      if (delchi2 > contourlvl && delchi2prev < contourlvl) {
        limit = ss2th13 + ((ss2th13 - sprev) / (delchi2 - delchi2prev)) * (contourlvl - delchi2);
        cout << delchi2prev << ", " << sprev << ", " << delchi2 << ", " << ss2th13 << ", " << limit << endl;
        break;
      }
      delchi2prev = delchi2;
      sprev = ss2th13;
    }
    lowlimit_invt.push_back(limit);

  }

  double * s_best_n = new double[nDeltaSteps + 1];
  double * s_best_i = new double[nDeltaSteps + 1];
  double * s_limit_n = new double[nDeltaSteps + 1];
  double * s_limit_i = new double[nDeltaSteps + 1];
  double * s_limit_nlow = new double[nDeltaSteps + 1];
  double * s_limit_ilow = new double[nDeltaSteps + 1];
  double * d = new double[nDeltaSteps + 1];
  for (i = 0; i < nDeltaSteps + 1; i++) {
    d[i] = deltapts.at(i);
    s_best_n[i] = bestfit_norm.at(i);
    s_best_i[i] = bestfit_invt.at(i);
    s_limit_n[i] = limit_norm.at(i);
    s_limit_i[i] = limit_invt.at(i);
    s_limit_nlow[i] = lowlimit_norm.at(i);
    s_limit_ilow[i] = lowlimit_invt.at(i);
  }

  TGraph * gn = new TGraph(nDeltaSteps + 1, s_best_n, d);
  gn->SetMarkerStyle(20);
  gn->SetTitle("");
  gn->GetYaxis()->SetTitle("#delta");
  gn->GetXaxis()->SetTitle("sin^{2}2#theta_{13}");
  gn->GetXaxis()->SetLimits(0, 0.6);
  gn->SetLineWidth(4);
  gn->SetMaximum(2);
  gn->SetName("BestFit_Normal");

  TGraph * gi = new TGraph(nDeltaSteps + 1, s_best_i, d);
  gi->SetMarkerStyle(20);
  gi->SetTitle("");
  gi->GetYaxis()->SetTitle("#delta");
  gi->GetXaxis()->SetTitle("sin^{2}2#theta_{13}");
  gi->GetXaxis()->SetLimits(0, 0.6);
  gi->SetLineWidth(4);
  gi->SetLineStyle(2);
  gi->SetMaximum(2);
  gi->SetName("BestFit_Inverted");

  TGraph * gn_limit = new TGraph(nDeltaSteps + 1, s_limit_n, d);
  gn_limit->SetMarkerStyle(20);
  gn_limit->SetTitle("");
  gn_limit->GetYaxis()->SetTitle("#delta");
  gn_limit->GetXaxis()->SetTitle("sin^{2}2#theta_{13}");
  gn_limit->GetXaxis()->SetLimits(0, 0.6);
  gn_limit->SetLineWidth(4);
  gn_limit->SetLineColor(kBlue);
  gn_limit->SetMarkerColor(kBlue);
  gn_limit->SetMaximum(2);
  gn_limit->SetName("UpperLimit_Normal");

  TGraph * gi_limit = new TGraph(nDeltaSteps + 1, s_limit_i, d);
  gi_limit->SetMarkerStyle(20);
  gi_limit->SetTitle("");
  gi_limit->GetYaxis()->SetTitle("#delta");
  gi_limit->GetXaxis()->SetTitle("sin^{2}2#theta_{13}");
  gi_limit->GetXaxis()->SetLimits(0, 0.6);
  gi_limit->SetLineWidth(4);
  gi_limit->SetLineColor(kRed);
  gi_limit->SetMarkerColor(kRed);
  gi_limit->SetMaximum(2);
  gi_limit->SetName("UpperLimit_Inverted");

  TGraph * gnl_limit = new TGraph(nDeltaSteps + 1, s_limit_nlow, d);
  gnl_limit->SetMarkerStyle(20);
  gnl_limit->SetTitle("");
  gnl_limit->GetYaxis()->SetTitle("#delta");
  gnl_limit->GetXaxis()->SetTitle("sin^{2}2#theta_{13}");
  gnl_limit->GetXaxis()->SetLimits(0, 0.6);
  gnl_limit->SetLineWidth(4);
  gnl_limit->SetLineColor(kBlue);
  gnl_limit->SetMarkerColor(kBlue);
  gnl_limit->SetMaximum(2);
  gnl_limit->SetName("LowerLimit_Normal");

  TGraph * gil_limit = new TGraph(nDeltaSteps + 1, s_limit_ilow, d);
  gil_limit->SetMarkerStyle(20);
  gil_limit->SetTitle("");
  gil_limit->GetYaxis()->SetTitle("#delta");
  gil_limit->GetXaxis()->SetTitle("sin^{2}2#theta_{13}");
  gil_limit->GetXaxis()->SetLimits(0, 0.6);
  gil_limit->SetLineWidth(4);
  gil_limit->SetLineColor(kRed);
  gil_limit->SetMarkerColor(kRed);
  gil_limit->SetMaximum(2);
  gil_limit->SetName("LowerLimit_Inverted");


  if (cl == 0) {
    cout << "90% ";
  }
  if (cl == 1) {
    cout << "68% ";
  }
  cout << "confidence level limit = " << limit_norm.at(0) << ", " << limit_invt.at(0) << endl;

  TFile * fout = new TFile(gSystem->ExpandPathName(outFileName.c_str()), "RECREATE");
  gn->Write();
  gi->Write();
  gn_limit->Write();
  gi_limit->Write();
  gnl_limit->Write();
  gil_limit->Write();
  fout->Write();
  fout->Close();

  delete [] s_best_n;
  delete [] s_best_i;
  delete [] s_limit_n;
  delete [] s_limit_i;
  delete [] s_limit_nlow;
  delete [] s_limit_ilow;
  delete [] d;

  return;
}
double NueFit2D::GetLikelihood(double t12, double t23, double t13, double dm2_32, double dm2_21, double delta) {
  if (NObs == 0) {
    cout << "NObs not set.  Quitting..." << endl;
    return 0;
  }
  if (FitMethod == 1 && (FracErr_Bkgd == 0 || FracErr_Sig == 0)) {
    cout << "FracErr_Bkgd and FracErr_Sig need to be set for ScaledChi2.  Quitting..." << endl;
    return 0;
  }
  if (Extrap.size() == 0) {
    cout << "No Extrapolate2D input.  Quitting..." << endl;
    return 0;
  }

  if (FitMethod == 3) {
    DefineStdDlnLMinuit();
  }
  if (FitMethod == 4) {
    DefineBinDlnLMinuit();
  }

  Bkgd = (TH1D *) NObs->Clone("Bkgd");
  Bkgd->Reset();
  Sig = (TH1D *) NObs->Clone("Sig");
  Sig->Reset();
  TH1D * NExp = (TH1D *) NObs->Clone("NExp");
  NExp->Reset();

  unsigned int ie;
  for (ie = 0; ie < Extrap.size(); ie++) {
    Extrap[ie]->GetPrediction();
    Extrap[ie]->SetOscPar(OscPar::kTh12, t12);
    Extrap[ie]->SetOscPar(OscPar::kTh23, t23);
    Extrap[ie]->SetOscPar(OscPar::kTh13, t13);
    Extrap[ie]->SetOscPar(OscPar::kDeltaM23, dm2_32);
    Extrap[ie]->SetOscPar(OscPar::kDeltaM12, dm2_21);
    Extrap[ie]->SetOscPar(OscPar::kDelta, delta);
    Extrap[ie]->OscillatePrediction();
  }

  Int_t i;

  if (ErrorMatrix == 0) {
    ErrorMatrix = new TH2D("ErrorMatrix", "", nBins, -0.5, nBins - 0.5, nBins, -0.5, nBins - 0.5);
  }
  if (InvErrorMatrix == 0) {
    InvErrorMatrix = new TH2D("InvErrorMatrix", "", nBins, -0.5, nBins - 0.5, nBins, -0.5, nBins - 0.5);
  }

  Double_t ss2th13 = 0;
  Double_t delchi2 = 0;
  Double_t best = 0;
  Double_t Th13increment = 0;
  if (nSinSq2Th13Steps > 0) {
    Th13increment = (SinSq2Th13High - SinSq2Th13Low) / (nSinSq2Th13Steps);
  }
  double chi2[3], val[3];
  TGraph * g3;
  TF1 * fit;
  double minchi2;

  best = -1.;
  minchi2 = 100;
  for (i = 0; i < 3; i++) {
    chi2[i] = -1.;
    val[i] = -1.;
  }
  for (i = 0; i < nSinSq2Th13Steps + 1; i++) {
    chi2[0] = chi2[1];
    chi2[1] = chi2[2];
    val[0] = val[1];
    val[1] = val[2];
    ss2th13 = i * Th13increment + SinSq2Th13Low;
    for (ie = 0; ie < Extrap.size(); ie++) {
      Extrap[ie]->SetSinSq2Th13(ss2th13);
      Extrap[ie]->OscillatePrediction();
    }
    Bkgd->Reset();
    Sig->Reset();
    NExp->Reset();
    for (ie = 0; ie < Extrap.size(); ie++) {
      Bkgd->Add(Extrap[ie]->Pred_TotalBkgd_VsBinNumber);
      Sig->Add(Extrap[ie]->Pred_Signal_VsBinNumber);
    }
    NExp->Add(Bkgd);
    NExp->Add(Sig);

    chi2[2] = 1e10;
    if (FitMethod == 0) {
      chi2[2] = PoissonChi2(NExp);
    } else if (FitMethod == 1) {
      chi2[2] = ScaledChi2(Bkgd, Sig);
    } else if (FitMethod == 2) {
      chi2[2] = StandardChi2(NExp);
    } else if (FitMethod == 3) {
      // Likelihood: "Standard" (N syst, N nuisance)
      // Calculate the likelihood (x2 for chi)
      chi2[2] = StandardLikelihood();
    } else if (FitMethod == 4) {
      // Likelihood: Bin by Bin Calculation of Systematics
      // Calculate the likelihood (x2 for chi)
      chi2[2] = BinLikelihood();
    } else {
      cout << "Error in GetLikelihood(): Unknown 'FitMethod'." << endl;
    }

    val[2] = ss2th13;

    if (i < 2) {
      continue;
    }

    if (i < 3 && chi2[2] > chi2[1] && chi2[1] > chi2[0]) { // first three points are increasing, first point is minimum.
      best = val[0];
      minchi2 = chi2[0];
      cout << "minimum at 1st point: " << minchi2 << " at " << best << endl;
      break;
    }

    if (chi2[2] > chi2[1] && chi2[0] > chi2[1]) { // found minimum
      g3 = new TGraph(3, val, chi2);
      fit = new TF1("pol2", "pol2");
      g3->Fit(fit, "Q"); // fit to second order polynominal
      if (fit->GetParameter(2) > 0) { // if the x^2 term is nonzero
        best = -fit->GetParameter(1) / (2 * fit->GetParameter(2)); // the location of the minimum is -p1/(2*p2)
        minchi2 = fit->GetParameter(0) + fit->GetParameter(1) * best + fit->GetParameter(2) * best * best;
        cout << "minimum with fit: " << minchi2 << " at " << best << endl;
      } else { // if the x^2 term is zero, then just use the minimum you got by scanning
        best = val[1];
        minchi2 = chi2[1];
        cout << "minimum with scan: " << minchi2 << " at " << best << endl;
      }
      break;
    }
  }

  for (ie = 0; ie < Extrap.size(); ie++) {
    Extrap[ie]->SetOscPar(OscPar::kTh13, t13);
    Extrap[ie]->OscillatePrediction();
  }
  Bkgd->Reset();
  Sig->Reset();
  NExp->Reset();
  for (ie = 0; ie < Extrap.size(); ie++) {
    Bkgd->Add(Extrap[ie]->Pred_TotalBkgd_VsBinNumber);
    Sig->Add(Extrap[ie]->Pred_Signal_VsBinNumber);
  }
  NExp->Add(Bkgd);
  NExp->Add(Sig);

  delchi2 = 1e10;
  if (FitMethod == 0) {
    delchi2 = PoissonChi2(NExp) - minchi2;
  } else if (FitMethod == 1) {
    delchi2 = ScaledChi2(Bkgd, Sig) - minchi2;
  } else if (FitMethod == 2) {
    delchi2 = StandardChi2(NExp) - minchi2;
  } else if (FitMethod == 3) {
    // Likelihood: "Standard" (N syst, N nuisance)
    // Calculate the likelihood (x2 for chi)
    delchi2 = StandardLikelihood() - minchi2;
  } else if (FitMethod == 4) {
    // Likelihood: Bin by Bin Calculation of Systematics
    // Calculate the likelihood (x2 for chi)
    delchi2 = BinLikelihood() - minchi2;
  } else {
    cout << "Error in GetLikelihood(): Unknown 'FitMethod'." << endl;
  }

  cout << "delchi2 = " << delchi2 << endl;

  return delchi2;
}
double NueFit2D::GetMinLikelihood(double delta, bool normalhier) {
  if (NObs == 0) {
    cout << "NObs not set.  Quitting..." << endl;
    return 0;
  }
  if (FitMethod == 1 && (FracErr_Bkgd == 0 || FracErr_Sig == 0)) {
    cout << "FracErr_Bkgd and FracErr_Sig need to be set for ScaledChi2.  Quitting..." << endl;
    return 0;
  }
  if (Extrap.size() == 0) {
    cout << "No Extrapolate2D input.  Quitting..." << endl;
    return 0;
  }

  if (FitMethod == 3) {
    DefineStdDlnLMinuit();
  }
  if (FitMethod == 4) {
    DefineBinDlnLMinuit();
  }
  if (FitMethod == 3 || FitMethod == 4) {
    if (ErrCalc != 0) {
      ErrCalc->SetUseGrid(false);
    }
  }

  Bkgd = (TH1D *) NObs->Clone("Bkgd");
  Bkgd->Reset();
  Sig = (TH1D *) NObs->Clone("Sig");
  Sig->Reset();
  TH1D * NExp = (TH1D *) NObs->Clone("NExp");
  NExp->Reset();

  unsigned int ie;
  for (ie = 0; ie < Extrap.size(); ie++) {
    Extrap[ie]->GetPrediction();
    if (normalhier) {
      Extrap[ie]->SetOscPar(OscPar::kTh12, grid_n_th12);
      Extrap[ie]->SetOscPar(OscPar::kTh23, grid_n_th23);
      Extrap[ie]->SetOscPar(OscPar::kDeltaM12, grid_n_dm2_21);
      Extrap[ie]->SetOscPar(OscPar::kDeltaM23, grid_n_dm2_32);
    } else {
      Extrap[ie]->SetOscPar(OscPar::kTh12, grid_i_th12);
      Extrap[ie]->SetOscPar(OscPar::kTh23, grid_i_th23);
      Extrap[ie]->SetOscPar(OscPar::kDeltaM12, grid_i_dm2_21);
      Extrap[ie]->SetOscPar(OscPar::kDeltaM23, grid_i_dm2_32);
    }
    Extrap[ie]->SetOscPar(OscPar::kDelta, delta);
    Extrap[ie]->OscillatePrediction();
  }

  Int_t i;

  if (ErrorMatrix == 0) {
    ErrorMatrix = new TH2D("ErrorMatrix", "", nBins, -0.5, nBins - 0.5, nBins, -0.5, nBins - 0.5);
  }
  if (InvErrorMatrix == 0) {
    InvErrorMatrix = new TH2D("InvErrorMatrix", "", nBins, -0.5, nBins - 0.5, nBins, -0.5, nBins - 0.5);
  }

  Double_t ss2th13 = 0;
  TF1 * fit;
  double minchi2;

  const int ns = 30;
  double max = 0.6;
  double width = max / ns;
  double th13[ns], c2[ns];
  for (i = 0; i < ns; i++) {
    ss2th13 = i * width;
    th13[i] = ss2th13;
    for (ie = 0; ie < Extrap.size(); ie++) {
      Extrap[ie]->SetSinSq2Th13(ss2th13);
      Extrap[ie]->OscillatePrediction();
    }
    Bkgd->Reset();
    Sig->Reset();
    NExp->Reset();
    for (ie = 0; ie < Extrap.size(); ie++) {
      Bkgd->Add(Extrap[ie]->Pred_TotalBkgd_VsBinNumber);
      Sig->Add(Extrap[ie]->Pred_Signal_VsBinNumber);
    }
    NExp->Add(Bkgd);
    NExp->Add(Sig);

    c2[i] = 1e10;
    if (FitMethod == 0) {
      c2[i] = PoissonChi2(NExp);
    } else if (FitMethod == 1) {
      c2[i] = ScaledChi2(Bkgd, Sig);
    } else if (FitMethod == 2) {
      c2[i] = StandardChi2(NExp);
    } else if (FitMethod == 3) {
      // Likelihood: "Standard" (N syst, N nuisance)
      // Calculate the likelihood (x2 for chi)
      c2[i] = StandardLikelihood();
    } else if (FitMethod == 4) {
      // Likelihood: Bin by Bin Calculation of Systematics
      // Calculate the likelihood (x2 for chi)
      c2[i] = BinLikelihood();
    } else {
      cout << "Error in GetMinLikelihood(): Unknown 'FitMethod'." << endl;
    }
  }

  TGraph * g = new TGraph(ns, th13, c2);
  fit = new TF1("pol6", "pol6");
  g->Fit(fit, "Q"); // fit to second order polynominal
  minchi2 = fit->GetMinimum(0, max);

  return minchi2;
}
double NueFit2D::GetMinLikelihood_Delta(bool normalhier) {
  if (NObs == 0) {
    cout << "NObs not set.  Quitting..." << endl;
    return 0;
  }
  if (FitMethod == 1 && (FracErr_Bkgd == 0 || FracErr_Sig == 0)) {
    cout << "FracErr_Bkgd and FracErr_Sig need to be set for ScaledChi2.  Quitting..." << endl;
    return 0;
  }
  if (Extrap.size() == 0) {
    cout << "No Extrapolate2D input.  Quitting..." << endl;
    return 0;
  }

  if (FitMethod == 3) {
    DefineStdDlnLMinuit();
  }
  if (FitMethod == 4) {
    DefineBinDlnLMinuit();
  }
  if (FitMethod == 3 || FitMethod == 4) {
    if (ErrCalc != 0) {
      ErrCalc->SetUseGrid(false);
    }
  }

  Bkgd = (TH1D *) NObs->Clone("Bkgd");
  Bkgd->Reset();
  Sig = (TH1D *) NObs->Clone("Sig");
  Sig->Reset();
  TH1D * NExp = (TH1D *) NObs->Clone("NExp");
  NExp->Reset();

  unsigned int ie;
  for (ie = 0; ie < Extrap.size(); ie++) {
    Extrap[ie]->GetPrediction();
    if (normalhier) {
      Extrap[ie]->SetOscPar(OscPar::kTh13, grid_n_th13);
      Extrap[ie]->SetOscPar(OscPar::kTh12, grid_n_th12);
      Extrap[ie]->SetOscPar(OscPar::kTh23, grid_n_th23);
      Extrap[ie]->SetOscPar(OscPar::kDeltaM12, grid_n_dm2_21);
      Extrap[ie]->SetOscPar(OscPar::kDeltaM23, grid_n_dm2_32);
    } else {
      Extrap[ie]->SetOscPar(OscPar::kTh13, grid_i_th13);
      Extrap[ie]->SetOscPar(OscPar::kTh12, grid_i_th12);
      Extrap[ie]->SetOscPar(OscPar::kTh23, grid_i_th23);
      Extrap[ie]->SetOscPar(OscPar::kDeltaM12, grid_i_dm2_21);
      Extrap[ie]->SetOscPar(OscPar::kDeltaM23, grid_i_dm2_32);
    }
    Extrap[ie]->OscillatePrediction();
  }

  Int_t i;

  if (ErrorMatrix == 0) {
    ErrorMatrix = new TH2D("ErrorMatrix", "", nBins, -0.5, nBins - 0.5, nBins, -0.5, nBins - 0.5);
  }
  if (InvErrorMatrix == 0) {
    InvErrorMatrix = new TH2D("InvErrorMatrix", "", nBins, -0.5, nBins - 0.5, nBins, -0.5, nBins - 0.5);
  }

  Double_t delta = 0;
  //   TF1* fit;
  double minchi2, mindelta;
  Double_t increment = 0;
  if (nDeltaSteps > 0) {
    increment = 2 * TMath::Pi() / (nDeltaSteps);
  }
  const int nd = nDeltaSteps + 1;
  double * d = new double[nDeltaSteps + 1];
  double * c2 = new double[nDeltaSteps + 1];
  for (i = 0; i < nd; i++) {
    delta = i * increment;
    d[i] = delta;
    for (ie = 0; ie < Extrap.size(); ie++) {
      Extrap[ie]->SetDeltaCP(delta);
      Extrap[ie]->OscillatePrediction();
    }
    Bkgd->Reset();
    Sig->Reset();
    NExp->Reset();
    for (ie = 0; ie < Extrap.size(); ie++) {
      Bkgd->Add(Extrap[ie]->Pred_TotalBkgd_VsBinNumber);
      Sig->Add(Extrap[ie]->Pred_Signal_VsBinNumber);
    }
    NExp->Add(Bkgd);
    NExp->Add(Sig);

    c2[i] = 1e10;
    if (FitMethod == 0) {
      c2[i] = PoissonChi2(NExp);
    } else if (FitMethod == 1) {
      c2[i] = ScaledChi2(Bkgd, Sig);
    } else if (FitMethod == 2) {
      c2[i] = StandardChi2(NExp);
    } else if (FitMethod == 3) {
      // Likelihood: "Standard" (N syst, N nuisance)
      // Calculate the likelihood (x2 for chi)
      c2[i] = StandardLikelihood();
    } else if (FitMethod == 4) {
      // Likelihood: Bin by Bin Calculation of Systematics
      // Calculate the likelihood (x2 for chi)
      c2[i] = BinLikelihood();
    } else {
      cout << "Error in GetMinLikelihood_Delta(): Unknown 'FitMethod'." << endl;
    }
  }

  // can we actually fit this to a nice function?
  //   TGraph *g = new TGraph(ns,th13,c2);
  //   fit = new TF1("pol6", "pol6");
  //   g->Fit(fit, "Q");//fit to second order polynominal
  //   minchi2 = fit->GetMinimum(0,max);

  minchi2 = 1000;
  mindelta = 0;
  for (i = 0; i < nd; i++) {
    if (c2[i] < minchi2) {
      minchi2 = c2[i];
      mindelta = d[i];
    }
  }

  cout << "min delta = " << mindelta << " , min chi2 = " << minchi2 << endl;

  delete [] d;
  delete [] c2;

  return minchi2;
}
void NueFit2D::RunMultiBinPseudoExpts_MHDeltaFit(bool Print) {
  if (FitMethod == 1 && (FracErr_Bkgd == 0 || FracErr_Sig == 0)) {
    cout << "FracErr_Bkgd and FracErr_Sig need to be set for ScaledChi2.  Quitting..." << endl;
    return;
  }
  if (Extrap.size() == 0) {
    cout << "No Extrapolate2D input.  Quitting..." << endl;
    return;
  }
  for (unsigned int ie = 0; ie < Extrap.size(); ie++) {
    Extrap[ie]->GetPrediction();
  }
  if (ErrCalc == 0) {
    cout << "Need to set ErrorCalc object!  Quitting..." << endl;
    return;
  }
  if (FitMethod == 3) {
    DefineStdDlnLMinuit();
  }
  if (FitMethod == 4) {
    DefineBinDlnLMinuit();
  }

  nBins = Extrap[0]->Pred_TotalBkgd_VsBinNumber->GetNbinsX();
  if (ErrorMatrix == 0) {
    ErrorMatrix = new TH2D("ErrorMatrix", "", nBins, -0.5, nBins - 0.5, nBins, -0.5, nBins - 0.5);
  }
  if (InvErrorMatrix == 0) {
    InvErrorMatrix = new TH2D("InvErrorMatrix", "", nBins, -0.5, nBins - 0.5, nBins, -0.5, nBins - 0.5);
  }

  TH2D * Error4Expts = new TH2D("Error4Expts", "", nBins, -0.5, nBins - 0.5, nBins, -0.5, nBins - 0.5);

  ReadGridFiles();

  if (nPts_Normal == 0 || nPts_Inverted == 0) {
    return;
  }

  gRandom->SetSeed(0);

  int i, u;
  unsigned int j, k;
  TH1D * nexp_bkgd_n = new TH1D("nexp_bkgd_n", "", nBins, -0.5, nBins - 0.5);
  TH1D * nexp_signal_n = new TH1D("nexp_signal_n", "", nBins, -0.5, nBins - 0.5);
  TH1D * nexp_n = new TH1D("nexp_n", "", nBins, -0.5, nBins - 0.5);

  TH1D * nexp_bkgd_i = new TH1D("nexp_bkgd_i", "", nBins, -0.5, nBins - 0.5);
  TH1D * nexp_signal_i = new TH1D("nexp_signal_i", "", nBins, -0.5, nBins - 0.5);
  TH1D * nexp_i = new TH1D("nexp_i", "", nBins, -0.5, nBins - 0.5);

  double chi2_norm, chi2_invt, chi2min_norm, chi2min_invt;
  double delchi2_norm, delchi2_invt, delchi2_bestmh;

  TH1D * chi2hist_deltafit = new TH1D("chi2hist_deltafit", "", 110000, -10, 100);
  TH1D * chi2hist_mhfit = new TH1D("chi2hist_mhfit", "", 110000, -10, 100);
  double ele;
  int noff;

  vector< vector<double> > nc_n, numucc_n, bnuecc_n, nutaucc_n, sig_n;
  vector< vector<double> > nc_i, numucc_i, bnuecc_i, nutaucc_i, sig_i;
  for (j = 0; j < nBins; j++) {
    nc_n.push_back(vector<double>() );
    numucc_n.push_back(vector<double>() );
    bnuecc_n.push_back(vector<double>() );
    nutaucc_n.push_back(vector<double>() );
    sig_n.push_back(vector<double>() );

    nc_i.push_back(vector<double>() );
    numucc_i.push_back(vector<double>() );
    bnuecc_i.push_back(vector<double>() );
    nutaucc_i.push_back(vector<double>() );
    sig_i.push_back(vector<double>() );

    for (k = 0; k < Extrap.size(); k++) {
      nc_n[j].push_back(0);
      numucc_n[j].push_back(0);
      bnuecc_n[j].push_back(0);
      nutaucc_n[j].push_back(0);
      sig_n[j].push_back(0);

      nc_i[j].push_back(0);
      numucc_i[j].push_back(0);
      bnuecc_i[j].push_back(0);
      nutaucc_i[j].push_back(0);
      sig_i[j].push_back(0);
    }
  }

  Bkgd = (TH1D *) nexp_n->Clone("Bkgd");
  Bkgd->Reset();
  Sig = (TH1D *) nexp_n->Clone("Sig");
  Sig->Reset();

  ofstream myfile;
  string file, ofile;

  // normal hierarchy

  TFile * f = new TFile(gSystem->ExpandPathName(outFileName.c_str()), "RECREATE");

  if (Print) {
    ofile = gSystem->ExpandPathName(outFileName.c_str());
    file = ofile.substr(0, ofile.length() - 5) + "_Normal.dat";
    myfile.open(gSystem->ExpandPathName(file.c_str()));
  }

  for (i = 0; i < nPts_Normal; i++) {
    cout << "point " << (i + 1) << "/" << nPts_Normal << " (normal hierarchy)" << endl;

    nexp_bkgd_n->Reset();
    nexp_signal_n->Reset();
    nexp_n->Reset();

    nexp_bkgd_i->Reset();
    nexp_signal_i->Reset();
    nexp_i->Reset();

    for (j = 0; j < nBins; j++) {
      GridTree_Normal[j]->GetEntry(i);
      nexp_bkgd_n->SetBinContent(j + 1, grid_background * GridScale_Normal);
      nexp_signal_n->SetBinContent(j + 1, grid_signal * GridScale_Normal);

      GridTree_Inverted[j]->GetEntry(i);
      nexp_bkgd_i->SetBinContent(j + 1, grid_background * GridScale_Inverted);
      nexp_signal_i->SetBinContent(j + 1, grid_signal * GridScale_Inverted);

      for (k = 0; k < Extrap.size(); k++) {
        GridTree_2_Normal[j][k]->GetEntry(i);
        nc_n[j][k] = grid_nc * GridScale_Normal;
        numucc_n[j][k] = grid_numucc * GridScale_Normal;
        bnuecc_n[j][k] = grid_bnuecc * GridScale_Normal;
        nutaucc_n[j][k] = grid_nutaucc * GridScale_Normal;
        sig_n[j][k] = grid_nue * GridScale_Normal;

        GridTree_2_Inverted[j][k]->GetEntry(i);
        nc_i[j][k] = grid_nc * GridScale_Inverted;
        numucc_i[j][k] = grid_numucc * GridScale_Inverted;
        bnuecc_i[j][k] = grid_bnuecc * GridScale_Inverted;
        nutaucc_i[j][k] = grid_nutaucc * GridScale_Inverted;
        sig_i[j][k] = grid_nue * GridScale_Inverted;
      }
    }
    nexp_n->Add(nexp_bkgd_n, nexp_signal_n, 1, 1);
    nexp_i->Add(nexp_bkgd_i, nexp_signal_i, 1, 1);

    // want to generate pseudos for normal hierarchy
    ErrCalc->SetGridPred(nBins, nc_n, numucc_n, bnuecc_n, nutaucc_n, sig_n);
    Error4Expts->Reset();
    ErrCalc->SetUseGrid(true);
    ErrCalc->CalculateSystErrorMatrix();
    Error4Expts->Add(ErrCalc->CovMatrix);
    ErrCalc->CalculateHOOError();
    Error4Expts->Add(ErrCalc->CovMatrix_Decomp);
    if (IncludeOscParErrs) {
      noff = 0;
      for (j = 0; j < nBins; j++) {
        GridTree_Normal[j]->GetEntry(i);
        ele = Error4Expts->GetBinContent(j + 1, j + 1);
        ele += (grid_bin_oscparerr[j] * grid_bin_oscparerr[j] * nexp_n->GetBinContent(j + 1) * nexp_n->GetBinContent(j + 1));
        Error4Expts->SetBinContent(j + 1, j + 1, ele);

        for (k = 0; k < nBins; k++) {
          if (k > j) {
            ele = Error4Expts->GetBinContent(j + 1, k + 1);
            ele += (grid_bin_oscparerr[j] * grid_bin_oscparerr[k] * nexp_n->GetBinContent(j + 1) * nexp_n->GetBinContent(k + 1));
            Error4Expts->SetBinContent(j + 1, k + 1, ele);

            ele = Error4Expts->GetBinContent(k + 1, j + 1);
            ele += (grid_bin_oscparerr[j] * grid_bin_oscparerr[k] * nexp_n->GetBinContent(j + 1) * nexp_n->GetBinContent(k + 1));
            Error4Expts->SetBinContent(k + 1, j + 1, ele);

            noff++;
          }
        }
      }
    }

    chi2hist_deltafit->Reset();
    chi2hist_deltafit->SetName(Form("Chi2Hist_DeltaFit_Normal_%i", i));

    chi2hist_mhfit->Reset();
    chi2hist_mhfit->SetName(Form("Chi2Hist_MHFit_Normal_%i", i));

    for (u = 0; u < NumExpts; u++) {
      cout << "expt " << (u + 1) << "/" << NumExpts << endl;

      GenerateOneCorrelatedExp(nexp_n, Error4Expts);
      if (Print) {
        myfile << grid_sinsq2th13 << " " << grid_delta << " ";
        for (j = 0; j < nBins; j++) {
          myfile << NObs->GetBinContent(j + 1) << " ";
        }
      }

      chi2min_norm = GetMinLikelihood_Delta(true);
      chi2min_invt = GetMinLikelihood_Delta(false);

      ErrCalc->SetGridPred(nBins, nc_n, numucc_n, bnuecc_n, nutaucc_n, sig_n);
      ErrCalc->SetUseGrid(true);
      chi2_norm = 1e10;
      if (FitMethod == 0) {
        chi2_norm = PoissonChi2(nexp_n);
      } else if (FitMethod == 1) {
        chi2_norm = ScaledChi2(nexp_bkgd_n, nexp_signal_n);
      } else if (FitMethod == 2) {
        chi2_norm = StandardChi2(nexp_n);
      } else if (FitMethod == 3) {
        // Likelihood: "Standard" (N syst, N nuisance)
        // Calculate the likelihood (x2 for chi)
        Bkgd->Reset();
        Bkgd->Add(nexp_bkgd_n);
        Sig->Reset();
        Sig->Add(nexp_signal_n);
        chi2_norm = StandardLikelihood();
      } else if (FitMethod == 4) {
        // Likelihood: Bin by Bin Calculation of Systematics
        // Calculate the likelihood (x2 for chi)
        Bkgd->Reset();
        Bkgd->Add(nexp_bkgd_n);
        Sig->Reset();
        Sig->Add(nexp_signal_n);
        chi2_norm = BinLikelihood();
      } else {
        cout << "Error in RunMultiBinPseudoExpts_MHDeltaFit(): Unknown 'FitMethod'." << endl;
      }

      ErrCalc->SetGridPred(nBins, nc_i, numucc_i, bnuecc_i, nutaucc_i, sig_i);
      ErrCalc->SetUseGrid(true);
      chi2_invt = 1e10;
      if (FitMethod == 0) {
        chi2_invt = PoissonChi2(nexp_i);
      } else if (FitMethod == 1) {
        chi2_invt = ScaledChi2(nexp_bkgd_i, nexp_signal_i);
      } else if (FitMethod == 2) {
        chi2_invt = StandardChi2(nexp_i);
      } else if (FitMethod == 3) {
        // Likelihood: "Standard" (N syst, N nuisance)
        // Calculate the likelihood (x2 for chi)
        Bkgd->Reset();
        Bkgd->Add(nexp_bkgd_i);
        Sig->Reset();
        Sig->Add(nexp_signal_i);
        chi2_invt = StandardLikelihood();
      } else if (FitMethod == 4) {
        // Likelihood: Bin by Bin Calculation of Systematics
        // Calculate the likelihood (x2 for chi)
        Bkgd->Reset();
        Bkgd->Add(nexp_bkgd_i);
        Sig->Reset();
        Sig->Add(nexp_signal_i);
        chi2_invt = BinLikelihood();
      } else {
        cout << "Error in RunMultiBinPseudoExpts_MHDeltaFit(): Unknown 'FitMethod'." << endl;
      }

      if (Print) {
        myfile << chi2_norm << " " << chi2_invt << " " << chi2min_norm << " " << chi2min_invt << endl;
      }

      delchi2_norm = chi2_norm - chi2min_norm;
      delchi2_invt = chi2_invt - chi2min_invt;
      delchi2_bestmh = delchi2_norm;
      if (delchi2_invt < delchi2_bestmh) {
        delchi2_bestmh = delchi2_invt;
      }

      chi2hist_deltafit->Fill(delchi2_norm);
      chi2hist_mhfit->Fill(delchi2_norm - delchi2_bestmh);
    }
    f->cd();
    chi2hist_deltafit->Write();
    chi2hist_mhfit->Write();
    f->Close();

    f = new TFile(gSystem->ExpandPathName(outFileName.c_str()), "UPDATE");
  }

  if (Print) {
    myfile.close();
  }

  for (j = 0; j < nBins; j++) {
    GridTree_Normal[j]->Write();

    for (k = 0; k < Extrap.size(); k++) {
      GridTree_2_Normal[j][k]->Write();
    }
  }

  f->Close();

  // inverted hierarchy

  f = new TFile(gSystem->ExpandPathName(outFileName.c_str()), "UPDATE");

  if (Print) {
    ofile = gSystem->ExpandPathName(outFileName.c_str());
    file = ofile.substr(0, ofile.length() - 5) + "_Inverted.dat";
    myfile.open(gSystem->ExpandPathName(file.c_str()));
  }

  for (i = 0; i < nPts_Inverted; i++) {
    cout << "point " << (i + 1) << "/" << nPts_Inverted << " (inverted hierarchy)" << endl;

    nexp_bkgd_n->Reset();
    nexp_signal_n->Reset();
    nexp_n->Reset();

    nexp_bkgd_i->Reset();
    nexp_signal_i->Reset();
    nexp_i->Reset();

    for (j = 0; j < nBins; j++) {
      GridTree_Normal[j]->GetEntry(i);
      nexp_bkgd_n->SetBinContent(j + 1, grid_background * GridScale_Normal);
      nexp_signal_n->SetBinContent(j + 1, grid_signal * GridScale_Normal);

      GridTree_Inverted[j]->GetEntry(i);
      nexp_bkgd_i->SetBinContent(j + 1, grid_background * GridScale_Inverted);
      nexp_signal_i->SetBinContent(j + 1, grid_signal * GridScale_Inverted);

      for (k = 0; k < Extrap.size(); k++) {
        GridTree_2_Normal[j][k]->GetEntry(i);
        nc_n[j][k] = grid_nc * GridScale_Normal;
        numucc_n[j][k] = grid_numucc * GridScale_Normal;
        bnuecc_n[j][k] = grid_bnuecc * GridScale_Normal;
        nutaucc_n[j][k] = grid_nutaucc * GridScale_Normal;
        sig_n[j][k] = grid_nue * GridScale_Normal;

        GridTree_2_Inverted[j][k]->GetEntry(i);
        nc_i[j][k] = grid_nc * GridScale_Inverted;
        numucc_i[j][k] = grid_numucc * GridScale_Inverted;
        bnuecc_i[j][k] = grid_bnuecc * GridScale_Inverted;
        nutaucc_i[j][k] = grid_nutaucc * GridScale_Inverted;
        sig_i[j][k] = grid_nue * GridScale_Inverted;
      }
    }
    nexp_n->Add(nexp_bkgd_n, nexp_signal_n, 1, 1);
    nexp_i->Add(nexp_bkgd_i, nexp_signal_i, 1, 1);

    // want to generate pseudos for inverted hierarchy
    ErrCalc->SetGridPred(nBins, nc_i, numucc_i, bnuecc_i, nutaucc_i, sig_i);
    Error4Expts->Reset();
    ErrCalc->SetUseGrid(true);
    ErrCalc->CalculateSystErrorMatrix();
    Error4Expts->Add(ErrCalc->CovMatrix);
    ErrCalc->CalculateHOOError();
    Error4Expts->Add(ErrCalc->CovMatrix_Decomp);
    if (IncludeOscParErrs) {
      noff = 0;
      for (j = 0; j < nBins; j++) {
        GridTree_Inverted[j]->GetEntry(i);
        ele = Error4Expts->GetBinContent(j + 1, j + 1);
        ele += (grid_bin_oscparerr[j] * grid_bin_oscparerr[j] * nexp_i->GetBinContent(j + 1) * nexp_i->GetBinContent(j + 1));
        Error4Expts->SetBinContent(j + 1, j + 1, ele);

        for (k = 0; k < nBins; k++) {
          if (k > j) {
            ele = Error4Expts->GetBinContent(j + 1, k + 1);
            ele += (grid_bin_oscparerr[j] * grid_bin_oscparerr[k] * nexp_i->GetBinContent(j + 1) * nexp_i->GetBinContent(k + 1));
            Error4Expts->SetBinContent(j + 1, k + 1, ele);

            ele = Error4Expts->GetBinContent(k + 1, j + 1);
            ele += (grid_bin_oscparerr[j] * grid_bin_oscparerr[k] * nexp_i->GetBinContent(j + 1) * nexp_i->GetBinContent(k + 1));
            Error4Expts->SetBinContent(k + 1, j + 1, ele);

            noff++;
          }
        }
      }
    }

    chi2hist_deltafit->Reset();
    chi2hist_deltafit->SetName(Form("Chi2Hist_DeltaFit_Inverted_%i", i));

    chi2hist_mhfit->Reset();
    chi2hist_mhfit->SetName(Form("Chi2Hist_MHFit_Inverted_%i", i));

    for (u = 0; u < NumExpts; u++) {
      cout << "expt " << (u + 1) << "/" << NumExpts << endl;

      GenerateOneCorrelatedExp(nexp_i, Error4Expts);
      if (Print) {
        myfile << grid_sinsq2th13 << " " << grid_delta << " ";
        for (j = 0; j < nBins; j++) {
          myfile << NObs->GetBinContent(j + 1) << " ";
        }
      }

      chi2min_norm = GetMinLikelihood_Delta(true);
      chi2min_invt = GetMinLikelihood_Delta(false);

      ErrCalc->SetGridPred(nBins, nc_n, numucc_n, bnuecc_n, nutaucc_n, sig_n);
      ErrCalc->SetUseGrid(true);
      chi2_norm = 1e10;
      if (FitMethod == 0) {
        chi2_norm = PoissonChi2(nexp_n);
      } else if (FitMethod == 1) {
        chi2_norm = ScaledChi2(nexp_bkgd_n, nexp_signal_n);
      } else if (FitMethod == 2) {
        chi2_norm = StandardChi2(nexp_n);
      } else if (FitMethod == 3) {
        // Likelihood: "Standard" (N syst, N nuisance)
        // Calculate the likelihood (x2 for chi)
        Bkgd->Reset();
        Bkgd->Add(nexp_bkgd_n);
        Sig->Reset();
        Sig->Add(nexp_signal_n);
        chi2_norm = StandardLikelihood();
      } else if (FitMethod == 4) {
        // Likelihood: Bin by Bin Calculation of Systematics
        // Calculate the likelihood (x2 for chi)
        Bkgd->Reset();
        Bkgd->Add(nexp_bkgd_n);
        Sig->Reset();
        Sig->Add(nexp_signal_n);
        chi2_norm = BinLikelihood();
      } else {
        cout << "Error in RunMultiBinPseudoExpts_MHDeltaFit(): Unknown 'FitMethod'." << endl;
      }

      ErrCalc->SetGridPred(nBins, nc_i, numucc_i, bnuecc_i, nutaucc_i, sig_i);
      ErrCalc->SetUseGrid(true);
      chi2_invt = 1e10;
      if (FitMethod == 0) {
        chi2_invt = PoissonChi2(nexp_i);
      } else if (FitMethod == 1) {
        chi2_invt = ScaledChi2(nexp_bkgd_i, nexp_signal_i);
      } else if (FitMethod == 2) {
        chi2_invt = StandardChi2(nexp_i);
      } else if (FitMethod == 3) {
        // Likelihood: "Standard" (N syst, N nuisance)
        // Calculate the likelihood (x2 for chi)
        Bkgd->Reset();
        Bkgd->Add(nexp_bkgd_i);
        Sig->Reset();
        Sig->Add(nexp_signal_i);
        chi2_invt = StandardLikelihood();
      } else if (FitMethod == 4) {
        // Likelihood: Bin by Bin Calculation of Systematics
        // Calculate the likelihood (x2 for chi)
        Bkgd->Reset();
        Bkgd->Add(nexp_bkgd_i);
        Sig->Reset();
        Sig->Add(nexp_signal_i);
        chi2_invt = BinLikelihood();
      } else {
        cout << "Error in RunMultiBinPseudoExpts_MHDeltaFit(): Unknown 'FitMethod'." << endl;
      }

      if (Print) {
        myfile << chi2_norm << " " << chi2_invt << " " << chi2min_norm << " " << chi2min_invt << endl;
      }

      delchi2_norm = chi2_norm - chi2min_norm;
      delchi2_invt = chi2_invt - chi2min_invt;
      delchi2_bestmh = delchi2_invt;
      if (delchi2_norm < delchi2_bestmh) {
        delchi2_bestmh = delchi2_norm;
      }

      chi2hist_deltafit->Fill(delchi2_invt);
      chi2hist_mhfit->Fill(delchi2_invt - delchi2_bestmh);
    }
    f->cd();
    chi2hist_deltafit->Write();
    chi2hist_mhfit->Write();
    f->Close();

    f = new TFile(gSystem->ExpandPathName(outFileName.c_str()), "UPDATE");
  }

  if (Print) {
    myfile.close();
  }

  for (j = 0; j < nBins; j++) {
    GridTree_Inverted[j]->Write();
    for (k = 0; k < Extrap.size(); k++) {
      GridTree_2_Inverted[j][k]->Write();
    }
  }

  f->Close();

  return;
}
void NueFit2D::RunMultiBinFC_MHDeltaFit() {
  if (NObs == 0) {
    cout << "NObs not set.  Quitting..." << endl;
    return;
  }
  if (Extrap.size() == 0) {
    cout << "No Extrapolate2D input.  Quitting..." << endl;
    return;
  }
  for (unsigned int ie = 0; ie < Extrap.size(); ie++) {
    Extrap[ie]->GetPrediction();
  }
  if (FitMethod == 1 && (FracErr_Bkgd == 0 || FracErr_Sig == 0)) {
    cout << "FracErr_Bkgd and FracErr_Sig need to be set for ScaledChi2.  Quitting..." << endl;
    return;
  }
  if (ErrCalc == 0) {
    cout << "No ErrorCalc object set!  Quitting..." << endl;
    return;
  }
  if (FitMethod == 3) {
    DefineStdDlnLMinuit();
  }
  if (FitMethod == 4) {
    DefineBinDlnLMinuit();
  }

  nBins = NObs->GetNbinsX();
  if (ErrorMatrix == 0) {
    ErrorMatrix = new TH2D("ErrorMatrix", "", nBins, -0.5, nBins - 0.5, nBins, -0.5, nBins - 0.5);
  }
  if (InvErrorMatrix == 0) {
    InvErrorMatrix = new TH2D("InvErrorMatrix", "", nBins, -0.5, nBins - 0.5, nBins, -0.5, nBins - 0.5);
  }

  ReadGridFiles();
  SetupChi2Hists();

  if (nPts_Normal == 0 || nPts_Inverted == 0) {
    return;
  }

  TH1D * nexp_bkgd = new TH1D("nexp_bkgd", "", nBins, -0.5, nBins - 0.5);
  TH1D * nexp_signal = new TH1D("nexp_signal", "", nBins, -0.5, nBins - 0.5);
  TH1D * nexp = new TH1D("nexp", "", nBins, -0.5, nBins - 0.5);

  int i;
  unsigned int j, k;

  if (gSystem->AccessPathName(gSystem->ExpandPathName(PseudoExpFile.c_str()))) {
    cout << "Pseudo-experiment file doesn't exist." << endl;
    return;
  }

  TFile * f = new TFile(gSystem->ExpandPathName(PseudoExpFile.c_str()), "READ");

  double * chi2data_NH = new double[nPts_Normal];
  double * chi2data_IH = new double[nPts_Inverted];
  double * chi2min_NH = new double[nPts_Normal];
  double * chi2min_IH = new double[nPts_Inverted];

  int chi2databin;
  TH1D * chi2hist = (TH1D *) f->Get("Chi2Hist_MHFit_Normal_0");
  double chi2binwidth = chi2hist->GetBinWidth(1);
  double chi2start = chi2hist->GetXaxis()->GetBinLowEdge(1);
  double frac;
  delete chi2hist;

  double * deltapts = new double[nPts_Normal];

  vector< vector<double> > nc, numucc, bnuecc, nutaucc, sig;
  for (j = 0; j < nBins; j++) {
    nc.push_back(vector<double>() );
    numucc.push_back(vector<double>() );
    bnuecc.push_back(vector<double>() );
    nutaucc.push_back(vector<double>() );
    sig.push_back(vector<double>() );
    for (k = 0; k < Extrap.size(); k++) {
      nc[j].push_back(0);
      numucc[j].push_back(0);
      bnuecc[j].push_back(0);
      nutaucc[j].push_back(0);
      sig[j].push_back(0);
    }
  }

  Bkgd = (TH1D *) NObs->Clone("Bkgd");
  Bkgd->Reset();
  Sig = (TH1D *) NObs->Clone("Sig");
  Sig->Reset();

  // normal hierarchy

  for (i = 0; i < nPts_Normal; i++) {
    if (i % 100 == 0) {
      cout << 100. * i / nPts_Normal << "% complete for normal hierarchy" << endl;
    }

    nexp_bkgd->Reset();
    nexp_signal->Reset();
    nexp->Reset();

    for (j = 0; j < nBins; j++) {
      GridTree_Normal[j]->GetEntry(i);
      nexp_bkgd->SetBinContent(j + 1, grid_background * GridScale_Normal);
      nexp_signal->SetBinContent(j + 1, grid_signal * GridScale_Normal);

      for (k = 0; k < Extrap.size(); k++) {
        GridTree_2_Normal[j][k]->GetEntry(i);
        nc[j][k] = grid_nc * GridScale_Normal;
        numucc[j][k] = grid_numucc * GridScale_Normal;
        bnuecc[j][k] = grid_bnuecc * GridScale_Normal;
        nutaucc[j][k] = grid_nutaucc * GridScale_Normal;
        sig[j][k] = grid_nue * GridScale_Normal;
      }
    }
    nexp->Add(nexp_bkgd, nexp_signal, 1, 1);
    ErrCalc->SetGridPred(nBins, nc, numucc, bnuecc, nutaucc, sig);

    chi2min_NH[i] = GetMinLikelihood_Delta(true);

    ErrCalc->SetUseGrid(true); // use grid predictions set up above
    chi2data_NH[i] = 1e10;
    if (FitMethod == 0) {
      chi2data_NH[i] = PoissonChi2(nexp);
    } else if (FitMethod == 1) {
      chi2data_NH[i] = ScaledChi2(Bkgd, Sig);
    } else if (FitMethod == 2) {
      chi2data_NH[i] = StandardChi2(nexp);
    } else if (FitMethod == 3) {
      // Likelihood: "Standard" (N syst, N nuisance)
      // Calculate the likelihood (x2 for chi)
      Bkgd->Reset();
      Bkgd->Add(nexp_bkgd);
      Sig->Reset();
      Sig->Add(nexp_signal);
      chi2data_NH[i] = StandardLikelihood();
    } else if (FitMethod == 4) {
      // Likelihood: Bin by Bin Calculation of Systematics
      // Calculate the likelihood (x2 for chi)
      Bkgd->Reset();
      Bkgd->Add(nexp_bkgd);
      Sig->Reset();
      Sig->Add(nexp_signal);
      chi2data_NH[i] = BinLikelihood();
    } else {
      cout << "Error in RunMultiBinFC_MHDeltaFit(): Unknown 'FitMethod'." << endl;
    }

    deltapts[i] = grid_delta / TMath::Pi();
  }

  // inverted hierarchy

  for (i = 0; i < nPts_Inverted; i++) {
    if (i % 100 == 0) {
      cout << 100. * i / nPts_Inverted << "% complete for inverted hierarchy" << endl;
    }

    nexp_bkgd->Reset();
    nexp_signal->Reset();
    nexp->Reset();

    for (j = 0; j < nBins; j++) {
      GridTree_Inverted[j]->GetEntry(i);
      nexp_bkgd->SetBinContent(j + 1, grid_background * GridScale_Inverted);
      nexp_signal->SetBinContent(j + 1, grid_signal * GridScale_Inverted);

      for (k = 0; k < Extrap.size(); k++) {
        GridTree_2_Inverted[j][k]->GetEntry(i);
        nc[j][k] = grid_nc * GridScale_Inverted;
        numucc[j][k] = grid_numucc * GridScale_Inverted;
        bnuecc[j][k] = grid_bnuecc * GridScale_Inverted;
        nutaucc[j][k] = grid_nutaucc * GridScale_Inverted;
        sig[j][k] = grid_nue * GridScale_Inverted;
      }
    }
    nexp->Add(nexp_bkgd, nexp_signal, 1, 1);
    ErrCalc->SetGridPred(nBins, nc, numucc, bnuecc, nutaucc, sig);

    chi2min_IH[i] = GetMinLikelihood_Delta(false);

    ErrCalc->SetUseGrid(true); // use grid predictions set up above
    chi2data_IH[i] = 1e10;
    if (FitMethod == 0) {
      chi2data_IH[i] = PoissonChi2(nexp);
    } else if (FitMethod == 1) {
      chi2data_IH[i] = ScaledChi2(Bkgd, Sig);
    } else if (FitMethod == 2) {
      chi2data_IH[i] = StandardChi2(nexp);
    } else if (FitMethod == 3) {
      // Likelihood: "Standard" (N syst, N nuisance)
      // Calculate the likelihood (x2 for chi)
      Bkgd->Reset();
      Bkgd->Add(nexp_bkgd);
      Sig->Reset();
      Sig->Add(nexp_signal);
      chi2data_IH[i] = StandardLikelihood();
    } else if (FitMethod == 4) {
      // Likelihood: Bin by Bin Calculation of Systematics
      // Calculate the likelihood (x2 for chi)
      Bkgd->Reset();
      Bkgd->Add(nexp_bkgd);
      Sig->Reset();
      Sig->Add(nexp_signal);
      chi2data_IH[i] = BinLikelihood();
    } else {
      cout << "Error in RunMultiBinFC_MHDeltaFit(): Unknown 'FitMethod'." << endl;
    }
  }

  if (nPts_Normal != nPts_Inverted) {
    cout << "different number of points" << endl;
  }

  double * deltafit_nh = new double[nPts_Normal];
  double * deltafit_ih = new double[nPts_Inverted];
  double * mhfit_nh = new double[nPts_Normal];
  double * mhfit_ih = new double[nPts_Inverted];
  double delchi2, delchi2_norm, delchi2_invt, delchi2_bestmh;

  // delta fit
  for (i = 0; i < nPts_Normal; i++) {
    // normal
    delchi2 = chi2data_NH[i] - chi2min_NH[i];
    chi2databin = int((delchi2 - chi2start) / chi2binwidth) + 1;

    TH1D * chi2hist = (TH1D *) f->Get(Form("Chi2Hist_DeltaFit_Normal_%i", i));

    if (chi2hist->Integral() < 1) {
      cout << "Warning, chi2hist is empty." << endl;
      frac = 0;
    } else {
      frac = chi2hist->Integral(1, chi2databin - 1) / chi2hist->Integral();
    }
    if (chi2databin == 1) {
      frac = 0;
    }

    deltafit_nh[i] = frac;

    // inverted
    delchi2 = chi2data_IH[i] - chi2min_IH[i];
    chi2databin = int((delchi2 - chi2start) / chi2binwidth) + 1;

    chi2hist = (TH1D *) f->Get(Form("Chi2Hist_DeltaFit_Inverted_%i", i));

    if (chi2hist->Integral() < 1) {
      cout << "Warning, chi2hist is empty." << endl;
      frac = 0;
    } else {
      frac = chi2hist->Integral(1, chi2databin - 1) / chi2hist->Integral();
    }
    if (chi2databin == 1) {
      frac = 0;
    }

    deltafit_ih[i] = frac;

    delete chi2hist;
  }

  // mh fit
  for (i = 0; i < nPts_Normal; i++) {
    delchi2_norm = chi2data_NH[i] - chi2min_NH[i];
    delchi2_invt = chi2data_IH[i] - chi2min_IH[i];
    delchi2_bestmh = delchi2_invt;
    if (delchi2_norm < delchi2_bestmh) {
      delchi2_bestmh = delchi2_norm;
    }
    cout << "NH: chi2min = " << chi2min_NH[i] << ", IH: chi2min = " << chi2min_IH[i] << endl;

    delchi2 = delchi2_norm - delchi2_bestmh;
    chi2databin = int((delchi2 - chi2start) / chi2binwidth) + 1;

    TH1D * chi2hist = (TH1D *) f->Get(Form("Chi2Hist_MHFit_Normal_%i", i));

    if (chi2hist->Integral() < 1) {
      cout << "Warning, chi2hist is empty." << endl;
      frac = 0;
    } else {
      frac = chi2hist->Integral(1, chi2databin - 1) / chi2hist->Integral();
    }
    if (chi2databin == 1) {
      frac = 0;
    }

    mhfit_nh[i] = frac;

    delchi2 = delchi2_invt - delchi2_bestmh;
    chi2databin = int((delchi2 - chi2start) / chi2binwidth) + 1;

    chi2hist = (TH1D *) f->Get(Form("Chi2Hist_MHFit_Inverted_%i", i));

    if (chi2hist->Integral() < 1) {
      cout << "Warning, chi2hist is empty." << endl;
      frac = 0;
    } else {
      frac = chi2hist->Integral(1, chi2databin - 1) / chi2hist->Integral();
    }
    if (chi2databin == 1) {
      frac = 0;
    }

    mhfit_ih[i] = frac;

    delete chi2hist;
  }

  f->Close();

  TGraph * g_delta_NH = new TGraph(nPts_Normal, deltapts, deltafit_nh);
  g_delta_NH->SetName("delta_NH");

  TGraph * g_delta_IH = new TGraph(nPts_Normal, deltapts, deltafit_ih);
  g_delta_IH->SetName("delta_IH");

  TGraph * g_MH_NH = new TGraph(nPts_Normal, deltapts, mhfit_nh);
  g_MH_NH->SetName("mh_NH");

  TGraph * g_MH_IH = new TGraph(nPts_Normal, deltapts, mhfit_ih);
  g_MH_IH->SetName("mh_IH");

  TFile * fout = new TFile(gSystem->ExpandPathName(outFileName.c_str()), "RECREATE");
  g_delta_NH->Write();
  g_delta_IH->Write();
  g_MH_NH->Write();
  g_MH_IH->Write();
  fout->Close();

  delete [] chi2data_NH;
  delete [] chi2data_IH;
  delete [] chi2min_NH;
  delete [] chi2min_IH;
  delete [] deltafit_nh;
  delete [] deltafit_ih;
  delete [] mhfit_nh;
  delete [] mhfit_ih;

  return;
}

// //////////////////////////////////////////////////////////////////////////////////////////////
// ////////////////////////////////Adam's Sterile Land Begins Here///////////////////////////////
// //////////////////////////////////////////////////////////////////////////////////////////////

void NueFit2D::RunSterileContour(int cl) {
  if (NObs == 0) {
    cout << "NObs not set.  Quitting..." << endl;
    return;
  }
  if (FitMethod == 1 && (FracErr_Bkgd == 0 || FracErr_Sig == 0)) {
    cout << "FracErr_Bkgd and FracErr_Sig need to be set for ScaledChi2.  Quitting..." << endl;
    return;
  }
  if (Extrap.size() == 0) {
    cout << "No Extrapolate2D input.  Quitting..." << endl;
    return;
  }

  if (FitMethod == 3) {
    DefineStdDlnLMinuit();
  }
  if (FitMethod == 4) {
    DefineBinDlnLMinuit();
  }

  Bkgd = (TH1D *) NObs->Clone("Bkgd");
  Bkgd->Reset();
  Sig = (TH1D *) NObs->Clone("Sig");
  Sig->Reset();
  NExp = (TH1D *) NObs->Clone("NExp");
  NExp->Reset();


  double sinsqth24 = TMath::Sin(Theta24) * TMath::Sin(Theta24);
  cout << Theta24 << " " << sinsqth24 << endl;

  unsigned int ie;
  for (ie = 0; ie < Extrap.size(); ie++) {
    Extrap[ie]->GetPrediction();
  }

  if (ErrorMatrix == 0) {
    ErrorMatrix = new TH2D("ErrorMatrix", "", nBins, -0.5, nBins - 0.5, nBins, -0.5, nBins - 0.5);
  }
  if (InvErrorMatrix == 0) {
    InvErrorMatrix = new TH2D("InvErrorMatrix", "", nBins, -0.5, nBins - 0.5, nBins, -0.5, nBins - 0.5);
  }

  Int_t i, j, dm;
  Double_t delta;
  Double_t Theta14;
  Double_t bth14;
  Double_t dm4[37];
  for (dm = 0; dm < 9; dm++) {
    dm4[dm] = .01 + dm * .01;
  }
  for (dm = 9; dm < 18; dm++) {
    dm4[dm] = .1 + (dm - 9) * .1;
  }
  for (dm = 18; dm < 27; dm++) {
    dm4[dm] = 1 + (dm - 18);
  }
  for (dm = 27; dm < 36; dm++) {
    dm4[dm] = 10 + 10 * (dm - 27);
  }
  dm4[36] = 100;

  Double_t ss2th14 = 0;
  Double_t Th14increment = 0;
  if (nSinSq2Th14Steps > 0) {
    Th14increment = (SinSq2Th14High - SinSq2Th14Low) / (nSinSq2Th14Steps);
  }

  double chi2[3], val[3];
  TGraph * g3;
  TF1 * fit;
  double best;
  double minchi2;

  double limit;
  double delchi2;
  double sprev, delchi2prev;
  double contourlvl = 0;
  if (cl == 0) { // 90% CL
    contourlvl = 2.71;
  } else if (cl == 1) { // 68.3% cL
    contourlvl = 1.0;
  } else {
    cout << "Error in RunDeltaChi2Contour(): Input value should be 0 or 1 for 90% or 68.3%.  Quitting..." << endl;
    return;
  }

  cout << "Seeking ";
  if (cl == 0) {
    cout << "90% ";
  } else { cout << "68% "; }
  cout << " CL upper limit" << endl;

  vector<double> deltapts;
  vector<double> bestfit_norm;
  vector<double> bestfit_invt;
  vector<double> limit_norm;
  vector<double> limit_invt;


  for (j = 0; j < 19; j++) {
    delta = dm4[j];
    for (ie = 0; ie < Extrap.size(); ie++) {
      Extrap[ie]->SetDm41(delta);

    }
    deltapts.push_back(delta);

    best = -1.;
    minchi2 = 100;
    for (i = 0; i < 3; i++) {
      chi2[i] = -1.;
      val[i] = -1.;
    }
    for (i = 0; i < nSinSq2Th14Steps + 1; i++) {
      chi2[0] = chi2[1];
      chi2[1] = chi2[2];
      val[0] = val[1];
      val[1] = val[2];
      ss2th14 = i * Th14increment + SinSq2Th14Low;
      for (ie = 0; ie < Extrap.size(); ie++) {
        Theta14 = TMath::ASin(TMath::Sqrt(ss2th14)) / 2.0;
        // sinsqth24 = (TMath::Tan(Theta14))*(TMath::Tan(Theta14));
        Extrap[ie]->SetSinSq2Th14(ss2th14);
        Extrap[ie]->SetSinSqTh24(sinsqth24);
        Extrap[ie]->OscillatePrediction();
      }
      Bkgd->Reset();
      Sig->Reset();
      NExp->Reset();
      for (ie = 0; ie < Extrap.size(); ie++) {
        Bkgd->Add(Extrap[ie]->Pred_TotalBkgd_VsBinNumber);
        Sig->Add(Extrap[ie]->Pred_Signal_VsBinNumber);
      }
      NExp->Add(Bkgd);
      NExp->Add(Sig);

      chi2[2] = 1e10;
      if (FitMethod == 0) {
        chi2[2] = PoissonChi2(NExp);
      } else if (FitMethod == 1) {
        chi2[2] = ScaledChi2(Bkgd, Sig);
      } else if (FitMethod == 2) {
        chi2[2] = StandardChi2(NExp);
      } else if (FitMethod == 3) {
        // Likelihood: "Standard" (N syst, N nuisance)
        // Calculate the likelihood (x2 for chi)
        chi2[2] = StandardLikelihood();
      } else if (FitMethod == 4) {
        // Likelihood: Bin by Bin Calculation of Systematics
        // Calculate the likelihood (x2 for chi)
        chi2[2] = BinLikelihood();
      } else {
        cout << "Error in RunDeltaChi2Contour(): Unknown 'FitMethod'." << endl;
      }

      val[2] = ss2th14;

      if (i < 2) {
        continue;
      }

      if (i < 3 && chi2[2] > chi2[1] && chi2[1] > chi2[0]) { // first three points are increasing, first point is minimum.
        best = val[0];
        minchi2 = chi2[0];
        cout << "minimum at 1st point: " << minchi2 << " at " << best << endl;
        bth14 = TMath::ASin(TMath::Sqrt(best)) / 2.0;
        // sinsqth24 = (TMath::Tan(bth14))*(TMath::Tan(bth14));
        break;
      }

      if (chi2[2] > chi2[1] && chi2[0] > chi2[1]) { // found minimum
        g3 = new TGraph(3, val, chi2);
        fit = new TF1("pol2", "pol2");
        g3->Fit(fit, "Q"); // fit to second order polynominal
        if (fit->GetParameter(2) > 0) { // if the x^2 term is nonzero
          best = -fit->GetParameter(1) / (2 * fit->GetParameter(2)); // the location of the minimum is -p1/(2*p2)
          minchi2 = fit->GetParameter(0) + fit->GetParameter(1) * best + fit->GetParameter(2) * best * best;
          cout << "minimum with fit: " << minchi2 << " at " << best << endl;
        } else { // if the x^2 term is zero, then just use the minimum you got by scanning
          best = val[1];
          minchi2 = chi2[1];
          cout << "minimum with scan: " << minchi2 << " at " << best << endl;
        }
        bth14 = TMath::ASin(TMath::Sqrt(best)) / 2.0;
        // sinsqth24 = (TMath::Tan(bth14))*(TMath::Tan(bth14));
        break;
      }
    }
    cout << sinsqth24 << endl;
    bestfit_norm.push_back(best);

    limit = 10000.;
    delchi2prev = 1000;
    sprev = 0;
    for (i = 0; i < nSinSq2Th14Steps + 1; i++) {
      ss2th14 = i * Th14increment + SinSq2Th14Low;
      for (ie = 0; ie < Extrap.size(); ie++) {
        Theta14 = TMath::ASin(TMath::Sqrt(ss2th14)) / 2.0;
        // sinsqth24 = (TMath::Tan(Theta14))*(TMath::Tan(Theta14));
        Extrap[ie]->SetSinSq2Th14(ss2th14);
        Extrap[ie]->SetSinSqTh24(sinsqth24);
        Extrap[ie]->OscillatePrediction();
      }
      Bkgd->Reset();
      Sig->Reset();
      NExp->Reset();
      for (ie = 0; ie < Extrap.size(); ie++) {
        Bkgd->Add(Extrap[ie]->Pred_TotalBkgd_VsBinNumber);
        Sig->Add(Extrap[ie]->Pred_Signal_VsBinNumber);
      }
      NExp->Add(Bkgd);
      NExp->Add(Sig);

      delchi2 = 1e10;
      if (FitMethod == 0) {
        delchi2 = PoissonChi2(NExp) - minchi2;
      } else if (FitMethod == 1) {
        delchi2 = ScaledChi2(Bkgd, Sig) - minchi2;
      } else if (FitMethod == 2) {
        delchi2 = StandardChi2(NExp) - minchi2;
      } else if (FitMethod == 3) {
        // Likelihood: "Standard" (N syst, N nuisance)
        // Calculate the likelihood (x2 for chi)
        delchi2 = StandardLikelihood() - minchi2;
      } else if (FitMethod == 4) {
        // Likelihood: Bin by Bin Calculation of Systematics
        // Calculate the likelihood (x2 for chi)
        delchi2 = BinLikelihood() - minchi2;
      } else {
        cout << "Error in RunDeltaChi2Contour(): Unknown 'FitMethod'." << endl;
      }

      if (i == 1) {
        continue;
      }

      if (delchi2 > contourlvl && delchi2prev < contourlvl) {
        limit = (ss2th14 + ((ss2th14 - sprev) / (delchi2 - delchi2prev)) * (contourlvl - delchi2));
        cout << delchi2prev << ", " << sprev << ", " << delchi2 << ", " << ss2th14 << ", " << limit << endl;
        bth14 = TMath::ASin(TMath::Sqrt(limit)) / 2.0;
        // sinsqth24 = (TMath::Tan(bth14))*(TMath::Tan(bth14));
        break;
      }
      delchi2prev = delchi2;
      sprev = ss2th14;
    }
    cout << sinsqth24 << endl;
    limit_norm.push_back(limit);
  }

  int nDeltaSteps = 18;
  double * s_sinsqth24 = new double[nDeltaSteps + 1];
  double * s_best_n = new double[nDeltaSteps + 1];
  double * s_limit_n = new double[nDeltaSteps + 1];
  double * s_LSND = new double[nDeltaSteps + 1];
  double * d = new double[nDeltaSteps + 1];
  for (i = 0; i < nDeltaSteps + 1; i++) {
    d[i] = deltapts.at(i);
    s_sinsqth24[i] = sinsqth24;
    s_best_n[i] = bestfit_norm.at(i);
    s_limit_n[i] = limit_norm.at(i);
    s_LSND[i] = sinsqth24 * limit_norm.at(i);

  }

  TGraph * gn = new TGraph(nDeltaSteps + 1, s_best_n, d);
  gn->SetMarkerStyle(20);
  gn->SetTitle("");
  gn->GetYaxis()->SetTitle("#Deltam^{2}_{41}");
  gn->GetYaxis()->SetRangeUser(0.01, 1);
  gn->SetMinimum(0.01);
  gn->SetMaximum(1);
  gn->GetXaxis()->SetTitle("sin^{2}2#theta_{14}");
  gn->GetXaxis()->SetLimits(0, 1);
  gn->SetLineWidth(4);
  gn->SetName("bfn");

  TGraph * gl = new TGraph(nDeltaSteps + 1, s_LSND, d);
  gl->SetMarkerStyle(20);
  gl->SetTitle("");
  gl->GetYaxis()->SetTitle("#Deltam^{2}_{41}");
  gl->GetYaxis()->SetRangeUser(0.01, 1);
  gl->SetMinimum(0.01);
  gl->SetMaximum(1);
  gl->GetXaxis()->SetTitle("sin^{2}#theta_{24}sin^{2}2#theta_{14}");
  gl->GetXaxis()->SetLimits(0, 1);
  gl->SetLineWidth(4);
  gl->SetName("lsnd");


  TGraph * gi = new TGraph(nDeltaSteps + 1, s_sinsqth24, d);
  gi->SetMarkerStyle(20);
  gi->SetTitle("");
  gi->GetYaxis()->SetTitle("#Deltam^{2}_{41}");
  gi->GetYaxis()->SetRangeUser(0.01, 1);
  gi->SetMinimum(0.01);
  gi->SetMaximum(1);
  gi->GetXaxis()->SetTitle("sin^{2}#theta_{24}");
  gi->GetXaxis()->SetLimits(0, 1);
  gi->SetLineWidth(4);
  gi->SetLineStyle(2);
  gi->SetName("sinsqth24");


  TGraph * gn_limit = new TGraph(nDeltaSteps + 1, s_limit_n, d);
  gn_limit->SetMarkerStyle(20);
  gn_limit->SetTitle("");
  gn_limit->GetYaxis()->SetTitle("#Deltam^{2}_{41}");
  gn_limit->GetYaxis()->SetRangeUser(0.01, 1);
  gn_limit->SetMinimum(0.01);
  gn_limit->SetMaximum(1);
  gn_limit->GetXaxis()->SetTitle("sin^{2}2#theta_{14}");
  gn_limit->GetXaxis()->SetLimits(0, 1);
  gn_limit->SetLineWidth(4);
  gn_limit->SetLineColor(kBlue);
  gn_limit->SetMarkerColor(kBlue);
  gn_limit->SetName("lmn");


  if (cl == 0) {
    cout << "90% ";
  }
  if (cl == 1) {
    cout << "68% ";
  }
  cout << "confidence level limit = " << limit_norm.at(0) << endl; // ", "<<limit_invt.at(0)<<endl;

  TFile * fout = new TFile(gSystem->ExpandPathName(outFileName.c_str()), "RECREATE");
  gn->Write();
  gl->Write();
  gn_limit->Write();
  gi->Write();
  fout->Write();
  fout->Close();

  delete [] s_best_n;
  delete [] s_limit_n;
  delete [] s_sinsqth24;
  delete [] s_LSND;
  delete [] d;

  return;
}


void NueFit2D::Run2DSterileSlice() {
  // Given a value of Dm41, produces a 2D chisquared graph of both sin(th14)^2 vs sin(th24)^2 and |Ue4|^2 vs |Umu4|^2



  if (NObs == 0) {
    cout << "NObs not set.  Quitting..." << endl;
    return;
  }
  if (FitMethod == 1 && (FracErr_Bkgd == 0 || FracErr_Sig == 0)) {
    cout << "FracErr_Bkgd and FracErr_Sig need to be set for ScaledChi2.  Quitting..." << endl;
    return;
  }
  if (Extrap.size() == 0) {
    cout << "No Extrapolate2D input.  Quitting..." << endl;
    return;
  }

  if (FitMethod == 3) {
    DefineStdDlnLMinuit();
  }
  if (FitMethod == 4) {
    DefineBinDlnLMinuit();
  }

  Bkgd = (TH1D *) NObs->Clone("Bkgd");
  Bkgd->Reset();
  Sig = (TH1D *) NObs->Clone("Sig");
  Sig->Reset();
  NExp = (TH1D *) NObs->Clone("NExp");
  NExp->Reset();


  unsigned int ie;
  for (ie = 0; ie < Extrap.size(); ie++) {
    Extrap[ie]->GetPrediction();
  }

  if (ErrorMatrix == 0) {
    ErrorMatrix = new TH2D("ErrorMatrix", "", nBins, -0.5, nBins - 0.5, nBins, -0.5, nBins - 0.5);
  }
  if (InvErrorMatrix == 0) {
    InvErrorMatrix = new TH2D("InvErrorMatrix", "", nBins, -0.5, nBins - 0.5, nBins, -0.5, nBins - 0.5);
  }

  Double_t Theta14;
  // Double_t Theta24;
  // Double_t D41;

  Int_t i, j, k, s;
  Int_t idx = 0; // Array indexer

  int arraysize = (nSinSqTh14Steps + 1) * (nSinSqTh24Steps + 1);
  cout << arraysize << endl;

  double Ue42[10201];
  double Um42[10201];
  double Ssqth24[10201];
  double chi[10201];
  double dchi[10201];

  double minchi2;
  double mc2;

  Double_t ssth14 = 0;
  Double_t csth14 = 0;
  Double_t Th14increment = 0;
  if (nSinSqTh14Steps > 0) {
    Th14increment = (SinSqTh14High - SinSqTh14Low) / (nSinSqTh14Steps);
  }

  Double_t ssth24 = 0;
  Double_t Th24increment = 0;
  if (nSinSqTh24Steps > 0) {
    Th24increment = (SinSqTh24High - SinSqTh24Low) / (nSinSqTh24Steps);
  }


  for (i = 0; i < nSinSqTh14Steps + 1; i++) {

    ssth14 = i * Th14increment + SinSqTh14Low;
    Theta14 = TMath::ASin(TMath::Sqrt(ssth14));
    csth14 = TMath::Cos(Theta14) * TMath::Cos(Theta14);


    for (ie = 0; ie < Extrap.size(); ie++) {
      Extrap[ie]->SetSinSqTh14(ssth14);
    }

    minchi2 = 100;
    mc2 = 1000;

    for (j = 0; j < nSinSqTh24Steps + 1; j++) {

      ssth24 = j * Th24increment + SinSqTh24Low;

      for (ie = 0; ie < Extrap.size(); ie++) {
        Extrap[ie]->SetSinSqTh24(ssth24);
        Extrap[ie]->OscillatePrediction();
      }

      Bkgd->Reset();
      Sig->Reset();
      NExp->Reset();

      for (ie = 0; ie < Extrap.size(); ie++) {
        Bkgd->Add(Extrap[ie]->Pred_TotalBkgd_VsBinNumber);
        Sig->Add(Extrap[ie]->Pred_Signal_VsBinNumber);
      }
      NExp->Add(Bkgd);
      NExp->Add(Sig);

      // delchi2 = 1e10;
      if (FitMethod == 0) {
        Ue42[idx] = ssth14;
        Um42[idx] = csth14 * ssth24;
        Ssqth24[idx] = ssth24;
        chi[idx] = PoissonChi2(NExp);
      } else if (FitMethod == 1) {
        Ue42[idx] = ssth14;
        Um42[idx] = csth14 * ssth24;
        Ssqth24[idx] = ssth24;
        chi[idx] = ScaledChi2(Bkgd, Sig);
      } else if (FitMethod == 2) {
        Ue42[idx] = ssth14;
        Um42[idx] = csth14 * ssth24;
        Ssqth24[idx] = ssth24;
        chi[idx] = StandardChi2(NExp);
      } else if (FitMethod == 3) {
        // Likelihood: "Standard" (N syst, N nuisance)
        // Calculate the likelihood (x2 for chi)
        Ue42[idx] = ssth14;
        Um42[idx] = csth14 * ssth24;
        Ssqth24[idx] = ssth24;
        chi[idx] = StandardLikelihood();
      } else if (FitMethod == 4) {
        // Likelihood: Bin by Bin Calculation of Systematics
        // Calculate the likelihood (x2 for chi)
        Ue42[idx] = ssth14;
        Um42[idx] = csth14 * ssth24;
        Ssqth24[idx] = ssth24;
        chi[idx] = BinLikelihood();
      } else {
        cout << "Error in Run2DSterileSlice(): Unknown 'FitMethod'." << endl;
      }

      idx++;

    } // TH24 loop


  } // TH14 loop

  double Ue42best[1];
  double Um42best[1];
  double Ssqth24best[1];

  for (k = 0; k < arraysize; k++) {
    if (chi[k] < mc2) {
      mc2 = chi[k];
      Ue42best[0] = Ue42[k];
      Um42best[0] = Um42[k];
      Ssqth24best[0] = Ssqth24[k];
    }
  }

  for (s = 0; s < arraysize; s++) {
    dchi[s] = chi[s] - mc2;
    if (abs(dchi[s]) < 0.00001) {
      cout << "I have a zero!" << endl;
    }

  }




  // make a TGraph
  TFile * w = new TFile(gSystem->ExpandPathName(outFileName.c_str()), "RECREATE");
  Int_t d = 1;

  TGraph * bpss = new TGraph(d, Ue42best, Ssqth24best);
  bpss->SetName("bpss");
  bpss->Write();

  TGraph * bpUU = new TGraph(d, Ue42best, Um42best);
  bpUU->SetName("bpUU");
  bpUU->Write();

  TGraph2D * ss = new TGraph2D(arraysize, Ue42, Ssqth24, chi);
  ss->SetName("ss");
  ss->GetXaxis()->SetTitle("sin^{2}#theta_{14}");
  ss->GetYaxis()->SetTitle("sin^{2}#theta_{24}");
  ss->GetZaxis()->SetTitle("-2lnL");
  ss->Write();

  TGraph2D * dss = new TGraph2D(arraysize, Ue42, Ssqth24, dchi);
  dss->SetName("dss");
  dss->GetXaxis()->SetTitle("sin^{2}#theta_{14}");
  dss->GetYaxis()->SetTitle("sin^{2}#theta_{24}");
  dss->GetZaxis()->SetTitle("-2#DeltalnL");
  dss->Write();

  TGraph2D * UU = new TGraph2D(arraysize, Ue42, Um42, chi);
  UU->SetName("UU");
  UU->GetXaxis()->SetTitle("sin^{2}#theta_{14}");
  UU->GetYaxis()->SetTitle("cos^{2}#theta_{14}sin^{2}#theta_{24}");
  UU->GetZaxis()->SetTitle("-2lnL");
  UU->Write();

  TGraph2D * dUU = new TGraph2D(arraysize, Ue42, Um42, dchi);
  dUU->SetName("dUU");
  dUU->GetXaxis()->SetTitle("sin^{2}#theta_{14}");
  dUU->GetYaxis()->SetTitle("cos^{2}#theta_{14}sin^{2}#theta_{24}");
  dUU->GetZaxis()->SetTitle("-2#DeltalnL");
  dUU->Write();

  w->Close();
  return;


} // 2D SterileSlice


/* Dung Phan (NueSterileSearch)
 */

void NueFit2D::RunSterileFit() {
  if (NObs == 0) {
    cout << "NObs not set.  Quitting..." << endl;
    return;
  }

  if (FitMethod == 1 && (FracErr_Bkgd == 0 || FracErr_Sig == 0)) {
    cout << "FracErr_Bkgd and FracErr_Sig need to be set for ScaledChi2.  Quitting..." << endl;
    return;
  }

  if (Extrap.size() == 0) {
    cout << "No Extrapolate2D input.  Quitting..." << endl;
    return;
  }

  if (FitMethod == 4) {
    DefineBinDlnLMinuit();
  }

  TH1D * nexp_bkgd = new TH1D("nexp_bkgd", "", nBins, -0.5, nBins - 0.5);
  TH1D * nexp_signal = new TH1D("nexp_signal", "", nBins, -0.5, nBins - 0.5);
  TH1D * nexp = new TH1D("nexp", "", nBins, -0.5, nBins - 0.5);


  TH1D * nexp_bkgd_min = new TH1D("nexp_bkgd_min", "", nBins, -0.5, nBins - 0.5);
  TH1D * nexp_signal_min = new TH1D("nexp_signal_min", "", nBins, -0.5, nBins - 0.5);
  TH1D * nexp_min = new TH1D("nexp_min", "", nBins, -0.5, nBins - 0.5);


  for (unsigned int ie = 0; ie < Extrap.size(); ie++) {
    Extrap[ie]->GetPrediction();
  }

  if (ErrorMatrix == 0) {
    ErrorMatrix = new TH2D("ErrorMatrix", "", nBins, -0.5, nBins - 0.5, nBins, -0.5, nBins - 0.5);
  }

  if (InvErrorMatrix == 0) {
    InvErrorMatrix = new TH2D("InvErrorMatrix", "", nBins, -0.5, nBins - 0.5, nBins, -0.5, nBins - 0.5);
  }

  unsigned int Nf = 1000; // number of fake experiments
  double * dchi = new double[Nf];
  TH1D * err = new TH1D("err", "", 1, -0.5, 0.5);            err->SetBinContent(1, 0.05);
  TH1D * sigerr = new TH1D("sigerr", "", 1, -0.5, 0.5);      sigerr->SetBinContent(1, 0.05);
  SetFracBkgdError(err);
  SetFracSigError(sigerr);


  gRandom->SetSeed(0);
  std::ofstream ofs;
  ostringstream strTh14, strTh24;
  strTh14 << SinSqTh14Low;
  strTh24 << SinSqTh24Low;
  // string filestr = "dm_" + strTh14.str() + "_" + strTh24.str();
  ofs.open("test.txt", std::ofstream::out | std::ofstream::app);
  for (unsigned int nf = 0; nf < Nf; nf++) {
    for (unsigned int ie = 0; ie < Extrap.size(); ie++) {
      Extrap[ie]->SetPrintResult();
      Extrap[ie]->SetSinSqTh14(SinSqTh14Low);
      Extrap[ie]->SetSinSqTh24(SinSqTh24Low);
      Extrap[ie]->SetOscPar(OscPar::kTh34, gRandom->Uniform(0., 0.5 * 3.1415));
      // cout << "Random: " << gRandom->Uniform(0., 0.5*3.1415) << endl;
      Extrap[ie]->SetOscPar(OscPar::kDelta, gRandom->Uniform(0., 2 * 3.1415));
      Extrap[ie]->SetOscPar(OscPar::kDelta14, gRandom->Uniform(0., 2 * 3.1415));
      Extrap[ie]->SetOscPar(OscPar::kDelta24, 0);
      Extrap[ie]->OscillatePrediction();
    }
    nexp_bkgd->Reset();
    nexp_signal->Reset();
    nexp->Reset();


    for (unsigned int ie = 0; ie < Extrap.size(); ie++) {
      nexp_bkgd->Add(Extrap[ie]->Pred_TotalBkgd_VsBinNumber);
      nexp_signal->Add(Extrap[ie]->Pred_Signal_VsBinNumber);
    }
    nexp->Add(nexp_bkgd);
    nexp->Add(nexp_signal);

    // Background only
    for (unsigned int ie = 0; ie < Extrap.size(); ie++) {
      Extrap[ie]->SetPrintResult();
      Extrap[ie]->SetSinSqTh14(0);
      Extrap[ie]->SetSinSqTh24(0);
      Extrap[ie]->SetOscPar(OscPar::kTh34, 0);
      // cout << "Random: " << gRandom->Uniform(0., 0.5*3.1415) << endl;
      Extrap[ie]->SetOscPar(OscPar::kDelta, 0);
      Extrap[ie]->SetOscPar(OscPar::kDelta14, 0);
      Extrap[ie]->SetOscPar(OscPar::kDelta24, 0);
      Extrap[ie]->OscillatePrediction();
    }
    nexp_bkgd_min->Reset();
    nexp_signal_min->Reset();
    nexp_min->Reset();

    if ((FracErr_Bkgd != 0) && (FracErr_Sig != 0)) {
      GenerateOneExperiment(nexp_bkgd, nexp_signal);
      dchi[nf] = BinLikelihood();
    }
    cout << "dchi[" << nf << "]: " << dchi[nf] << endl;
    ofs << dchi[nf] << endl;
  }
  ofs.close();


  return;
}


































/* Dung Code
*/

void NueFit2D::RunMultiBinPseudoExptsSterileFit(bool Print) {
  if (FitMethod == 1 && (FracErr_Bkgd == 0 || FracErr_Sig == 0)) {
    cout << "FracErr_Bkgd and FracErr_Sig need to be set for ScaledChi2.  Quitting..." << endl;
    return;
  }
  if (Extrap.size() == 0) {
    cout << "No Extrapolate2D input.  Quitting..." << endl;
    return;
  }
  for (unsigned int ie = 0; ie < Extrap.size(); ie++) {
    Extrap[ie]->GetPrediction();
  }
  if (ErrCalc == 0) {
    cout << "Need to set ErrorCalc object!  Quitting..." << endl;
    return;
  }
  if (FitMethod == 3) {
    DefineStdDlnLMinuit();
  }
  if (FitMethod == 4) {
    DefineBinDlnLMinuit();
  }

  nBins = Extrap[0]->Pred_TotalBkgd_VsBinNumber->GetNbinsX();
  if (ErrorMatrix == 0) {
    ErrorMatrix = new TH2D("ErrorMatrix", "", nBins, -0.5, nBins - 0.5, nBins, -0.5, nBins - 0.5);
  }
  if (InvErrorMatrix == 0) {
    InvErrorMatrix = new TH2D("InvErrorMatrix", "", nBins, -0.5, nBins - 0.5, nBins, -0.5, nBins - 0.5);
  }

  TH2D * Error4Expts = new TH2D("Error4Expts", "", nBins, -0.5, nBins - 0.5, nBins, -0.5, nBins - 0.5);

  ReadGridFilesSterileFit();

  if (nPts_Normal == 0) {
    return;
  }

  gRandom->SetSeed(0);

  TH1D * nexp_bkgd = new TH1D("nexp_bkgd", "", nBins, -0.5, nBins - 0.5);
  TH1D * nexp_signal = new TH1D("nexp_signal", "", nBins, -0.5, nBins - 0.5);
  TH1D * nexp = new TH1D("nexp", "", nBins, -0.5, nBins - 0.5);
  double delchi2, chi2min;
  TH1D * chi2hist = new TH1D("chi2hist", "", 110000, -10, 100);
  double ele;
  int noff;

  vector< vector<double> > nc, numucc, bnuecc, nutaucc, sig;
  for (unsigned int j = 0; j < nBins; j++) {
    nc.push_back(vector<double>() );
    numucc.push_back(vector<double>() );
    bnuecc.push_back(vector<double>() );
    nutaucc.push_back(vector<double>() );
    sig.push_back(vector<double>() );
    for (unsigned int k = 0; k < Extrap.size(); k++) {
      nc[j].push_back(0);
      numucc[j].push_back(0);
      bnuecc[j].push_back(0);
      nutaucc[j].push_back(0);
      sig[j].push_back(0);
    }
  }

  Bkgd = (TH1D *) nexp->Clone("Bkgd");
  Bkgd->Reset();
  Sig = (TH1D *) nexp->Clone("Sig");
  Sig->Reset();

  ofstream myfile;
  string file, ofile;

  // normal hierarchy
  TFile * f = new TFile(gSystem->ExpandPathName(outFileName.c_str()), "RECREATE");
  if (Print) {
    ofile = gSystem->ExpandPathName(outFileName.c_str());
    file = ofile.substr(0, ofile.length() - 5) + "_Normal.dat";
    myfile.open(gSystem->ExpandPathName(file.c_str()));
  }

  for (unsigned int i = 0; i < 1; i++) { // Testing with only one point
  //for (unsigned int i = 0; i < nPts_Normal; i++) {
    cout << "point " << (i + 1) << "/" << nPts_Normal << " (normal hierarchy)" << endl;

    nexp_bkgd->Reset();
    nexp_signal->Reset();
    nexp->Reset();

    for (unsigned int j = 0; j < nBins; j++) {
      GridTree_Normal[j]->GetEntry(i);
      //nexp_bkgd->SetBinContent(j + 1, grid_background * GridScale_Normal);
      nexp_bkgd->SetBinContent(j + 1, grid_background);
      nexp_signal->SetBinContent(j + 1, grid_signal);

      for (unsigned int k = 0; k < Extrap.size(); k++) {
        GridTree_2_Normal[j][k]->GetEntry(i);
        nc[j][k] = grid_nc;
        numucc[j][k] = grid_numucc;
        bnuecc[j][k] = grid_bnuecc;
        nutaucc[j][k] = grid_nutaucc;
        sig[j][k] = grid_nue;
      }
    }
    nexp->Add(nexp_bkgd, nexp_signal, 1, 1);
    ErrCalc->SetGridPred(nBins, nc, numucc, bnuecc, nutaucc, sig);

    Error4Expts->Reset();
    ErrCalc->SetUseGrid(true);
    ErrCalc->CalculateSystErrorMatrix();
    Error4Expts->Add(ErrCalc->CovMatrix);
    ErrCalc->CalculateHOOError();
    Error4Expts->Add(ErrCalc->CovMatrix_Decomp);
    if (IncludeOscParErrs) {
      noff = 0;
      for (unsigned int j = 0; j < nBins; j++) {
        ele = Error4Expts->GetBinContent(j + 1, j + 1);
        ele += (grid_bin_oscparerr[j] * grid_bin_oscparerr[j] * nexp->GetBinContent(j + 1) * nexp->GetBinContent(j + 1));
        Error4Expts->SetBinContent(j + 1, j + 1, ele);

        for (unsigned int k = 0; k < nBins; k++) {
          if (k > j) {
            ele = Error4Expts->GetBinContent(j + 1, k + 1);
            ele += (grid_bin_oscparerr[j] * grid_bin_oscparerr[k] * nexp->GetBinContent(j + 1) * nexp->GetBinContent(k + 1));
            Error4Expts->SetBinContent(j + 1, k + 1, ele);

            ele = Error4Expts->GetBinContent(k + 1, j + 1);
            ele += (grid_bin_oscparerr[j] * grid_bin_oscparerr[k] * nexp->GetBinContent(j + 1) * nexp->GetBinContent(k + 1));
            Error4Expts->SetBinContent(k + 1, j + 1, ele);

            noff++;
          }
        }
      }
    }

    chi2hist->Reset();
    chi2hist->SetName(Form("Chi2Hist_Normal_%i", i));

    for (unsigned int u = 0; u < NumExpts; u++) {
      cout << "expt " << (u + 1) << "/" << NumExpts << endl;

      GenerateOneCorrelatedExp(nexp, Error4Expts);
      if (Print) {
        myfile << grid_sinsqth14 << " " << grid_sinsqth24 << " ";
        for (unsigned int j = 0; j < nBins; j++) {
          myfile << NObs->GetBinContent(j + 1) << " ";
        }
        myfile << endl;
      }

      vector<double> SSS;
      SSS.push_back(grid_sinsqth14);
      SSS.push_back(grid_sinsqth24);
      chi2min = DoGlobalMinSearchSterileFit(SSS);

      ErrCalc->SetUseGrid(true); // will use the grid predictions set above
      delchi2 = 1e10;
      if (FitMethod == 0) {
        delchi2 = PoissonChi2(nexp) - chi2min;
      } else if (FitMethod == 1) {
        delchi2 = ScaledChi2(nexp_bkgd, nexp_signal) - chi2min;
      } else if (FitMethod == 2) {
        delchi2 = StandardChi2(nexp) - chi2min;
      } else if (FitMethod == 3) {
        // Likelihood: "Standard" (N syst, N nuisance)
        // Calculate the likelihood (x2 for chi)
        Bkgd->Reset();
        Bkgd->Add(nexp_bkgd);
        Sig->Reset();
        Sig->Add(nexp_signal);
        delchi2 = StandardLikelihood() - chi2min;
      } else if (FitMethod == 4) {
        // Likelihood: Bin by Bin Calculation of Systematics
        // Calculate the likelihood (x2 for chi)
        Bkgd->Reset();
        Bkgd->Add(nexp_bkgd);
        Sig->Reset();
        Sig->Add(nexp_signal);
        delchi2 = BinLikelihood() - chi2min;
	cout << "Delta Chisquare min pseudo-exp " << u << " of point " << i << " is " << chi2min << endl;
	cout << "Delta Chisquare for pseudo-exp " << u << " of point " << i << " is " << delchi2 << endl;
      } else {
        cout << "Error in RunMultiBinPseudoExpts(): Unknown 'FitMethod'." << endl;
      }
      chi2hist->Fill(delchi2);

    }
    f->cd();
    chi2hist->Write();
    f->Close();

    f = new TFile(gSystem->ExpandPathName(outFileName.c_str()), "UPDATE");
  }

  if (Print) {
    myfile.close();
  }

  for (unsigned int j = 0; j < nBins; j++) {
    GridTree_Normal[j]->Write();

    for (unsigned int k = 0; k < Extrap.size(); k++) {
      GridTree_2_Normal[j][k]->Write();
    }
  }

  f->Close();

  return;
}


void NueFit2D::ReadGridFilesSterileFit() {
  double fp;
  TTree * temp;

  grid_bin_oscparerr.clear();
  for (unsigned int i = 0; i < nBins; i++) {
    grid_bin_oscparerr.push_back(0);
  }

  TFile * fnorm;
  GridTree_Normal.clear();
  GridTree_2_Normal.clear();
  nPts_Normal = 0;
  if (!gSystem->AccessPathName(gSystem->ExpandPathName(GridFileName_Normal.c_str()))) { // if file exists
    fnorm = new TFile(gSystem->ExpandPathName(GridFileName_Normal.c_str()), "READ");
    for (unsigned int i = 0; i < nBins; i++) {
      temp = (TTree *) fnorm->Get(Form("Bin_%i", i));
      temp->SetName(Form("Bin_%i_Normal", i));
      temp->SetBranchAddress("Background", &grid_background);
      temp->SetBranchAddress("Signal", &grid_signal);
      temp->SetBranchAddress("SinSq2Th14", &grid_sinsqth14);
      temp->SetBranchAddress("SinSq2Th24", &grid_sinsqth24);
      temp->SetBranchAddress("Dmsq41", &grid_dmsq41);
      temp->SetBranchAddress("Delta13", &grid_delta13);
      temp->SetBranchAddress("Delta14", &grid_delta14);
      temp->SetBranchAddress("Delta24", &grid_delta24);
      temp->SetBranchAddress("Theta34", &grid_th34);
      if (IncludeOscParErrs) {
        temp->SetBranchAddress("DNExp_DOscPars", &grid_bin_oscparerr[i]);
      }
      GridTree_Normal.push_back(temp);

      GridTree_2_Normal.push_back(vector<TTree *>() );

      for (unsigned int j = 0; j < Extrap.size(); j++) {
        temp = (TTree *) fnorm->Get(Form("Bin_%i_Run_%i", i, j));
        temp->SetName(Form("Bin_%i_Run_%i_Normal", i, j));
        temp->SetBranchAddress("NC", &grid_nc);
        temp->SetBranchAddress("NuMuCC", &grid_numucc);
        temp->SetBranchAddress("BNueCC", &grid_bnuecc);
        temp->SetBranchAddress("NuTauCC", &grid_nutaucc);
        temp->SetBranchAddress("Signal", &grid_nue);
        temp->SetBranchAddress("SinSq2Th14", &grid_sinsqth14);
        temp->SetBranchAddress("SinSq2Th24", &grid_sinsqth24);
        temp->SetBranchAddress("Dmsq41", &grid_dmsq41);
        temp->SetBranchAddress("Delta13", &grid_delta13);
        temp->SetBranchAddress("Delta14", &grid_delta14);
        temp->SetBranchAddress("Delta24", &grid_delta24);
        temp->SetBranchAddress("Theta34", &grid_th34);
        GridTree_2_Normal[i].push_back(temp);
      }
    }
    nPts_Normal = GridTree_Normal[0]->GetEntries();

    paramtree_Normal = (TTree *) fnorm->Get("paramtree");
    paramtree_Normal->SetName("paramtree_Normal");
    paramtree_Normal->SetBranchAddress("farPOT", &fp);
    paramtree_Normal->SetBranchAddress("Theta12", &grid_n_th12);
    paramtree_Normal->SetBranchAddress("Theta13", &grid_n_th13);
    paramtree_Normal->SetBranchAddress("Theta23", &grid_n_th23);
    paramtree_Normal->SetBranchAddress("DeltaMSq23", &grid_n_dm2_32);
    paramtree_Normal->SetBranchAddress("DeltaMSq12", &grid_n_dm2_21);
    paramtree_Normal->SetBranchAddress("DeltaMSq41", &grid_dmsq41);
    paramtree_Normal->SetBranchAddress("SinSq2Th14", &grid_sinsqth14);
    paramtree_Normal->SetBranchAddress("SinSq2Th24", &grid_sinsqth24);
    paramtree_Normal->GetEntry(0);
  } else {
    cout << "Grid file (normal hierarchy) doesn't exist." << endl;
    return;
  }

  return;
}

/*
double NueFit2D::GetMinLikelihoodSterileFit() {
  if (NObs == 0) {
    cout << "NObs not set.  Quitting..." << endl;
    return 0;
  }
  if (FitMethod == 1 && (FracErr_Bkgd == 0 || FracErr_Sig == 0)) {
    cout << "FracErr_Bkgd and FracErr_Sig need to be set for ScaledChi2.  Quitting..." << endl;
    return 0;
  }
  if (Extrap.size() == 0) {
    cout << "No Extrapolate2D input.  Quitting..." << endl;
    return 0;
  }

  if (FitMethod == 3) {
    DefineStdDlnLMinuit();
  }
  if (FitMethod == 4) {
    DefineBinDlnLMinuit();
  }
  if (FitMethod == 3 || FitMethod == 4) {
    if (ErrCalc != 0) {
      ErrCalc->SetUseGrid(false);
    }
  }

  Bkgd = (TH1D *) NObs->Clone("Bkgd");
  Bkgd->Reset();
  Sig = (TH1D *) NObs->Clone("Sig");
  Sig->Reset();
  TH1D * NExp = (TH1D *) NObs->Clone("NExp");
  NExp->Reset();

  for (unsigned int ie = 0; ie < Extrap.size(); ie++) {
    Extrap[ie]->GetPrediction();
    Extrap[ie]->SetOscPar(OscPar::kTh12, grid_n_th12);
    Extrap[ie]->SetOscPar(OscPar::kTh13, grid_n_th13);
    Extrap[ie]->SetOscPar(OscPar::kTh23, grid_n_th23);
    Extrap[ie]->SetOscPar(OscPar::kDeltaM12, grid_n_dm2_21);
    Extrap[ie]->SetOscPar(OscPar::kDeltaM23, grid_n_dm2_32);
    Extrap[ie]->SetOscPar(OscPar::kDm41, grid_dmsq41);
    Extrap[ie]->SetOscPar(OscPar::kTh34, grid_th34);
    Extrap[ie]->SetOscPar(OscPar::kDelta, grid_delta13);
    Extrap[ie]->SetOscPar(OscPar::kDelta14, grid_delta14);
    Extrap[ie]->SetOscPar(OscPar::kDelta24, grid_delta24);
    Extrap[ie]->OscillatePrediction();
  }

  if (ErrorMatrix == 0) {
    ErrorMatrix = new TH2D("ErrorMatrix", "", nBins, -0.5, nBins - 0.5, nBins, -0.5, nBins - 0.5);
  }
  if (InvErrorMatrix == 0) {
    InvErrorMatrix = new TH2D("InvErrorMatrix", "", nBins, -0.5, nBins - 0.5, nBins, -0.5, nBins - 0.5);
  }

  double theLikelihood = 1e10;
  double minSinSqTh14;
  double minSinSqTh24;

  for (unsigned int i = 0; i < 11; i++) {
    double gridSinSTH14 = i / 10;
    for (unsigned int j = 0; j < 11; j++) {
      double gridSinSTH24 = j / 10;
      for (unsigned int ie = 0; ie < Extrap.size(); ie++) {
        Extrap[ie]->SetOscPar(OscPar::kTh14, TMath::ASin(TMath::Sqrt(gridSinSTH14)));
        Extrap[ie]->SetOscPar(OscPar::kTh24, TMath::ASin(TMath::Sqrt(gridSinSTH24)));
        Extrap[ie]->OscillatePrediction();
      }
      Bkgd->Reset();
      Sig->Reset();
      NExp->Reset();
      for (unsigned int ie = 0; ie < Extrap.size(); ie++) {
        Bkgd->Add(Extrap[ie]->Pred_TotalBkgd_VsBinNumber);
        Sig->Add(Extrap[ie]->Pred_Signal_VsBinNumber);
      }
      NExp->Add(Bkgd);
      NExp->Add(Sig);

      double c2;
      switch (FitMethod) {
        case 0: {
          c2 = PoissonChi2(NExp);
          break;
        }
        case 1: {
          c2 = ScaledChi2(Bkgd, Sig);
          break;
        }
        case 2: {
          c2 = StandardChi2(NExp);
          break;
        }
        case 3: {
          c2 = StandardLikelihood();
          break;
        }
        case 4: {
          c2 = BinLikelihood();
          break;
        }
        default: {
          cout << "Error in GetMinLikelihoodSterileFit(): Unknown 'FitMethod'. BinLikelihood() is used." << endl;
          c2 = BinLikelihood();
          break;
        }
      }

      if (c2 < theLikelihood) {
        theLikelihood = c2;
        minSinSqTh14 = gridSinSTH14;
        minSinSqTh24 = gridSinSTH24;
      }
    }
  }


  for (int i = -5; i < 6; i++) {
    double gridSinSTH14 = grid_sinsqth14 + i * 0.01;
    for (int j = -5; j < 6; j++) {
      double gridSinSTH24 = grid_sinsqth24 + j * 0.01;
      for (unsigned int ie = 0; ie < Extrap.size(); ie++) {
        Extrap[ie]->SetOscPar(OscPar::kTh14, TMath::ASin(TMath::Sqrt(gridSinSTH14)));
        Extrap[ie]->SetOscPar(OscPar::kTh24, TMath::ASin(TMath::Sqrt(gridSinSTH24)));
        Extrap[ie]->OscillatePrediction();
      }
      Bkgd->Reset();
      Sig->Reset();
      NExp->Reset();
      for (unsigned int ie = 0; ie < Extrap.size(); ie++) {
        Bkgd->Add(Extrap[ie]->Pred_TotalBkgd_VsBinNumber);
        Sig->Add(Extrap[ie]->Pred_Signal_VsBinNumber);
      }
      NExp->Add(Bkgd);
      NExp->Add(Sig);

      double c2;
      switch (FitMethod) {
        case 0: {
          c2 = PoissonChi2(NExp);
          break;
        }
        case 1: {
          c2 = ScaledChi2(Bkgd, Sig);
          break;
        }
        case 2: {
          c2 = StandardChi2(NExp);
          break;
        }
        case 3: {
          c2 = StandardLikelihood();
          break;
        }
        case 4: {
          c2 = BinLikelihood();
          break;
        }
        default: {
          cout << "Error in GetMinLikelihoodSterileFit(): Unknown 'FitMethod'. BinLikelihood() is used." << endl;
          c2 = BinLikelihood();
          break;
        }
      }

      if (c2 < theLikelihood) {
        theLikelihood = c2;
        minSinSqTh14 = gridSinSTH14;
        minSinSqTh24 = gridSinSTH24;
      }
    }
  }
  cout << "The min likelihood is " << theLikelihood << " at " << minSinSqTh14 << " - " << minSinSqTh24 << endl;

  return theLikelihood;
}
*/
