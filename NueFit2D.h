#ifndef NueFit2D_h
#define NueFit2D_h

#include "TROOT.h"
#include "NueAna/MultiBinAna/Extrapolate2D.h"
#include "NueAna/MultiBinAna/ErrorCalc.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TMath.h"
#include "TCanvas.h"
#include <vector>
#include "TMatrixD.h"
#include <string>
#include "TMinuit.h"
#include "TObject.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include <iostream>
#include "TSystem.h"
#include "TRandom.h"
#include "TF1.h"
#include "TMatrixDEigen.h"

using namespace std;

double Poisson(double mean, int n, bool numOnly = false);
double Gaussian(double x, double mean, double sig);

class NueFit2D: public TObject {
  public:
    NueFit2D();
    virtual ~NueFit2D();

    void RunScaledChi2Sensitivity();
    void RunStandardChi2Sensitivity();
    void RunFCAnalytical();
    void RunPseudoExperiments();

    virtual void RunMultiBinPseudoExpts(bool Print=true);
    void RunFCTraditional();
    void RunMultiBinFC();
    virtual double GetSensitivityAt(double delta=0, bool normalhier=true);//returns sin^2 2th13 90% CL upper limit for given values of delta and the hierarchy
    virtual void RunDeltaChi2Contour(int cl=0);//0 = 90%, 1 = 68.4%
    virtual void RunSterileContour(int cl=0); // 0 = 90%, 1 = 68.4%
    virtual void Run2DSterileSlice();

    virtual double GetLikelihood(double t12=0.6,double t23=0.785398,double t13=0,double dm2_32=2.32e-3,double dm2_21=7.59e-5,double delta=0);//gets the likelihood value for NObs and prediction at these osc par values
    virtual void RunDataGrid(string filenorm,string fileinvt);

    virtual void AddExtrap(Extrapolate2D* E);
    void SetNObs(TH1D *n);//x axis is 2D bin number, content is number of events in that bin

    void SetFracBkgdError(TH1D *n);//x axis is 2D bin number, content is fractional error on number of events
    void SetFracSigError(TH1D *n);//x axis is 2D bin number, content is fractional error on number of events

    void SetSystErrorMatrix(TH2D *n);//SYSTEMATIC ONLY error matrix where nij is covariance of 2D bin number i and 2D bin number j
    void SetIncludeOscParErrs() { IncludeOscParErrs = true; };
    void SetOutputFile(string s="SensOut.root");
    void SetNExperiments(int n = 2000) { NumExpts = n; };
    void FormatChi2Hists(int nx=504,double xlow=-0.0005,double xhigh=0.5035,int ny=101,double ylow=-0.01,double yhigh=2.01);
    void SetPseudoExperimentInputFile(string s="PseudoExp.root") { PseudoExpFile = s; };
    void SetGridFiles(string snorm="Grid_1stAna_Norm.root",string sinvt="Grid_1stAna_Invt.root");
    void SetGridNorm(double n=0) { GridNorm = n; };//in case you want a different normalization than the grid files were made with

    void SetNDeltaSteps(int n=20) { nDeltaSteps = n; };
    void SetNSinSq2Th13Steps(int n=3000) { nSinSq2Th13Steps = n; };
    void SetNSinSq2Th14Steps(int n=3000) { nSinSq2Th14Steps = n; };
    void SetDeltaRange(double l=0, double h=2);//in units of pi
    void SetSinSq2Th13Range(double l=0.05, double h = 0.35);
    void SetSinSq2Th14Range(double l=0.001, double h = 0.301);

    void SetNSinSqTh14Steps(int n=100) { nSinSqTh14Steps = n; };
    void SetNSinSqTh24Steps(int n=100) { nSinSqTh24Steps = n; };
    void SetSinSqTh14Range(double l=0.0, double h = 1.0 );
    void SetSinSqTh24Range(double l=0.0, double h = 1.0 );

    void SetNepsetauSteps(int n=3000) { nepsetauSteps = n; };
    void SetepsetauRange(double l=-4.0, double h=4.0);

    void SetFitMethod(int m=0) { FitMethod = m; };//0 = poisson, 1 = scaled chi2, 2 = standard chi2, 3 = standard (nuisance param for each syst) likelihood, 4 = bin by bin likelihood

    void AddError(ErrorCalc *Err); //Set the errorcalc object.
    void SetFracError(vector<TH1D*> bkglist, vector<TH1D*> siglist);//This will go directly to Error object

    double StdLikeComparison(vector<double> npar);
    double BinLikeComparison(vector<double> npar);

    virtual void RunMultiBinPseudoExpts_MHDeltaFit(bool Print=true);
    virtual void RunMultiBinFC_MHDeltaFit();



    //Delta Fit Additions
    Double_t Theta12;
    Double_t Theta12Unc;
    Double_t Theta23;
    Double_t Theta23Unc;
    Double_t Theta13;
    Double_t Theta13Unc;
    Double_t DeltaMSq32;
    Double_t DeltaMSq32Unc;
    Double_t DeltaMSq21;
    Double_t DeltaMSq21Unc;
    Double_t delta;
    //NSI
    double Eps_ee;
    double Eps_emu;
    double Eps_etau;
    double Eps_etau_unc;
    double Eps_mumu;
    double Eps_mutau;
    double Eps_tautau;
    double Delta_emu;
    double Delta_etau;
    double Delta_etau_unc;
    double Delta_mutau;
    //Sterile
    double Theta24;
    double Theta14;
    double DeltaMSq41;

    void DefTheta13(double t13, double t13unc){
      Theta13 = t13;
      Theta13Unc = t13unc;
      return;
    }

    void DefTheta23(double t23, double t23unc){
      Theta23 = t23;
      Theta23Unc = t23unc;
      return;
    }

    void DefTheta12(double t12, double t12unc){
      Theta12 = t12;
      Theta12Unc = t12unc;
      return;
    }

    void DefDMS32(double dms32, double dms32unc){
      DeltaMSq32 = dms32;
      DeltaMSq32Unc = dms32unc;
      return;
    }

    void DefDMS21(double dms21, double dms21unc){
      DeltaMSq21 = dms21;
      DeltaMSq21Unc = dms21unc;
      return;
    }

    void DefDeltaCP(double deltaCP){
      delta = deltaCP;
      return;
    }

    void DefTheta24(double t24){
      Theta24 = t24;
      return;
    }
    void DefTheta14(double t14){
      Theta14 = t14;
      return;
    }
    void DefDMS41(double dms41){
      DeltaMSq41 = dms41;
      return;
    }


    // NSI


    void DefEps_ee(double eps_ee){
      Eps_ee = eps_ee;
      return;
    }

    void DefEps_emu(double eps_emu){
      Eps_emu = eps_emu;
      return;
    }

    void DefEps_etau(double eps_etau, double eps_etau_unc){
      Eps_etau = eps_etau;
      Eps_etau_unc = eps_etau_unc;
      return;
    }

    void DefDelta_etau(double delta_etau, double delta_etau_unc){
      Delta_etau = delta_etau;
      Delta_etau_unc = delta_etau_unc;
      return;
    }


  protected:

    double DoStdMinParam();
    double StandardLikelihood();
    virtual void DefineStdDlnLMinuit();

    double DoBinMinParam();
    double BinLikelihood();
    virtual void DefineBinDlnLMinuit();

    virtual void CalculateDlnLMatrix(TH2D *SystMatrix, TH2D *HOOHEMatrix);

    ErrorCalc *ErrCalc;
    std::vector<TH1D*> FracErr_Bkgd_List;
    std::vector<TH1D*> FracErr_Sig_List;
    TMinuit *minuit;

    virtual void ReadGridFiles();
    void SetupChi2Hists();

    void GenerateOneExperiment(TH1D *nexp_bkgd,TH1D *nexp_signal,TH1D *dnexp_oscpar=0);
    void GenerateOneCorrelatedExp(TH1D *nexp,TH2D *err);//generates a pseudo-experiment including systematic correlations based on the error matrix

    double ScaledChi2(TH1D *nexp_bkgd,TH1D *nexp_signal);
    double PoissonChi2(TH1D *nexp);
    double StandardChi2(TH1D *nexp);
    void CalculateErrorMatrixInv(TH1D *nexp);

    virtual double GetMinLikelihood(double delta=0,bool normalhier=true);//used only within FC code - gets the best fit likelihood value for NObs; other parameters fixed at values in grid files and given delta

    virtual double GetMinLikelihood_Delta(bool normalhier=true);//used only within FC code - gets the best fit likelihood value for NObs; other parameters fixed at values in grid files (including theta13).  Fit is along delta.

    //Josh's functions for analytical FC
    double EvaluateOmega(double signal, double background);
    double CalculateRank(int n, double s, double b, double errBg, double k, double errK);
    bool FindBestFit(int n,  double s, double b, double errBg, double k, double errK, double* res);

    string outFileName;

    vector<Extrapolate2D*> Extrap;
    TH1D *NObs;
    TH1D *FracErr_Bkgd;
    TH1D *FracErr_Sig;
    TH1D *Bkgd;
    TH1D *Sig;
    TH1D *NExp;

    TH2D *ErrorMatrix;
    TH2D *InvErrorMatrix;
    TH2D *ExternalErrorMatrix;

    unsigned int nBins;

    TH2D *Chi2_Normal;
    TH2D *Chi2_Inverted;

    string GridFileName_Normal;
    string GridFileName_Inverted;

    vector<TTree*> GridTree_Normal;
    vector<TTree*> GridTree_Inverted;
    vector< vector<TTree*> > GridTree_2_Normal;
    vector< vector<TTree*> > GridTree_2_Inverted;
    TTree *paramtree_Normal,*paramtree_Inverted;

    double grid_background,grid_signal,grid_delta,grid_sinsq2th13,grid_oscparerr;
    double grid_nc,grid_numucc,grid_bnuecc,grid_nutaucc,grid_nue;
    vector<double> grid_bin_oscparerr;
    double grid_n_th12,grid_n_th23,grid_n_dm2_32,grid_n_dm2_21,grid_n_th13;
    double grid_i_th12,grid_i_th23,grid_i_dm2_32,grid_i_dm2_21,grid_i_th13;

    int nPts_Normal;
    int nPts_Inverted;

    int NumExpts;

    string PseudoExpFile;

    double GridNorm;
    double GridScale_Normal;
    double GridScale_Inverted;

    bool IncludeOscParErrs;

    int nSinSq2Th13_Chi2Hist;
    int nDelta_Chi2Hist;
    double SinSq2Th13Low_Chi2Hist,SinSq2Th13High_Chi2Hist;
    double DeltaLow_Chi2Hist,DeltaHigh_Chi2Hist;

    int nDeltaSteps;
    int nSinSq2Th13Steps;
    int nSinSq2Th14Steps;
    int nSinSqTh14Steps;
    int nSinSqTh24Steps;
    double DeltaLow, DeltaHigh;
    double SinSq2Th13Low, SinSq2Th13High;
    double SinSq2Th14Low, SinSq2Th14High;
    double SinSqTh14Low, SinSqTh14High;
    double SinSqTh24Low, SinSqTh24High;

    int nepsetauSteps;
    double epsetauLow, epsetauHigh;

    int FitMethod;

    /** Dung Phan
    *** Sterile package
    **/
public:
    void ReadGridFilesSterileFit();
    void RunMultiBinPseudoExptsSterileFit(bool);
    void RunSterileFit(); // Obsolete (don't use except when you're Dung)
    double GetMinLikelihoodSterileFit();
    double DoGlobalMinSearchSterileFit(vector<double>);
    double StdLikeComparisonForGlobalMinSterileFit(vector<double>);

protected:
    double grid_sinsqth14, grid_sinsqth24, grid_dmsq41, grid_delta13, grid_delta14, grid_delta24, grid_th34;

    ClassDef(NueFit2D,1)
};


#endif
