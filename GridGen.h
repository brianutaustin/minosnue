#ifndef GridGen_h
#define GridGen_h

#include "TROOT.h"
#include "NueAna/MultiBinAna/Extrapolate2D.h"
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "OscProb/OscCalc.h"
#include "NueAna/Extrapolation/Background.h"
#include "TRandom.h"

using namespace std;

class GridGen{
  public:
    GridGen();
    virtual ~GridGen();
    
    void SetOutputFile(string s="GridOut.root") { outFileName = s; };
    void AddExtrap(Extrapolate2D* E);
    
    void SetNDeltaSteps(int n=20) { nDeltaSteps = n; };
    void SetNSinSq2Th13Steps(int n=3000) { nSinSq2Th13Steps = n; };
    void SetDeltaRange(double l=0, double h=2);//in units of pi
    void SetSinSq2Th13Range(double l=0.05, double h = 0.35);
    void SetNormalHierarchy(Bool_t normal = true) { NormalHier = normal; };
    
    void SetTheta12(double val = 0.59365, double errup = 0.041, double errdn = 0.041);
    void SetTheta23(double val = 0.78540, double errup = 0.122, double errdn = 0.122);
    void SetAbsValDeltaMSq23(double val = 2.32e-3, double errup = 0.11e-3, double errdn = 0.11e-3);
    void SetDeltaMSq12(double val = 8.0e-5, double errup = 0.6e-5, double errdn = 0.6e-5);

    void SetTheta13(double val = 0.1594, double errup = 0.0110, double errdn = 0.0110);
    
    void FreezeTheta23(Bool_t flag) { FrozenTheta23 = flag;}

    void SetNExperiments(int n=2000) { NumExpts = n; };
    
    void Run();
    void RunWithOscParErrs(string s = "OscParErrDistributions.root");
    virtual void RunMultiBinOscParErrs(string s = "OscParErrDistributions.root");
    virtual void RunMultiBin_VaryTheta13(string s = "OscParErrDistributions.root");
    
    // Dung's code
    virtual void SetNStepSinSqTh14SterileFit(int n = 100);
    virtual void SetNStepSinSqTh24SterileFit(int n = 100);
    virtual void SetSinSqTh14SterileFit(double);
    virtual void SetSinSqTh24SterileFit(double);
    virtual void RunMultiBinOscParErrsSterileFit(string s = "OscParErrDistributions.root", double dm41 = 0.1);
    
  protected:
    
    double AsymGaus(Double_t sm, Double_t sp);
    double DrawTheta23(double dtheta);
    
    vector<Extrapolate2D*> Extrap;
    
    Int_t nDeltaSteps;
    Int_t nSinSq2Th13Steps;
    
    Double_t DeltaLow;
    Double_t DeltaHigh;
    Double_t SinSq2Th13Low;
    Double_t SinSq2Th13High;

    Bool_t FrozenTheta23;
    
    Bool_t NormalHier;
    
    string outFileName;
    
    //these values will remain fixed
    Double_t Theta12;
    Double_t Theta23;
    Double_t DeltaMSq23;
    Double_t DeltaMSq12;
    Double_t dTheta12_up,dTheta12_dn;
    Double_t dTheta23_up,dTheta23_dn;
    Double_t dDeltaMSq23_up,dDeltaMSq23_dn;
    Double_t dDeltaMSq12_up,dDeltaMSq12_dn;
    
    Double_t Theta13;
    Double_t dTheta13_up,dTheta13_dn;
    
    Int_t NumExpts;

private:
    // Dung's code
    Int_t NStepSinSqTh14;
    Int_t NStepSinSqTh24;
    double SinSqTh14GridValue;
    double SinSqTh24GridValue;
};

#endif
