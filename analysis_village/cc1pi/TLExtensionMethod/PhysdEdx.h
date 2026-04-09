#ifndef PhysdEdx_H
#define PhysdEdx_H

#include <map>
#include "Math/VavilovAccurate.h"
#include "TF1.h"
#include "TF1Convolution.h"
#include "TSpline.h"
#include <iostream>
#include <cmath>
#include <algorithm>


using namespace std;

class TSpline3;

class PhysdEdx :  public TObject{

public:

  PhysdEdx();
  PhysdEdx(int pdg);
  ~PhysdEdx();

  void SetPdgCode(int pdg);
  int GetPdgCode(){ return pdgcode;};

  double Landau_xi(double KE, double pitch);
  double Get_Wmax(double KE);
  double meandEdx(double KE);
  double MPVdEdx(double KE, double pitch);
  double IntegratedEdx(double KE0, double KE1, int n = 10000);
  double RangeFromKE(double KE);
  double RangeFromKESpline(double KE);
  double KEFromRangeSpline(double range);
  double KEAtLength(double KE0, double tracklength);
  double KEtoMomentum(double KE);
  double MomentumtoKE(double momentum);
  void CreateSplineAtKE(int iKE);

  double dEdx_PDF(double KE, double pitch, double dEdx);
  double dEdx_PDF_max(double KE, double pitch, double dEdx);
  double dEdx_Gaus_Sigma(double KE, double pitch);
  double dEdx_PDF_w_convolution(double KE, double pitch, double dEdx, int target_plane, int PDG);
  double dEdx_PDF_max_w_convolution(double KE, double pitch, double dEdx, int target_plane, int PDG);
  double dEdx_PDF_w_convolution_f1(double KE, double rr, double dEdx, double pitch, vector<TF1*> tf1_vec);
  double dEdx_PDF_max_w_convolution_f1(double KE, double rr, double pitch, vector<TF1*> tf1_vec);
  static double dEdx_PDF_function(double *x, double *par);
  static vector<map<int, vector<double>>> pdg_plane_map;

private:

  int pdgcode;
  double mass;
  int charge;

  TSpline3 *sp_KE_range;
  TSpline3 *sp_range_KE;

  map<int, TSpline3*> spmap;

  double densityEffect(double beta, double gamma);

  double betaGamma(double KE);

  void CreateSplines(int np = 1000, double minke = .01, double maxke = 2e5);

  // == Bethe-Bloch parameters, https://indico.fnal.gov/event/14933/contributions/28526/attachments/17961/22583/Final_SIST_Paper.pdf
  const double rho = 1.39; // [g/cm3], density of LAr
  const double K = 0.307075; // [MeV cm2 / mol]
  const double Z = 18.; // atomic number of Ar
  const double A = 39.948; // [g / mol], atomic mass of Ar
  const double I = 197.0e-6; // [MeV], mean excitation energy, JINST 19 (2024) 01, P01009
  const double me = 0.511; // [Mev], mass of electron
  // == Parameters for the density correction
  const double density_C = 5.2146;
  const double density_y0 = 0.2;
  const double density_y1 = 3.0;
  const double density_a = 0.19559;
  const double density_k = 3.0;

public:
  ClassDef(PhysdEdx,1)  //
};




#endif