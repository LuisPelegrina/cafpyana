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
  void dEdx_PDF_and_max_w_convolution(double KE, double pitch, double dEdx,
                                               int target_plane, int PDG, const string &mode,
                                               double &out_likelihood, double &out_likelihood_max);
      
  double dEdx_PDF_w_convolution_f1(double KE, double rr, double dEdx, double pitch, vector<TF1*> tf1_vec);
  double dEdx_PDF_max_w_convolution_f1(double KE, double rr, double pitch, vector<TF1*> tf1_vec);
  static double dEdx_PDF_function(double *x, double *par);
  static vector<map<int, vector<double>>> pdg_plane_map;
  static vector<map<int, vector<double>>> pdg_plane_shift_map;
  static vector<map<int, vector<double>>> pdg_plane_map_data;
  static vector<map<int, vector<double>>> pdg_plane_shift_map_data;
     static vector<map<int, vector<double>>> pdg_plane_mpv_map;
    static vector<map<int, vector<double>>> pdg_plane_gsigma_map;
    static vector<map<int, vector<double>>> pdg_plane_width_map;
    static vector<map<int, vector<double>>> pdg_plane_mpv_map_data;
    static vector<map<int, vector<double>>> pdg_plane_gsigma_map_data;
    static vector<map<int, vector<double>>> pdg_plane_width_map_data;
// Low Residual Range maps (3-entry vectors: {mpv, gsigma, width})
  static vector<map<int, vector<double>>> pdg_plane_low_rr_3p5;
  static vector<map<int, vector<double>>> pdg_plane_low_rr_3p5_data;
  static vector<map<int, vector<double>>> pdg_plane_low_rr_2p5;
  static vector<map<int, vector<double>>> pdg_plane_low_rr_2p5_data;

/*
// PhysdEdx.h

// Flat replacement for vector<vector<map<int, vector<double>>>>.
// Key: (x_range_idx, plane_idx, pdg) -> {p0, p1, p2}
using XRangePlanePdgKey = std::tuple<int, int, int>;

static std::map<XRangePlanePdgKey, std::array<double, 3>> pdg_plane_map_w_x_ranges;
static std::map<XRangePlanePdgKey, std::array<double, 3>> pdg_plane_shift_map_w_x_ranges;

// Small accessor so call sites don't need to spell out std::get<> everywhere.
// Returns {-1,-1,-1} if the (x_range, plane, pdg) combo hasn't been filled in yet.
static const std::array<double, 3>& GetPlaneMapEntry(
    const std::map<XRangePlanePdgKey, std::array<double, 3>> &m,
    int x_range_idx, int plane_idx, int pdg) {
  static const std::array<double, 3> not_found = {-1, -1, -1};
  auto it = m.find({x_range_idx, plane_idx, pdg});
  return it != m.end() ? it->second : not_found;
}
*/
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