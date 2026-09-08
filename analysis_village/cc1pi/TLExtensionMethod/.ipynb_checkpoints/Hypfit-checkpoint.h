//
// Created by Luis Pelegrina Gutiérrez on 24/7/25.
//
#ifndef Hypfit_h
#define Hypfit_h

#include <iostream>
#include <vector>

#include "PhysdEdx.h"

using namespace std;

struct Likelihood_point {
  double rr;                  // some variable
  double Likelihood;          // likelihood value
  int num_hits_not_used;      // number of not used hits
};
/*
// Hypfit.h (or wherever shared constants live)
// x-ranges in cm, matching the python analysis: [0,50), [50,100), [100,150), [150,200)
static const int N_X_RANGES = 4;
static const std::vector<std::pair<double, double>> x_ranges = {
  {0.0,   50.0},
  {50.0,  100.0},
  {100.0, 150.0},
  {150.0, 200.0}
};

// Returns the index into x_ranges for a given x position, or -1 if out of range.
// Uses |x|, consistent with the python-side filtering (0 < |x| < bound).
int GetXRangeIndex(double x) {
  double abs_x = fabs(x);
  for (size_t i = 0; i < x_ranges.size(); i++) {
    if (abs_x >= x_ranges[i].first && abs_x < x_ranges[i].second) {
      return static_cast<int>(i);
    }
  }
  return -1; // outside all defined x-ranges
}
*/

class Hypfit {


public:
  Hypfit();

    
  virtual ~Hypfit();

  std::map< int, PhysdEdx* > map_PhysdEdx;
  double Gaussian(const vector<double> & dEdx, const vector<double> & ResRange, int PID);
  double Likelihood(const vector<double> & dEdx, const vector<double> & ResRange, const vector<double> & pitch, const vector<int> bad_hits, int PID);
  double NormLikelihood(const vector<double> & dEdx, const vector<double> & ResRange, const vector<double> & pitch, const vector<int> bad_hits, bool include_small_likelihood, int PID);
  double NormLikelihood_w_convolution(const vector<double> & dEdx, const vector<double> & ResRange, const vector<double> & pitch, const vector<int> bad_hits, bool include_small_likelihood, int PID, vector<TF1*> tf1_vec);
  double NormLikelihood_w_convolution_each_time(const vector<double> & dEdx, const vector<double> & ResRange,
                                               const vector<double> & pitch, const vector<int> bad_hits,
                                               bool include_small_likelihood, int PID, int target_plane,
                                               const string &mode = "");
  vector<Likelihood_point> NormLikelihoodPairVector_w_convolution(const vector<double> & dEdx, const vector<double> & ResRange, const vector<double> & pitch, const vector<int> bad_hits, bool include_small_likelihood, int PID, vector<TF1*> tf1_vec);


  map<int, vector<TF1*>> get_conv_function_map(int target_pdg, string mode, int max_rr, bool use_data_map = false);
  std::vector<int> get_hits_to_ignore( vector<double> rr, vector<double> dEdx, string mode);
  double GetTLExtensionP(int target_PDG, vector<double> this_rr_vec, vector<double> this_dEdx_vec, vector<double> this_pitch_vec, int best_plane, string cleaning_method);
    
  double robust_max_x(TF1 *f, double xmin, double xmax, int n_scan = 2000);
  map<int, vector<TF1*>> get_langau_map(int target_pdg, int max_rr, bool use_data_map = false);
 static double langau_cpp(double *x, double *p);

/*
// == Header additions (Hypfit.h)
    std::vector<double> get_pitch_bin_edges() const;
    int get_pitch_bin_index(double pitch, const std::vector<double> &edges) const;
    std::map<int, std::vector<std::vector<TF1*>>> get_conv_function_map_w_pitch(int target_pdg, std::string mode, int max_rr);
    double NormLikelihood_w_convolution_w_pitch(const std::vector<double> & dEdx, const std::vector<double> & ResRange,
        const std::vector<double> & pitch, const std::vector<int> bad_hits, bool include_small_likelihood,
        int PID, std::vector<std::vector<TF1*>> tf1_map);

// Hypfit.h
double NormLikelihood_w_convolution_w_x_ranges(
    const vector<double> & dEdx,
    const vector<double> & ResRange,
    const vector<double> & pitch,
    const vector<double> & Xpos,
    const vector<int> bad_hits,
    bool include_small_likelihood,
    int PID,
    const vector<vector<TF1*>> & tf1_vec_by_xrange);
// Hypfit.h
vector<map<int, vector<TF1*>>> get_conv_function_map_w_x_ranges(
    int target_pdg, string mode, int max_rr);
    */

private:

  // == Tunable parameters
  int N_skip = 3;
  double max_additional_res_length_pion = 250.; // == [cm]
  double max_additional_res_length_proton = 120.; // == [cm]
  double res_length_step_pion = 1.0; // == [cm]
  double res_length_step_proton = 0.2; // == [cm]
  double dEdx_truncate_upper_pion = 5.;
  double dEdx_truncate_bellow_pion = 0.5;
  double dEdx_truncate_upper_proton = 20.;
  double dEdx_truncate_bellow_proton = 0.2;
  double pion_mass = 0.13957039;
  double muon_mass = 0.1056583755;
  double proton_mass = 0.93827208943;

  std::map<int, std::vector<TF1*>> conv_tf1_map_muon;
  std::map<int, std::vector<TF1*>> conv_tf1_map_pion;
  std::map<int, std::vector<TF1*>> conv_tf1_map_proton;

};

#endif
