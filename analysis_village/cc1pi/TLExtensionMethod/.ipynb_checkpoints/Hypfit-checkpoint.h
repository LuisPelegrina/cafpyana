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


class Hypfit {


public:
  Hypfit();

    
  virtual ~Hypfit();

  std::map< int, PhysdEdx* > map_PhysdEdx;
  double Gaussian(const vector<double> & dEdx, const vector<double> & ResRange, int PID);
  double Likelihood(const vector<double> & dEdx, const vector<double> & ResRange, const vector<double> & pitch, const vector<int> bad_hits, int PID);
  double NormLikelihood(const vector<double> & dEdx, const vector<double> & ResRange, const vector<double> & pitch, const vector<int> bad_hits, bool include_small_likelihood, int PID);
  double NormLikelihood_w_tf1_vec(const vector<double> & dEdx, const vector<double> & ResRange, const vector<double> & pitch, const vector<int> bad_hits, bool include_small_likelihood, int PID, vector<TF1*> tf1_vec);
  vector<Likelihood_point> NormLikelihoodPairVector_w_tf1_vec(const vector<double> & dEdx, const vector<double> & ResRange, const vector<double> & pitch, const vector<int> bad_hits, bool include_small_likelihood, int PID, vector<TF1*> tf1_vec);


  map<int, vector<TF1*>> get_conv_function_map(int target_pdg, string mode, int max_rr, bool use_data_map = false);
  std::vector<int> get_hits_to_ignore( vector<double> rr, vector<double> dEdx, string mode);
 double GetTLExtensionP(int target_PDG, std::vector<double> this_rr_vec, std::vector<double> this_dEdx_vec, std::vector<double> this_pitch_vec, int best_plane, std::string cleaning_method, std::map<int, std::vector<TF1*>> tf1_map);
    
  double robust_max_x(TF1 *f, double xmin, double xmax, int n_scan = 2000);
  map<int, vector<TF1*>> get_langau_map(int target_pdg, int max_rr, bool use_data_map = false);
 static double langau_cpp(double *x, double *p);


private:

  // == Tunable parameters
  int N_skip = 3;
  long unsigned int N_minimum_hits = 10;
  double max_additional_res_length_pion = 250.; // == [cm]
  double max_additional_res_length_proton = 120.; // == [cm]
  double res_length_step_pion = 1.0; // == [cm]
  double res_length_step_proton = 0.2; // == [cm]
  double dEdx_truncate_upper_harsh = 5.;
  double dEdx_truncate_bellow_harsh = 0.5;
  double dEdx_truncate_upper = 10.;
  double dEdx_truncate_bellow = 0.5;
  double pion_mass = 0.13957039;
  double muon_mass = 0.1056583755;
  double proton_mass = 0.93827208943;
};

#endif
