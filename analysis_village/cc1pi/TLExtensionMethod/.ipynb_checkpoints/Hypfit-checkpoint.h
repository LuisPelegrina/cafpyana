//
// Created by Luis Pelegrina Gutiérrez on 24/7/25.
//
#ifndef Hypfit_h
#define Hypfit_h

#include <iostream>
#include <vector>

#include "PhysdEdx.h"

using namespace std;

class Hypfit {

public:
  Hypfit();

  virtual ~Hypfit();

  std::map< int, PhysdEdx* > map_PhysdEdx;
  double Gaussian(const vector<double> & dEdx, const vector<double> & ResRange, int PID);
  double Likelihood(const vector<double> & dEdx, const vector<double> & ResRange, const vector<double> & pitch, const vector<int> bad_hits, int PID);
  double NormLikelihood(const vector<double> & dEdx, const vector<double> & ResRange, const vector<double> & pitch, const vector<int> bad_hits, bool include_small_likelihood, int PID);
  double NormLikelihood_w_convolution(const vector<double> & dEdx, const vector<double> & ResRange, const vector<double> & pitch, const vector<int> bad_hits, bool include_small_likelihood, int PID, vector<TF1*> tf1_vec);
  double NormLikelihood_w_convolution_each_time(const vector<double> & dEdx, const vector<double> & ResRange, const vector<double> & pitch, const vector<int> bad_hits, bool include_small_likelihood, int PID, int target_plane);

  map<int, vector<TF1*>>  get_conv_function_map(int target_pdg, string mode, int max_rr);
  std::vector<int> get_hits_to_ignore( vector<double> rr, vector<double> dEdx, string mode);
  double GetTLExtensionP(int target_PDG, vector<double> this_rr_vec, vector<double> this_dEdx_vec, vector<double> this_pitch_vec, int best_plane, string cleaning_method);

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
};

#endif
