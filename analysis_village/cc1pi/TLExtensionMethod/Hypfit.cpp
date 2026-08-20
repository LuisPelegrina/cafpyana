//
// Created by Luis Pelegrina Gutiérrez on 24/7/25.
//

#include "Hypfit.h"

Hypfit::Hypfit(){
  map_PhysdEdx[13] = new PhysdEdx(13);
  map_PhysdEdx[211] = new PhysdEdx(211);
  map_PhysdEdx[321] = new PhysdEdx(321);
  map_PhysdEdx[2212] = new PhysdEdx(2212);

  std::cout << "Hypfit: Pre-calculating convolution maps..." << std::endl;
  conv_tf1_map_muon   = get_conv_function_map(13,   "none", 2000);
  conv_tf1_map_pion   = get_conv_function_map(211,  "none", 2000);
  conv_tf1_map_proton = get_conv_function_map(2212, "none", 2000);
}

Hypfit::~Hypfit(){
  for(auto const& [pdg, phys] : map_PhysdEdx) {
    delete phys;
  }
}
double Hypfit::Gaussian(const vector<double> & dEdx, const vector<double> & ResRange, int PID){

  // == PID input : mass hypothesis, valid only for muons, charged pions, and protons
  if(!(PID == 13 || PID == 2212 || PID == 211)){
    return -9999.;
  }
  // == Tunable parameters
  double min_additional_res_length = 0.;
  double max_additional_res_length = max_additional_res_length_pion;
  double res_length_step = res_length_step_pion;
  double dEdx_truncate_upper = dEdx_truncate_upper_pion;
  double dEdx_truncate_bellow = dEdx_truncate_bellow_pion;
  if(PID == 2212){
    max_additional_res_length = max_additional_res_length_proton;
    dEdx_truncate_upper = dEdx_truncate_upper_proton;
    dEdx_truncate_bellow = dEdx_truncate_bellow_proton;
    res_length_step = res_length_step_proton;
  }
  int res_length_trial = (max_additional_res_length - min_additional_res_length) / res_length_step;

  // == Initialize
  double best_additional_res_length = -0.1;
  double best_chi2 = 99999.;

  int this_N_calo = dEdx.size();
  
  if(this_N_calo <= 10){
    return -8888.; // == Too small number of hits
  }
  int i_bestfit = -1;
  int this_N_hits = this_N_calo;

  // == Fit
  for(int i = 0; i < res_length_trial; i++){
    double this_additional_res_length = min_additional_res_length + (i + 0.) * res_length_step;
    double this_chi2 = 0.;
    for(int j = N_skip; j < this_N_hits - N_skip; j++){ // == Do not use first and last N_skip hits
      double this_res_length = ResRange.at(j) + this_additional_res_length;
      double this_KE = map_PhysdEdx[PID]->KEFromRangeSpline(this_res_length);
      double dEdx_theory = map_PhysdEdx[PID]->meandEdx(this_KE);
      double dEdx_measured = dEdx.at(j);
      if(dEdx_measured < dEdx_truncate_bellow || dEdx_measured > dEdx_truncate_upper) continue; // == Truncate
      // == Gaussian approx.
      //double dEdx_theory_err = dEdx_theory * 0.02;
      this_chi2 += pow(dEdx_measured - dEdx_theory, 2);
    }
    this_chi2 = this_chi2 / (this_N_hits + 0.); // == chi2 / n.d.f
    if(this_chi2 < best_chi2){
      best_chi2 = this_chi2;
      best_additional_res_length = this_additional_res_length;
      i_bestfit = i;
    }
  }


  double original_res_length = ResRange.at(this_N_calo - 1); // == [cm]
  double best_total_res_length = best_additional_res_length + original_res_length;

  // == Define fitting failed cases
  if(i_bestfit == res_length_trial - 1){ // == Fit failed : no mimumum
    return -7777.;
  }
  else if(best_chi2 > 99990.){ // == Fit failed : best_chi2 > 99990."
    return -6666.;
  }
  else if(best_chi2 < 1.0e-11){ // == Fit failed : best_chi2 < 1.0e-11
    return -5555.;
  }

  double best_KE = map_PhysdEdx[PID] -> KEFromRangeSpline(best_total_res_length);
  return best_KE;
}

double Hypfit::Likelihood(const vector<double> & dEdx, const vector<double> & ResRange, const vector<double> & pitch, const vector<int> bad_hits, int PID){

  // == PID input : mass hypothesis, valid only for muons, charged pions, and protons
  if(!(PID == 13 || PID == 2212 || PID == 211)){
    return -9999.;
  }
  // == Tunable parameters
  double min_additional_res_length = 0.;
  double max_additional_res_length = max_additional_res_length_pion;
  double res_length_step = res_length_step_pion;
  double dEdx_truncate_upper = dEdx_truncate_upper_pion;
  double dEdx_truncate_bellow = dEdx_truncate_bellow_pion;
  if(PID == 2212){
    max_additional_res_length = max_additional_res_length_proton;
    dEdx_truncate_upper = dEdx_truncate_upper_proton;
    dEdx_truncate_bellow = dEdx_truncate_bellow_proton;
    res_length_step = res_length_step_proton;
  }
  int res_length_trial = (max_additional_res_length - min_additional_res_length) / res_length_step;

  // == Initialize
  double best_additional_res_length = -0.1;
  double best_m2lnL = 99999.;

  int this_N_calo = dEdx.size();
  if(this_N_calo-1 <= 0){
    return -2222.; // ==No hits
  }
  if(this_N_calo - bad_hits.size() <= 4) {
    return -8888.; // == Too small number of hits
  }
  int i_bestfit = -1;
  int this_N_hits = this_N_calo;

  // == Fit
  for(int i = 0; i < res_length_trial; i++){
    double this_additional_res_length = min_additional_res_length + (i + 0.) * res_length_step;
    double this_m2lnL = 0.;
    for(int j = 0; j < this_N_hits; j++){ // == Do not use first and last N_skip this
      if (std::find(bad_hits.begin(), bad_hits.end(), j) != bad_hits.end()) continue;

      double this_res_length = ResRange.at(j) + this_additional_res_length;
      double this_KE = map_PhysdEdx[PID]->KEFromRangeSpline(this_res_length);
      double dEdx_measured = dEdx.at(j);
      if(dEdx_measured < dEdx_truncate_bellow || dEdx_measured > dEdx_truncate_upper) continue; // == Truncate
      // == Likelihood
      double this_pitch = pitch.at(j);
      double this_likelihood = map_PhysdEdx[PID] -> dEdx_PDF(this_KE, this_pitch, dEdx_measured);
      if(this_likelihood > 1e-6) this_m2lnL += (-2.0) * log(this_likelihood);
    }
    if(this_m2lnL < best_m2lnL){
      best_m2lnL = this_m2lnL;
      best_additional_res_length = this_additional_res_length;
      i_bestfit = i;
    }
  }

  // == Result
  double original_res_length = ResRange.at(this_N_calo - 1); // == [cm]
  double best_total_res_length = best_additional_res_length + original_res_length; // == [cm]
  // == Define fitting failed cases
  if(i_bestfit == res_length_trial - 1){ // == Fit failed : no mimumum
    return -7777.;
  }
  else if(best_m2lnL > 99990.){ // == Fit failed : best likelihood > 99990.
    return -6666.;
  }
  else if(fabs(best_m2lnL) < 1.0e-11){ // == Fit failed : |best_m2lnL| < 1.0e-11, no valid likelihood value
    return -5555.;
  }

  // == Return
  double best_KE = map_PhysdEdx[PID] -> KEFromRangeSpline(best_total_res_length);
  return best_KE;
}

double Hypfit::NormLikelihood(const vector<double> & dEdx, const vector<double> & ResRange, const vector<double> & pitch, const vector<int> bad_hits, bool include_small_likelihood, int PID){

  // == PID input : mass hypothesis, valid only for muons, charged pions, and protons
  if(!(PID == 13 || PID == 2212 || PID == 211)){
    return -9999.;
  }
  // == Tunable parameters
  double min_additional_res_length = 0.;
  double max_additional_res_length = max_additional_res_length_pion;
  double res_length_step = res_length_step_pion;
  if(PID == 2212){
    max_additional_res_length = max_additional_res_length_proton;
    res_length_step = res_length_step_proton;
  }
  int res_length_trial = (max_additional_res_length - min_additional_res_length) / res_length_step;

  // == Initialize
  double best_additional_res_length = -0.1;
  double best_m2lnL = 99999.;

  int this_N_calo = dEdx.size();
  if(this_N_calo-1 <= 0){
    return -2222.; // ==No hits
  }
  if(this_N_calo - bad_hits.size() <= 4) {
    return -8888.; // == Too small number of hits
  }
  int i_bestfit = -1;
  int this_N_hits = this_N_calo;

  // == Fit
  for(int i = 0; i < res_length_trial; i++) {

    double this_additional_res_length = min_additional_res_length + (i + 0.) * res_length_step;
    double this_m2lnL = 0.;
    int num_invalid_hits = 0;
    int num_valid_hits = 0;

    for (int j = 0; j < this_N_hits; j++) { // == Do not use first and last N_skip this

      if (std::find(bad_hits.begin(), bad_hits.end(), j) != bad_hits.end()) continue;
      // == Do not use first and last N_skip this
      double this_res_length = ResRange.at(j) + this_additional_res_length;
      double this_KE = map_PhysdEdx[PID]->KEFromRangeSpline(this_res_length);
      double dEdx_measured = dEdx.at(j);
      // == Likelihood
      double this_pitch = pitch.at(j);
      double this_likelihood = map_PhysdEdx[PID]->dEdx_PDF(this_KE, this_pitch, dEdx_measured);
      double this_likelihood_max = map_PhysdEdx[PID]->dEdx_PDF_max(this_KE, this_pitch, dEdx_measured);

      if (this_likelihood > 1e-6) {
        num_valid_hits++;
        this_m2lnL += (2.0) * (log(this_likelihood_max) - log(this_likelihood));
      } else {
        if (include_small_likelihood) this_m2lnL += (2.0) * (log(this_likelihood_max) - log(1e-6));
        num_invalid_hits++;
      }
    }

    if (num_valid_hits > 0) {
      if (include_small_likelihood) {
        this_m2lnL = this_m2lnL / (num_valid_hits + num_invalid_hits);
      } else {
        this_m2lnL = this_m2lnL / num_valid_hits;
      }
    }


    if(this_m2lnL < best_m2lnL){
      best_m2lnL = this_m2lnL;
      best_additional_res_length = this_additional_res_length;
      i_bestfit = i;
    }
  }
  // == Result
  double original_res_length = ResRange.at(this_N_calo - 1); // == [cm]
  double best_total_res_length = best_additional_res_length + original_res_length; // == [cm]
  // == Define fitting failed cases
  if(i_bestfit == res_length_trial - 1){ // == Fit failed : no mimumum
    return -7777.;
  }
  else if(best_m2lnL > 99990.){ // == Fit failed : best likelihood > 99990.
    return -6666.;
  }
  else if(fabs(best_m2lnL) < 1.0e-11){ // == Fit failed : |best_m2lnL| < 1.0e-11, no valid likelihood value
    return -5555.;
  }

  // == Return
  double best_KE = map_PhysdEdx[PID] -> KEFromRangeSpline(best_total_res_length);
  return best_KE;

}

double Hypfit::NormLikelihood_w_convolution(const vector<double> & dEdx, const vector<double> & ResRange, const vector<double> & pitch, const vector<int> bad_hits, bool include_small_likelihood, int PID, vector<TF1*> tf1_vec) {

  // == PID input : mass hypothesis, valid only for muons, charged pions, and protons
  if(!(PID == 13 || PID == 2212 || PID == 211)){
    return -9999.;
  }
  // == Tunable parameters
  double min_additional_res_length = 0.;
  double max_additional_res_length = max_additional_res_length_pion;
  double res_length_step = res_length_step_pion;
  if(PID == 2212){
    max_additional_res_length = max_additional_res_length_proton;
    res_length_step = res_length_step_proton;
  }
  int res_length_trial = (max_additional_res_length - min_additional_res_length) / res_length_step;

  // == Initialize
  double best_additional_res_length = -0.1;
  double best_m2lnL = 99999.;


  int this_N_calo = dEdx.size();
  if(this_N_calo-1 <= 0){
    return -2222.; // ==No hits
  }
  if(this_N_calo - bad_hits.size() <= 4) {
    return -8888.; // == Too small number of hits
  }
  int i_bestfit = -1;
  int this_N_hits = this_N_calo;

  // == Fit

  for(int i = 0; i < res_length_trial; i++) {

    double this_additional_res_length = min_additional_res_length + (i + 0.) * res_length_step;
    double this_m2lnL = 0.;
    int num_invalid_hits = 0;
    int num_valid_hits = 0;

    for (int j = 0; j < this_N_hits; j++) { // == Do not use first and last N_skip this

      if (std::find(bad_hits.begin(), bad_hits.end(), j) != bad_hits.end()) continue;
      // == Do not use first and last N_skip this
      double this_res_length = ResRange.at(j) + this_additional_res_length;
      double this_KE = map_PhysdEdx[PID]->KEFromRangeSpline(this_res_length);
      double dEdx_measured = dEdx.at(j);
      // == Likelihood
      double this_pitch = pitch.at(j);
      double this_likelihood = map_PhysdEdx[PID]->dEdx_PDF_w_convolution_f1(this_KE, this_res_length, dEdx_measured, this_pitch, tf1_vec);
      double this_likelihood_max = map_PhysdEdx[PID]->dEdx_PDF_max_w_convolution_f1(this_KE, this_res_length, this_pitch, tf1_vec);
      //double this_likelihood = map_PhysdEdx[PID]->dEdx_PDF_w_convolution(this_KE, this_pitch, dEdx_measured);
      //double this_likelihood_max = map_PhysdEdx[PID]->dEdx_PDF_max_w_convolution(this_KE, this_pitch, dEdx_measured);

      if (this_likelihood > 1e-6) {
        num_valid_hits++;
        this_m2lnL += (2.0) * (log(this_likelihood_max) - log(this_likelihood));
      } else {
        if (include_small_likelihood) this_m2lnL += (2.0) * (log(this_likelihood_max) - log(1e-6));
        num_invalid_hits++;
      }
    }

    if (num_valid_hits > 0) {
      if (include_small_likelihood) {
        this_m2lnL = this_m2lnL / (num_valid_hits + num_invalid_hits);
      } else {
        this_m2lnL = this_m2lnL / num_valid_hits;
      }
    }

    if(this_m2lnL < best_m2lnL){
      best_m2lnL = this_m2lnL;
      best_additional_res_length = this_additional_res_length;
      i_bestfit = i;
    }

  }

  // == Result
  double original_res_length = ResRange.at(this_N_calo - 1); // == [cm]

  double best_total_res_length = best_additional_res_length + original_res_length; // == [cm]
  // == Define fitting failed cases
  if(i_bestfit == res_length_trial - 1){ // == Fit failed : no mimumum
    return -7777.;
  }
  else if(best_m2lnL > 99990.){ // == Fit failed : best likelihood > 99990.
    return -6666.;
  }
  else if(fabs(best_m2lnL) < 1.0e-11){ // == Fit failed : |best_m2lnL| < 1.0e-11, no valid likelihood value
    return -5555.;
  } else if(std::isnan(best_total_res_length)) {
    return -3333.;
  }


  // == Return
  double best_KE = map_PhysdEdx[PID] -> KEFromRangeSpline(best_total_res_length);
  return best_KE;

}

vector<Likelihood_point> Hypfit::NormLikelihoodPairVector_w_convolution(const vector<double> & dEdx, const vector<double> & ResRange, const vector<double> & pitch, const vector<int> bad_hits, bool include_small_likelihood, int PID, vector<TF1*> tf1_vec) {

  vector<Likelihood_point> likelihood_point_vector;

  // == PID input : mass hypothesis, valid only for muons, charged pions, and protons
  if(!(PID == 13 || PID == 2212 || PID == 211)){
    likelihood_point_vector = {{-9999, -9999, -9999}};
    return likelihood_point_vector;
  }

  // == Tunable parameters
  double min_additional_res_length = 0.;
  double max_additional_res_length = max_additional_res_length_pion;
  double res_length_step = res_length_step_pion;
  if(PID == 2212){
    max_additional_res_length = max_additional_res_length_proton;
    res_length_step = res_length_step_proton;
  }
  int res_length_trial = (max_additional_res_length - min_additional_res_length) / res_length_step;

  int this_N_calo = dEdx.size();
  if(this_N_calo - 1 <= 0){
    likelihood_point_vector = {{-2222, -2222, -2222}};
    return likelihood_point_vector; // == No hits
  }
  if(this_N_calo - bad_hits.size() <= 4) {
    likelihood_point_vector = {{-8888, -8888, -8888}};
    return likelihood_point_vector; // == Too small number of hits
  }
  int this_N_hits = this_N_calo;

  // == Fit
  for(int i = 0; i < res_length_trial; i++) {

    double this_additional_res_length = min_additional_res_length + (i + 0.) * res_length_step;
    double this_m2lnL = 0.;
    int num_invalid_hits = 0;
    int num_valid_hits = 0;

    for (int j = 0; j < this_N_hits; j++) {

      if (std::find(bad_hits.begin(), bad_hits.end(), j) != bad_hits.end()) continue;

      double this_res_length = ResRange.at(j) + this_additional_res_length;
      double this_KE = map_PhysdEdx[PID]->KEFromRangeSpline(this_res_length);
      double dEdx_measured = dEdx.at(j);
      double this_pitch = pitch.at(j);

      // == Corrected PDF calls matching NormLikelihood_w_convolution
      double this_likelihood = map_PhysdEdx[PID]->dEdx_PDF_w_convolution_f1(this_KE, this_res_length, dEdx_measured, this_pitch, tf1_vec);
      double this_likelihood_max = map_PhysdEdx[PID]->dEdx_PDF_max_w_convolution_f1(this_KE, this_res_length, this_pitch, tf1_vec);

      if (this_likelihood > 1e-6) {
        num_valid_hits++;
        this_m2lnL += (2.0) * (log(this_likelihood_max) - log(this_likelihood));
      } else {
        if (include_small_likelihood) this_m2lnL += (2.0) * (log(this_likelihood_max) - log(1e-6));
        num_invalid_hits++;
      }
    }

    if (num_valid_hits > 0) {
      if (include_small_likelihood) {
        this_m2lnL = this_m2lnL / (num_valid_hits + num_invalid_hits);
      } else {
        this_m2lnL = this_m2lnL / num_valid_hits;
      }
    }

    Likelihood_point lp = {this_additional_res_length + ResRange.at(this_N_calo - 1), this_m2lnL, num_invalid_hits};
    likelihood_point_vector.push_back(lp);
  }

  // == Return
  return likelihood_point_vector;
}




double Hypfit::NormLikelihood_w_convolution_each_time(const vector<double> & dEdx, const vector<double> & ResRange, const vector<double> & pitch, const vector<int> bad_hits, bool include_small_likelihood, int PID, int target_plane) {

  // == PID input : mass hypothesis, valid only for muons, charged pions, and protons
  if(!(PID == 13 || PID == 2212 || PID == 211)){
    return -9999.;
  }
  // == Tunable parameters
  double min_additional_res_length = 0.;
  double max_additional_res_length = max_additional_res_length_pion;
  double res_length_step = res_length_step_pion;
  if(PID == 2212){
    max_additional_res_length = max_additional_res_length_proton;
    res_length_step = res_length_step_proton;
  }
  int res_length_trial = (max_additional_res_length - min_additional_res_length) / res_length_step;

  // == Initialize
  double best_additional_res_length = -0.1;
  double best_m2lnL = 99999.;


  int this_N_calo = dEdx.size();
  if(this_N_calo-1 <= 0){
    return -2222.; // ==No hits
  }
  if(this_N_calo - bad_hits.size() <= 4) {
    return -8888.; // == Too small number of hits
  }
  int i_bestfit = -1;
  int this_N_hits = this_N_calo;

  // == Fit

  for(int i = 0; i < res_length_trial; i++) {

    double this_additional_res_length = min_additional_res_length + (i + 0.) * res_length_step;
    double this_m2lnL = 0.;
    int num_invalid_hits = 0;
    int num_valid_hits = 0;

    for (int j = 0; j < this_N_hits; j++) { // == Do not use first and last N_skip this

      if (std::find(bad_hits.begin(), bad_hits.end(), j) != bad_hits.end()) continue;
      // == Do not use first and last N_skip this
      double this_res_length = ResRange.at(j) + this_additional_res_length;
      double this_KE = map_PhysdEdx[PID]->KEFromRangeSpline(this_res_length);
      double dEdx_measured = dEdx.at(j);
      // == Likelihood
      double this_pitch = pitch.at(j);
      double this_likelihood = map_PhysdEdx[PID]->dEdx_PDF_w_convolution(this_KE, this_pitch, dEdx_measured, target_plane, PID);
      double this_likelihood_max = map_PhysdEdx[PID]->dEdx_PDF_max_w_convolution(this_KE, this_pitch, dEdx_measured, target_plane, PID);

      if (this_likelihood > 1e-6) {
        num_valid_hits++;
        this_m2lnL += (2.0) * (log(this_likelihood_max) - log(this_likelihood));
      } else {
        if (include_small_likelihood) this_m2lnL += (2.0) * (log(this_likelihood_max) - log(1e-6));
        num_invalid_hits++;
      }
    }

    if (num_valid_hits > 0) {
      if (include_small_likelihood) {
        this_m2lnL = this_m2lnL / (num_valid_hits + num_invalid_hits);
      } else {
        this_m2lnL = this_m2lnL / num_valid_hits;
      }
    }

    if(this_m2lnL < best_m2lnL){
      best_m2lnL = this_m2lnL;
      best_additional_res_length = this_additional_res_length;
      i_bestfit = i;
    }

  }

  // == Result
  double original_res_length = ResRange.at(this_N_calo - 1); // == [cm]

  double best_total_res_length = best_additional_res_length + original_res_length; // == [cm]
  // == Define fitting failed cases
  if(i_bestfit == res_length_trial - 1){ // == Fit failed : no mimumum
    return -7777.;
  }
  else if(best_m2lnL > 99990.){ // == Fit failed : best likelihood > 99990.
    return -6666.;
  }
  else if(fabs(best_m2lnL) < 1.0e-11){ // == Fit failed : |best_m2lnL| < 1.0e-11, no valid likelihood value
    return -5555.;
  } else if(std::isnan(best_total_res_length)) {
    return -3333.;
  }

  // == Return
  double best_KE = map_PhysdEdx[PID] -> KEFromRangeSpline(best_total_res_length);
  return best_KE;

}


std::vector<int> Hypfit::get_hits_to_ignore( vector<double> rr, vector<double> dEdx, string mode){


  std::vector<int> hits_id_to_ignore;

  if(rr.size() <= 6) {
    for(size_t i_rr = 0; i_rr < rr.size(); i_rr++) {
      hits_id_to_ignore.push_back(i_rr); // == Truncate
    }
    return hits_id_to_ignore;
  }

  double dEdx_truncate_upper = 20.;
  double dEdx_truncate_bellow = 0.5;
  if (mode.find("truncate_harsh") != std::string::npos) {
    dEdx_truncate_upper = 5.;
    dEdx_truncate_bellow = 0.5;
  }
  for(size_t i_rr = 0; i_rr < rr.size(); i_rr++) {
    if(dEdx.at(i_rr) < dEdx_truncate_bellow || dEdx.at(i_rr) > dEdx_truncate_upper)  {
      hits_id_to_ignore.push_back(i_rr); // == Truncate
      continue;
    }

    if (mode.find("include_1cm") == std::string::npos)  {
      if(rr.at(i_rr) < 1) {
        hits_id_to_ignore.push_back(i_rr);
        continue;
      }
    }

    if (mode.find("ignore_3cm") != std::string::npos) {
      if(rr.at(i_rr) < 3) {
        hits_id_to_ignore.push_back(i_rr);
        continue;
      }
    }

    if((int)i_rr < N_skip)  {
      hits_id_to_ignore.push_back(i_rr);
      continue;
    }

    if(i_rr >= rr.size() - N_skip)  {
      hits_id_to_ignore.push_back(i_rr);
      continue;
    }
  }
  if (mode.find("skip_only") != std::string::npos) return hits_id_to_ignore;

  TGraph *gr_dEdx_corrected = new TGraph();
  for(size_t i_rr = 0; i_rr < rr.size(); i_rr++) {
    if(std::find(hits_id_to_ignore.begin(), hits_id_to_ignore.end(), i_rr) != hits_id_to_ignore.end()) continue;
    gr_dEdx_corrected->SetPoint(gr_dEdx_corrected->GetN(), rr.at(i_rr), dEdx.at(i_rr));
  }

// Suppose you already have TGraph* gr
  TF1* fexp = new TF1("fexp", "[0] + [1]*exp(-[2]*x)", 0, 500); // adjust range
  fexp->SetParNames("A", "B");

  // Set initial guesses
  fexp->SetParameters(1.0, 0.1);

  //Enforce positivity
  fexp->SetParLimits(0, 0, 1e6); // A > 0
  fexp->SetParLimits(1, 0, 1e6); // A > 0
  fexp->SetParLimits(2, 0, 1e3); // B > 0

  gr_dEdx_corrected->Fit(fexp, "0Q"); // "R" restricts to given function range
  double sum =0;
  int cont = 0;
  for(size_t i_rr = 5; i_rr < rr.size()-5; i_rr++) {
    if(std::find(hits_id_to_ignore.begin(), hits_id_to_ignore.end(), i_rr) == hits_id_to_ignore.end()) {
      sum+= dEdx.at(i_rr)  - fexp->Eval(rr.at(i_rr));
      cont++;
    }
  }
  double mean = sum / cont;


// --- Standard deviation (population version)
  double sq_sum = 0.0;
  for(size_t i_rr =  5; i_rr < rr.size()-5; i_rr++) {
    if(std::find(hits_id_to_ignore.begin(), hits_id_to_ignore.end(), i_rr) == hits_id_to_ignore.end()) {
      double val =  dEdx.at(i_rr) - fexp->Eval(rr.at(i_rr));
      sq_sum += (val - mean) * (val - mean);
    }
  }
  double stdev = std::sqrt(sq_sum / cont);

  for(size_t i_rr =  15; i_rr < rr.size(); i_rr++) {
    if(std::find(hits_id_to_ignore.begin(), hits_id_to_ignore.end(), i_rr) == hits_id_to_ignore.end()) {
      if(abs(dEdx.at(i_rr) - fexp->Eval(rr.at(i_rr)) - mean) > 1.5*stdev)  {
        hits_id_to_ignore.push_back(i_rr);
      }
    }
  }


  delete gr_dEdx_corrected;
  delete fexp;

  return hits_id_to_ignore;
};



double Hypfit::robust_max_x(TF1 *f, double xmin, double xmax, int n_scan) {
  double step = (xmax - xmin) / (n_scan - 1);
  std::vector<double> ys(n_scan);
  int best_i = 0;
  double best_y = f->Eval(xmin);
  ys[0] = best_y;
  for (int i = 1; i < n_scan; i++) {
    double y = f->Eval(xmin + i * step);
    ys[i] = y;
    if (y > best_y) { best_y = y; best_i = i; }
  }
  double best_x = xmin + best_i * step;
  if (best_i <= 0 || best_i >= n_scan - 1) return best_x;  // edge, skip refinement
  double y0 = ys[best_i - 1], y1 = ys[best_i], y2 = ys[best_i + 1];
  double denom = (y0 - 2.0 * y1 + y2);
  if (denom == 0) return best_x;
  double delta = 0.5 * (y0 - y2) / denom;  // vertex offset, in units of `step`
  return best_x + delta * step;
}
    
map<int, vector<TF1*>> Hypfit::get_conv_function_map(int target_pdg, string mode, int max_rr) {
  map<int, vector<TF1*>> conv_tf1_map;

  double mass = pion_mass;
  if (target_pdg == 13)   mass = muon_mass;
  if (target_pdg == 2212) mass = proton_mass;

  for (int i_rr = 0; i_rr < max_rr; i_rr++) {
    double pitch = 0.32;
    double rr = (2 * i_rr + 1) / 2.;
    double this_KE = map_PhysdEdx[target_pdg]->KEFromRangeSpline(rr);
    double gamma = (this_KE / mass) + 1.0;
    double beta  = TMath::Sqrt(1 - (1.0 / (gamma * gamma)));
    double this_xi   = map_PhysdEdx[target_pdg]->Landau_xi(this_KE, pitch);
    double this_Wmax = map_PhysdEdx[target_pdg]->Get_Wmax(this_KE);
    double this_kappa = this_xi / this_Wmax;
    double this_dEdx_BB = map_PhysdEdx[target_pdg]->meandEdx(this_KE);
    double par[5] = {this_kappa, beta * beta, this_xi, this_dEdx_BB, pitch};

    TF1 *PDF = new TF1("", PhysdEdx::dEdx_PDF_function, -10., 20., 5);
    PDF->SetParameters(par[0], par[1], par[2], par[3], par[4]);
    double PDF_max = PDF->GetMaximumX();

    for (int i_p = 0; i_p < 3; i_p++) {
      double sigma = PhysdEdx::pdg_plane_map[i_p][target_pdg][0] +
        PhysdEdx::pdg_plane_map[i_p][target_pdg][1] * pow(PDF_max, PhysdEdx::pdg_plane_map[i_p][target_pdg][2]);

      TF1 *f_gaus = new TF1("f_gaus", "gaus", -10, 10);
      f_gaus->SetParameters(1.0, 0.0, sigma);

      TF1Convolution *fconv = new TF1Convolution(PDF, f_gaus, true);
      fconv->SetRange(0, 20);
      fconv->SetNofPointsFFT(10000);
      TF1 *f0 = new TF1("f0", fconv, 0, 20, 0);

      double f0_max_x = robust_max_x(f0, 0., 20., 2000.);
      double dx = PDF_max - f0_max_x;
      TF1 *f_shifted = new TF1("f_shifted",
        [f0, dx](double *x, double *p){ return f0->Eval(x[0] - dx); },
        0, 20, 0);

      conv_tf1_map[i_p].push_back(mode.find("shift") != std::string::npos ? f_shifted : f0);
    }
  }
  return conv_tf1_map;
}


double Hypfit::GetTLExtensionP(int target_PDG, vector<double> this_rr_vec, vector<double> this_dEdx_vec, vector<double> this_pitch_vec, int best_plane, string cleaning_method) {

  // 1. Select the correct pre-calculated map
  std::map<int, vector<TF1*>>* selected_map;
  if (target_PDG == 13)        selected_map = &conv_tf1_map_muon;
  else if (target_PDG == 211)   selected_map = &conv_tf1_map_pion;
  else if (target_PDG == 2212)  selected_map = &conv_tf1_map_proton;
  else return -1.0; // Unsupported PDG

  // 2. Get the vector for the specific plane
  vector<TF1*> this_conv_tf1_vec = (*selected_map)[best_plane];

  // 3. Clean hits
  vector<int> bad_hits = get_hits_to_ignore(this_rr_vec, this_dEdx_vec, cleaning_method);

  // 4. Calculate
  double reco_KE = NormLikelihood_w_convolution(this_dEdx_vec, this_rr_vec, this_pitch_vec, bad_hits,
                                                true,
                                                target_PDG, this_conv_tf1_vec);

  if(reco_KE < 0) return reco_KE;
  double reco_P = map_PhysdEdx[target_PDG]->KEtoMomentum(reco_KE)/1000.0;

  return reco_P;
}