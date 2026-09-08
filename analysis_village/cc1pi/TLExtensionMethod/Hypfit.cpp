
#include "Hypfit.h"

Hypfit::Hypfit(){
  map_PhysdEdx[13] = new PhysdEdx(13);
  map_PhysdEdx[211] = new PhysdEdx(211);
  map_PhysdEdx[321] = new PhysdEdx(321);
  map_PhysdEdx[2212] = new PhysdEdx(2212);

  std::cout << "Hypfit: Pre-calculating convolution maps..." << std::endl;
  conv_tf1_map_muon   = get_conv_function_map(13,   "none", 750);
  conv_tf1_map_pion   = get_conv_function_map(211,  "none", 750);
  conv_tf1_map_proton = get_conv_function_map(2212, "none", 750);
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





double Hypfit::NormLikelihood_w_convolution_each_time(const vector<double> & dEdx, const vector<double> & ResRange, const vector<double> & pitch, const vector<int> bad_hits, bool include_small_likelihood, int PID, int target_plane, const string &mode) {
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

  double best_additional_res_length = -0.1;
  double best_m2lnL = 99999.;

  int this_N_calo = dEdx.size();
  if(this_N_calo-1 <= 0){
    return -2222.;
  }
  if(this_N_calo - bad_hits.size() <= 4) {
    return -8888.;
  }
  int i_bestfit = -1;
  int this_N_hits = this_N_calo;

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

      double this_likelihood, this_likelihood_max;
      map_PhysdEdx[PID]->dEdx_PDF_and_max_w_convolution(
          this_KE, this_pitch, dEdx_measured, target_plane, PID, mode,
          this_likelihood, this_likelihood_max);

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

  double original_res_length = ResRange.at(this_N_calo - 1);
  double best_total_res_length = best_additional_res_length + original_res_length;

  if(i_bestfit == res_length_trial - 1){
    return -7777.;
  }
  else if(best_m2lnL > 99990.){
    return -6666.;
  }
  else if(fabs(best_m2lnL) < 1.0e-11){
    return -5555.;
  } else if(std::isnan(best_total_res_length)) {
    return -3333.;
  }

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

    if (mode.find("include_2cm") == std::string::npos)  {
      if(rr.at(i_rr) < 2) {
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

/*
map<int, vector<TF1*>> Hypfit::get_conv_function_map(int target_pdg, string mode, int max_rr) {
  map<int, vector<TF1*>> conv_tf1_map;
  double pitch = 0.32;

  double mass = pion_mass;
  if (target_pdg == 13)   mass = muon_mass;
  if (target_pdg == 2212) mass = proton_mass;

  for (int i_rr = 0; i_rr < max_rr; i_rr++) {
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
    double PDF_max = robust_max_x(PDF, 0., 10., 2000);

    for (int i_p = 0; i_p < 3; i_p++) {
      double sigma = PhysdEdx::pdg_plane_map[i_p][target_pdg][0] +
        PhysdEdx::pdg_plane_map[i_p][target_pdg][1] * pow(PDF_max, PhysdEdx::pdg_plane_map[i_p][target_pdg][2]);

      TF1 *f_gaus = new TF1("f_gaus", "gaus", -10, 10);
      f_gaus->SetParameters(1.0, 0.0, sigma);

      TF1Convolution *fconv = new TF1Convolution(PDF, f_gaus, true);
      fconv->SetRange(0, 20);
      fconv->SetNofPointsFFT(10000);

      TF1 *f0 = new TF1("f0", fconv, 0, 20, 0);
      double f0_max_x = robust_max_x(f0, 0., 10., 2000);
      double dx = PDF_max - f0_max_x;
      TF1 *f_shifted = new TF1("f_shifted",
        [f0, dx](double *x, double *p){ return f0->Eval(x[0] - dx); },
        0, 20, 0);

      conv_tf1_map[i_p].push_back(mode.find("shift") != std::string::npos ? f_shifted : f0);
    }
  }
  return conv_tf1_map;
}
*/


map<int, vector<TF1*>> Hypfit::get_conv_function_map(int target_pdg, string mode, int max_rr, bool use_data_map) {
  map<int, vector<TF1*>> conv_tf1_map;
  double pitch = 0.32;
  double mass = pion_mass;
  if (target_pdg == 13)   mass = muon_mass;
  if (target_pdg == 2212) mass = proton_mass;
  bool apply_shift = mode.find("shift") != std::string::npos;

  // Select MC vs. data-derived sigma/shift maps once, up front.
  const vector<std::map<int, vector<double>>> &sigma_map =
      use_data_map ? PhysdEdx::pdg_plane_map_data : PhysdEdx::pdg_plane_map;
  const vector<std::map<int, vector<double>>> &shift_map =
      use_data_map ? PhysdEdx::pdg_plane_shift_map_data : PhysdEdx::pdg_plane_shift_map;

  for (int i_rr = 0; i_rr < max_rr; i_rr++) {
    double rr = (2 * i_rr + 1) / 2.;
    double this_KE = map_PhysdEdx[target_pdg]->KEFromRangeSpline(rr);
    double gamma = (this_KE / mass) + 1.0;
    double beta  = TMath::Sqrt(1 - (1.0 / (gamma * gamma)));
    double this_xi   = map_PhysdEdx[target_pdg]->Landau_xi(this_KE, pitch);
    double this_Wmax = map_PhysdEdx[target_pdg]->Get_Wmax(this_KE);
    double this_kappa = this_xi / this_Wmax;
    double this_dEdx_BB = map_PhysdEdx[target_pdg]->meandEdx(this_KE);
    double par[5] = {this_kappa, beta * beta, this_xi, this_dEdx_BB, pitch};

    // Stack allocation: TF1Convolution clones these internally, so local instances auto-destruct
    TF1 PDF("", PhysdEdx::dEdx_PDF_function, -10., 20., 5);
    PDF.AddToGlobalList(kFALSE);
    PDF.SetParameters(par[0], par[1], par[2], par[3], par[4]);
    double PDF_max = robust_max_x(&PDF, 0., 10., 2000);

    for (int i_p = 0; i_p < 3; i_p++) {
      double sigma = sigma_map[i_p].at(target_pdg)[0] +
        sigma_map[i_p].at(target_pdg)[1] * pow(PDF_max, sigma_map[i_p].at(target_pdg)[2]);

      // Stack allocation for Gaussian resolution
      TF1 f_gaus("", "gaus", -10, 10);
      f_gaus.AddToGlobalList(kFALSE);
      f_gaus.SetParameters(1.0, 0.0, sigma);

      // TF1Convolution takes owned clones of PDF & f_gaus
      TF1Convolution *fconv = new TF1Convolution(&PDF, &f_gaus, true);
      fconv->SetRange(0, 20);
      fconv->SetNofPointsFFT(2000);

      // TF1 constructed with fconv takes ownership of fconv
      TF1 *f0 = new TF1("f0", fconv, 0, 10, 0);
      f0->AddToGlobalList(kFALSE);

      if (!apply_shift) {
        conv_tf1_map[i_p].push_back(f0);
        continue;
      }

      const vector<double> &sp = shift_map[i_p].at(target_pdg);
      bool has_shift_fit = !(sp[0] == -1 && sp[1] == -1 && sp[2] == -1);
      double dx = has_shift_fit ? (sp[0] + sp[1] * TMath::Exp(-rr / sp[2])) : 0.0;

      TF1 *f_shifted = new TF1("f_shifted",
        [f0, dx](double *x, double *p){ return f0->Eval(x[0] - dx); },
        0, 20, 0);
      f_shifted->AddToGlobalList(kFALSE);

      conv_tf1_map[i_p].push_back(f_shifted);
    }
  }
  return conv_tf1_map;
}


/*
double exp_decay_cpp(double x, double a, double b, double c) {
    return a * TMath::Exp(-b * x) + c;
}

double Hypfit::langau_cpp(double *x, double *p) {
    const double width = p[0];
    const double mpv   = p[1];
    const double area  = p[2];
    const double sigG  = p[3];
    if (width <= 0 || sigG <= 0) return 0.0;
    const double invsq2pi = 0.398942280401;
    const double np = 500.0;   // integration slices, same as the macro
    const double sc = 5.0;     // convolution extends +/- sc * sigG
    const double xx = x[0];
    double xlow = xx - sc * sigG;
    double xupp = xx + sc * sigG;
    double step = (xupp - xlow) / np;
    double sum = 0.0;
    for (double i = 1.0; i <= np / 2.0; i += 1.0) {
        double t1 = xlow + (i - 0.5) * step;
        double fland1 = TMath::Landau(t1, mpv, width) / width;
        sum += fland1 * TMath::Gaus(xx, t1, sigG);
        double t2 = xupp - (i - 0.5) * step;
        double fland2 = TMath::Landau(t2, mpv, width) / width;
        sum += fland2 * TMath::Gaus(xx, t2, sigG);
    }
    return area * step * sum * invsq2pi / sigG;
}
 
map<int, vector<TF1*>> Hypfit::get_langau_map(int target_pdg, int max_rr, bool use_data_map) {
  map<int, vector<TF1*>> langau_tf1_map;
  // Select MC vs. data-derived width/mpv/gsigma maps up front
  const vector<std::map<int, vector<double>>> &width_map =
      use_data_map ? PhysdEdx::pdg_plane_width_map_data : PhysdEdx::pdg_plane_width_map;
  const vector<std::map<int, vector<double>>> &mpv_map =
      use_data_map ? PhysdEdx::pdg_plane_mpv_map_data : PhysdEdx::pdg_plane_mpv_map;
  const vector<std::map<int, vector<double>>> &gsigma_map =
      use_data_map ? PhysdEdx::pdg_plane_gsigma_map_data : PhysdEdx::pdg_plane_gsigma_map;
  // Select low-RR maps up front
  const vector<std::map<int, vector<double>>> &low_rr_map_2p5 =
      use_data_map ? PhysdEdx::pdg_plane_low_rr_2p5_data : PhysdEdx::pdg_plane_low_rr_2p5;
  const vector<std::map<int, vector<double>>> &low_rr_map_3p5 =
      use_data_map ? PhysdEdx::pdg_plane_low_rr_3p5_data : PhysdEdx::pdg_plane_low_rr_3p5;
  const vector<std::map<int, vector<double>>> &low_rr_map_4p5 =
      use_data_map ? PhysdEdx::pdg_plane_low_rr_4p5_data : PhysdEdx::pdg_plane_low_rr_4p5;

  const int nsamp = 2000;  // matches SetNofPointsFFT(2000) in get_conv_function_map

  for (int i_rr = 0; i_rr < max_rr; i_rr++) {
    double rr = (2 * i_rr + 1) / 2.;
    for (int i_p = 0; i_p < 3; i_p++) {
      double this_width = -1.0;
      double this_mpv = -1.0;
      double this_gsigma = -1.0;
      bool has_low_rr_val = false;
       // Check for RR = 2.5 cm override
      if (std::abs(rr - 2.5) < 1e-4) {
        auto it = low_rr_map_2p5[i_p].find(target_pdg);
        if (it != low_rr_map_2p5[i_p].end() && it->second.size() >= 3 && it->second[0] != -1) {
          this_mpv    = it->second[0];
          this_gsigma = it->second[1];
          this_width  = it->second[2];
          has_low_rr_val = true;
        }
      } 
        
      // Check for RR = 3.5 cm override
      if (std::abs(rr - 3.5) < 1e-4) {
        auto it = low_rr_map_3p5[i_p].find(target_pdg);
        if (it != low_rr_map_3p5[i_p].end() && it->second.size() >= 3 && it->second[0] != -1) {
          this_mpv    = it->second[0];
          this_gsigma = it->second[1];
          this_width  = it->second[2];
          has_low_rr_val = true;
        }
      } 
      // Check for RR = 4.5 cm override
      else if (std::abs(rr - 4.5) < 1e-4) {
        auto it = low_rr_map_4p5[i_p].find(target_pdg);
        if (it != low_rr_map_4p5[i_p].end() && it->second.size() >= 3 && it->second[0] != -1) {
          this_mpv    = it->second[0];
          this_gsigma = it->second[1];
          this_width  = it->second[2];
          has_low_rr_val = true;
        }
      }
      // Fallback to exponential decay calculation if no low-RR override is available
      if (!has_low_rr_val) {
        const vector<double> &wp = width_map[i_p].at(target_pdg);
        const vector<double> &mp = mpv_map[i_p].at(target_pdg);
        const vector<double> &gp = gsigma_map[i_p].at(target_pdg);
        bool has_fit =
            !(wp[0] == -1 && wp[1] == -1 && wp[2] == -1) &&
            !(mp[0] == -1 && mp[1] == -1 && mp[2] == -1) &&
            !(gp[0] == -1 && gp[1] == -1 && gp[2] == -1);
        if (!has_fit) {
          // No fit for this (plane, pdg) -- push nullptr to preserve positions
          langau_tf1_map[i_p].push_back(nullptr);
          continue;
        }
        this_width  = exp_decay_cpp(rr, wp[0], wp[1], wp[2]);
        this_mpv    = exp_decay_cpp(rr, mp[0], mp[1], mp[2]);
        this_gsigma = exp_decay_cpp(rr, gp[0], gp[1], gp[2]);
      }

      // Build the raw langau TF1 as before, but only to sample it -- never
      // evaluated directly by consumers, so its per-call integration cost
      // is paid once here rather than once per hit downstream.
      TF1 f_langau_raw("", Hypfit::langau_cpp, 0., 20., 4);
      f_langau_raw.SetParameters(this_width, this_mpv, 1.0, this_gsigma);

      std::vector<double> xs(nsamp), ys(nsamp);
      for (int k = 0; k < nsamp; k++) {
        double xx = 20.0 * k / (nsamp - 1);
        xs[k] = xx;
        ys[k] = f_langau_raw.Eval(xx);
      }
      TGraph *g_langau = new TGraph(nsamp, xs.data(), ys.data());
      g_langau->SetBit(kCanDelete, kFALSE);

      TF1 *f_langau = new TF1("f_langau",
        [g_langau](double *x, double *p){ return g_langau->Eval(x[0]); },
        0, 20, 0);
      f_langau->AddToGlobalList(kFALSE);
      langau_tf1_map[i_p].push_back(f_langau);
    }
  }
  return langau_tf1_map;
}
*/

double power_law_cpp(double x, double a, double b, double c) {
    return a * std::pow(x, -b) + c;
}

double Hypfit::langau_cpp(double *x, double *p) {
    const double width = p[0];
    const double mpv   = p[1];
    const double area  = p[2];
    const double sigG  = p[3];
    if (width <= 0 || sigG <= 0) return 0.0;
    const double invsq2pi = 0.398942280401;
    const double np = 200.0;   // 200 integration steps is more than sufficient for high precision
    const double sc = 5.0;     // convolution extends +/- sc * sigG
    const double xx = x[0];
    double xlow = xx - sc * sigG;
    double xupp = xx + sc * sigG;
    double step = (xupp - xlow) / np;
    double sum = 0.0;
    for (double i = 1.0; i <= np / 2.0; i += 1.0) {
        double t1 = xlow + (i - 0.5) * step;
        double fland1 = TMath::Landau(t1, mpv, width) / width;
        sum += fland1 * TMath::Gaus(xx, t1, sigG);
        double t2 = xupp - (i - 0.5) * step;
        double fland2 = TMath::Landau(t2, mpv, width) / width;
        sum += fland2 * TMath::Gaus(xx, t2, sigG);
    }
    return area * step * sum * invsq2pi / sigG;
}

map<int, vector<TF1*>> Hypfit::get_langau_map(int target_pdg, int max_rr, bool use_data_map) {
  map<int, vector<TF1*>> langau_tf1_map;

  // Select MC vs. data-derived width/mpv/gsigma maps up front
  const vector<std::map<int, vector<double>>> &width_map =
      use_data_map ? PhysdEdx::pdg_plane_width_map_data : PhysdEdx::pdg_plane_width_map;
  const vector<std::map<int, vector<double>>> &mpv_map =
      use_data_map ? PhysdEdx::pdg_plane_mpv_map_data : PhysdEdx::pdg_plane_mpv_map;
  const vector<std::map<int, vector<double>>> &gsigma_map =
      use_data_map ? PhysdEdx::pdg_plane_gsigma_map_data : PhysdEdx::pdg_plane_gsigma_map;

  // Select low-RR maps up front
  const vector<std::map<int, vector<double>>> &low_rr_map_2p5 =
      use_data_map ? PhysdEdx::pdg_plane_low_rr_2p5_data : PhysdEdx::pdg_plane_low_rr_2p5;
  const vector<std::map<int, vector<double>>> &low_rr_map_3p5 =
      use_data_map ? PhysdEdx::pdg_plane_low_rr_3p5_data : PhysdEdx::pdg_plane_low_rr_3p5;

  constexpr int nsamp = 1000;    // 1000 points gives 0.02 MeV/cm resolution on [0, 20]
  constexpr double xmax = 20.0;
  constexpr double inv_dx = (nsamp - 1) / xmax;

  for (int i_rr = 0; i_rr < max_rr; i_rr++) {
    double rr = (2 * i_rr + 1) / 2.;

    for (int i_p = 0; i_p < 3; i_p++) {
      double this_width  = -1.0;
      double this_mpv    = -1.0;
      double this_gsigma = -1.0;
      bool has_low_rr_val = false;

      // Check for RR = 2.5 cm override
      if (std::abs(rr - 2.5) < 1e-4) {
        auto it = low_rr_map_2p5[i_p].find(target_pdg);
        if (it != low_rr_map_2p5[i_p].end() && it->second.size() >= 3 && it->second[0] != -1) {
          this_mpv       = it->second[0];
          this_gsigma    = it->second[1];
          this_width     = it->second[2];
          has_low_rr_val = true;
        }
      } 
      // Check for RR = 3.5 cm override
      else if (std::abs(rr - 3.5) < 1e-4) {
        auto it = low_rr_map_3p5[i_p].find(target_pdg);
        if (it != low_rr_map_3p5[i_p].end() && it->second.size() >= 3 && it->second[0] != -1) {
          this_mpv       = it->second[0];
          this_gsigma    = it->second[1];
          this_width     = it->second[2];
          has_low_rr_val = true;
        }
      } 


      // Fallback to Power-Law calculation if no low-RR override is available
      if (!has_low_rr_val) {
        const vector<double> &wp = width_map[i_p].at(target_pdg);
        const vector<double> &mp = mpv_map[i_p].at(target_pdg);
        const vector<double> &gp = gsigma_map[i_p].at(target_pdg);

        bool has_fit = !(wp[0] == -1 && wp[1] == -1 && wp[2] == -1) &&
                       !(mp[0] == -1 && mp[1] == -1 && mp[2] == -1) &&
                       !(gp[0] == -1 && gp[1] == -1 && gp[2] == -1);

        if (!has_fit) {
          langau_tf1_map[i_p].push_back(nullptr);
          continue;
        }

        this_width  = power_law_cpp(rr, wp[0], wp[1], wp[2]);
        this_mpv    = power_law_cpp(rr, mp[0], mp[1], mp[2]);
        this_gsigma = power_law_cpp(rr, gp[0], gp[1], gp[2]);
      }

      // 1. Precompute lookup array
      TF1 f_langau_raw("", Hypfit::langau_cpp, 0., xmax, 4);
      f_langau_raw.SetParameters(this_width, this_mpv, 1.0 /* area */, this_gsigma);

      std::vector<double> ys(nsamp);
      for (int k = 0; k < nsamp; k++) {
        double xx = xmax * k / (nsamp - 1);
        ys[k] = f_langau_raw.Eval(xx);
      }

      // 2. Direct O(1) index calculation inside lambda (No TGraph, No Binary Search)
      TF1 *f_langau = new TF1(
          "f_langau",
          [ys](double *x, double *p) -> double {
            double xx = x[0];
            if (xx <= 0.0) return ys[0];
            if (xx >= 20.0) return ys.back();

            double pos = xx * inv_dx;
            int idx = static_cast<int>(pos);
            double frac = pos - idx;
            return ys[idx] + frac * (ys[idx + 1] - ys[idx]);
          },
          0, xmax, 0);

      f_langau->AddToGlobalList(kFALSE);
      langau_tf1_map[i_p].push_back(f_langau);
    }
  }
  return langau_tf1_map;
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


// == Hypfit.cxx
/*
// == Asymmetric pitch binning informed by:
// ==   (1) the hit-pitch distribution, sharply peaked near the geometric
// ==       floor (~0.3 cm) with a long, sparse tail to ~2 cm
// ==   (2) the theory MPV-vs-pitch curve, which varies fastest below
// ==       ~0.5 cm and flattens out towards 2 cm
// == Lower edge fixed at 0.30 cm (minimum physical wire pitch).
vector<double> Hypfit::get_pitch_bin_edges() const {
  static const vector<double> edges = {
    0.30, 0.315, 0.325, 0.34, 0.36, 0.38, 0.40, 0.43, 0.45, 0.48,
    0.51, 0.55, 0.60, 0.65, 0.72, 0.80, 0.90, 1.00,
    1.15, 1.30, 1.45, 1.60, 1.80, 2.00
  };
  return edges;
}
// == Returns the bin index for a given pitch value.
// == Clamps to the first bin if pitch < edges.front() and to the
// == last bin if pitch > edges.back() (i.e. defaults to the last bin
// == when pitch exceeds the max bin boundary, as requested).
int Hypfit::get_pitch_bin_index(double pitch, const vector<double> &edges) const {
  int n_bins = (int)edges.size() - 1;
  if (pitch <= edges.front()) return 0;
  if (pitch >= edges.back())  return n_bins - 1;
  for (int i_bin = 0; i_bin < n_bins; i_bin++) {
    if (pitch >= edges[i_bin] && pitch < edges[i_bin + 1]) return i_bin;
  }
  return n_bins - 1; // == fallback, should not be reached
}

map<int, vector<vector<TF1*>>> Hypfit::get_conv_function_map_w_pitch(int target_pdg, string mode, int max_rr) {
  map<int, vector<vector<TF1*>>> conv_tf1_map; // == [plane][pitch_bin][rr]
  vector<double> pitch_edges = get_pitch_bin_edges();
  int n_pitch_bins = (int)pitch_edges.size() - 1;
  bool use_shifted = (mode.find("shift") != std::string::npos);

  double mass = pion_mass;
  if (target_pdg == 13)   mass = muon_mass;
  if (target_pdg == 2212) mass = proton_mass;

  for (int plane = 0; plane < 3; plane++)
    conv_tf1_map[plane] = vector<vector<TF1*>>(n_pitch_bins);

  for (int i_pitch = 0; i_pitch < n_pitch_bins; i_pitch++) {
    double this_pitch = 0.5 * (pitch_edges[i_pitch] + pitch_edges[i_pitch + 1]);

    for (int i_rr = 0; i_rr < max_rr; i_rr++) {
      double rr = (2 * i_rr + 1) / 2.;
      double this_KE = map_PhysdEdx[target_pdg]->KEFromRangeSpline(rr);
      double gamma = (this_KE / mass) + 1.0;
      double beta  = TMath::Sqrt(1 - (1.0 / (gamma * gamma)));
      double this_xi     = map_PhysdEdx[target_pdg]->Landau_xi(this_KE, this_pitch);
      double this_Wmax   = map_PhysdEdx[target_pdg]->Get_Wmax(this_KE);
      double this_kappa  = this_xi / this_Wmax;
      double this_dEdx_BB = map_PhysdEdx[target_pdg]->meandEdx(this_KE);
      double par[5] = {this_kappa, beta * beta, this_xi, this_dEdx_BB, this_pitch};

      TF1 *PDF = new TF1("", PhysdEdx::dEdx_PDF_function, -10., 20., 5);
      PDF->AddToGlobalList(kFALSE); // == keep gROOT's function list from growing
      PDF->SetParameters(par[0], par[1], par[2], par[3], par[4]);
      double PDF_max = robust_max_x(PDF, 0., 10., 2000);

      for (int plane = 0; plane < 3; plane++) {
        double sigma = PhysdEdx::pdg_plane_map[plane][target_pdg][0] +
          PhysdEdx::pdg_plane_map[plane][target_pdg][1] * pow(PDF_max, PhysdEdx::pdg_plane_map[plane][target_pdg][2]);

        TF1 *f_gaus = new TF1("f_gaus", "gaus", -10, 10);
        f_gaus->AddToGlobalList(kFALSE);
        f_gaus->SetParameters(1.0, 0.0, sigma);

        TF1Convolution *fconv = new TF1Convolution(PDF, f_gaus, true);
        fconv->SetRange(0, 20);
        fconv->SetNofPointsFFT(1000);
        TF1 *f0 = new TF1("f0", fconv, 0, 20, 0);
        f0->AddToGlobalList(kFALSE);
        double f0_max_x = robust_max_x(f0, 0., 10., 2000);
        double dx = PDF_max - f0_max_x;
        TF1 *f_shifted = new TF1("f_shifted",
          [f0, dx](double *x, double *p){ return f0->Eval(x[0] - dx); },
          0, 20, 0);
        f_shifted->AddToGlobalList(kFALSE);

        conv_tf1_map[plane][i_pitch].push_back(use_shifted ? f_shifted : f0);
      } // == plane
    } // == rr
  } // == pitch

  return conv_tf1_map;
}

double Hypfit::NormLikelihood_w_convolution_w_pitch(const vector<double> & dEdx, const vector<double> & ResRange,
    const vector<double> & pitch, const vector<int> bad_hits, bool include_small_likelihood,
    int PID, vector<vector<TF1*>> tf1_map) {
  // == tf1_map : [pitch_bin][rr], i.e. get_conv_function_map_w_pitch(...)[plane]
  if(!(PID == 13 || PID == 2212 || PID == 211)){
    return -9999.;
  }
  double min_additional_res_length = 0.;
  double max_additional_res_length = max_additional_res_length_pion;
  double res_length_step = res_length_step_pion;
  if(PID == 2212){
    max_additional_res_length = max_additional_res_length_proton;
    res_length_step = res_length_step_proton;
  }
  int res_length_trial = (max_additional_res_length - min_additional_res_length) / res_length_step;

  double best_additional_res_length = -0.1;
  double best_m2lnL = 99999.;
  int this_N_calo = dEdx.size();
  if(this_N_calo-1 <= 0){
    return -2222.;
  }
  if(this_N_calo - bad_hits.size() <= 4) {
    return -8888.;
  }
  int i_bestfit = -1;
  int this_N_hits = this_N_calo;
  vector<double> pitch_edges = get_pitch_bin_edges();

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

      // == Defaults to the last pitch bin if this_pitch exceeds the max boundary
      int this_pitch_bin = get_pitch_bin_index(this_pitch, pitch_edges);
      const vector<TF1*> & tf1_vec = tf1_map.at(this_pitch_bin);

      double this_likelihood     = map_PhysdEdx[PID]->dEdx_PDF_w_convolution_f1(this_KE, this_res_length, dEdx_measured, this_pitch, tf1_vec);
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
    if(this_m2lnL < best_m2lnL){
      best_m2lnL = this_m2lnL;
      best_additional_res_length = this_additional_res_length;
      i_bestfit = i;
    }
  }

  double original_res_length = ResRange.at(this_N_calo - 1);
  double best_total_res_length = best_additional_res_length + original_res_length;
  if(i_bestfit == res_length_trial - 1){
    return -7777.;
  }
  else if(best_m2lnL > 99990.){
    return -6666.;
  }
  else if(fabs(best_m2lnL) < 1.0e-11){
    return -5555.;
  } else if(std::isnan(best_total_res_length)) {
    return -3333.;
  }
  double best_KE = map_PhysdEdx[PID] -> KEFromRangeSpline(best_total_res_length);
  return best_KE;
}






// Hypfit.cpp
vector<map<int, vector<TF1*>>> Hypfit::get_conv_function_map_w_x_ranges(
    int target_pdg, string mode, int max_rr) {

  vector<map<int, vector<TF1*>>> conv_tf1_map_w_x_ranges(x_ranges.size());

  double pitch = 0.32;
  double mass = pion_mass;
  if (target_pdg == 13)   mass = muon_mass;
  if (target_pdg == 2212) mass = proton_mass;
  bool apply_shift = mode.find("shift") != std::string::npos;

  for (size_t i_xr = 0; i_xr < x_ranges.size(); i_xr++) {
    map<int, vector<TF1*>> &conv_tf1_map = conv_tf1_map_w_x_ranges[i_xr];

    for (int i_rr = 0; i_rr < max_rr; i_rr++) {
      double rr = (2 * i_rr + 1) / 2.;
      double this_KE = map_PhysdEdx[target_pdg]->KEFromRangeSpline(rr);
      double gamma = (this_KE / mass) + 1.0;
      double beta  = TMath::Sqrt(1 - (1.0 / (gamma * gamma)));
      double this_xi   = map_PhysdEdx[target_pdg]->Landau_xi(this_KE, pitch);
      double this_Wmax = map_PhysdEdx[target_pdg]->Get_Wmax(this_KE);
      double this_kappa = this_xi / this_Wmax;
      double this_dEdx_BB = map_PhysdEdx[target_pdg]->meandEdx(this_KE);
      double par[5] = {this_kappa, beta * beta, this_xi, this_dEdx_BB, pitch};

      // Stack allocation: TF1Convolution clones these internally, so local instances auto-destruct
      TF1 PDF("", PhysdEdx::dEdx_PDF_function, -10., 20., 5);
      PDF.AddToGlobalList(kFALSE);
      PDF.SetParameters(par[0], par[1], par[2], par[3], par[4]);
      double PDF_max = robust_max_x(&PDF, 0., 10., 2000);
        
      for (int i_p = 0; i_p < 3; i_p++) {
        const auto &sigma_pars = PhysdEdx::GetPlaneMapEntry(
            PhysdEdx::pdg_plane_map_w_x_ranges, i_xr, i_p, target_pdg);
        double sigma = sigma_pars[0] + sigma_pars[1] * pow(PDF_max, sigma_pars[2]);

        // Stack allocation for Gaussian resolution
        TF1 f_gaus("", "gaus", -10, 10);
        f_gaus.AddToGlobalList(kFALSE);
        f_gaus.SetParameters(1.0, 0.0, sigma);

        // TF1Convolution takes owned clones of PDF & f_gaus
        TF1Convolution *fconv = new TF1Convolution(&PDF, &f_gaus, true);
        fconv->SetRange(0, 20);
        fconv->SetNofPointsFFT(1000);

        // TF1 constructed with fconv takes ownership of fconv
        TF1 *f0 = new TF1("f0", fconv, 0, 10, 0);
        f0->AddToGlobalList(kFALSE);

        if (!apply_shift) {
          conv_tf1_map[i_p].push_back(f0);
          continue;
        }

        const auto &sp = PhysdEdx::GetPlaneMapEntry(
    PhysdEdx::pdg_plane_shift_map_w_x_ranges, i_xr, i_p, target_pdg);
          
        bool has_shift_fit = !(sp[0] == -1 && sp[1] == -1 && sp[2] == -1);
        double dx = has_shift_fit ? (sp[0] + sp[1] * TMath::Exp(-rr / sp[2])) : 0.0;

        TF1 *f_shifted = new TF1("f_shifted",
          [f0, dx](double *x, double *p){ return f0->Eval(x[0] - dx); },
          0, 20, 0);
        f_shifted->AddToGlobalList(kFALSE);
        conv_tf1_map[i_p].push_back(f_shifted);
      }
    }
  }
  return conv_tf1_map_w_x_ranges;
}



// Hypfit.cpp
double Hypfit::NormLikelihood_w_convolution_w_x_ranges(
    const vector<double> & dEdx,
    const vector<double> & ResRange,
    const vector<double> & pitch,
    const vector<double> & Xpos,
    const vector<int> bad_hits,
    bool include_small_likelihood,
    int PID,
    const vector<vector<TF1*>> & tf1_vec_by_xrange) {

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

      // == Find which x-range this hit falls into, and use the matching conv functions
      int xrange_idx = GetXRangeIndex(Xpos.at(j));
      if (xrange_idx < 0 || xrange_idx >= static_cast<int>(tf1_vec_by_xrange.size())) {
        // Hit falls outside all defined x-ranges: treat like an invalid-likelihood hit
        num_invalid_hits++;
        continue;
      }
      const vector<TF1*> & tf1_vec = tf1_vec_by_xrange[xrange_idx];

      // == Do not use first and last N_skip this
      double this_res_length = ResRange.at(j) + this_additional_res_length;
      double this_KE = map_PhysdEdx[PID]->KEFromRangeSpline(this_res_length);
      double dEdx_measured = dEdx.at(j);
      // == Likelihood
      double this_pitch = pitch.at(j);
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
*/