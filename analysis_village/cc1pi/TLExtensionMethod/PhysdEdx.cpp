#include "PhysdEdx.h"

ClassImp(PhysdEdx)  // link with ROOT dictionary

ROOT::Math::VavilovAccurate vav;

using namespace std;

PhysdEdx::PhysdEdx()
  : pdgcode(0)
  , mass(0)
  , charge(0)
  , sp_KE_range(0)
  , sp_range_KE(0){
}

PhysdEdx::PhysdEdx(int pdg)
  : pdgcode(0)
  , mass(0)
  , charge(0)
  , sp_KE_range(0)
  , sp_range_KE(0){
  SetPdgCode(pdg);
}

void PhysdEdx::SetPdgCode(int pdg){

  pdgcode = pdg;

  if (abs(pdgcode) == 13){//muon
    mass = 105.6583755;
    charge = 1;
  }
  else if (abs(pdgcode) == 211){//pion
    mass = 139.57039;
    charge = 1;
  }
  else if (abs(pdgcode) == 321){//kaon
    mass = 493.677;
    charge = 1;
  }
  else if (pdgcode == 2212){//proton
    mass = 938.27208816;
    charge = 1;
  }
  else{
    cout<<"Unknown pdg code "<<pdgcode<<endl;
    exit(1);
  }

  CreateSplines();
}

double PhysdEdx::densityEffect(double beta, double gamma){
  // == Estimate the density correction
  double density_y = TMath::Log10(beta * gamma);
  double ln10 = TMath::Log(10);
  double this_delta = 0.;
  if(density_y > density_y1){
    this_delta = 2.0 * ln10 * density_y - density_C;
  }
  else if (density_y < density_y0){
    this_delta = 0.;
  }
  else{
    this_delta = 2.0 * ln10 * density_y - density_C + density_a * pow(density_y1 - density_y, density_k);
  }

  return this_delta;
}

double PhysdEdx::betaGamma(double KE){
  double gamma, beta;
  gamma = (KE + mass) / mass;
  beta = sqrt( 1 - 1/pow(gamma,2));

  return beta*gamma;
}

double PhysdEdx::Landau_xi(double KE, double pitch){
  double gamma = (KE/mass)+1.0;
  double beta = TMath::Sqrt(1-(1.0/(gamma*gamma)));
  if (beta < 1e-6) beta = 1e-6; // Safety floor
  double xi = rho * pitch * 0.5 * K * (Z / A) * pow(1. / beta, 2);
  return xi;
}

double PhysdEdx::Get_Wmax(double KE){
  double gamma = (KE/mass)+1.0;
  double beta = TMath::Sqrt(1-(1.0/(gamma*gamma)));
  double Wmax = (2.0 * me * pow(beta * gamma, 2)) / (1.0 + 2.0 * me * (gamma / mass) + pow((me / mass),2));

  return Wmax;
}

double PhysdEdx::meandEdx(double KE){

  double gamma = (KE + mass) / mass;
  double beta = sqrt( 1 - 1/pow(gamma,2));
  double wmax = Get_Wmax(KE);
  double dEdX = (rho*K*Z*pow(charge,2))/(A*pow(beta,2))*(0.5*log(2*me*pow(gamma,2)*pow(beta,2)*wmax/pow(I,2)) - pow(beta,2) - densityEffect( beta, gamma )/2 );

  return dEdX;
}

double PhysdEdx::MPVdEdx(double KE, double pitch){

  //KE is kinetic energy in MeV
  //pitch is in cm
  double gamma = (KE + mass) / mass;
  double beta = sqrt( 1 - 1/pow(gamma,2));

  double xi = Landau_xi(KE, pitch);

  double eloss_mpv = xi*(log( 2*me*pow(gamma,2)*pow(beta,2) / I ) + log( xi / I ) + 0.2 - pow(beta,2) - densityEffect( beta, gamma ) )/pitch;

  return eloss_mpv;
}

double PhysdEdx::IntegratedEdx(double KE0, double KE1, int n){

  if (KE0>KE1) swap(KE0, KE1);

  double step = (KE1-KE0)/n;

  double area = 0;

  for (int i = 0; i<n; ++i){
    double dEdx = meandEdx(KE0 + (i+0.5)*step);
    if (dEdx)
      area += 1/dEdx*step;
  }
  return area;
}

double PhysdEdx::RangeFromKE(double KE){

  return IntegratedEdx(0, KE);
}

void PhysdEdx::CreateSplines(int np, double minke, double maxke){

  if (sp_KE_range) delete sp_KE_range;
  if (sp_range_KE) delete sp_range_KE;

  for (const auto & x : spmap){
    if (x.second) delete x.second;
  }
  spmap.clear();

  double *KE = new double[np];
  double *Range = new double[np];

  for (int i = 0; i<np; ++i){
    double ke = pow(10, log10(minke)+i*log10(maxke/minke)/np);
    KE[i] = ke;
    Range[i] = RangeFromKE(ke);
  }

  sp_KE_range = new TSpline3("sp_KE_range", KE, Range, np, "b2e2", 0, 0);
  sp_range_KE = new TSpline3("sp_range_KE", Range, KE, np, "b2e2", 0, 0);

  delete[] KE;
  delete[] Range;

  cout<<"Done creating splines for particle with pdgcode "<<pdgcode<<endl;
}

double PhysdEdx::RangeFromKESpline(double KE){
  if (!sp_KE_range){
    cout<<"Spline does not exist."<<endl;
    exit(1);
  }
  return sp_KE_range->Eval(KE);
}

double PhysdEdx::KEFromRangeSpline(double range){
  if (!sp_range_KE){
    cout<<"Spline does not exit."<<endl;
    exit(1);
  }
  return sp_range_KE->Eval(range);
}

double PhysdEdx::KEAtLength(double KE0, double tracklength){

  int iKE = int(KE0);

  if (spmap.find(iKE)==spmap.end()){
    CreateSplineAtKE(iKE);
  }
  double deltaE = spmap[iKE]->Eval(tracklength);

  if (deltaE < 0) return 0;//cout<<"Negative delta E: "<<deltaE<<endl;
  if (KE0 - deltaE < 0) return 0;//cout<<"Negative KE: "<<KE0 - deltaE<<endl;

  return KE0 - deltaE;
}

double PhysdEdx::KEtoMomentum(double KE){
  return sqrt(pow(KE, 2) + 2.0 * KE * mass);
}

double PhysdEdx::MomentumtoKE(double momentum){
  return sqrt(pow(momentum, 2) + pow(mass, 2)) - mass;
}

void PhysdEdx::CreateSplineAtKE(int iKE){

  double KE0 = iKE;

  // Sample every 10 MeV
  int np = int(KE0/10);
  double *deltaE;
  double *trklength;
  if (np>1){
    deltaE = new double[np];
    trklength = new double[np];
    for (int i = 0; i<np; ++i){
      double KE = KE0 - i*10;
      deltaE[i] = KE0 - KE;
      trklength[i] = IntegratedEdx(KE, KE0);
    }
  }
  else{
    cout<<"KE too low: "<<iKE<<endl;
    np = 2;
    deltaE = new double[np];
    trklength = new double[np];
    deltaE[0] = 0;
    trklength[0] = 0;
    deltaE[1] = KE0;
    trklength[1] = RangeFromKE(KE0);
  }

  spmap[iKE] = new TSpline3(Form("KE %d",iKE), trklength, deltaE, np, "b2e2", 0, 0);
  delete[] trklength;
  delete[] deltaE;
}

double PhysdEdx::dEdx_PDF_function(double *x, double *par){
  // == par[5] = {kappa, beta^2, xi, <dE/dx>BB, width}

  // 1. Extreme Value Guard
  if (par[0] < 1e-5 || par[1] < 1e-5 || std::isnan(par[0]) || std::isinf(par[0])) {
    return 0.0; 
  }
    
  double a = par[2] / par[4];
  double b = (0.422784 + par[1] + log(par[0])) * par[2] / par[4] + par[3];
  double y = (x[0] - b) / a;

  double this_vav = 0.;

  if(par[0] < 0.01){ // == Landau
    this_vav = TMath::Landau(y);
    this_vav =  this_vav / a;
  }
  else if(par[0] > 10.){ // == Gaussian
    double mu = vav.Mean(par[0], par[1]);
    double sigma = sqrt(vav.Variance(par[0], par[1]));
    this_vav =  TMath::Gaus(y, mu, sigma);
  }
  else{ // == Vavilov
    this_vav =  vav.Pdf(y, par[0], par[1]);
    this_vav =  this_vav / a;
  }
    
  return this_vav;
}

double PhysdEdx::dEdx_PDF(double KE, double pitch, double dEdx){
  if (KE < 0.0001) KE = 0.0001; // Don't allow KE to go to zero
    
  double gamma = (KE/mass)+1.0;
  double beta = TMath::Sqrt(1-(1.0/(gamma*gamma)));
  double this_xi = Landau_xi(KE, pitch);
  double this_Wmax = Get_Wmax(KE);
  double this_kappa = this_xi / this_Wmax;
  double this_dEdx_BB = meandEdx(KE);
  double par[5] = {this_kappa, beta * beta, this_xi, this_dEdx_BB, pitch};

  TF1 *PDF = new TF1("", PhysdEdx::dEdx_PDF_function, -10., 20., 5); 
  PDF -> SetParameters(par[0], par[1], par[2], par[3], par[4]);

  double out = PDF -> Eval(dEdx);
  delete PDF;
  return out;
}

void PhysdEdx::dEdx_PDF_and_max_w_convolution(double KE, double pitch, double dEdx,
                                               int target_plane, int PDG, const string &mode,
                                               double &out_likelihood, double &out_likelihood_max){
  double gamma = (KE/mass)+1.0;
  double beta = TMath::Sqrt(1-(1.0/(gamma*gamma)));
  double this_xi = Landau_xi(KE, pitch);
  double this_Wmax = Get_Wmax(KE);
  double this_kappa = this_xi / this_Wmax;
  double this_dEdx_BB = meandEdx(KE);
  double par[5] = {this_kappa, beta * beta, this_xi, this_dEdx_BB, pitch};

  TF1 *PDF = new TF1("", PhysdEdx::dEdx_PDF_function, -10., 20., 5);
  PDF -> SetParameters(par[0], par[1], par[2], par[3], par[4]);
  double PDF_max = PDF->GetMaximumX();

  double sigma = pdg_plane_map[target_plane][PDG][0] +
                     pdg_plane_map[target_plane][PDG][1] * pow(PDF_max, pdg_plane_map[target_plane][PDG][2]);

  TF1 *f_gaus = new TF1("", "gaus", -10, 10);
  f_gaus->SetParameters(1.0, 0.0, sigma); // norm, mean, sigma

  // Create the convolution object -- built exactly once per hit now,
  // shared between the likelihood-value and peak-height computations.
  TF1Convolution *fconv = new TF1Convolution(PDF, f_gaus, true);
  fconv->SetRange(0, 20);
  fconv->SetNofPointsFFT(10000); // resolution of FFT
  TF1 *f0 = new TF1("f", fconv, 0, 20, 0);

  // Peak height is invariant under a pure horizontal shift, so f0's max
  // IS the shifted function's max too -- no need to build f_shifted at all.
  out_likelihood_max = f0->GetMaximum();

  if (mode.find("shift") != std::string::npos) {
    double f0_max_x = f0->GetMaximumX();
    double dx = PDF_max - f0_max_x;
    out_likelihood = f0->Eval(dEdx - dx); // algebraically identical to f_shifted->Eval(dEdx)
  } else {
    out_likelihood = f0->Eval(dEdx);
  }

  delete PDF;
  delete f0;
  delete f_gaus;
  delete fconv;
}


double PhysdEdx::dEdx_PDF_w_tf1_vec(double KE, double rr, double dEdx, double pitch, vector<TF1*> tf1_vec){
  if(rr < tf1_vec.size()) return tf1_vec[(int)rr]->Eval(dEdx);
  else return dEdx_PDF(KE, pitch, dEdx);
}

double PhysdEdx::dEdx_PDF_max_w_tf1_vec(double KE, double rr, double pitch, vector<TF1*> tf1_vec){
  if(rr < tf1_vec.size()) return tf1_vec[(int)rr]->GetMaximum();
  else return dEdx_PDF_max(KE, pitch, 0);
}

double PhysdEdx::dEdx_PDF_max(double KE, double pitch, double dEdx){

  double gamma = (KE/mass)+1.0;
  double beta = TMath::Sqrt(1-(1.0/(gamma*gamma)));
  double this_xi = Landau_xi(KE, pitch);
  double this_Wmax = Get_Wmax(KE);
  double this_kappa = this_xi / this_Wmax;
  double this_dEdx_BB = meandEdx(KE);
  double par[5] = {this_kappa, beta * beta, this_xi, this_dEdx_BB, pitch};

  TF1 *PDF = new TF1("", PhysdEdx::dEdx_PDF_function, -10., 20., 5);
  PDF -> SetParameters(par[0], par[1], par[2], par[3], par[4]);

  double out = PDF ->GetMaximum();
  delete PDF;
  return out;
}

double PhysdEdx::dEdx_Gaus_Sigma(double KE, double pitch){

  double gamma = (KE/mass)+1.0;
  double beta = TMath::Sqrt(1-(1.0/(gamma*gamma)));
  double this_xi = Landau_xi(KE, pitch);
  double this_Wmax = Get_Wmax(KE);
  double this_kappa = this_xi / this_Wmax;

  double sigma = sqrt(vav.Variance(this_kappa, beta * beta));

  return sigma;
}

PhysdEdx::~PhysdEdx(){

}

vector<std::map<int, vector<double>>> PhysdEdx::pdg_plane_map = {
  // Plane 0
  {
    {13, {0.15414618, 0.0073370046, 3.0989568}},
    {2212, {-1, -1, -1}},
    {211, {0.19965545, 0.002953214, 3.3149893}}
  },
  // Plane 1
  {
    {13, {0.21482707, 0.0094663776, 2.9415465}},
    {2212, {-1, -1, -1}},
    {211, {0.25435054, 0.0048757179, 3.0440928}}
  },
  // Plane 2
  {
    {13, {0.052186123, 0.0055492335, 3.3267398}},
    {2212, {-1, -1, -1}},
    {211, {0.079418371, 0.0031402731, 3.4434963}}
  }
};

vector<std::map<int, vector<double>>> PhysdEdx::pdg_plane_shift_map = {
  // Plane 0
  {
    {13, {0.031548719, -0.35704881, 6.3348627}},
    {2212, {-1, -1, -1}},
    {211, {-0.0013209031, -0.94333909, 3.3408714}}
  },
  // Plane 1
  {
    {13, {0.018715757, -0.44487653, 5.8912638}},
    {2212, {-1, -1, -1}},
    {211, {-0.012385964, -0.84908051, 3.618248}}
  },
  // Plane 2
  {
    {13, {0.027761254, -0.19153523, 6.2212352}},
    {2212, {-1, -1, -1}},
    {211, {0.020253528, -0.79405887, 3.2851615}}
  }
};

vector<std::map<int, vector<double>>> PhysdEdx::pdg_plane_map_data = {
  // Plane 0
  {
    {13, {0.17635079, 0.0070148444, 3.2658645}},
    {2212, {-1, -1, -1}},
    {211, {0.21591654, 0.0039822186, 3.2254805}}
  },
  // Plane 1
  {
    {13, {0.24245998, 0.0095666313, 3.0463578}},
    {2212, {-1, -1, -1}},
    {211, {0.27999877, 0.0051358339, 3.106631}}
  },
  // Plane 2
  {
    {13, {0.0022153164, 0.012065375, 2.8947866}},
    {2212, {-1, -1, -1}},
    {211, {0.045043346, 0.005089308, 3.2042905}}
  }
};

vector<std::map<int, vector<double>>> PhysdEdx::pdg_plane_shift_map_data = {
  // Plane 0
  {
    {13, {0.070567998, -0.17169503, 9.444994}},
    {2212, {-1, -1, -1}},
    {211, {0.033207138, -0.55974003, 3.700684}}
  },
  // Plane 1
  {
    {13, {0.030349388, -0.25544324, 10.484708}},
    {2212, {-1, -1, -1}},
    {211, {-0.0055525167, -0.48336116, 4.584093}}
  },
  // Plane 2
  {
    {13, {0.0010592869, -0.16968473, 12.269416}},
    {2212, {-1, -1, -1}},
    {211, {-0.028664558, -0.40390918, 5.0316504}}
  }
};




/*
Langau fits
*/

vector<std::map<int, vector<double>>> PhysdEdx::pdg_plane_mpv_map = {
  // Plane 0
  {
    {13, {-1, -1, -1}},
    {2212, {-1, -1, -1}},
    {211, {6.159844, 0.58102941, 1.1984958}}
  },
  // Plane 1
  {
    {13, {-1, -1, -1}},
    {2212, {-1, -1, -1}},
    {211, {6.0261733, 0.57380522, 1.1851848}}
  },
  // Plane 2
  {
    {13, {-1, -1, -1}},
    {2212, {-1, -1, -1}},
    {211, {6.1511992, 0.5632623, 1.1547571}}
  }
};

vector<std::map<int, vector<double>>> PhysdEdx::pdg_plane_gsigma_map = {
  // Plane 0
  {
    {13, {-1, -1, -1}},
    {2212, {-1, -1, -1}},
    {211, {1.2230179, 0.75143221, 0.14647473}}
  },
  // Plane 1
  {
    {13, {-1, -1, -1}},
    {2212, {-1, -1, -1}},
    {211, {0.94553995, 0.4224341, 0.079181631}}
  },
  // Plane 2
  {
    {13, {-1, -1, -1}},
    {2212, {-1, -1, -1}},
    {211, {1.7481279, 1.0718866, 0.094531357}}
  }
};

vector<std::map<int, vector<double>>> PhysdEdx::pdg_plane_width_map = {
  // Plane 0
  {
    {13, {-1, -1, -1}},
    {2212, {-1, -1, -1}},
    {211, {0.054083349, 0.63642048, 0.096753482}}
  },
  // Plane 1
  {
    {13, {-1, -1, -1}},
    {2212, {-1, -1, -1}},
    {211, {6.3736161, 3.202817, 0.11485853}}
  },
  // Plane 2
  {
    {13, {-1, -1, -1}},
    {2212, {-1, -1, -1}},
    {211, {19.899956, 0.0010827431, -19.718831}}
  }
};

vector<std::map<int, vector<double>>> PhysdEdx::pdg_plane_mpv_map_data = {
  // Plane 0
  {
    {13, {-1, -1, -1}},
    {2212, {-1, -1, -1}},
    {211, {6.5151238, 0.62522296, 1.3099476}}
  },
  // Plane 1
  {
    {13, {-1, -1, -1}},
    {2212, {-1, -1, -1}},
    {211, {6.3491059, 0.62106716, 1.276556}}
  },
  // Plane 2
  {
    {13, {-1, -1, -1}},
    {2212, {-1, -1, -1}},
    {211, {6.3065685, 0.60197032, 1.1963313}}
  }
};

vector<std::map<int, vector<double>>> PhysdEdx::pdg_plane_gsigma_map_data = {
  // Plane 0
  {
    {13, {-1, -1, -1}},
    {2212, {-1, -1, -1}},
    {211, {1.333518, 0.70926261, 0.11899788}}
  },
  // Plane 1
  {
    {13, {-1, -1, -1}},
    {2212, {-1, -1, -1}},
    {211, {1.5296072, 0.82589198, 0.21487886}}
  },
  // Plane 2
  {
    {13, {-1, -1, -1}},
    {2212, {-1, -1, -1}},
    {211, {2.6312587, 1.3023461, 0.074483758}}
  }
};

vector<std::map<int, vector<double>>> PhysdEdx::pdg_plane_width_map_data = {
  // Plane 0
  {
    {13, {-1, -1, -1}},
    {2212, {-1, -1, -1}},
    {211, {375.81207, 5.3262913, 0.12765877}}
  },
  // Plane 1
  {
    {13, {-1, -1, -1}},
    {2212, {-1, -1, -1}},
    {211, {0.15128567, 0.77813038, 0.13661804}}
  },
  // Plane 2
  {
    {13, {-1, -1, -1}},
    {2212, {-1, -1, -1}},
    {211, {0.61247895, 0.058084994, -0.38261576}}
  }
};


vector<std::map<int, vector<double>>> PhysdEdx::pdg_plane_low_rr_2p5 = {
  // Plane 0
  {
    {13, {-1, -1, -1}},
    {2212, {-1, -1, -1}},
    {211, {4.4694652, 0.8211343, 0.18246969}}
  },
  // Plane 1
  {
    {13, {-1, -1, -1}},
    {2212, {-1, -1, -1}},
    {211, {4.4521248, 0.86622439, 0.2265413}}
  },
  // Plane 2
  {
    {13, {-1, -1, -1}},
    {2212, {-1, -1, -1}},
    {211, {4.6600885, 0.79475692, 0.18291234}}
  }
};


vector<std::map<int, vector<double>>> PhysdEdx::pdg_plane_low_rr_2p5_data = {
  // Plane 0
  {
    {13, {-1, -1, -1}},
    {2212, {-1, -1, -1}},
    {211, {4.6324529, 0.64602767, 0.42133296}}
  },
  // Plane 1
  {
    {13, {-1, -1, -1}},
    {2212, {-1, -1, -1}},
    {211, {4.5292116, 0.80924871, 0.36116977}}
  },
  // Plane 2
  {
    {13, {-1, -1, -1}},
    {2212, {-1, -1, -1}},
    {211, {4.6965529, 0.78220324, 0.24326973}}
  }
};


vector<std::map<int, vector<double>>> PhysdEdx::pdg_plane_low_rr_3p5 = {
  // Plane 0
  {
    {13, {-1, -1, -1}},
    {2212, {-1, -1, -1}},
    {211, {4.0185155, 0.63042954, 0.14723379}}
  },
  // Plane 1
  {
    {13, {-1, -1, -1}},
    {2212, {-1, -1, -1}},
    {211, {4.0001484, 0.67703628, 0.18890132}}
  },
  // Plane 2
  {
    {13, {-1, -1, -1}},
    {2212, {-1, -1, -1}},
    {211, {4.1330978, 0.5850288, 0.14605994}}
  }
};


vector<std::map<int, vector<double>>> PhysdEdx::pdg_plane_low_rr_3p5_data = {
  // Plane 0
  {
    {13, {-1, -1, -1}},
    {2212, {-1, -1, -1}},
    {211, {4.1286335, 0.69667365, 0.19698588}}
  },
  // Plane 1
  {
    {13, {-1, -1, -1}},
    {2212, {-1, -1, -1}},
    {211, {4.0934128, 0.96187172, 0.08740676}}
  },
  // Plane 2
  {
    {13, {-1, -1, -1}},
    {2212, {-1, -1, -1}},
    {211, {4.1313592, 0.61029554, 0.15711856}}
  }
};
