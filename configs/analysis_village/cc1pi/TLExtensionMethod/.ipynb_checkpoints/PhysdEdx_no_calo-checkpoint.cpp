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
    {13, {0.14239127, 0.007574283, 3.0387078}},
    {2212, {-1, -1, -1}},
    {211, {0.19342207, 0.0026053117, 3.3601791}}
  },
  // Plane 1
  {
    {13, {0.19833565, 0.0098315159, 2.8896978}},
    {2212, {-1, -1, -1}},
    {211, {0.2444265, 0.0046725243, 3.0363741}}
  },
  // Plane 2
  {
    {13, {0.046899734, 0.0051396115, 3.3410388}},
    {2212, {-1, -1, -1}},
    {211, {0.072227891, 0.0030891832, 3.4204308}}
  }
};

vector<std::map<int, vector<double>>> PhysdEdx::pdg_plane_shift_map = {
  // Plane 0
  {
    {13, {0.014364484, -0.39140505, 6.8488855}},
    {2212, {-1, -1, -1}},
    {211, {-0.016998084, -0.93561805, 3.6427571}}
  },
  // Plane 1
  {
    {13, {0.0064494297, -0.47022613, 6.4193064}},
    {2212, {-1, -1, -1}},
    {211, {-0.022293774, -0.87420722, 3.8857692}}
  },
  // Plane 2
  {
    {13, {0.0056085459, -0.23560982, 6.7604309}},
    {2212, {-1, -1, -1}},
    {211, {0.0015916355, -0.84947135, 3.4513382}}
  }
};

vector<std::map<int, vector<double>>> PhysdEdx::pdg_plane_map_data = {
  // Plane 0
  {
    {13, {0.1431099, 0.005111082, 3.293785}},
    {2212, {-1, -1, -1}},
    {211, {0.17006244, 0.0048448422, 2.9276276}}
  },
  // Plane 1
  {
    {13, {0.21409109, 0.0059681254, 3.2682334}},
    {2212, {-1, -1, -1}},
    {211, {0.25001509, 0.0032153943, 3.2902294}}
  },
  // Plane 2
  {
    {13, {2.3283896e-18, 0.0080824503, 3.1037726}},
    {2212, {-1, -1, -1}},
    {211, {0.03632486, 0.003460479, 3.367582}}
  }
};

vector<std::map<int, vector<double>>> PhysdEdx::pdg_plane_shift_map_data = {
  // Plane 0
  {
    {13, {-0.0084319398, -0.4555194, 10.058967}},
    {2212, {-1, -1, -1}},
    {211, {-0.057187877, -0.77567604, 5.2097971}}
  },
  // Plane 1
  {
    {13, {-0.012287711, -0.47588277, 10.737676}},
    {2212, {-1, -1, -1}},
    {211, {-0.061211541, -0.81367363, 4.9102044}}
  },
  // Plane 2
  {
    {13, {-0.0089006432, -0.32902245, 11.700444}},
    {2212, {-1, -1, -1}},
    {211, {-0.048480585, -0.62156304, 5.2789952}}
  }
};



/*
Langau fits
*/

vector<std::map<int, vector<double>>> PhysdEdx::pdg_plane_mpv_map = {
  // Plane 0
  {
    {13, {1.4156122, 5.9500048, -0.72348789}},
    {2212, {-1, -1, -1}},
    {211, {1.1736011, 6.0258736, -0.57156766}}
  },
  // Plane 1
  {
    {13, {1.3999468, 5.6704702, -0.7058217}},
    {2212, {-1, -1, -1}},
    {211, {1.136704, 5.8294632, -0.55197804}}
  },
  // Plane 2
  {
    {13, {1.3857295, 6.1615277, -0.72385948}},
    {2212, {-1, -1, -1}},
    {211, {1.1178919, 5.9916339, -0.54968754}}
  }
};

vector<std::map<int, vector<double>>> PhysdEdx::pdg_plane_gsigma_map = {
  // Plane 0
  {
    {13, {0.20078661, 3.688381, -1.4388254}},
    {2212, {-1, -1, -1}},
    {211, {0.2251219, 1.9870644, -1.1977204}}
  },
  // Plane 1
  {
    {13, {0.2453847, 3.8854823, -1.4750141}},
    {2212, {-1, -1, -1}},
    {211, {0.25873027, 1.6295292, -1.0735347}}
  },
  // Plane 2
  {
    {13, {0.1179547, 4.0659829, -1.5163307}},
    {2212, {-1, -1, -1}},
    {211, {0.10745817, 2.3081002, -1.2135255}}
  }
};

vector<std::map<int, vector<double>>> PhysdEdx::pdg_plane_width_map = {
  // Plane 0
  {
    {13, {4.6434675e-15, 0.1720967, -0.19891404}},
    {2212, {-1, -1, -1}},
    {211, {8.0565571e-13, 0.17584904, -0.18717068}}
  },
  // Plane 1
  {
    {13, {1.4924704e-12, 0.22980875, -0.2133934}},
    {2212, {-1, -1, -1}},
    {211, {8.1544952e-12, 0.22031106, -0.1932778}}
  },
  // Plane 2
  {
    {13, {6.5186121e-20, 0.17681463, -0.22122481}},
    {2212, {-1, -1, -1}},
    {211, {6.8510805e-20, 0.20742332, -0.19934856}}
  }
};

vector<std::map<int, vector<double>>> PhysdEdx::pdg_plane_mpv_map_data = {
  // Plane 0
  {
    {13, {1.4352118, 5.6577727, -0.74919611}},
    {2212, {-1, -1, -1}},
    {211, {1.2589061, 5.9093785, -0.61336339}}
  },
  // Plane 1
  {
    {13, {1.4363751, 5.3695765, -0.74037308}},
    {2212, {-1, -1, -1}},
    {211, {1.1633165, 5.5418691, -0.55931939}}
  },
  // Plane 2
  {
    {13, {1.4164432, 5.9848807, -0.76243861}},
    {2212, {-1, -1, -1}},
    {211, {1.1874204, 5.9029821, -0.58732478}}
  }
};

vector<std::map<int, vector<double>>> PhysdEdx::pdg_plane_gsigma_map_data = {
  // Plane 0
  {
    {13, {0.19240686, 4.3896431, -1.5699705}},
    {2212, {-1, -1, -1}},
    {211, {0.19283491, 1.8886534, -1.1190206}}
  },
  // Plane 1
  {
    {13, {0.24281805, 2.7702142, -1.3590671}},
    {2212, {-1, -1, -1}},
    {211, {0.25233707, 2.2740085, -1.2776945}}
  },
  // Plane 2
  {
    {13, {0.081280963, 4.5588244, -1.513978}},
    {2212, {-1, -1, -1}},
    {211, {0.078152964, 3.7170766, -1.462581}}
  }
};

vector<std::map<int, vector<double>>> PhysdEdx::pdg_plane_width_map_data = {
  // Plane 0
  {
    {13, {5.3770581e-14, 0.20448145, -0.22212097}},
    {2212, {-1, -1, -1}},
    {211, {5.5658115e-13, 0.18493867, -0.15854188}}
  },
  // Plane 1
  {
    {13, {0.088907007, 0.44777931, -0.81891915}},
    {2212, {-1, -1, -1}},
    {211, {9.736046e-12, 0.25730013, -0.18880002}}
  },
  // Plane 2
  {
    {13, {2.8408533e-16, 0.20151147, -0.23611108}},
    {2212, {-1, -1, -1}},
    {211, {2.3469286e-14, 0.25362471, -0.23129325}}
  }
};



vector<std::map<int, vector<double>>> PhysdEdx::pdg_plane_low_rr_2p5 = {
  // Plane 0
  {
    {13, {4.2257398, 0.96958629, 0.015484037}},
    {2212, {-1, -1, -1}},
    {211, {4.615885, 0.99091159, 0.010629962}}
  },
  // Plane 1
  {
    {13, {4.1952815, 1.0587858, 0.014189238}},
    {2212, {-1, -1, -1}},
    {211, {4.5402558, 1.0324138, 0.043352271}}
  },
  // Plane 2
  {
    {13, {4.4049141, 0.84803046, 0.012295129}},
    {2212, {-1, -1, -1}},
    {211, {4.7558236, 1.0080562, 0.014716591}}
  }
};


vector<std::map<int, vector<double>>> PhysdEdx::pdg_plane_low_rr_2p5_data = {
  // Plane 0
  {
    {13, {4.0132275, 0.92002705, 0.0081794528}},
    {2212, {-1, -1, -1}},
    {211, {4.2391205, 0.71906254, 0.34027046}}
  },
  // Plane 1
  {
    {13, {4.0261344, 1.0182415, 0.018494527}},
    {2212, {-1, -1, -1}},
    {211, {4.2962684, 0.95673289, 0.18302125}}
  },
  // Plane 2
  {
    {13, {4.2322827, 0.90949315, 0.014490555}},
    {2212, {-1, -1, -1}},
    {211, {4.7228627, 0.98399852, 0.015731804}}
  }
};


vector<std::map<int, vector<double>>> PhysdEdx::pdg_plane_low_rr_3p5 = {
  // Plane 0
  {
    {13, {3.7835452, 0.765144, 0.010325646}},
    {2212, {-1, -1, -1}},
    {211, {4.0315088, 0.7140519, 0.066182227}}
  },
  // Plane 1
  {
    {13, {3.726832, 0.83477749, 0.018212562}},
    {2212, {-1, -1, -1}},
    {211, {3.9848542, 0.76736432, 0.10720551}}
  },
  // Plane 2
  {
    {13, {3.8056138, 0.63369108, 0.049400759}},
    {2212, {-1, -1, -1}},
    {211, {4.1694973, 0.71645858, 0.037324229}}
  }
};


vector<std::map<int, vector<double>>> PhysdEdx::pdg_plane_low_rr_3p5_data = {
  // Plane 0
  {
    {13, {3.599905, 0.726092, 0.017741719}},
    {2212, {-1, -1, -1}},
    {211, {3.8513415, 0.63485499, 0.15212235}}
  },
  // Plane 1
  {
    {13, {3.3924895, 0.70248802, 0.20516684}},
    {2212, {-1, -1, -1}},
    {211, {3.8257499, 0.68359835, 0.20339456}}
  },
  // Plane 2
  {
    {13, {3.6358288, 0.65731862, 0.061828874}},
    {2212, {-1, -1, -1}},
    {211, {4.0076445, 0.65047411, 0.081336994}}
  }
};
