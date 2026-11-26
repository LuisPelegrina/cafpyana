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
    {13, {0.1563192, 0.0053635199, 3.2887171}},
    {2212, {-1, -1, -1}},
    {211, {0.19108573, 0.0031730093, 3.2175612}}
  },
  // Plane 1
  {
    {13, {0.18917943, 0.015541558, 2.533299}},
    {2212, {-1, -1, -1}},
    {211, {0.23898418, 0.0063814324, 2.8174429}}
  },
  // Plane 2
  {
    {13, {0.054760258, 0.0037983289, 3.5680351}},
    {2212, {-1, -1, -1}},
    {211, {0.071101581, 0.0030216455, 3.4161522}}
  }
};

vector<std::map<int, vector<double>>> PhysdEdx::pdg_plane_shift_map = {
  // Plane 0
  {
    {13, {0.018389603, -1.5478875, -1.264582, -3.7331641e-05}},
    {2212, {-1, -1, -1, -1}},
    {211, {-0.010712217, -3.0269304, -1.627109, -5.2597545e-06}}
  },
  // Plane 1
  {
    {13, {0.005025401, -1.7790517, -1.2881737, 1.4252063e-05}},
    {2212, {-1, -1, -1, -1}},
    {211, {-0.018775045, -2.472112, -1.486998, 3.0301201e-05}}
  },
  // Plane 2
  {
    {13, {0.025472829, -0.7934845, -1.0925058, -0.00021812304}},
    {2212, {-1, -1, -1, -1}},
    {211, {0.0075837051, -2.6959241, -1.6983712, -0.00014603333}}
  }
};

vector<std::map<int, vector<double>>> PhysdEdx::pdg_plane_map_data = {
  // Plane 0
  {
    {13, {0.1535849, 0.0037832005, 3.5183933}},
    {2212, {-1, -1, -1}},
    {211, {0.17565064, 0.0043624964, 2.985106}}
  },
  // Plane 1
  {
    {13, {0.19921491, 0.011694506, 2.7464344}},
    {2212, {-1, -1, -1}},
    {211, {0.23079055, 0.0069077225, 2.7732983}}
  },
  // Plane 2
  {
    {13, {0.00052705518, 0.008279583, 3.0732207}},
    {2212, {-1, -1, -1}},
    {211, {0.03397799, 0.0036560656, 3.3141642}}
  }
};

vector<std::map<int, vector<double>>> PhysdEdx::pdg_plane_shift_map_data = {
  // Plane 0
  {
    {13, {0.040501002, -1.2706342, -0.79439389, -5.8164754e-05}},
    {2212, {-1, -1, -1, -1}},
    {211, {0.01828599, -1.5077292, -0.87517218, -5.95002e-06}}
  },
  // Plane 1
  {
    {13, {0.020719559, -1.3860764, -0.82679197, 4.075841e-05}},
    {2212, {-1, -1, -1, -1}},
    {211, {0.0016511046, -1.5158901, -0.91001961, 8.3061458e-05}}
  },
  // Plane 2
  {
    {13, {0.063682755, -0.84302777, -0.63105158, -0.00019602276}},
    {2212, {-1, -1, -1, -1}},
    {211, {0.028419688, -0.90752796, -0.71735372, -0.00010630219}}
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
    {211, {1.1792513, 6.0797908, -0.57742801}}
  },
  // Plane 1
  {
    {13, {1.3999468, 5.6704702, -0.7058217}},
    {2212, {-1, -1, -1}},
    {211, {1.1457378, 5.937448, -0.56372017}}
  },
  // Plane 2
  {
    {13, {1.3857295, 6.1615277, -0.72385948}},
    {2212, {-1, -1, -1}},
    {211, {1.149363, 6.1300518, -0.56860394}}
  }
};

vector<std::map<int, vector<double>>> PhysdEdx::pdg_plane_gsigma_map = {
  // Plane 0
  {
    {13, {0.20078661, 3.688381, -1.4388254}},
    {2212, {-1, -1, -1}},
    {211, {0.23813997, 1.9347357, -1.2014479}}
  },
  // Plane 1
  {
    {13, {0.2453847, 3.8854823, -1.4750141}},
    {2212, {-1, -1, -1}},
    {211, {0.28181119, 1.88798, -1.1721308}}
  },
  // Plane 2
  {
    {13, {0.1179547, 4.0659829, -1.5163307}},
    {2212, {-1, -1, -1}},
    {211, {0.1154325, 2.2336085, -1.2115678}}
  }
};

vector<std::map<int, vector<double>>> PhysdEdx::pdg_plane_width_map = {
  // Plane 0
  {
    {13, {4.6434675e-15, 0.1720967, -0.19891404}},
    {2212, {-1, -1, -1}},
    {211, {8.8805024e-12, 0.17651874, -0.22314605}}
  },
  // Plane 1
  {
    {13, {1.4924704e-12, 0.22980875, -0.2133934}},
    {2212, {-1, -1, -1}},
    {211, {3.3666314e-12, 0.21705877, -0.22658084}}
  },
  // Plane 2
  {
    {13, {6.5186121e-20, 0.17681463, -0.22122481}},
    {2212, {-1, -1, -1}},
    {211, {5.0611498e-12, 0.20117205, -0.21128822}}
  }
};

vector<std::map<int, vector<double>>> PhysdEdx::pdg_plane_mpv_map_data = {
  // Plane 0
  {
    {13, {1.4352118, 5.6577727, -0.74919611}},
    {2212, {-1, -1, -1}},
    {211, {1.2416689, 5.793851, -0.60225503}}
  },
  // Plane 1
  {
    {13, {1.4363751, 5.3695765, -0.74037308}},
    {2212, {-1, -1, -1}},
    {211, {1.2361791, 5.7618752, -0.60114316}}
  },
  // Plane 2
  {
    {13, {1.4164432, 5.9848807, -0.76243861}},
    {2212, {-1, -1, -1}},
    {211, {1.207361, 6.0650699, -0.60536701}}
  }
};

vector<std::map<int, vector<double>>> PhysdEdx::pdg_plane_gsigma_map_data = {
  // Plane 0
  {
    {13, {0.19240686, 4.3896431, -1.5699705}},
    {2212, {-1, -1, -1}},
    {211, {0.19877704, 1.4483317, -0.99185276}}
  },
  // Plane 1
  {
    {13, {0.24281805, 2.7702142, -1.3590671}},
    {2212, {-1, -1, -1}},
    {211, {0.28881181, 4.6425465, -1.7153678}}
  },
  // Plane 2
  {
    {13, {0.081280963, 4.5588244, -1.513978}},
    {2212, {-1, -1, -1}},
    {211, {0.081001739, 3.0206604, -1.3553818}}
  }
};

vector<std::map<int, vector<double>>> PhysdEdx::pdg_plane_width_map_data = {
  // Plane 0
  {
    {13, {5.3770581e-14, 0.20448145, -0.22212097}},
    {2212, {-1, -1, -1}},
    {211, {3.5606185e-12, 0.17823138, -0.17530255}}
  },
  // Plane 1
  {
    {13, {0.088907007, 0.44777931, -0.81891915}},
    {2212, {-1, -1, -1}},
    {211, {1.8635815e-15, 0.66964202, -0.62312629}}
  },
  // Plane 2
  {
    {13, {2.8408533e-16, 0.20151147, -0.23611108}},
    {2212, {-1, -1, -1}},
    {211, {3.428302e-13, 0.23903187, -0.23196764}}
  }
};

vector<std::map<int, vector<double>>> PhysdEdx::pdg_plane_low_rr_2p5 = {
  // Plane 0
  {
    {13, {4.2257398, 0.96958629, 0.015484037}},
    {2212, {-1, -1, -1}},
    {211, {4.6541363, 0.97122055, 0.0091236178}}
  },
  // Plane 1
  {
    {13, {4.1952815, 1.0587858, 0.014189238}},
    {2212, {-1, -1, -1}},
    {211, {4.6440426, 1.0363343, 0.012575878}}
  },
  // Plane 2
  {
    {13, {4.4049141, 0.84803046, 0.012295129}},
    {2212, {-1, -1, -1}},
    {211, {4.8430167, 0.96179222, 0.0064191566}}
  }
};


vector<std::map<int, vector<double>>> PhysdEdx::pdg_plane_low_rr_2p5_data = {
  // Plane 0
  {
    {13, {4.0132275, 0.92002705, 0.0081794528}},
    {2212, {-1, -1, -1}},
    {211, {4.2487909, 0.66536962, 0.36733542}}
  },
  // Plane 1
  {
    {13, {4.0261344, 1.0182415, 0.018494527}},
    {2212, {-1, -1, -1}},
    {211, {4.296531, 0.91656055, 0.1864586}}
  },
  // Plane 2
  {
    {13, {4.2322827, 0.90949315, 0.014490555}},
    {2212, {-1, -1, -1}},
    {211, {4.6217074, 0.90453154, 0.08768914}}
  }
};


vector<std::map<int, vector<double>>> PhysdEdx::pdg_plane_low_rr_3p5 = {
  // Plane 0
  {
    {13, {3.7835452, 0.765144, 0.010325646}},
    {2212, {-1, -1, -1}},
    {211, {4.0485334, 0.71118671, 0.062844022}}
  },
  // Plane 1
  {
    {13, {3.726832, 0.83477749, 0.018212562}},
    {2212, {-1, -1, -1}},
    {211, {4.0480196, 0.79187788, 0.065620585}}
  },
  // Plane 2
  {
    {13, {3.8056138, 0.63369108, 0.049400759}},
    {2212, {-1, -1, -1}},
    {211, {4.1604983, 0.67811126, 0.058778597}}
  }
};


vector<std::map<int, vector<double>>> PhysdEdx::pdg_plane_low_rr_3p5_data = {
  // Plane 0
  {
    {13, {3.599905, 0.726092, 0.017741719}},
    {2212, {-1, -1, -1}},
    {211, {3.8711376, 0.6445285, 0.14584675}}
  },
  // Plane 1
  {
    {13, {3.3924895, 0.70248802, 0.20516684}},
    {2212, {-1, -1, -1}},
    {211, {3.8187993, 0.68023025, 0.20621984}}
  },
  // Plane 2
  {
    {13, {3.6358288, 0.65731862, 0.061828874}},
    {2212, {-1, -1, -1}},
    {211, {4.0291312, 0.63458447, 0.087469887}}
  }
};