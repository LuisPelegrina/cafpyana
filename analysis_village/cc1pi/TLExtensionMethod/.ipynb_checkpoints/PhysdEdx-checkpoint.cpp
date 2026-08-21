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

double PhysdEdx::dEdx_PDF_w_convolution(double KE, double pitch, double dEdx, int target_plane, int PDG){

  double gamma = (KE/mass)+1.0;
  double beta = TMath::Sqrt(1-(1.0/(gamma*gamma)));
  double this_xi = Landau_xi(KE, pitch);
  double this_Wmax = Get_Wmax(KE);
  double this_kappa = this_xi / this_Wmax;
  double this_dEdx_BB = meandEdx(KE);
  double par[5] = {this_kappa, beta * beta, this_xi, this_dEdx_BB, pitch};

  TF1 *PDF = new TF1("", PhysdEdx::dEdx_PDF_function, -10., 20., 5);
  PDF -> SetParameters(par[0], par[1], par[2], par[3], par[4]);

  TF1 *f_gaus = new TF1("", "gaus", -10, 10);
  double PDF_max = PDF->GetMaximumX();
  
  double sigma = pdg_plane_map[target_plane][PDG][0] +
                     pdg_plane_map[target_plane][PDG][1] * pow(PDF_max, pdg_plane_map[target_plane][PDG][2]);

  f_gaus->SetParameters(1.0, 0.0, sigma); // norm, mean, sigma

// Create the convolution object
  TF1Convolution *fconv = new TF1Convolution(PDF, f_gaus, true);
  fconv->SetRange(0, 20);
  fconv->SetNofPointsFFT(1000); // resolution of FFT

  double out = (*fconv)(&dEdx, nullptr);
  delete PDF;
  delete f_gaus;
  delete fconv;
  return out;
}

double PhysdEdx::dEdx_PDF_max_w_convolution(double KE, double pitch, double dEdx, int target_plane, int PDG){

  double gamma = (KE/mass)+1.0;
  double beta = TMath::Sqrt(1-(1.0/(gamma*gamma)));
  double this_xi = Landau_xi(KE, pitch);
  double this_Wmax = Get_Wmax(KE);
  double this_kappa = this_xi / this_Wmax;
  double this_dEdx_BB = meandEdx(KE);
  double par[5] = {this_kappa, beta * beta, this_xi, this_dEdx_BB, pitch};
    
  TF1 *PDF = new TF1("", PhysdEdx::dEdx_PDF_function, -10., 20., 5);
  PDF -> SetParameters(par[0], par[1], par[2], par[3], par[4]);

  TF1 *f_gaus = new TF1("", "gaus", -10, 10);
  double PDF_max = PDF->GetMaximumX();

  double sigma = pdg_plane_map[target_plane][PDG][0] +
                     pdg_plane_map[target_plane][PDG][1] * pow(PDF_max, pdg_plane_map[target_plane][PDG][2]);
    
  f_gaus->SetParameters(1.0, 0.0, sigma); // norm, mean, sigma

// Create the convolution object
  TF1Convolution *fconv = new TF1Convolution(PDF, f_gaus, true);
  fconv->SetRange(0, 20);
  fconv->SetNofPointsFFT(1000); // resolution of FFT

  TF1 *f = new TF1("f", fconv, 0, 20, 0);

  double out = f->GetMaximum();
  delete PDF;
  delete f;
  delete f_gaus;
  delete fconv;
  return out;
}

double PhysdEdx::dEdx_PDF_w_convolution_f1(double KE, double rr, double dEdx, double pitch, vector<TF1*> tf1_vec){
  if(rr < tf1_vec.size()) return tf1_vec[(int)rr]->Eval(dEdx);
  else return dEdx_PDF(KE, pitch, dEdx);
}

double PhysdEdx::dEdx_PDF_max_w_convolution_f1(double KE, double rr, double pitch, vector<TF1*> tf1_vec){
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


vector< std::map<int, vector<double>>> PhysdEdx::pdg_plane_map = {
  //Plane 0 functions
  {
    {13, {0.091143,0.0328459, 2.41535416}},
    {2212, {-1,-1,-1}},
    {211, {0.16513332, 0.01404078, 2.35118304}}
  },
  //Plane 1 functions
  {
      {13, {0.0842804,0.0632091,1.66645}},
      {2212, {-1,-1,-1}},
      {211, {0.16284117, 0.03529789, 1.7890453}}
  },
  //Plane 2 functions
  {
      {13, {0.074171,0.0094173,2.74065}},
      {2212, {-1,-1,-1}},
      {211, {0.08217781, 0.00918892, 2.66855711}}
  }
};