//This program solves the Lippman-Schwinger Equation 
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <armadillo>
//#include <limits>
#include <algorithm>



struct MeshPoint{
	double k;//position
	double w;//weight

  MeshPoint(double ik, double iw) : k(ik), w(iw) {}
};
struct LowEnergyConstants{
	double c0;
	double c2;
	double c4;
	double c4p;

  LowEnergyConstants(double ic0, double ic2, double ic4, double ic4p) : c0(ic0), c2(ic2), c4(ic4), c4p(ic4p) {}
};


//The function GaussLegendreQuadrature() takes the lower and upper limits of 
//integration x1, x2, calculates and return the abcissas in x[0,...,n - 1] 
//and the weights in w[0,...,n - 1]
//of length n of the Gauss--Legendre n--point quadrature formulae.
void GaussLegendreQuadrature(double x1, double x2, double x[], double w[], int n) {
   int          m,j,i;
   double       z1,z,xm,xl,pp,p3,p2,p1;
   const double PI = 3.14159265359;
   const double ZERO = 1e-10; double       *x_low, *x_high, *w_low, *w_high;

   m  = (n + 1)/2;                             // roots are symmetric in the interval
   xm = 0.5 * (x2 + x1);
   xl = 0.5 * (x2 - x1);

   x_low  = x;                                       // pointer initialization
   x_high = x + n - 1;
   w_low  = w;
   w_high = w + n - 1;

   for(i = 1; i <= m; i++) {                             // loops over desired roots
      z = cos(PI * (i - 0.25)/(n + 0.5));
			//
			//  Starting with the above approximation to the ith root
			//  we enter the main loop of refinement bt Newtons method.
			// 
			do {
				p1 =1.0;
				p2 =0.0;

				// loop up recurrence relation to get the
				// Legendre polynomial evaluated at x

				for(j = 1; j <= n; j++) {
					p3 = p2;
					p2 = p1;
					p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3)/j;
				}

				//
				// p1 is now the desired Legrendre polynomial. Next compute
				// ppp its derivative by standard relation involving also p2,
				// polynomial of one lower order.
				//

				pp = n * (z * p1 - p2)/(z * z - 1.0);
				z1 = z;
				z  = z1 - p1/pp;                   // Newton's method
			} while(fabs(z - z1) > ZERO);

			//
			//  Scale the root to the desired interval and put in its symmetric
			//  counterpart. Compute the weight and its symmetric counterpart
			// 
			*(x_low++)  = xm - xl * z;
			*(x_high--) = xm + xl * z;
			*w_low      = 2.0 * xl/((1.0 - z * z) * pp * pp);
			*(w_high--) = *(w_low++);
	 }
} // End_ function GaussLegendreQuadrature()

//  Set up mesh for integration using Gaussian Quadrature with Legendre
//polynomials.
//  Input: number of mesh points desired (n)
//  Output: vector of mesh points
std::vector<MeshPoint> SetupMesh(double n){
  std::vector<MeshPoint> mesh;

	std::vector<double> x;
	std::vector<double> w;

	x.resize(n);
	w.resize(n);
  GaussLegendreQuadrature(-1,1, &x[0], &w[0], n);
  //remap weights from interval [-1,1] to [0,infinity]
  //Note C can be either 1 (corresponds to k in 1/fm) or 
  //it can be 200 (hbarc = 197 MeVfm)
  double t;//cache value for use in both remappings

  const double C = 197;//determines units of k 
  const double PI = 3.14159265359;
  for (int i = 0; i < n; i++){
    t = PI/4. * (1+x[i]);
    x[i] = C*tan(t);
    w[i] = C*(PI/4.)*w[i]/pow(cos(t),2.);
    mesh.push_back(MeshPoint(x[i],w[i]));
  }

  if (mesh.size() != n){
    std::cout << "Mesh size ("<<mesh.size()<<") seems inconsistent with n ("<<n<<").\n";
  }

  return mesh;
}

double CalculatePotential(const double kp, const double k, const std::string &pot_type, LowEnergyConstants lec, double cutoff){
  double potential = 0;
  const double PI = 3.14159265359;

  if (pot_type == "nopieft"){
//    const double CUTOFF = 138;//MeV

    double regulator_p = exp(-pow(kp,4.)/pow(cutoff,4.));//for kprime
    double regulator = exp(-pow(k,4.)/pow(cutoff,4.));//for k
    
    potential += lec.c0;
    potential += lec.c2*(pow(k,2.)+pow(kp,2.));
    potential += lec.c4*(pow(k,4.)+pow(kp,4.));
    potential += lec.c4p*pow(k*kp,2.);
    
    return regulator_p * potential * regulator;
  }
  else if (pot_type == "pieft"){
    const double A = 1;
    const double MU = 138; //Pion Mass, MeV
    const double Va = -10.463;//MeV
    double regulator_p = exp(-pow(kp,4.)/pow(cutoff,4.));//for kprime
    double regulator = exp(-pow(k,4.)/pow(cutoff,4.));//for k
    potential += lec.c0;
    potential += lec.c2*(pow(k,2.)+pow(kp,2.));
    potential += lec.c4*(pow(k,4.)+pow(kp,4.));
    potential += lec.c4p*pow(k*kp,2.);
    double kpplusk2 = pow(kp+k,2.);
    double kpminusk2 = pow(kp-k,2.);
    potential += Va/(4.*MU*kp*k) * log( (kpplusk2+pow(A*MU,2.))/(kpminusk2+pow(A*MU,2.)) );
   
    return regulator_p * potential * regulator;
  }
  else{
    std::cout << "Invalid potential type! Must be EFT for this many parameters!\n";
    return sqrt(-1);
  }
}


double CalculatePotential(const double kp, const double k, const std::string &pot_type){
  double potential = 0;
  const double PI = 3.14159265359;

  if (pot_type == "yukawa"){
    const double A = 1;
    const double B = 4;
    const double C = 7;
    const double MU = 138; //Pion Mass, MeV

    const double Va = -10.463;//MeV
    const double Vb = -1650.6;//MeV
    const double Vc = 6484.3;//MeV

    double kpplusk2 = pow(kp+k,2.);
    double kpminusk2 = pow(kp-k,2.);
    potential += Va/(4.*MU*kp*k) * log( (kpplusk2+pow(A*MU,2.))/(kpminusk2+pow(A*MU,2.)) );
    potential += Vb/(4.*MU*kp*k) * log( (kpplusk2+pow(B*MU,2.))/(kpminusk2+pow(B*MU,2.)) );
    potential += Vc/(4.*MU*kp*k) * log( (kpplusk2+pow(C*MU,2.))/(kpminusk2+pow(C*MU,2.)) );
  }

  else if (pot_type == "well"){
    const double V0 = 50; // MeV
    const double R = 0.01;//1 / MeV
    if (k == kp){
      potential =V0/(2.*k*kp) * (sin((k+kp)*R)/(k+kp) - R); 
    }
    else{
      potential = V0/(2.*k*kp) * (sin((k+kp)*R)/(k+kp) - sin((k-kp)*R)/(k-kp));
    }
  }
  
  else{
    std::cout << "Invalid potential type!\n";
    return sqrt(-1);
  }

  return potential;
}

//Ensure mesh does not contain values close to k0 value, which would cause
//singular values to appear
const int VerifyMesh(double k0, const std::vector<MeshPoint> &mesh){
  const double THRESHOLD = 0.1;
  for (auto &mp : mesh){
    if(fabs(mp.k - k0) < THRESHOLD){
      std::cout << "Difference (" << fabs(mp.k-k0) << ") between mesh point (" << mp.k 
                <<") and k0 ("<<k0<<") is below threshold (" << THRESHOLD <<")!\n";
      return -1;
    }
  }
  return 0;
}

//Build potential matrix based on "pot_type", which can either be the
//parametrized Yukawa type potential (if pot_type is "yukawa"), or a square
//well (if pot_type is "well")

void BuildPotentialMatrix(arma::mat &potential, const std::vector<MeshPoint> &mesh, const double k0, const std::string &pot_type, LowEnergyConstants lec, double cutoff){
  const size_t n = mesh.size();
  for (int i = 0; i < n+1; i++){
    for (int j = i; j < n+1; j++){
      if (i == n && j ==n){
        potential(n, n) = CalculatePotential(k0, k0, pot_type, lec, cutoff); 
      }
      else if (j == n){
        potential(i, n) = CalculatePotential(mesh.at(i).k, k0, pot_type, lec, cutoff); 
      }
      else{
        potential(i,j) = CalculatePotential(mesh.at(i).k, mesh.at(j).k, pot_type, lec, cutoff);
      }
    }
  }
  potential = arma::symmatu(potential);
}

void BuildPotentialMatrix(arma::mat &potential, const std::vector<MeshPoint> &mesh, const double k0, const std::string &pot_type){
  const size_t n = mesh.size();
  for (int i = 0; i < n+1; i++){
    for (int j = i; j < n+1; j++){
      if (i == n && j ==n){
        potential(n, n) = CalculatePotential(k0, k0, pot_type); 
      }
      else if (j == n){
        potential(i, n) = CalculatePotential(mesh.at(i).k, k0, pot_type); 
      }
      else{
        potential(i,j) = CalculatePotential(mesh.at(i).k, mesh.at(j).k, pot_type);
      }
    }
  }
  potential = arma::symmatu(potential);
}

//Build A matrix defined in notes
void BuildAMatrix(arma::mat &A, const arma::mat &potential, const std::vector<MeshPoint> &mesh, const double k0){
  const double PI = 3.14159265359;
  const double M = 938.;//MeV
  const size_t n = mesh.size();
  double u_N_sum = 0;
  for (int i = 0; i < n; i++){
    u_N_sum += mesh.at(i).w*pow(k0,2.)/(pow(k0,2.)-pow(mesh.at(i).k,2.));
  }
  double u_N = -2.*M/PI*u_N_sum;

  for (int i = 0; i < n+1; i++){
    for (int j = 0; j < n+1; j++){
      if (i == j){
        A(i,j) += 1;
      }

      if (j == n){
        A(i,n) -= potential(i,j)*u_N;
      }
      else{
        double uj = 2.*M/PI * mesh.at(j).w*pow(mesh.at(j).k,2.)/(pow(k0,2.)-pow(mesh.at(j).k,2.));
        A(i,j) -= potential(i,j)*uj;
      }
    }//loop over j
  }//loop over i
}

int main(int argc, char **argv){

  std::string USAGE("LSESolver [number of mesh points] [value of E in Lab units] [well,yukawa,nopieft] [order for eft] [verbose, 0 or 1]\n");
  bool verbose;
  if (argc < 4) {
    std::cout << USAGE;
    return -1;
  }
  if (argc == 6){
    verbose = true;
  }
  if (argc == 7){
    verbose = std::stoi(argv[6]);
  }

  int n = std::stoi(argv[1]);
  const double M = 938.;//nucleon mass in MeV
  const double MPI = 938.;//pion mass in MeV 
  const double HBARC = 197.;//MeV*fm
  //Note: input E should be in  Lab Frame
  double E_lab = std::stod(argv[2]);
  //Formula taken from Scott to convert E_lab to k in center-of-mass
  double k0 = sqrt(E_lab/83.) * HBARC; //Multiplied by 197 to convert from 1/fm to MeV
  double E_cm =  E_lab/(83.*M)*pow(HBARC,2.);
  std::string pot_type(argv[3]);
  //lo,nlo,nnlo
  std::string order(argv[4]);
  double cutoff = std::stod(argv[5]);

  std::transform(pot_type.begin(), pot_type.end(), pot_type.begin(), ::tolower);


  //Create mesh for integration
  std::vector<MeshPoint> mesh = SetupMesh(n);


  //The trick for solving the principal value problem requires k0 to not be in
  //mesh.
  if (VerifyMesh(k0, mesh) != 0){
    std::cout << "Failed to verify mesh\n";
    return -1;
  }

  arma::mat potential(n+1, n+1, arma::fill::zeros);
  arma::mat A(n+1, n+1, arma::fill::zeros);

  const double PI = 3.14159265359;

  if (pot_type == "nopieft" || pot_type == "pieft"){

    //Need to calculate expected phase shift for fitting purposes 
    BuildPotentialMatrix(potential, mesh, k0, "yukawa");
    BuildAMatrix(A, potential, mesh, k0);
    arma::mat r_matrix = inv(A)*potential;
    double yukawa_phase_shift = atan(-M*k0*r_matrix(n,n))*180./PI;

    r_matrix.zeros();
    A.zeros();
    potential.zeros();

    std::vector<double> c0_vals;
    std::vector<double> c2_vals;
    std::vector<double> c4_vals;
    std::vector<double> c4p_vals;

    if (order == "lo"){
      //parameters for breakdown  = 138
      if (pot_type == "nopieft"){
        c0_vals  = {-1.48091e-05};
      }
     
     // paramaters for cutoff = 983
//    c0_vals  = {-2.08e-06};

      //parameters for breakdown = 983 pioneft
      if (pot_type == "pieft"){
        c0_vals = {-1.65e-06};
      }
      c2_vals  = {0};
      c4_vals  = {0};
      c4p_vals = {0};
    }

    if (order == "nlo"){

      //parameters for breakdown  = 983
//    c0_vals  = {2560};
//    c2_vals  = {3e-07};

      //parameters for breakdown  = 138
      if (pot_type == "nopieft"){
      c0_vals  = {1380};
      c2_vals  = {2.85e-05};
      }

      //paramets for breakdown = 69
//    c0_vals  = {550};
//    c2_vals  = {0.0001};
//    Parameaters for breakdown = 35
//    c0_vals  = {0.1};
//    c2_vals  = {0.13};
//    Parameters for breakdown = 2000
//    c0_vals  = {10000};
//    c2_vals  = {1e-07};

      //parameters for breakdown = 983 pioneft
      if (pot_type == "pieft"){
        c0_vals  = {2.4e5};
        c2_vals  = {2.75e-06};
      }
      c4_vals  = {0};
      c4p_vals = {0};
    }//nlo

    if (order == "nnlo"){
      //parameters for breakdown  = 138
      if (pot_type == "nopieft"){
        c0_vals  = {6e-08};
        c2_vals  = {1.2e-06};
        c4_vals  = {6000};
        c4p_vals = {0.0011};
      }
      //parameters for breakdown = 983
//    c0_vals  = {9e-06};
//    c2_vals  = {3e-05};
//    c4_vals  = {1e5};
//    c4p_vals = {0};

      //parameters for breakdown = 983 pioneft
      if (pot_type == "pieft"){
        c0_vals  = {5e-08};
        c2_vals  = {1.05e-06};
        c4_vals  = {3.054e6};
        c4p_vals = {0};
      }
    }

    std::vector<LowEnergyConstants> lec_search;//low energy constants
    std::vector<double> phase_shifts;

    std::vector<LowEnergyConstants> lecs;
    for (auto &c0 : c0_vals){
      for (auto &c2 : c2_vals){
        for (auto &c4 : c4_vals){
          for (auto &c4p : c4p_vals){
            lecs.push_back( LowEnergyConstants(c0, c2, c4, c4p) );
          }
        }
      }
    }

    double cur_best_eft_phase_shift = 10000;
    LowEnergyConstants best_lecs(0,0,0,0);

    for (auto &lec : lecs){
      BuildPotentialMatrix(potential, mesh, k0, pot_type, lec, cutoff);
      BuildAMatrix(A, potential, mesh, k0);
      r_matrix = inv(A)*potential;
      double phase_shift = atan(-M*k0*r_matrix(n,n))*180./PI;

      if (abs(phase_shift - yukawa_phase_shift) < abs(cur_best_eft_phase_shift-yukawa_phase_shift)){
        best_lecs.c0 = lec.c0;
        best_lecs.c2 = lec.c2;
        best_lecs.c4 = lec.c4;
        best_lecs.c4p = lec.c4p;
        cur_best_eft_phase_shift = phase_shift;
      }
      std::cout << E_lab << " "  << lec.c0 << " " << lec.c2 << " " << lec.c4 
                << " " << lec.c4p << " " << phase_shift << "\n";
      r_matrix.zeros();
      A.zeros();
      potential.zeros();
    }//loop over LECs
  }//potential == eft

  else{
    BuildPotentialMatrix(potential, mesh, k0, pot_type);
    BuildAMatrix(A, potential, mesh, k0);
    arma::mat r_matrix = inv(A)*potential;

    double phase_shift = atan(-M*k0*r_matrix(n,n))*180./PI;

    if (pot_type == "well" && phase_shift < 0){
      phase_shift += 180;
    }
    if (verbose){
      std::cout << "Phase shift: " << phase_shift << " degrees\n";
      if (pot_type == "well"){
        double V0 = 50.;//MeV
        double R = 0.01;//1/MeV
        double expected_phase_shift = (atan(sqrt(E_cm/(E_cm+V0)) * tan(R*sqrt(M*(E_cm+V0)))) - R*sqrt(M*E_cm))*180./PI;
        if (expected_phase_shift < 0){
          expected_phase_shift += 180;
        }
        std::cout << "Expected spherical well phase shift: " << expected_phase_shift << std::endl;
      }//if well
    }//if verbose
    else{
      std::cout << E_lab << " " << phase_shift << "\n";
    }//not verbose
  }//pot_type != eft
}
