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

//NOTE: From IntegrationExample by Morten/Scott; Want to rewrite this
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

//Ensure mesh does not contain values close to k0 value 
const int VerifyMesh(double k0, const std::vector<MeshPoint> &mesh){
  const double THRESHOLD = 0.1;
  for (auto &mp : mesh){
    if(abs(mp.k - k0) < THRESHOLD){
      std::cout << "Difference between mesh point (" << mp.k 
                <<") and k0 ("<<k0<<") is below threshold!\n";
      return -1;
    }
  }
  return 0;
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

  std::string USAGE("LSESolver [number of mesh points] [value of E] [well or yukawa]\n");
  if (argc < 4) {
    std::cout << USAGE;
    return -1;
  }
  int n = std::stoi(argv[1]);
  double E = std::stol(argv[2]);
  std::string pot_type(argv[3]);

  std::transform(pot_type.begin(), pot_type.end(), pot_type.begin(), ::tolower);


  //Create mesh for integration
  std::vector<MeshPoint> mesh = SetupMesh(n);

  const double M = 938.;//MeV
  double k0 = sqrt(M*E);

  if (VerifyMesh(k0, mesh) != 0){
    std::cout << "Failed to verify mesh\n";
    return -1;
  }

  arma::mat potential(n+1, n+1, arma::fill::zeros);
  arma::mat A(n+1, n+1, arma::fill::zeros);
  BuildPotentialMatrix(potential, mesh, k0, pot_type);
  BuildAMatrix(A, potential, mesh, k0);
  arma::mat r_matrix = inv(A)*potential;

  double phase_shift = atan(-M*k0*r_matrix(n,n));
  const double PI = 3.14159265359;
  std::cout << "Phase shift: " << phase_shift*180./PI << " degrees\n";

  if (pot_type == "well"){
    double V0 = 50.;//MeV
    double R = 0.01;//1/MeV
    double expected_phase_shift = atan(sqrt(E/(E+V0)) * tan(R*sqrt(M*(E+V0)))) - R*sqrt(M*E);
    std::cout << "Expected spherical well phase shift: " << expected_phase_shift*180./PI << std::endl;
  }
}
