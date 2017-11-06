//This program solves the Lippman-Schwinger Equation 
#include <iostream>
#include <string>
#include <vector>
#include <cmath>



struct MeshPoint{
	double k;//position
	double w;//weight
}

//NOTE: From IntegrationExample by Morten/Scott; Want to rewrite this
//The function GaussLegendreQuadrature() takes the lower and upper limits of 
//integration x1, x2, calculates and return the abcissas in x[0,...,n - 1] 
//and the weights in w[0,...,n - 1]
//of length n of the Gauss--Legendre n--point quadrature formulae.
 
void GaussLegendreQuadrature(double x1, double x2, double x[], double w[], int n) {
   int          m,j,i;
   double       z1,z,xm,xl,pp,p3,p2,p1;
   const double PI = 3.14159265359;
   double       *x_low, *x_high, *w_low, *w_high;

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
  GaussLegendreQuadrature(-1,1, &x[0], &w[0], n)
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

double CalculatePotential(double kp, double k){
  double potential = 0;
  const double Va = -10.463;//MeV
  const double Vb = -1650.6;//MeV
  const double Vc = 6484.3;//MeV

  const double A = 1;
  const double B = 4;
  const double C = 7;

  double kpplusk2 = pow(kp+k,2.);
  double kpminusk2 = pow(kp-k,2.);
  potential += Va/(4.*kp*k) * log((kpplusk2+pow(A,2.))/(kpminusk2+pow(A,2.)));
  potential += Vb/(4.*kp*k) * log((kpplusk2+pow(B,2.))/(kpminusk2+pow(B,2.)));
  potential += Vc/(4.*kp*k) * log((kpplusk2+pow(C,2.))/(kpminusk2+pow(C,2.)));

  return potential;
}
int main(int argc, char **argv){

    std::string USAGE("LSESolver [number of mesh points]\n");
	if (argc < 2) {
        std::cout << USAGE;
		return -1;
    }
    int n = std::strtoi(argv[2]);

    //Create mesh for integration
    std::vector<MeshPoint> mesh = SetupMesh(n);

    arma::mat potential(n+1, n+1);
    arma::mat A(n+1, n+1);

    const double K0 = 0;

    
    const double M = 938.//MeV
    double u_N_sum = 0;
    for (int i = 0; i < n+1; i++){
        u_N_sum += mesh.at(i).w*pow(K0,2.)/(pow(K0,2.)-pow(mesh.at(i).k,2.));
    }
    double u_N = -2.*PI*u_N_sum*M;
    potential(0,0) = potential(n+1,n+1) = CalculatePotential(k0,k0);
    for (int i = 0; i < n+1; i++){
        for (int j = i; j < n+1; j++){
            if (i == n && j ==n){
                potential(n, n) = CalculatePotential(K0, K0); 
            }
            else if (j == n){
                potential(i, n) = CalculatePotential(mesh.at(i).k, K0); 
            }
            else{
                potential(i,j) = CalculatePotential(mesh.at(i).k,
                                                    mesh.at(j).k);
            }
            if (i == j){
                A(i,j) += 1;
            }

            if (j == n){
                A(i,n) -= potential(i,j)*u_N;
            }
            else{
                uj = 2.*M/PI * mesh.at(j).w*pow(mesh.at(j).k,2.)/(pow(K0,2.)-pow(mesh.at(j).k,2.));
                A(i,j) -= potential(i,j)*uj
            }
        }


    }

    

}
