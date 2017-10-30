//This program solves the Lippman-Schwinger Equation 
#include <iostream>
#include <string>


struct MeshPoint{
	double x;//position
	double w;//weight
}

//NOTE: From IntegrationExample by Morten/Scott; Want to rewrite this
//The function GaussLegendreQuadrature() takes the lower and upper limits of 
//integration x1, x2, calculates and return the abcissas in x[0,...,n - 1] 
//and the weights in w[0,...,n - 1]
//of length n of the Gauss--Legendre n--point quadrature formulae.
 
void GaussLegendreQuadrature(double x1, double x2, double x[], double w[], int n) {
   int         m,j,i;
   double      z1,z,xm,xl,pp,p3,p2,p1;
   double      const  pi = 3.14159265359;
   double      *x_low, *x_high, *w_low, *w_high;

   m  = (n + 1)/2;                             // roots are symmetric in the interval
   xm = 0.5 * (x2 + x1);
   xl = 0.5 * (x2 - x1);

   x_low  = x;                                       // pointer initialization
   x_high = x + n - 1;
   w_low  = w;
   w_high = w + n - 1;

   for(i = 1; i <= m; i++) {                             // loops over desired roots
      z = cos(pi * (i - 0.25)/(n + 0.5));
			//
			//  Starting with the above approximation to the ith root
			//  we enter the main loop of refinement bt Newtons method.
			// 
			do {
				p1 =1.0;
				p2 =0.0;

				//
				// loop up recurrence relation to get the
				// Legendre polynomial evaluated at x
				//

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

//Input: integration limits (a,b) and number of mesh points desired (n)
//Returns vector of mesh points
std::vector<MeshPoint> setup_mesh(double a, double b, double n){
  std::vector<MeshPoint> mesh;

	std::vector<double> x;
	std::vector<double> w;

	x.resize(n);
	w.resize(n);
  GaussLegendreQuadrature(a,b, &x[0], &w[0], n)

	for (int i = 0; i < n; i++){
		
	}


  return mesh;
}

int main(int argc, char **argv){

	std::string USAGE("lsesolver [number of mesh points]\n");
	if (argc < 2) {
		std::cout << USAGE;
		return -1;
	}
  std::vector<MeshPoint> mesh = setup_mesh(-1, 1, n);
  //setup mesh points and weights


       
}
