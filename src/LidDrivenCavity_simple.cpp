//============================================================================
// Name        : LidDrivenCavity_simple.cpp
// Author      : Karla Kauffmann
// Version     :
// Copyright   : 
// Description : LidDrivenCavity
//============================================================================

#include <iostream>
#include <fstream>
#include <ctime>
#include <cmath>
#include <chrono>
// This library is used to maximize the output precision
#include <limits>

using namespace std;
using namespace chrono;

int main (void)
{
   /* Display numbers with maximum precision possible */
	cout.precision(numeric_limits<double>::digits10 + 1);


   /* Variable for tracking execution time */
	high_resolution_clock::time_point t1 = high_resolution_clock::now();

	const int Nx = 128;    			/* number of cells for the CV computation in the x and y direction */
	const int Ny = Nx;				/* not including ghost cells */

	double Lx = 1.0; 			/* length of the cavity in the x and y direction */
	double Ly = 1.0;
	double p_fix = 0.0; 		/* fixed lower west corner pressure */
	double density = 1.0; 		/* density of the fluid */
	double Reynolds_num = 100.0; /* Reynolds number */
	double Gamma_constant = 1/Reynolds_num; /* constant diffusion coefficient */

	double urfu = 0.5; 			/* under-relaxation factors for velocity u, v and pressure. Velocity factors cannot be 0.0 									since they go in the denominator of a_p in the momentum equations */
	double urfv = 0.5;
	double urfp = 0.5;

	double residual_threshold = 1.0E-5; 	/* the value of the residuals to end the SIMPLE iterations */
	double residual_threshold_GS = 1.0E-4; 	/* the value of the residuals to end the Gauss-Seidel iterations */
	int i_iter = 1; 			/* current number of SIMPLE iteration */
	int i_iter_GS = 1; 			/* current number of Gauss-Seidel iteration */
	int MAX_ITER = 10000; 		/* set the maximum number of iterations to store in the residual vector */
	int MAX_ITER_GS = 1; 	/* set the maximum number of iterations allowed for the Gauss-Seidel iteration */
	double first_iter=1/1.0E-7;
	double max_residual = 1.0;
	double max_residual_GS= 1.0;
	double residual_GS;
	double residual_simple;

	double u_north = 1.0; 		/* Velocities at boundaries */
	double u_south = 0.0;
	double v_west = 0.0;
	double v_east = 0.0;

	double deltax = Lx/Nx;		/* size of cell */
	double deltay = Ly/Ny;

	double F_w_u, F_e_u, F_s_u, F_n_u; /* coefficients for the momentum equations of u and v  */
	double D_w_u, D_e_u, D_s_u, D_n_u;
	double a_w_u[Nx+1][Ny+2], a_e_u[Nx+1][Ny+2], a_s_u[Nx+1][Ny+2], a_n_u[Nx+1][Ny+2];
	double S_p_u, S_u;
	double a_p_u[Nx+1][Ny+2];
	double delta_F_u;

	double F_w_v, F_e_v, F_s_v, F_n_v;
	double D_w_v, D_e_v, D_s_v, D_n_v;
	double a_w_v[Nx+2][Ny+1], a_e_v[Nx+2][Ny+1], a_s_v[Nx+2][Ny+1], a_n_v[Nx+2][Ny+1];
	double S_p_v, S_v;
	double a_p_v[Nx+2][Ny+1];
	double delta_F_v;

	double a_w_p[Nx+2][Ny+2], a_e_p[Nx+2][Ny+2], a_s_p[Nx+2][Ny+2], a_n_p[Nx+2][Ny+2], a_p_p[Nx+2][Ny+2], b_p[Nx+2][Ny+2];

	//============================================================================
	//
	// Initialization
	//
	//============================================================================

	double u[Nx+1][Ny+2] = {0};		/* velocity in the x direction */
	double u_old[Nx+1][Ny+2] = {0};
	double u_star[Nx+1][Ny+2] = {0};
	double u_prime[Nx+1][Ny+2] = {0};
	double uc[Nx+1][Ny+2] = {0};

	double v[Nx+2][Ny+1] = {0};		/* velocity in the y direction */
	double v_old[Nx+2][Ny+1] = {0};
	double v_star[Nx+2][Ny+1] = {0};
	double v_prime[Nx+2][Ny+1] = {0};
	double vc[Nx+2][Ny+1] = {0};

	double p[Nx+2][Ny+2] = {0};		/* pressure */
	double p_old[Nx+2][Ny+2] = {0};
	double p_prime[Nx+2][Ny+2] = {0};
	double pc[Nx+2][Ny+2] = {0};

	double du[Nx+1][Ny+2] = {0};		/* parameter d for the pressure correction equation for u and v velocities */
	double dv[Nx+2][Ny+1] = {0};

	//============================================================================
	//
	// Start of SIMPLE iteration
	//
	//============================================================================

//    while ((max_residual > residual_threshold) && (i_iter < (MAX_ITER + 1))) {
	while (i_iter < (MAX_ITER + 1)) {

		// Restarts maximum SIMPLE residual to zero at each iteration
		max_residual = 0.0;

    	//============================================================================
    	//
    	// Boundary conditions for ghost cells. For ghost cells that are not exactly
    	// at the boundary we use linear interpolation.
    	//
    	//============================================================================

    	for (int i=0; i < (Nx+1); i++) {
    		u_old[i][0] = 2.0*u_south - u_old[i][1];
    		u_old[i][Ny+1] = 2.0*u_north - u_old[i][Ny];
    	}

    	for (int j=0; j < (Ny+1); j++) {
    		v_old[0][j] = 2.0*v_west - v_old[1][j];
    		v_old[Nx+1][j] = 2.0*v_east - v_old[Nx][j];
    	}

    	//============================================================================
    	//
    	// We obtain the guessed velocities (u_star*) by solving the system of
    	// momentum equations.
    	//
    	//============================================================================

    	// We calculate the coefficients for the momentum equation
    	for(int i=1;i<Nx;i++){
    		for(int j=1;j<(Ny+1);j++){
				/* Equations for u_velocity for north (n), south (s), east (e) and west (w) */

				F_w_u = 0.5 * density * ( u_old[i][j] + u_old[i-1][j]);
				F_e_u = 0.5 * density * ( u_old[i+1][j] + u_old[i][j]);
				F_s_u = 0.5 * density * (v_old[i+1][j-1] + v_old[i][j-1]);
				F_n_u = 0.5 * density * (v_old[i][j] + v_old[i+1][j]);
				D_w_u = Gamma_constant/deltax;
				D_e_u = Gamma_constant/deltax;
				D_s_u = Gamma_constant/deltay;
				D_n_u = Gamma_constant/deltay;
				a_w_u[i][j] = D_w_u + 0.5 * F_w_u;
				a_e_u[i][j] = D_e_u - 0.5 * F_e_u;
				a_s_u[i][j] = D_s_u + 0.5 * F_s_u;
				a_n_u[i][j] = D_n_u - 0.5 * F_n_u;
				S_p_u = 0.0;
				S_u = 0.0;
				delta_F_u = F_e_u + F_n_u - F_s_u - F_w_u;
				a_p_u[i][j] = a_w_u[i][j] + a_e_u[i][j] + a_s_u[i][j] + a_n_u[i][j] + delta_F_u - S_p_u;
				a_p_u[i][j] = a_p_u[i][j]/urfu; /* underrelaxation */
 				du[i][j] = deltax*deltay/a_p_u[i][j];
    		}
    	}


    	// u_star = u_old
    	for(int i=0;i<(Nx+1);i++){
    		for(int j=0;j<(Ny+2);j++){
    			u_star[i][j] = u_old[i][j];
    		}
    	}

    	// Reset values for GS iteration
    	max_residual_GS= 1.0;
    	i_iter_GS = 1;

    	// solve u_star using Gauss-Seidel
    	while ((max_residual_GS > residual_threshold_GS) && (i_iter_GS < (MAX_ITER_GS + 1))) {
        	// Restarts maximum residual to zero at each iteration
        	max_residual_GS = 0.0;

			for(int i=1;i<Nx;i++){

	        	residual_GS = 0.0;

				for(int j=1;j<(Ny+1);j++){

					residual_GS = abs(a_p_u[i][j]*u_star[i][j] - a_n_u[i][j]*u_star[i][j+1] - a_s_u[i][j]*u_star[i][j-1] - a_e_u[i][j] * u_star[i+1][j] - a_w_u[i][j] * u_star[i-1][j] - ( p_old[i-1][j] - p_old[i][j] ) * deltax * deltay -  S_u - (1 - urfu) * a_p_u[i][j] * u_old[i][j]);

					u_star[i][j] = (a_s_u[i][j] * u_star[i][j-1] + a_n_u[i][j] * u_star[i][j+1] + a_e_u[i][j] * u_star[i+1][j] + a_w_u[i][j] * u_star[i-1][j] - ( p_old[i+1][j] - p_old[i][j] ) * deltax * deltay +  S_u + (1 - urfu) * a_p_u[i][j] * u_old[i][j])/a_p_u[i][j];

					// This variable stores the maximum residual of the vector residuals[i], by comparing each element
			        // with itself. If it is higher, it stores the value.
					max_residual_GS = max(residual_GS,max_residual_GS);
				}
			}

	    	if (i_iter_GS == 1){
	    		// We save the first iteration value of the residual
	    		//to use as the normalization factor for the the next residuals.
	    		first_iter = max_residual_GS;

	    		// To avoid dividing by zero
	    		if (first_iter == 0.0){
	    			first_iter = 1/residual_threshold_GS;
	    		}
	    	}
	    	// Value of the residual at each iteration normalized to the value
	    	// of the first residual
	    	max_residual_GS = max_residual_GS / first_iter;

//	    	// We print the iteration number and the values of the solution in each iteration
//	    	cout << "GS iteration for u " << i_iter_GS << ": error " << max_residual_GS << endl;

			i_iter_GS++;
//	     	if (i_iter_GS == (MAX_ITER_GS + 1)) {
//	      		cout << "You have reached the maximum number of Gauss-Seidel iterations before converging to the set threshold" << endl;
//	     	}
    	}

    	//============================================================================
    	//
    	// We obtain the guessed velocities (v_star*) by solving the system of
    	// momentum equations.
    	//
    	//============================================================================

    	// We calculate the coefficients for the momentum equation
    	for(int i=1;i<(Nx+1);i++){
    		for(int j=1;j<Ny;j++){
				/* Equations for u_velocity for north (n), south (s), east (e) and west (w) */

				F_w_v = 0.5 * density * (u_old[i-1][j] + u_old[i-1][j+1]);
				F_e_v = 0.5 * density * (u_old[i][j] + u_old[i][j+1]);
				F_s_v = 0.5 * density * (v_old[i][j-1] + v_old[i][j]);
				F_n_v = 0.5 * density * (v_old[i][j] + v_old[i][j+1]);
				D_w_v = Gamma_constant/deltax;
				D_e_v = Gamma_constant/deltax;
				D_s_v = Gamma_constant/deltay;
				D_n_v = Gamma_constant/deltay;
				a_w_v[i][j] = D_w_v + 0.5 * F_w_v;
				a_e_v[i][j] = D_e_v - 0.5 * F_e_v;
				a_s_v[i][j] = D_s_v + 0.5 * F_s_v;
				a_n_v[i][j] = D_n_v - 0.5 * F_n_v;
				S_p_v = 0.0;
				S_v = 0.0;
				delta_F_v = F_e_v + F_n_v - F_s_v - F_w_v;
				a_p_v[i][j] = a_w_v[i][j] + a_e_v[i][j] + a_s_v[i][j] + a_n_v[i][j] + delta_F_v - S_p_v;
				a_p_v[i][j] = a_p_v[i][j]/urfv; /* underrelaxation */
 				dv[i][j] = deltax*deltay/a_p_v[i][j];

    		}
    	}


    	// v_star = v_old
    	for(int i=0;i<(Nx+2);i++){
    		for(int j=0;j<(Ny+1);j++){
    			v_star[i][j] = v_old[i][j];
    		}
    	}

       	// Reset values for GS iteration
    	max_residual_GS= 1.0;
    	i_iter_GS = 1;

    	// solve v_star using Gauss-Seidel
    	while ((max_residual_GS > residual_threshold_GS) && (i_iter_GS < (MAX_ITER_GS + 1))) {
        	// Restarts maximum residual to zero at each iteration
        	max_residual_GS = 0.0;

			for(int i=1;i<(Nx+1);i++){

	        	residual_GS = 0.0;

				for(int j=1;j<Ny;j++){

					residual_GS = abs(a_p_v[i][j]*v_star[i][j] - a_n_v[i][j]*v_star[i][j+1] - a_s_v[i][j]*v_star[i][j-1] - a_e_v[i][j] * v_star[i+1][j] - a_w_v[i][j] * v_star[i-1][j] - ( p_old[i][j-1] - p_old[i][j] ) * deltax * deltay -  S_v - (1 - urfv) * a_p_v[i][j] * v_old[i][j]);

					v_star[i][j] = (a_s_v[i][j] * v_star[i][j-1] + a_n_v[i][j] * v_star[i][j+1] + a_e_v[i][j] * v_star[i+1][j] + a_w_v[i][j] * v_star[i-1][j] - ( p_old[i][j+1] - p_old[i][j] ) * deltax * deltay +  S_v + (1 - urfv) * a_p_v[i][j] * v_old[i][j])/a_p_v[i][j];

					// This variable stores the maximum residual of the vector residuals[i], by comparing each element
			        // with itself. If it is higher, it stores the value.
					max_residual_GS = max(residual_GS,max_residual_GS);
				}
			}

	    	if (i_iter_GS == 1){
	    		// We save the first iteration value of the residual
	    		//to use as the normalization factor for the the next residuals.
	    		first_iter = max_residual_GS;

	    		// To avoid dividing by zero
	    		if (first_iter == 0.0){
	    			first_iter = 1/residual_threshold_GS;
	    		}
	    	}
	    	// Value of the residual at each iteration normalized to the value
	    	// of the first residual
	    	max_residual_GS = max_residual_GS / first_iter;

//	    	// We print the iteration number and the values of the solution in each iteration
//	    	cout << "GS iteration for v " << i_iter_GS << ": error " << max_residual_GS << endl;

			i_iter_GS++;
//	     	if (i_iter_GS == (MAX_ITER_GS + 1)) {
//	      		cout << "You have reached the maximum number of Gauss-Seidel iterations before converging to the set threshold" << endl;
//	     	}
    	}

    	//============================================================================
    	//
    	// We obtain the pressure correction (p_prime') by solving the system of
    	// equations.
    	//
    	//============================================================================

    	for(int i=1;i<(Nx+1);i++){
    		for(int j=1;j<(Ny+1);j++){
			   a_w_p[i][j] = density * du[i-1][j] * deltax*deltay;
			   a_e_p[i][j] = density * du[i][j] * deltax*deltay;
			   a_s_p[i][j] = density * dv[i][j-1] * deltax*deltay;
			   a_n_p[i][j] = density * dv[i][j] * deltax*deltay;
			   a_p_p[i][j] = a_w_p[i][j] + a_e_p[i][j] + a_s_p[i][j] + a_n_p[i][j];
			   b_p[i][j] = -(density * u_star[i][j] * deltax*deltay) + (density * u_star[(i-1)][j] * deltax*deltay) - (density * v_star[i][j] * deltax*deltay) + (density * v_star[i][j-1] * deltax*deltay);

    		}
    	}

    	// Setting the reference value of the pressure on the lower west corner
		a_w_p[1][1] = 0.0;
		a_e_p[1][1] = 0.0;
		a_s_p[1][1] = 0.0;
		a_n_p[1][1] = 0.0;
		a_p_p[1][1] = 1.0;
		b_p[1][1] = p_fix;

		// Initialize p_prime' to zero
		for(int i=0;i<(Nx+2);i++){
			for(int j=0;j<(Ny+2);j++){
				p_prime[i][j] = 0.0;
			}
		}

       	// Reset values for GS iteration
    	max_residual_GS= 1.0;
    	i_iter_GS = 1;

    	// solve v_star using Gauss-Seidel
    	while ((max_residual_GS > residual_threshold_GS) && (i_iter_GS < (MAX_ITER_GS + 1))) {
        	// Restarts maximum residual to zero at each iteration
        	max_residual_GS = 0.0;

			for(int i=1;i<(Nx+1);i++){

	        	residual_GS = 0.0;

				for(int j=1;j<(Ny+1);j++){

					residual_GS = abs(a_p_p[i][j]*p_prime[i][j] - a_n_p[i][j]*p_prime[i][j+1] - a_s_p[i][j]*p_prime[i][j-1] - a_e_p[i][j] * p_prime[i+1][j] - a_w_p[i][j] * p_prime[i-1][j] - b_p[i][j]);

					p_prime[i][j] = (a_s_p[i][j] * p_prime[i][j-1] + a_n_p[i][j] * p_prime[i][j+1] + a_e_p[i][j] * p_prime[i+1][j] + a_w_p[i][j] * p_prime[i-1][j] + b_p[i][j])/a_p_p[i][j];

					// This variable stores the maximum residual of the vector residuals[i], by comparing each element
			        // with itself. If it is higher, it stores the value.
					max_residual_GS = max(residual_GS,max_residual_GS);
				}
			}

	    	if (i_iter_GS == 1){
	    		// We save the first iteration value of the residual
	    		//to use as the normalization factor for the the next residuals.
	    		first_iter = max_residual_GS;

	    		// To avoid dividing by zero
	    		if (first_iter == 0.0){
	    			first_iter = 1/residual_threshold_GS;
	    		}
	    	}
	    	// Value of the residual at each iteration normalized to the value
	    	// of the first residual
	    	max_residual_GS = max_residual_GS / first_iter;

//	    	// We print the iteration number and the values of the solution in each iteration
//	    	cout << "GS iteration for p " << i_iter_GS << ": error " << max_residual_GS << endl;

			i_iter_GS++;
//	     	if (i_iter_GS == (MAX_ITER_GS + 1)) {
//	      		cout << "You have reached the maximum number of Gauss-Seidel iterations before converging to the set threshold" << endl;
//	     	}
    	}

    	//============================================================================
		//
		// Obtain u and v velocity corrections: u_prime' and v_prime'
		//
		//============================================================================

       	for(int i=1;i<Nx;i++){
        		for(int j=1;j<(Ny+1);j++){
        			u_prime[i][j] = du[i][j]*(p_prime[i][j]-p_prime[i+1][j]);
        		}
       	}

       	for(int i=1;i<(Nx+1);i++){
        		for(int j=1;j<Ny;j++){
        			v_prime[i][j] = dv[i][j]*(p_prime[i][j]-p_prime[i][j+1]);
        		}
       	}

    	//============================================================================
		//
		// Correct u, v and p.
		//
		//============================================================================

       	for(int i=1;i<(Nx+1);i++){
       		for(int j=1;j<(Ny+1);j++){
       			p[i][j] = p_old[i][j] + urfp*p_prime[i][j];
       		}
       	}

       	for(int i=1;i<Nx;i++){
       		for(int j=1;j<(Ny+1);j++){
       			u[i][j] = u_star[i][j] + u_prime[i][j];
       		}
       	}

       	for(int i=1;i<(Nx+1);i++){
       		for(int j=1;j<Ny;j++){
       			v[i][j] = v_star[i][j] + v_prime[i][j];
       		}
       	}

//       	for(int i=1;i<(Nx+1);i++){
//       		residual_simple = 0.0;
//       		for(int j=1;j<Ny;j++){
//       			residual_simple = abs(-(density * u_star[i][j] * deltax*deltay) + (density * u_star[(i-1)][j] * deltax*deltay) - (density * v_star[i][j] * deltax*deltay) + (density * v_star[i][j-1] * deltax*deltay));
//
//
//
//				// This variable stores the maximum residual of the vector residuals[i], by comparing each element
//		        // with itself. If it is higher, it stores the value.
//				max_residual = max(residual_simple,max_residual);
//       		}
//       	}
//
//    	if (i_iter == 1){
//    		// We save the first iteration value of the residual
//    		//to use as the normalization factor for the the next residuals.
//    		first_iter = max_residual;
//
//    		// To avoid dividing by zero
//    		if (first_iter == 0.0){
//    			first_iter = 1/residual_threshold;
//    		}
//    	}
//    	// Value of the residual at each iteration normalized to the value
//    	// of the first residual
//    	max_residual = max_residual / first_iter;

       	//============================================================================
		//
		// Write to screen
		//
		//============================================================================

    	// We print the iteration number and the values of the solution in each iteration
    	if((i_iter)%100 == 0){
    		cout << "SIMPLE iteration " << i_iter << ": error " << max_residual << endl;
    	}

       	//============================================================================
		//
		// u_old = u
       	// v_old = v
       	// p_old = p
		//
		//============================================================================

       	for(int i=0;i<(Nx+2);i++){
       		for(int j=0;j<(Ny+2);j++){
       			p_old[i][j] = p[i][j];
       		}
       	}

       	for(int i=0;i<(Nx+1);i++){
       		for(int j=0;j<(Ny+2);j++){
       			u_old[i][j] = u[i][j];
       		}
       	}

       	for(int i=0;i<(Nx+2);i++){
       		for(int j=0;j<(Ny+1);j++){
       			v_old[i][j] = v[i][j];
       		}
       	}

    	i_iter++;
     	if (i_iter == (MAX_ITER + 1)) {
      		cout << "You have reached the maximum number of SIMPLE iterations before converging to the set threshold" << endl;
     	}

    }



	// OUTPUT DATA

	for (int i=0; i<=Nx; i++)
	{
		for (int j=0; j<=Ny; j++)
		{
			uc[i][j] = 0.5*(u[i][j]+u[i][j+1]);
			vc[i][j] = 0.5*(v[i][j]+v[i+1][j]);
			pc[i][j] = 0.25*(p[i][j]+p[i+1][j]+p[i][j+1]+p[i+1][j+1]);
		}
	}

    /* write to a file the u_velocity along vertical line through center
     * write to a file the v_velocity along horizontal line through center*/
    ofstream myfile[2];

    myfile[0].open("UVP.txt");
    myfile[1].open("u_velocity_center.txt");
    if (myfile[0].is_open() && myfile[1].is_open())
    {
    	myfile[0] << "VARIABLES=\"X\",\"Y\",\"U\",\"V\",\"P\"\n";
    	myfile[0] << "ZONE  F=POINT\n";
    	myfile[0] << "I=" << Nx <<", J=" << Ny << "\n";
        for(int j=0; j<(Ny+1); j++) {
        	for(int i=0; i<(Nx+1); i++){

        		double xpos, ypos;
        		xpos = i*deltax;
        		ypos = j*deltay;

        		myfile[0] << xpos << "\t" << ypos << "\t" << uc[i][j] << "\t" << vc[i][j] << "\t" << pc[i][j] << "\n";

        	}
        }

        for ( int j = 1 ; j < Ny+1 ; j++ ){
        	double ypos;
			ypos = (double) j*deltay;

        	myfile[1] << (uc[Nx/2][j] + uc[(Nx/2)+1][j])/(2.) << "\t" << ypos << "\n";
        }

        myfile[0].close();
        myfile[1].close();
    }
    else cout << "Unable to open file";


  /* Ending the execution and printing the execution time */
  high_resolution_clock::time_point t2 = high_resolution_clock::now();

  duration<double> time_span = duration_cast<duration<double>>(t2 - t1);

  cout << endl;
  cout << "It took me " << time_span.count() << " seconds.";
  cout << endl;

  return 0;

}
