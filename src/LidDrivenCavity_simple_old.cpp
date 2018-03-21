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
#define grid 128

using namespace std;
using namespace std::chrono;

int main (void)
{
   /* Display numbers with maximum precision possible */
	cout.precision(numeric_limits<double>::digits10 + 1);


   /* Variable for tracking execution time */
	high_resolution_clock::time_point t1 = high_resolution_clock::now();

	double u[grid][grid+1], un[grid][grid+1], uc[grid][grid];
	double v[grid+1][grid], vn[grid+1][grid], vc[grid][grid];
	double p[grid+1][grid+1], pn[grid+1][grid+1], pc[grid][grid];
	double m[grid+1][grid+1];
	int i, j, step;
	double dx, dy, dt, delta, error, Re, residual_threshold;
	step =1;
	dx = 1.0/(grid-1);
	dy = 1.0/(grid-1);
	dt = 0.001;
	delta = 4.5;
	error = 1.0;
	Re = 100.0;
	residual_threshold = 1.0E-4;

	// Initializing u
		for (i=0; i<=(grid-1); i++)
		{
			for (j=0; j<=(grid); j++)
			{
				u[i][j] = 0.0;
				u[i][grid] = 1.0;
				u[i][grid-1] = 1.0;
			}
		}

	// Initializing v
		for (i=0; i<=(grid); i++)
		{
			for (j=0; j<=(grid-1); j++)
			{
				v[i][j] = 0.0;
			}
		}

	// Initializing p
		for (i=0; i<=(grid); i++)
		{
			for (j=0; j<=(grid); j++)
			{
				p[i][j] = 1.0;
			}
		}

	while (error > residual_threshold)
	{
		// Solve u-momentum equation
		for (i=1; i<=(grid-2); i++)
		{
			for (j=1; j<=(grid-1); j++)
			{
				un[i][j] = u[i][j] - dt*(  (u[i+1][j]*u[i+1][j]-u[i-1][j]*u[i-1][j])/2.0/dx
							+0.25*( (u[i][j]+u[i][j+1])*(v[i][j]+v[i+1][j])-(u[i][j]+u[i][j-1])*(v[i+1][j-1]+v[i][j-1]) )/dy  )
								- dt/dx*(p[i+1][j]-p[i][j])
									+ dt*1.0/Re*( (u[i+1][j]-2.0*u[i][j]+u[i-1][j])/dx/dx +(u[i][j+1]-2.0*u[i][j]+u[i][j-1])/dy/dy );
			}
		}

		// Boundary conditions
		for (j=1; j<=(grid-1); j++)
		{
			un[0][j] = 0.0;
			un[grid-1][j] = 0.0;
		}

		for (i=0; i<=(grid-1); i++)
		{
			un[i][0] = -un[i][1];
			un[i][grid] = 2 - un[i][grid-1];
		}


		// Solves v-momentum
		for (i=1; i<=(grid-1); i++)
		{
			for (j=1; j<=(grid-2); j++)
			{
				vn[i][j] = v[i][j] - dt* ( 0.25*( (u[i][j]+u[i][j+1])*(v[i][j]+v[i+1][j])-(u[i-1][j]+u[i-1][j+1])*(v[i][j]+v[i-1][j]) )/dx
							+(v[i][j+1]*v[i][j+1]-v[i][j-1]*v[i][j-1])/2.0/dy )
								- dt/dy*(p[i][j+1]-p[i][j])
									+ dt*1.0/Re*( (v[i+1][j]-2.0*v[i][j]+v[i-1][j])/dx/dx+(v[i][j+1]-2.0*v[i][j]+v[i][j-1])/dy/dy );
			}
		}

		// Boundary conditions
		for (j=1; j<=(grid-2); j++)
		{
			vn[0][j] = -vn[1][j];
			vn[grid][j] = -vn[grid-1][j];
		}

		for (i=0; i<=(grid); i++)
		{
			vn[i][0] = 0.0;
			vn[i][grid-1] = 0.0;
		}

		// Solves continuity equation
		for (i=1; i<=(grid-1); i++)
		{
			for (j=1; j<=(grid-1); j++)
			{
				pn[i][j] = p[i][j]-dt*delta*(  ( un[i][j]-un[i-1][j] )/dx + ( vn[i][j]-vn[i][j-1] ) /dy  );
			}
		}


		// Boundary conditions
		for (i=1; i<=(grid-1); i++)
		{
			pn[i][0] = pn[i][1];
			pn[i][grid] = pn[i][grid-1];
		}

		for (j=0; j<=(grid); j++)
		{
			pn[0][j] = pn[1][j];
			pn[grid][j] = pn[grid-1][j];
		}

		// Displaying error
		error = 0.0;

		for (i=1; i<=(grid-1); i++)
		{
			for (j=1; j<=(grid-1); j++)
			{
				m[i][j] = (  ( un[i][j]-un[i-1][j] )/dx + ( vn[i][j]-vn[i][j-1] )/dy  );
				error = error + fabs(m[i][j]);
			}
		}

		if (step%1000 ==1)
		{
			cout << "Error is " << error << " for the step " << step << "\n";
		}


		// Iterating u
		for (i=0; i<=(grid-1); i++)
		{
			for (j=0; j<=(grid); j++)
			{
				u[i][j] = un[i][j];
			}
		}

		// Iterating v
		for (i=0; i<=(grid); i++)
		{
			for (j=0; j<=(grid-1); j++)
			{
				v[i][j] = vn[i][j];
			}
		}

		// Iterating p
		for (i=0; i<=(grid); i++)
		{
			for (j=0; j<=(grid); j++)
			{
				p[i][j] = pn[i][j];
			}
		}

		step = step + 1;

	}

	for (i=0; i<=(grid-1); i++)
	{
		for (j=0; j<=(grid-1); j++)
		{
			uc[i][j] = 0.5*(u[i][j]+u[i][j+1]);
			vc[i][j] = 0.5*(v[i][j]+v[i+1][j]);
			pc[i][j] = 0.25*(p[i][j]+p[i+1][j]+p[i][j+1]+p[i+1][j+1]);
		}
	}



	// OUTPUT DATA

    /* write to a file the u_velocity along vertical line through center
     * write to a file the v_velocity along horizontal line through center*/
    ofstream myfile[2];

    myfile[0].open("UVP.txt");
    myfile[1].open("u_velocity_center.txt");
    if (myfile[0].is_open() && myfile[1].is_open())
    {
    	myfile[0] << "VARIABLES=\"X\",\"Y\",\"U\",\"V\",\"P\"\n";
    	myfile[0] << "ZONE  F=POINT\n";
    	myfile[0] << "I=" << grid <<", J=" << grid << "\n";
        for(int j=0; j<(grid); j++) {
        	for(int i=0; i<(grid); i++){

        		double xpos, ypos;
        		xpos = i*dx;
        		ypos = j*dy;

        		myfile[0] << xpos << "\t" << ypos << "\t" << uc[i][j] << "\t" << vc[i][j] << "\t" << pc[i][j] << "\n";

        	}
        }

        for ( j = 0 ; j < grid ; j++ ){
        	double ypos;
			ypos = (double) j*dy;

        	myfile[1] << (uc[grid/2][j] + uc[(grid/2)+1][j])/(2.) << "\t" << ypos << "\n";
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
