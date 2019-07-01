// Evaluation.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

#include<vector>

//#ifndef __cplusplus
#define max(a,b)    (((a) > (b)) ? (a) : (b))
#define min(a,b)    (((a) < (b)) ? (a) : (b))
//#endif  /* __cplusplus */

using namespace std;

#include <math.h>
#define PI  3.1415926535897932384626433832795

void SEC18_MaOP1(vector<double> &x, vector<double> &f, int nobj, int nvar);
void SEC18_MaOP2(vector<double> &x, vector<double> &f, int nobj, int nvar);
void SEC18_MaOP3(vector<double> &x, vector<double> &f, int nobj, int nvar);
void SEC18_MaOP4(vector<double> &x, vector<double> &f, int nobj, int nvar);
void SEC18_MaOP5(vector<double> &x, vector<double> &f, int nobj, int nvar);
void SEC18_MaOP6(vector<double> &x, vector<double> &f, int nobj, int nvar);
void SEC18_MaOP7(vector<double> &x, vector<double> &f, int nobj, int nvar);
void SEC18_MaOP8(vector<double> &x, vector<double> &f, int nobj, int nvar);
void SEC18_MaOP9(vector<double> &x, vector<double> &f, int nobj, int nvar);
void SEC18_MaOP10(vector<double> &x, vector<double> &f, int nobj, int nvar);


void SEC18_MaOP1(vector<double> &x, vector<double> &f, int nobj, int nvar)
{
	f = vector<double>(nobj, 0);
	// computing distance function
	double g = 0;
	for(int n=nobj; n<=nvar; n++)    // nobj-1 --> nobj
	{
		g +=(x[n-1] - 0.5)*(x[n-1] - 0.5) + 1 - cos(20*PI*(x[n-1]-0.5));  // n --> n-1
	}
	g = g/nvar;  // 100 -->10

	// computing position functions
	double prod_x = 1;
	for(int m=nobj; m>=1; m--)
	{
		int id = nobj-m + 1;                         // id starts with 1 and ends with
		if(m>1)
		{
			f[m-1] = (1 + g)*(1 - prod_x*(1 - x[id-1]));
	        prod_x    = prod_x*x[id-1];                      // x1 --- x1x2...x(m-1)
		}
		else   // the first objective function   f[0]
		{
		    f[m-1] = (1 + g)*(1 - prod_x);
		}

		f[m-1] = (0.1 + 10*m)*f[m-1];    // m --> m - 1
	}
}

void SEC18_MaOP2(vector<double> &x, vector<double> &f, int nobj, int nvar)
{
	f = vector<double>(nobj, 0);
	double g = 0;
	double tmp = 1;
	for(int i=0; i<nobj-1; i++)
	{
		tmp*=sin(0.5*PI*x[i]);
	}

	for(int n=nobj; n<=nvar; n++)
	{
		if((n%5)==0)
		{
			g = g + (x[n-1] - tmp)*(x[n-1] - tmp);
		}
		else
		{
			g = g + (x[n-1] - 0.5)*(x[n-1] - 0.5);   // n --> n-1   2018.5.8
		}
	}

	g = 200*g;

	double tmp2 = 1;

	for(int m = nobj; m>=1; m--)
	{
		int p = pow(2, (m%2)+1);
		if(m==nobj)
		{
			f[m-1] = (1 + g)*pow(sin(0.5*x[0]*PI), p);
		}
		else if(m<nobj&&m>=2)
		{
			tmp2   = tmp2*cos(0.5*PI*x[nobj-m-1]);
		    f[m-1] = (1 + g)*pow(tmp2*sin(0.5*PI*x[nobj-m]), p);
		}
		else
		{
			f[m-1] = (1 + g)*pow(tmp2*cos(0.5*PI*x[nobj-2]), p);
		}
	}
}


void SEC18_MaOP3(vector<double> &x, vector<double> &f, int nobj, int nvar)
{
    f = vector<double>(nobj, 0);
    double g = 0;
	double tmp = 1;
	for(int i=0; i<nobj-1; i++)
	{
		tmp*=sin(0.5*PI*x[i]);
	}
    for(int n=nobj; n<=nvar; n++)
	{
        if(n%5==0)
		{
            g = g + n*pow(abs(x[n-1] - tmp), 0.1);
		}
		else
		{
			g = g + n*pow(abs(x[n-1] - 0.5), 0.1);
	    }
	}

	double tmp2 = 1;

	for(int m = nobj; m>=1; m--)
	{
		if(m==nobj)
		{
			f[m-1] = (1 + g)*sin(0.5*x[0]*PI);
		}
		else if(m<nobj&&m>=2)
		{
			tmp2   = tmp2*cos(0.5*PI*x[nobj-m-1]);
		    f[m-1] = (1 + g)*tmp2*sin(0.5*PI*x[nobj-m]);
		}
		else
		{
			f[m-1] = (1 + g)*tmp2*cos(0.5*PI*x[nobj-2]);
		}
	}
}


void SEC18_MaOP4(vector<double> &x, vector<double> &f, int nobj, int nvar)
{
    f = vector<double>(nobj, 0);
    double g = 0;
	double tmp = 1;
	for(int i=0; i<nobj-1; i++)
	{
		tmp*=sin(0.5*PI*x[i]);
	}
    for(int n=nobj; n<=nvar; n++)
	{
		if(n%5==0)
		{
            g = g + 2*sin(x[0]*PI)*(abs(-0.9*pow(x[n-1] - tmp, 2)) + pow(abs(x[n-1] - tmp), 0.6));
		}
		else
		{
			g = g + 2*sin(x[0]*PI)*(abs(-0.9*pow(x[n-1] - 0.5, 2)) + pow(abs(x[n-1] - 0.5), 0.6));
		}
	}

	g = g*10;

	double tmp2 = 1;

	for(int m = nobj; m>=1; m--)
	{
		if(m==nobj)
		{
			f[m-1] = (1 + g)*sin(0.5*x[0]*PI);
		}
		else if(m<nobj&&m>=2)
		{
			tmp2 = tmp2*cos(0.5*PI*x[nobj-m-1]);
		    f[m-1] = (1 + g)*tmp2*sin(0.5*PI*x[nobj-m]);
		}
		else
		{
			f[m-1] = (1 + g)*tmp2*cos(0.5*PI*x[nobj-2]);
		}
	}
}


void SEC18_MaOP5(vector<double> &x, vector<double> &f, int nobj, int nvar)
{
	f = vector<double>(nobj, 0);
	vector<double> g = vector<double>(nobj, 0);
	for(int m=1; m<=nobj; m++)
	{
		if(m<=3)
		{
			double sum = 0;
			for(int j=3; j<=nvar; j++)
			{
				sum = sum + pow(x[j-1] - x[0]*x[1], 2);
			}
			g[m-1] = max(0, -1.4*cos(2*x[0]*PI)) + sum;
		}
		else
		{
			g[m-1] = exp(pow(x[m-1] - x[0]*x[1], 2)) - 1;  // order of pow and exp is not correct.
		}
		g[m-1] = 10*g[m-1];
	}
	double alpha1 = cos(0.5*PI*x[0])*cos(0.5*PI*x[1]);
	double alpha2 = cos(0.5*PI*x[0])*sin(0.5*PI*x[1]);
	double alpha3 = sin(0.5*PI*x[0]);
	f[0] = (1 + g[0])*alpha1;
	f[1] = 4*(1 + g[1])*alpha2;
	f[2] = (1 + g[2])*alpha3;
	for(int m = 4; m<=nobj; m++)
	{
		double ratio = 1.0*m/nobj;                                                               // m/nobj equals to zero
		f[m-1] = (1 + g[m-1])*(ratio*alpha1 + (1.0-ratio)*alpha2 + sin(0.5*m*PI/nobj)*alpha3);
	}
}

void SEC18_MaOP6(vector<double> &x, vector<double> &f, int nobj, int nvar)
{
	f = vector<double>(nobj, 0);
	vector<double> g = vector<double>(nobj, 0);
	for(int m=1; m<=nobj; m++)
	{
		if(m<=3)
		{
			double sum = 0;
			for(int j=3; j<=nvar; j++)
			{
				sum = sum + pow(x[j-1] - x[0]*x[1], 2);

			}
			g[m-1] = max(0, 1.4*sin(4*x[0]*PI)) + sum;
		}
		else
		{
			g[m-1] = exp(pow(x[m-1] - x[0]*x[1], 2)) - 1;
		}
		g[m-1] = g[m-1]*10;

	}

	double alpha1 = x[0]*x[1];
    double alpha2 = x[0]*(1-x[1]);
    double alpha3 = (1-x[0]);
    f[0] = (1 + g[0])*alpha1;
    f[1] = 2*(1 + g[1])*alpha2;
    f[2] = 6*(1 + g[2])*alpha3;
	for(int m = 4; m<=nobj; m++)
	{
		double ratio = 1.0*m/nobj;
		f[m-1] = (1 + g[m-1])*(ratio*alpha1 + (1.0-ratio)*alpha2 + sin(0.5*m*PI/nobj)*alpha3);
	}
}

void SEC18_MaOP7(vector<double> &x, vector<double> &f, int nobj, int nvar)
{
    f = vector<double>(nobj, 0);

    double g = 0;
	double tmp = 1;
	for(int i=0; i<nobj-1; i++)
	{
		tmp*=sin(0.5*PI*x[i]);
	}

    for(int n=nobj; n<=nvar; n++)
	{
        if(n%5==0)
		{
            g = g + pow(x[n-1] - tmp, 2);
		}
        else
		{
            g = g + pow(x[n-1] - 0.5, 2);
		}
	}

	g = 100*g;

	vector<double> alpha = vector<double>(nobj, 0);
    double tau = sqrt(2)/2;

    alpha[0]  = -pow(2*x[0]-1, 3) +1;
    int  T    = floor((nobj-1.0)/2.0);
    for(int i=1; i<=T; i++)
	{
        alpha[2*i-1] = x[0] + (2*x[i]-1)*tau + tau*pow(abs(2*x[i]-1), 0.5+x[0]);
        alpha[2*i]   = x[0] - (2*x[i]-2)*tau + tau*pow(abs(2*x[i]-1), 0.5+x[0]);
	}
	if(nobj%2==0)
	{
        alpha[nobj-1] = 1 - alpha[0];
    }

    for(int m=1; m<=nobj; m++)
	{
        f[m-1] = (1 + g)*alpha[m-1];
	}
}


void SEC18_MaOP8(vector<double> &x, vector<double> &f, int nobj, int nvar)
{
    f = vector<double>(nobj, 0);

    double g = 0;
	double tmp = 1;
	for(int i=0; i<nobj-1; i++)
	{
		tmp*=sin(0.5*PI*x[i]);
	}

    for(int n=nobj; n<=nvar; n++)
	{
        if(n%5==0)
		{
            g = g + pow(x[n-1] - tmp, 2);
		}
        else
		{
            g = g + pow(x[n-1] - 0.5, 2);
		}
	}

	g = g*100;

	vector<double> alpha = vector<double>(nobj, 0);
    double tau = sqrt(2)/2;

    alpha[0]  = -pow(2*x[0]-1, 3) + 1;
    int  T    = floor((nobj-1.0)/2.0);
    for(int i=1; i<=T; i++)
	{
        alpha[2*i-1] = x[0] + (2*x[i])*tau + tau*pow(abs(2*x[i]-1), 1-0.5*sin(4*PI*x[0]));
        alpha[2*i]   = x[0] - (2*x[i]-2)*tau + tau*pow(abs(2*x[i]-1), 1-0.5*sin(4*PI*x[0]));
	}

	if(nobj%2==0)
	{
        alpha[nobj-1] = 1 - alpha[0];
    }


    for(int m=1; m<=nobj; m++)
	{
        f[m-1] = (1 + g)*alpha[m-1];
	}
}

void SEC18_MaOP9(vector<double> &x, vector<double> &f, int nobj, int nvar)
{
    f = vector<double>(nobj, 0);

    double g = 0;
	double tmp = 1;
	for(int i=0; i<nobj-1; i++)
	{
		tmp*=sin(0.5*PI*x[i]);
	}

    for(int n=nobj; n<=nvar; n++)
	{
        if(n%5==0)
		{
            g = g + pow(x[n-1] - tmp, 2);
		}
        else
		{
            g = g + pow(x[n-1] - 0.5, 2);
		}
	}

	g = g*100;

	vector<double> alpha = vector<double>(nobj, 0);
    double tau = sqrt(2)/2;

    alpha[0]  = -pow(2*x[0]-1, 3) +1 ;
    int  T    = floor((nobj-1.0)/2.0);
    for(int i=1; i<=T; i++)
	{
		double z = 2*(2*x[i] - floor(2*x[i])) - 1;
        alpha[2*i-1] = x[0] + 2*x[i]*tau + tau*pow(abs(z), 0.5 + x[0]);
        alpha[2*i]   = x[0] - (2*x[i]-2)*tau + tau*pow(abs(z), 0.5 + x[0]);
	}

	if(nobj%2==0)
	{
        alpha[nobj-1] = 1 - alpha[0];
    }


    for(int m=1; m<=nobj; m++)
	{
        f[m-1] = (1 + g)*alpha[m-1];
	}
}

void SEC18_MaOP10(vector<double> &x, vector<double> &f, int nobj, int nvar)
{
    f = vector<double>(nobj, 0);

    double g = 0;
	double tmp = 1;
	for(int i=0; i<nobj-1; i++)
	{
		tmp*=sin(0.5*PI*x[i]);
	}

    for(int n=nobj; n<=nvar; n++)
	{
        if(n%5==0)
		{
            g = g + pow(x[n-1] - tmp, 2);
		}
        else
		{
            g = g + pow(x[n-1] - 0.5, 2);
		}
	}

	g = g*100;

	vector<double> alpha = vector<double>(nobj, 0);
    double tau = sqrt(2)/2, p, z;

    alpha[0]  = -pow(2*x[0]-1, 3) + 1;
    int  T    = floor((nobj-1.0)/2.0);
    for(int i=1; i<=T; i++)
	{
        z = 2*(2*x[i] - floor(2*x[i])) - 1;  // i+1 --> i  2018.5.8
        if(x[i]<0.5)
		{
            p = 0.5 + x[0];
		}
		else
		{
			p = 1.5 - x[0];
		}

		alpha[2*i-1] = x[0] + 2*x[i]*tau + tau*pow(abs(z), p);
        alpha[2*i]   = x[0] - (2*x[i]-2)*tau + tau*pow(abs(z), p);
	}

	if(nobj%2==0)
	{
        alpha[nobj-1] = 1 - alpha[0];
    }

    for(int m=1; m<=nobj; m++)
	{
        f[m-1] = (1 + g)*alpha[m-1];
	}

}