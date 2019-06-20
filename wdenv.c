#include <stdio.h>
#include <math.h>

/* This determines the smallest step in w-space */
/* if the wstep goes bellow this limit stop searching for better w */ 
#define MinWstep 1e-300

/* Prototype declarations of all functions used by the program */
void init();
int menu(double *RadiusStart,double *RadiusEnd,int *NRadiusSteps);
void findw(double radius);
double Pe(double a,double b,double w);
double Pd(double a,double b,double w);
double X(double a,double b,double w);
double RADIUS(double a,double w,double R);
double MASS(double a,double w);

/* Some constants used by the program */
const double c      =2.99791e8, // m/s
             pi     =3.14159,
             me     =9.10953e-31, // kg 
             mp     =1.67265e-27, // kg
             planck =6.6237e-34, // Js
             k      =1.38024e-23, // J/K
             G      =6.668e-11, // N m^2/kg^2
             Lsun   =3.86e26, // W
             Msun   =1.991e30, // kg
             Rsun   =6.960e8, // m
             mue    =2e0,
             mu     =1.3e0;

/* These variables hold the values of the model we are looking for: */
/* M       =  mass of the dwarf in solarmass. */
/* L =  luminosity in solarluminosity */
/* They are set to some default values. */

double M = 1.053,
       L = 3.02e-3;

/* pointer to file to hold output */
  FILE *fp;

/* name of file to use for output */
char outfile[40] = "wdenv.dat";
         
void main()
{
/* dummy variable to hold user responses to questions */
/* we are not really interested in their value */
  char cdummy;

/* RadiusStart,RadiusEnd hold minimum and maximum radii to make models for */
/* in solarradii  also set to default values*/
  double 
    RadiusStart = .001,
    RadiusEnd   = .011;

/* NRadiusSteps hold number of radii to make envelope for */
  int 
    NRadiusSteps = 21,
    i,idummy;

  do /* while (1) => repeat endlessly */
    {
    
/*  menu of choices */
      do 
	{
	  idummy = menu(&RadiusStart,&RadiusEnd,&NRadiusSteps);
	} while(idummy > 0);
      if (idummy == -1) exit(0);

/* seems to be needed to clear the keyboard buffer */
      scanf("%c",&cdummy);

      fp = fopen(outfile,"w+");
/* make header for table */
      printf("  Router     Rinner    MassInner   x\n");
      fprintf(fp,"Router,Rinner,MassInner,x\n");
      for (i=0;i<NRadiusSteps;i++)
	{
	  findw(RadiusStart+
		(double)i/(double)(NRadiusSteps-1)*(RadiusEnd-RadiusStart));
	}
      fclose(fp);
      printf("[Return] to continue\n");
      scanf("%c",&cdummy);
    } while(1);
}

/* Ask the user for model parameters: */
/* mass,luminosity */
/* Also let the user decide on radii to make envelope for: */
/* RadiusStart,RadiusEnd,NRadiusSteps */

/* input to init(): */
/*   *RadiusStart,*RadiusEnd,*NRadiusSteps */

/* output by init(): */
/*    -1 user chose stop */
/*     0 user chose start model */
/* other user changed some value */

/* global variables changable by init(): */
/* M,L */
int menu(double *RadiusStart,double *RadiusEnd,int *NRadiusSteps)
{
/* variable to hold user input on menu */
  int input;

/* variable to hold value entered by user */
  float fdummy;

  printf("**********************************************\n");
  printf("1) mass in solarmass: %f\n",M);
  printf("2) luminosity in solarluminosity: %f\n",L);
  printf("3) Smallest Radius: %f\n",*RadiusStart);
  printf("4) Largest Radius: %f\n",*RadiusEnd);
  printf("5) Number of Radii: %d\n",*NRadiusSteps);
  printf("6) Output filename: %s\n",outfile);
  printf("9) start model\n");
  printf("0) EXIT\n");
  printf("**********************************************\n");
  printf("Choice: ");
  scanf("%d",&input);  

  switch(input)
  {
  case 9: return(0);
  case 1:
    {
      printf("mass in solarmass:\n");
      scanf("%f",&fdummy);
      M = (double)fdummy;
      return(1);
    }
  case 2:
    {
      printf("luminosity in solarluminosity:\n");
      scanf("%f",&fdummy);
      L = (double)fdummy;
      return(2);
    }
  case 3:
    {
      printf("Smallest Radius:\n");
      scanf("%f",&fdummy);
      *RadiusStart = (double)fdummy;
      return(3);
    }
  case 4:
    {
      printf("Largest Radius:\n");
      scanf("%f",&fdummy);
      *RadiusEnd = (double)fdummy;
      return(4);
    }
  case 5:
    {
      printf("Number of Radii:\n");
      scanf("%d",&*NRadiusSteps);
      return(5);
    }
  case 6:
    {
      printf("output filename:\n");
      scanf("%s",&outfile);
      return(6);
    }
  case 0:
    {
      return(-1);
    }
  }
  return(0);
}

/* findw(R) */
/* procedure that takes a radius(R) and searches for the w */
/* that makes the electron pressure on the inner edge of the envelope */
/* equal to the degenerate pressure until the w differs less than */
/* MinWstep from the true w. */

/* input to find(radius): */
/*   outer radius of the dwarf */

void findw(double R)
{
/* w holds the value at which we evaluate Pe and Pd */
/* wstep is the amount we change w each loop to get a better w */
/* alpha holds value of opacity parameter */
/* radius,mass,x hold the radius,mass and x for the final w value */

  double w,wstep,alpha,beta,radius,mass,x;

  alpha = 6.2716e-3*pow(mu*pow(L,2e0)*R/pow(M,3e0),1e0/4e0);
  beta = 1.44029e-3*mu*M/R;

/* we assume right value of w is between 0 ans 2 */
/* set w=wstep=1 */
/* 1 make wstep half of wstep */
/* 2 if Pe(w) < Pd(w) then decrease w by wstep else increase by wstep  */
/* repeat 1 and 2 if wstep > MinWstep */

/* start in between 0 and 2 */
  w = wstep = 1e0;
 
  do
    {
      wstep /= 2e0;
      (Pd(alpha,beta,w) > Pe(alpha,beta,w)) ?  (w -= wstep) : (w += wstep);
    } while (wstep > MinWstep);

  radius = RADIUS(alpha,w,R);
  mass = MASS(alpha,w);
  x = X(alpha,beta,w);
  fprintf(fp,"%e %e %e %e\n",
	 R,radius,mass,x);
  printf("%f %e %e %e\n",
	 R,radius,mass,x);
}

/* Pe(a,b,w) */
/* Returns value of the electron pressure for given a,b,w */
double Pe(double a,double b,double w)
{
  return(124e0/15e0*pow(b,4e0)*mu/a/mue*pow(w,9e0)*pow(w+a,8e0));
}

/* Pd(a,b,w) */
/* Returns value of the degenerate pressure for given a,b,w */
double Pd(double a,double b,double w)
{
  double x;
  x = X(a,b,w);
  return(x*(2e0*x*x-3e0)*sqrt(x*x+1e0)+3e0*log(x+sqrt(x*x+1e0)));
}

/* MASS(a,w) */
/* returns interior-mass for a given a,w */
/* uses M the total starmass */
double MASS(double a,double w)
{
  const double sqrt2 = 1.41421356237309504880;
  double f,w4,w4w41;

  w4 = pow(w,4e0);
  w4w41 = w4/(w4+1);

  f =  195e0 / 32e0*w
     - 195e0 / 128e0 * sqrt2 
     * ( log( w*w + w*sqrt2 + 1e0 ) + 2e0 * atan((w*w+sqrt(1+w4)-1e0)/w/sqrt2))
     + (195e0/256e0*sqrt2-253e0/51e0*a)*log(w4+1)
     + 172e0/51e0*a*pow(w4w41,4e0)
     - (w/3e0 - 253e0/153e0*a)*pow(w4w41,3e0)
     - (13e0/24e0*w-253e0/102e0*a)*w4w41*w4w41
     - (117e0/96e0*w-253e0/51e0*a)*w4w41;

  return(M*(1e0-.0101243*pow(mu,4e0)*M*M*f/a));
}

/* RADIUS(a,w,R) */
/* returns radius for a given a,w */
/* uses R, the starradius */
double RADIUS(double a,double w,double R)
{
  return(R/(pow(w+a,3e0)*(w+19e0/51e0*a) - 19e0/51e0*pow(a,4e0) + 1e0));
}

/* X(a,b,w) */
/* returns degeneracy parameter for a given a,b,w */
double X(double a,double b,double w)
{
  return(pow(31e0*pi/60e0*mu/a/mue,1e0/3e0)*b*pow(w,7e0/3e0)*pow(w+a,2e0));
}

