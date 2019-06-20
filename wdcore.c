#include <stdio.h>
#include <math.h>

/* maximum number of array points */
#define Ncoremax 1000
/* This determines the size of the mass steps to try in the */
/* runge kutta integration. */
/* Each mass step is initially mass/InitialMSteps */
/* Could be much smaller then 200 but this high number gives nicer plots */
#define InitialMSteps 2e2

/* If we need a massstep smaller than MinimalMstep to get the */
/* required aquiracy stop the integration */
#define MinimalMStep 1e-20

/* MaxRelErrorQ determines how accurate we want each core-integration to */
/* approach the requested density at the outer edge of the core */ 
#define MaxRelErrorQ 1e-8 

/* MaxExponent gives upperlevel to values allowed to be fed into exp() */
#define MaxExponent 1e2

/* Prototype declarations of all functions used by the program */
void init();
int menu(double *rhoc);
int findrhoc(double rhocstart);
int intcore(double rhoc);
int rk5(int j);
double dsdm(double ls,double lq);
double dqdm(double lm,double ls,double lq);

/* Some constants used by the program */
const double c      =2.99791e8, /* m/s */
             pi     =3.14159,
             me     =9.10953e-31, /* kg */ 
             mp     =1.67265e-27, /* kg */
             planck =6.6237e-34, /* Js */
             k      =1.38024e-23, /* J/K */
             G      =6.668e-11, /* N m^2/kg^2 */
             Lsun   =3.86e26, /* W */
             Msun   =1.991e30, /* kg */
             Rsun   =6.960e8, /* m */
             mue    =2e0,
             mu     =1.3e0;

/* More constants used by the program. They are set in the function init */
double A,B,X,Y;

/* These variables hold the values of the model we are looking for: */
/* ta   =  absolute error tolerated in one integration step. */
/* tr   =  relative error tolerated in one integration step. */
/* MaxRelErrorM determines how accurate the final coremass should compare */
/* to the requested coremass */
/* mass =  requested mass of the core. */
/* x    =  requested value of the degeneracy parameter at the outside */
/*         of the core. */
/* rhocfinal = value of rhoc chosen to get satisfactory core model */
/* They are set to some default values. */
double ta = 0.03,
       tr = 1.2e-2,
       MaxRelErrorM=1e-6,
       mass = 7.795353, /* (= 1.053*Msun/Y) */
       x    = .193369,
       finalrhoc = 0e0;

/* These arrays store the integrated values of the mass,density and radius */
double m[Ncoremax],q[Ncoremax],s[Ncoremax];

/* These arrays store the integration-errors of density and radius */
double qerror[Ncoremax],serror[Ncoremax];

/* flag indicating that if set an exponent that was too big was */
/* used in one of the structure equations (dsdm,dqdm) */
int dfdmwarn;

/* name of file to use for output */
char outfile[40] = "wdcore.dat";

void main()
{
/* pointer to file to hold output */
  FILE *fp;

/* IOuter will receive index of the m[],q[] and s[] that contain the border */
/* of the model that has the right mass and outside density */ 
/* idummy is a dummy variable to hold the choice made in the menu */
  int iOuter,idummy,i;

/* dummy variable to hold user responses to questions */
/* we are not really interested in their value */
  char cdummy;

/* rhoc hold the initial guess for the central density */
  double rhoc = 1e9;

/*  just set some constant expressions, and set some variables to zero */
  init();
  
  do /* while (1) => repeat endlessly */
    {
    
/*  menu of choices */
      do 
	{idummy = menu(&rhoc);} while(idummy > 0);
      if (idummy == -1) exit(0);

/* seems to be needed to clear the keyboard buffer */
      scanf("%c",&cdummy);

/*  Make model with right mass and x by looking for the right rhoc */
      iOuter = findrhoc(rhoc);

/* output results */
      fp = fopen(outfile,"w+");
      fprintf(fp,"central density: %e\n",finalrhoc);
      fprintf(fp,"number,mass,radius,density,s-error,q-error\n");

      for(i=0;i<=iOuter;i++)
	{
	  fprintf(fp,"%d %e %e %e %e %e\n", /* write to the file */
		  i,
		  m[i]*Y/Msun,
		  exp(s[i])*X/Rsun,
		  exp(q[i])*B,
		  serror[i],
		  qerror[i]);
	}
         fclose(fp);
      printf("[Return] to continue\n");
      scanf("%c",&cdummy);
    } while(1);
}

/* Set the additional constants A,B,X,Y. Set dfdmwarn,m[],q[],s[] */
/* and qerror[],serror[]  to 0 */
void init()
{
  int j;
  A      =6.002608e21;
  B      =9.810486e8*mue;
  X      =4.448136e6/mue;
  Y      =1.075781e30/(mue*mue);
  dfdmwarn = 0;
  for (j=0;j<Ncoremax;j++) s[j]=q[j]=m[j]=qerror[j]=serror[j]=0e0;
}

/* Ask the user for model parameters: */
/* mass,x,ta,tr,MaxRelErrorM */
/* Also let the user decide on initial value of rhoc to try */

/* input to init(): */
/*   *rhoc  */

/* output by init(): */
/*    -1 user chose stop */
/*     0 user chose start model */
/* other user changed some value */

/* global variables changable by init(): */
/* mass,x,tr,ta,MaxRelErrorM,outfile */
int menu(double *rhoc)
{
/* variable to hold user input on menu */
  int input;

/* variable to hold value entered by user */
  float fdummy;

  printf("**********************************************\n");
  printf("1) coremass in solarmass: %f\n",mass*Y/Msun);
  printf("2) x at the outside of the core: %f\n",x);
  printf("3) initial value of central density(kg/m**3): %e\n",*rhoc);
  printf("4) maximum absolute error in integration step: %f\n",ta);
  printf("5) maximum relative error in integration step: %f\n",tr);
  printf("6) maximum relative error in mass: %e\n",MaxRelErrorM);
  printf("7) output filename: %s\n",outfile);
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
      printf("coremass in solarmass:\n");
      scanf("%f",&fdummy);
      mass = Msun*(double)fdummy/Y;
      return(1);
    }
  case 2:
    {
      printf("x at the outside of the core:\n");
      scanf("%f",&fdummy);
      x = (double)fdummy;
      return(2);
    }
  case 3:
    {
      printf("initial value of central density(kg/m**3):\n");
      scanf("%f",&fdummy);
      *rhoc = (double)fdummy;
      return(3);
    }
  case 4:
    {
      printf("maximum absolute error in integration step:\n");
      scanf("%f",&fdummy);
      ta = (double)fdummy;
      return(4);
    }
  case 5:
    {
      printf("maximum relative error in integration step:\n");
      scanf("%f",&fdummy);
      tr = (double)fdummy;
      return(5);
    }
  case 6:
    {
      printf("maximum relative error in mass:\n");
      scanf("%f",&fdummy);
      MaxRelErrorM = (double)fdummy;
      return(5);
    }
  case 7:
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

/* findrhoc(rhocstart) */
/* Function that takes a central density(rhoc) and calculates a model  */
/* coremass. Next it adjusts the rhoc, to make the coremass match better  */
/* with requested mass, until the difference between those is less than  */
/* the tolerated error MaxRelErrorM. */
/* input to findrhoc(): */
/*   rhocstart first guess for rhoc  */

/* output by findrhoc(): */
/*    index where the integrated values of m,q and s can be found */
/*    for the total core. */

int findrhoc(double rhocstart)
{
/* rhoc: central density to do integration for. */
/* rhoc0: previous value of rhoc used. */
/* rhocstep: amount to change rhoc to get to new rhoc, */
/*   so rhocstep=rhoc-rhoc0 */
/* mm: resulting coremass */
  double rhoc,rhoc0,rhocstep,mm;

/* index where the integrated values of m,q and s can be found */
  int iOuter;
  
  rhoc0 = 0e0;
  rhoc = rhocstart;
  rhocstep = rhoc;
  
  do /* first increase rhoc until mass is bigger than requested */
    {
      /* make a model with the current value of rhoc */
      iOuter = intcore(rhoc); 

      /* get the resulting mass */
      mm = m[iOuter];

      /* if we accidently run into a good solution now we can stop */
      if (fabs(mm/mass-1e0) < MaxRelErrorM) return(iOuter);

      /* determine better rhoc */ 
      if (mass > mm)
	{
	  /* far from the requested mass so probably rhoc is  */
	  /* also far off make a big step thus */
    	  if (mass/mm > 2e0)
	    {
    	      printf("nowhere near %e\n",mass/mm);
    	      rhoc0 = rhoc;

	      /* assume mass = constant*rhoc */
	      /* then real rhoc = mass/mm*rhoc  */
    	      rhoc = mass/mm*rhoc;
    	      rhocstep = rhoc - rhoc0;
    	    }
    	  else
    	    {
    	      rhoc0 = rhoc;
    	      rhoc += rhocstep;
    	    }
	}
    } while ( mass > mm ); 

/* we found a lowerlimit of the rhoc (rhoc0) and an upperlimit rhoc */
/* now make models in between until the difference between the model mass */
/* and the requested mass is smaller than MaxRelErrorM */
/* we do this with an halfstep refinement methode for rhoc: */
/* 1 make rhocstep half of previous rhocstep */
/* 2 if last modelmass was too low increase rhoc with rhocstep */
/* else decrease rhoc with rhocstep */
/* 3 calculate a modelcoremass with this new rhoc */
/* repeat 1,2 and 3 until error-condition satified */


  rhocstep /= 2e0;
  rhoc -= rhocstep;
  do
    {
      iOuter = intcore(rhoc);
      mm = m[iOuter];

      rhocstep /= 2e0;
      if ( mm > mass)
	{
	  rhoc -= rhocstep;
	}
      else
	{
	  rhoc += rhocstep;
	}
    } while (fabs(mm/mass-1e0) > MaxRelErrorM);
  return(iOuter);
}

/* intcore(lrhoc) */
/* This function will given a central density integrate the structure */
/* equations for the density and radius until the density is close to */
/* the requested density rho(x) */

/* input to intcore(): */
/*   lrhoc: central density to do the integration for  */

/* output by intcore(): */
/*    index where the integrated values of m,q and s can be found */
/*    for the total core. */

/* global variables changable by intcore(): */
/* m[],q[],s[] */

int intcore(double lrhoc)
{
  /* i: index to array element to do integration for */
  /* rk5failed: receives value from rk5 indicating accuracy was not reached */
  int i,rk5failed;

/* qend: density at the outside of the core */
/* mstep: stepsize we use when integrating from old value to new value */
  double qend,mstep;

  mstep = (mass/InitialMSteps);
  qend = 3e0*log(x);
  i = 1;

/*   we put in for the start of the integration : */
/*   m= 1/1000000 of Msun, q = ln(lrhoc/B) => r = (3/(4*pi) m/lrhoc)^1/3 */
/*   s = 1/3*ln(3/(4*pi)*m/lrhoc) - ln(X) */
  m[0] = 1e-6*Msun/Y;
  q[0] = log(lrhoc/B);
  s[0] = 1e0/3e0*log(3e0/(4e0*pi)*m[0]*Y/lrhoc)-log(X);


/* We set finalrhoc to lrhoc, if it does not get set again here */
/* it was indeed the final lrhoc ;-) */
  finalrhoc = lrhoc;

  do /* increase mass until density is lower then the outside density qend */
    {
      m[i] = m[i-1]+mstep;
      rk5failed = rk5(i-1);
      if (rk5failed) /* did not make the accuracy test so decrease stepsize*/
	{
	  mstep /= 2e0;

	  if (mstep < MinimalMStep) /* don't half mstep endlessly */
	    {
	      printf("The mass step has become too small\n");
	      exit(0); 
	    }
	}
      else  /* accurate enough next mass point */
	{
	  /* we accidently bumped into a good solution here,so stop */
	  if (fabs(q[i]/qend-1e0) < MaxRelErrorQ) return(i);
	  i++;
	  if (i >= Ncoremax)  /* more masspoints needed than allocated */
	    {
	      printf("Couldn't find solution, too few grid points\n");
	      exit(0);
	    }
	  mstep = (mass/InitialMSteps); /* make mstep big again */
	}
    } while (exp(q[i-1]) > exp(qend)); /* stepped outside of core */

/* we found the core edge between m[i-1] and m[i-2] */
/* now make models in between until the difference between the model density */
/* and the requested density (qend) is smaller than MaxRelErrorQ */
/* we do this with an halfstep refinement methode for m[i]: */
/* 1 make mstep half of previous mstep */
/* 2 if last density was too high increase m[i] with mstep */
/* else decrease m[i] with mstep */
/* 3 calculate a density with this new m[i] */
/* repeat 1,2 and 3 until error-condition satified */

  i--;
  mstep = m[i] - m[i-1];
  do /* do halfstep procedure until error-condition satified */
    {	  
      mstep /= 2e0;
      if (mstep < MinimalMStep)  /* don't half mstep endlessly */
	{
	  printf("The mass step has become too small\n");
	  exit(0);
	}
      if (exp(q[i]) > exp(qend)) 
	{
	  m[i] += mstep;
	}
      else
	{
	  m[i] -= mstep;
	}
      rk5failed = rk5(i-1);
      /* we decreased mstep but now error has become bigger */
      if (rk5failed) printf("something stinks");
    } while (fabs(q[i]/qend-1e0) > MaxRelErrorQ);
  printf("central density: %e\nmass:%f radius:%e density:%e\n",
	 lrhoc,m[i]*Y/Msun,exp(s[i])*X/Rsun,exp(q[i])*B);
  return(i);
}

/* rk5(j) */
/* This function will make 1 integration step for s and q */
/* from mass point m[j] to m[j+1] resulting in s[j+1] and q[j+1] */
/* it will also evaluate the error in these integrations */

/* input to rk5(): */
/*   j: index of m[],s[] and q[] values to start integration from */

/* output by rk5(): */
/*    boolean value indicating whether the requested accuracy prescription */
/*    has been violated */

/* global variables changable by intcore(): */
/* q[j+1],s[j+1],qerror[j+1],serror[j+1] */

int rk5(int j)
{
/* h: mass stepsize */
  double h,k[2][7],n[2][6];

/* clear dfdmwarn flag for new integration */
  dfdmwarn = 0;

  h = m[j+1]-m[j];

  k[0][0] = h * dsdm(s[j],q[j]);
  k[1][0] = h * dqdm(m[j],s[j],q[j]);
  n[0][0] = s[j] + 2e0/9e0*k[0][0];
  n[1][0] = q[j] + 2e0/9e0*k[1][0];

  k[0][1] = h * dsdm(n[0][0],n[1][0]);
  k[1][1] = h * dqdm(m[j]+2e0/9e0*h,n[0][0],n[1][0]);
  n[0][1] = s[j] + k[0][0]/12e0+k[0][1]/4e0;
  n[1][1] = q[j] + k[1][0]/12e0+k[1][1]/4e0;

  k[0][2] = h * dsdm(n[0][1],n[1][1]);
  k[1][2] = h * dqdm(m[j]+h/3e0,n[0][1],n[1][1]);
  n[0][2] = s[j] + k[0][0]/8e0+3e0/8e0*k[0][2];
  n[1][2] = q[j] + k[1][0]/8e0+3e0/8e0*k[1][2];

  k[0][3] = h * dsdm(n[0][2],n[1][2]);
  k[1][3] = h * dqdm(m[j]+h/2e0,n[0][2],n[1][2]);
  n[0][3] = s[j] +  53e0/125e0*k[0][0] -  27e0/25e0*k[0][1]
                 + 126e0/125e0*k[0][2] + 56e0/125e0*k[0][3];
  n[1][3] = q[j] +  53e0/125e0*k[1][0] -  27e0/25e0*k[1][1]
                 + 126e0/125e0*k[1][2] + 56e0/125e0*k[1][3];

  k[0][4] = h * dsdm(n[0][3],n[1][3]);
  k[1][4] = h * dqdm(m[j]+4e0/5e0*h,n[0][3],n[1][3]);
  n[0][4] = s[j] - 9e0/4e0*k[0][0] +  27e0/4e0*k[0][1]- 9e0/7e0*k[0][2] 
                 -   4e0*k[0][3] + 25e0/14e0*k[0][4];
  n[1][4] = q[j] - 9e0/4e0*k[1][0] +  27e0/4e0*k[1][1]- 9e0/7e0*k[1][2] 
                 -   4e0*k[1][3] + 25e0/14e0*k[1][4];

  k[0][5] = h * dsdm(n[0][4],n[1][4]);
  k[1][5] = h * dqdm(m[j]+h,n[0][4],n[1][4]);
  n[0][5] = s[j] + 19e0/24e0*k[0][0] -    9e0/4e0*k[0][1]+ 23e0/14e0*k[0][2] 
                 +   2e0/3e0*k[0][3] + 25e0/168e0*k[0][4];
  n[1][5] = q[j] + 19e0/24e0*k[1][0] -    9e0/4e0*k[1][1]+ 23e0/14e0*k[1][2] 
                 +   2e0/3e0*k[1][3] + 25e0/168e0*k[1][4];

  k[0][6] = h * dsdm(n[0][5],n[1][5]);
  k[1][6] = h * dqdm(m[j]+h,n[0][5],n[1][5]);

  s[j+1] = s[j] +    5e0/48e0*k[0][0] + 27e0/56e0*k[0][2]
                + 125e0/336e0*k[0][4] +  k[0][5]/24e0;
  q[j+1] = q[j] +    5e0/48e0*k[1][0] + 27e0/56e0*k[1][2]
                + 125e0/336e0*k[1][4] +  k[1][5]/24e0;


  serror[j+1] = fabs(     3e0/2e0*k[0][0] - 81e0/7e0*k[0][2] + 16e0*k[0][3] 
	      - 125e0/14e0*k[0][4] +      3e0*k[0][6])
            /(ta*fabs(h)+tr*fabs(k[0][0]));

  qerror[j+1] = fabs(     3e0/2e0*k[1][0] - 81e0/7e0*k[1][2] + 16e0*k[1][3] 
	      - 125e0/14e0*k[1][4] +      3e0*k[1][6])
            /(ta*fabs(h)+tr*fabs(k[1][0]));

/* if either error is too big return with flag set */
  return((qerror[j+1] > 1e0) || (serror[j+1] > 1e0) || dfdmwarn);

}

/* function dsdm(ls,lq) */
/* structure equation for s */
/* gives warning if the exponent is too big */

/* input to dsdm: */
/*   ls: value of s to evaluate dsdm for */
/*   lq: value of q to evaluate dsdm for */

/* output by dsdm(): */
/*   value of dsdm calculated */

/* global variables changable by dsdm(): */
/* dfdmwarn indicating whether the exponent has grown too big */

double dsdm(double ls,double lq)
{
  double te;
  te = -lq-3e0*ls;

/* if we want to feed to big a number to the exp function raise warnflag */
  if ( te > MaxExponent ) dfdmwarn = 1;
  return(exp(te));
}

/* function dqdm(lm,ls,lq) */
/* structure equation for q */
/* gives warning if the exponent is too big */

/* input to dqdm: */
/*   lm: value of m to evaluate dqdm for */
/*   ls: value of s to evaluate dqdm for */
/*   lq: value of q to evaluate dqdm for */
/* output by dqdm(): */
/*   value of dqdm calculated */
/* global variables changable by dqdm(): */
/* dfdmwarn indicating whether the exponent has grown too big */

double dqdm(double lm,double ls,double lq)
{
  double te1,te2;
  te1=-4e0*ls-5e0*lq/3e0;
  te2=2e0*lq/3e0;
/* if we want to feed to big a number to the exp function raise warnflag */
  if ( (te1 > MaxExponent) || (te2 > MaxExponent) ) dfdmwarn = 1;
  return(-1e0*lm*exp(te1)*sqrt(exp(te2)+1e0));
}
