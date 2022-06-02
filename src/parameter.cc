#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "parameter.hpp"
#include "alloy.hpp"
#include "wanglandau.hpp"
#include "main.hpp"


/*
NOTE - when adding a new input parameter, it is necessary to update all of 
the functions in this file
*/

//This function reads in the input parameters from the input file given at the command line
////////////////////////////////////////////////////////////////////////////
void ReadInput(const char* filename)
{

  char line[200];
  char pname[51];

  FILE* f = fopen(filename, "r");
  if (f == NULL) ErrorMsg(1, filename);

  while (fgets(line, sizeof(line), f) != NULL) {

    if ((line[0] != '#') && (line[0] != '\n')) {

      if (sscanf(line, "%50s", pname) == 1) {

//Number of Atoms
	if (!strcmp(pname, "N")) {
	  if (sscanf(line, "%*50s %d", &alloyState.N) != 1) ErrorMsg(2, pname);
	  if (alloyState.N < 1) ErrorMsg(3, pname);
	  continue;
	}
	
//Nonbonded Potential Selection
	if (!strcmp(pname, "NBINTERACTION")) {
	  if (sscanf(line, "%*50s %d", &alloyState.NBINTERACTION) != 1) ErrorMsg(2, pname);
	  if ( (alloyState.NBINTERACTION < 1) || (alloyState.NBINTERACTION > 4) ) ErrorMsg(3, pname);
	  continue;
	}

//Bin Width for Primary Sampling Direction	
	if (!strcmp(pname, "dWLD1")) {
	  if (sscanf(line, "%*50s %lg", &wlState.dWLD1) != 1) ErrorMsg(2, pname);
	  if (wlState.dWLD1 > 10.0 || (wlState.dWLD1 < 0.000001)) ErrorMsg(3, pname);  //Generic Boundary Limits
	  continue;
	}

//Maximum Boundary for Primary Sampling Direction	
	if (!strcmp(pname, "WLD1min")) {
	  if (sscanf(line, "%*50s %lg", &wlState.WLD1min) != 1) ErrorMsg(2, pname);
	  if (wlState.WLD1min > 20.0) ErrorMsg(3, pname);  
	  continue;
	}
	
//Minimum Boundary for Primary Sampling Direction	
	if (!strcmp(pname, "WLD1max")) {
	  if (sscanf(line, "%*50s %lg", &wlState.WLD1max) != 1) ErrorMsg(2, pname);
//	  if (WLD1max < -20.0) ErrorMsg(3, pname);  
	  continue;
	}

//Flatness Criteria
	if (!strcmp(pname, "Flatness")) {
	  if (sscanf(line, "%*50s %lg", &wlState.Flatness) != 1) ErrorMsg(2, pname);
	  if ( (wlState.Flatness <= 0.0) || (wlState.Flatness >= 1.0) ) ErrorMsg(3, pname);
	  continue;
	}

//Initial Modification Factor
	if (!strcmp(pname, "ModFactorInit")) {
	  if (sscanf(line, "%*50s %lg", &wlState.ModFactorInit) != 1) ErrorMsg(2, pname);
	  if (wlState.ModFactorInit <= 0.0) ErrorMsg(3, pname);
	  continue;
	}

//Modification Factor Iterator
	if (!strcmp(pname, "IterationFactor")) {
	  if (sscanf(line, "%*50s %lg", &wlState.IterationFactor) != 1) ErrorMsg(2, pname);
	  if (wlState.IterationFactor <= 0.0) ErrorMsg(3, pname);
	  continue;
	}

//Final Modification Factor	
	if (!strcmp(pname, "ModFactorFinal")) {
	  if (sscanf(line, "%*50s %lg", &wlState.ModFactorFinal) != 1) ErrorMsg(2, pname);
	  if (wlState.ModFactorFinal <= 0.0) ErrorMsg(3, pname);
	  continue;
	}

//Turns on the production run (0 = off, 1 = on)	
	if (!strcmp(pname, "ProductionBinSamps")) {
	  if (sscanf(line, "%*50s %d", &wlState.ProductionBinSamps) != 1) ErrorMsg(2, pname);
	  if ( wlState.ProductionBinSamps < 1 ) ErrorMsg(3, pname);
	  continue;
	}
		
//Initial Temperture used in Thermoqs()
	if (!strcmp(pname, "TTi")) {
	  if (sscanf(line, "%*50s %lg", &alloyState.TTi) != 1) ErrorMsg(2, pname);
	  if ( alloyState.TTi < 0 ) ErrorMsg(3, pname);
	  continue;
	}
	
//Final Temperture used in Thermoqs()
	if (!strcmp(pname, "TTf")) {
	  if (sscanf(line, "%*50s %lg", &alloyState.TTf) != 1) ErrorMsg(2, pname);
	  if ( alloyState.TTf < 0 ) ErrorMsg(3, pname);
	  continue;
	}
	
//Temperture Increment used in Thermoqs()
	if (!strcmp(pname, "dTT")) {
	  if (sscanf(line, "%*50s %lg", &alloyState.dTT) != 1) ErrorMsg(2, pname);
	  if ( alloyState.dTT < 0 ) ErrorMsg(3, pname);
	  continue;
	}

//Turns on Metropolis Sampling
	if (!strcmp(pname, "MetropolisSampling")) {
	  if (sscanf(line, "%*50s %d", &ptState.MetropolisSampling) != 1) ErrorMsg(2, pname);
	  if ( (ptState.MetropolisSampling < 0) || (ptState.MetropolisSampling > 1) ) ErrorMsg(3, pname);
	  continue;
	}
	
//Initial Temperture used in simple Seq() loop
	if (!strcmp(pname, "MTi")) {
	  if (sscanf(line, "%*50s %lg", &ptState.MTi) != 1) ErrorMsg(2, pname);
	  if ( ptState.MTi < 0 ) ErrorMsg(3, pname);
	  continue;
	}
	
//Final Temperture used in simple Seq() loop
	if (!strcmp(pname, "MTf")) {
	  if (sscanf(line, "%*50s %lg", &ptState.MTf) != 1) ErrorMsg(2, pname);
	  if ( ptState.MTf < 0 ) ErrorMsg(3, pname);
	  continue;
	}
	
//Temperture Increment used in simple Seq() loop
	if (!strcmp(pname, "MdT")) {
	  if (sscanf(line, "%*50s %lg", &ptState.MdT) != 1) ErrorMsg(2, pname);
	  if ( ptState.MdT < 0 ) ErrorMsg(3, pname);
	  continue;
	}

//Number of MC Samples taken in the Seq() function
	if (!strcmp(pname, "MSAMPS")) {
	  if (sscanf(line, "%*50s %lg", &ptState.MSAMPS) != 1) ErrorMsg(2, pname);
	  if ( ptState.MSAMPS < 0 ) ErrorMsg(3, pname);
	  continue;
	}
	
//Number of Samples separating data in Seq()
	if (!strcmp(pname, "MSEP")) {
	  if (sscanf(line, "%*50s %lg", &ptState.MSEP) != 1) ErrorMsg(2, pname);
	  if ( ptState.MSEP < 0 ) ErrorMsg(3, pname);
	  continue;
	}
	
//Number of Dropped Samples before each T run
	if (!strcmp(pname, "MDROP")) {
	  if (sscanf(line, "%*50s %lg", &ptState.MDROP) != 1) ErrorMsg(2, pname);
	  if ( ptState.MDROP < 0 ) ErrorMsg(3, pname);
	  continue;
	}

	ErrorMsg(4, pname);

      }

      else ErrorMsg(5, "");

    }

  }

  fclose(f);

}

//Writes the input parameters to see if they were read in correctly
////////////////////////////////////////////////////////////////////////////
void WriteInput()
{

  printf("\n### Program Parameters #################################\n\n");

  //if(ptState.MetropolisSampling == 0)	
  //{	
	printf("#  Number of Atoms  %d\n", alloyState.N);
	/*
	printf("#  Interaction Constant  %g\n", Jkb);
	printf("#  Number of Energy Bins  %d\n", EBINS);
	printf("#  Maximum Energy Bound  %g\n", WLEmax);
	printf("#  Minimum Energy Bound  %g\n", WLEmin);
	printf("#  Frontier Sampling  %d\n", FrontierSampling);
	printf("#  Maximum Possible Boost Factor  %g\n", MaxBoost);
	printf("#  Number of Sweeps in between Boosts  %d\n", BoostSweeps);
	printf("#  Top Portion of the DOS kept after Boosting  %g\n", Mountain);
	printf("#  Steps Between Smoothing  %d\n", SmoothSkip);
	printf("#  Flatness Criteria  %g\n", Flatness);
	printf("#  Initial Modification Factor  %g\n", ModFactorInit);
	printf("#  Modification Factor Divider  %g\n", IterationFactor);
	printf("#  Final Modification Factor  %g\n", ModFactorFinal);
	*/
  //};
  /*
  if(ptState.MetropolisSampling == 1)
  {
	printf("#  Metropolis Sampling  %d\n", ptState.MetropolisSampling);
	printf("#  Initial Temperature  %g\n", ptState.MTi);
	printf("#  Final Temperature  %g\n", ptState.MTf);
	printf("#  Temperature Increment  %g\n", ptState.MdT);
	printf("#  Number of Samples in Metropolis  %g\n", SeqSAMPS);
	printf("#  Number of Samples Separating Taken Data  %g\n", SeqSEP);
	printf("#  Number of Point Initially Dropped  %g\n", SeqDROP);
  };
  */
  printf("\n########################################################\n\n");

}


//Error Messages
////////////////////////////////////////////////////////////////////////////
void ErrorMsg(int message, const char* arg)
{
  switch (message) {
  case  0 : {
    printf("Syntax: error in or command line arguments [input-file]\n");
    break;
  }
  case  1 : {
    printf("Error: Input-file '%s' not found.\n", arg);
    break;
  }
  case  2 : {
    printf("Error: Parameter '%s' has invalid format.\n", arg);
    break;
  }
  case  3 : {
    printf("Error: Parameter '%s' is out of range.\n", arg);
    break;
  }
  case  4 : {
    printf("Error: Parameter '%s' is unknown.\n", arg);
    break;
  }
  
  case  5 : {
    printf("Error: Unreadable parameter.\n");
    break;
  }
  
  default : printf("Error.\n");
  }
  printf("Program aborted.\n");
  exit(1);
}


