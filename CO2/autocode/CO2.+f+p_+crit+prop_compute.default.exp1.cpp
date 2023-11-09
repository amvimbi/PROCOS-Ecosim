//==================================================================
// Code generated automatically
// Description: Experiment class body 
//==================================================================
#include "CO2.+f+p_+crit+prop_compute.default.exp1.h"

CO2__FP_CritProp_compute_default_exp1* CO2__FP_CritProp_compute_default_exp1::s_current= NULL;
CO2__FP_CritProp_compute_default_exp1::CO2__FP_CritProp_compute_default_exp1(const char* mgr,const char* installDir,bool dmode):CO2__FP_CritProp_compute_default(mgr,installDir,dmode)
{
m_infoLibraryName = "CO2";
m_infoGroupName = "";
m_infoExperimentName = "exp1";
m_infoExperimentFileName = "CO2.+f+p_+crit+prop_compute.default.exp1";
m_infoExperimentFileNameExtra = "exp1";
m_infoExperimentDate = "15/06/2023 14:58:53.041000";

m_perfFlag = false;

// Global variables
n_typ_exp=0;
// Experiment variables (initialisation)
nBounds = 4;// Number of boundary variables
if (m_boundaryBranch==NULL)
{
	m_boundaryBranch= new int[4];
	for (int ib=0; ib < 4 ; ib++)
		m_boundaryBranch[ib]= 0;
}
setDefaultReportSeparator("\t");
s_current= this;
}

  // init internal partition models if used
void CO2__FP_CritProp_compute_default_exp1::initInternalModels()
{
  CO2__FP_CritProp_compute_default::initInternalModels();
}


//Add the experiment variables to the symbols table
void CO2__FP_CritProp_compute_default_exp1::addExptVariables()
{
varHasEquationInBoundBlock("h");
varHasEquationInBoundBlock("rho");
varHasEquationInBoundBlock("u");
varHasEquationInBoundBlock("x");

}


//Run the experiment code
void CO2__FP_CritProp_compute_default_exp1::runExperiment()
{
	DEBUG_LEVEL = 1 ;
	IMETHOD = 11 ;
	REL_ERROR = 1e-06 ;
	ABS_ERROR = 1e-06 ;
	TOLERANCE = 1e-06 ;
	INIT_INTEG_STEP = -1. ;
	MAX_INTEG_STEP = -1. ;
	NSTEPS = 1 ;
	REPORT_MODE = 9 ;
	TIME = 0. ;
	TSTOP = 15. ;
	CINT = 0.1 ;
	INTEG() ;
}

/* Initialisation of defaults for global variables only used in experiment*/
void CO2__FP_CritProp_compute_default_exp1::initDefaultsLibraryGlobalsInExp()
{
}

/* Initialisation of variables in experiment*/
void CO2__FP_CritProp_compute_default_exp1::initDefaultsExp()
{
}

/* Initialisation of boundaries*/
void CO2__FP_CritProp_compute_default_exp1::evalBoundsExp(double TIME)
{
	if (m_boundsLaw)
		{(*m_boundsLaw)(TIME); return;}
	if (m_boundaryBranch[0] == 0)  
		unkR[2198] = 0. ;
	if (m_boundaryBranch[1] == 0)  
		unkR[2199] = 0. ;
	if (m_boundaryBranch[2] == 0)  
		unkR[2200] = 0. ;
	if (m_boundaryBranch[3] == 0)  
		unkR[2201] = 0. ;
}
/* Initialisation of delays in experiments*/
void CO2__FP_CritProp_compute_default_exp1::initDelaysExpt()
{
}

#ifndef SIMULA_USE_DECK_SYMBOLS
bool CO2__FP_CritProp_compute_default_exp1::gcs(unsigned int* size, unsigned int* chunkSize, const char*** chunkedStr)
{
return false;
}
#endif //SIMULA_USE_DECK_SYMBOLS

#ifndef SIMULA_USE_GRAPHICAL_MAIN
/* It creates a main program to begin the simulation */
#if SIMULA_PLATFORM == 1
int wmain( int argc, wchar_t * wargv[] )
{
	char **argv = NULL;
	if (!createUtf8argv(argc, wargv, argv))
	{
		argc = 0;
	}
#else
int main( int argc, char * argv[] )
{
#endif
	try
	{
#ifndef SIMULA_USE_DECK_SYMBOLS
		CO2__FP_CritProp_compute_default_exp1 ecomodel;
#else
		CO2__FP_CritProp_compute_default_exp1 ecomodel(0,0,true);
#endif //SIMULA_USE_DECK_SYMBOLS
		mainLoop(argc, argv, &ecomodel);
	}
	catch(...)
	{
		printf("Program stops due to abnormal condition\n");
	}
#if SIMULA_PLATFORM == 1
	releaseUtf8argv(argc, argv);
#endif
	return 0;	
}
#endif //SIMULA_USE_GRAPHICAL_MAIN

/* Function to create an instance of the experiment class */
SIMULA_EXPORT_C_EXP void *fcnExpCreate(const char *name=NULL, const char *dirInstall=NULL, bool d=false, int t=0)
{
	void *ee = 0;
	try
	{
		INTEG_simula *tmp = new CO2__FP_CritProp_compute_default_exp1(name,dirInstall,d);
		ee = (void*)tmp->createThinModel(t,(void(*)())fcnExpCreate);
	}
	catch(...)
	{
	}
	return ee;
}

/* Function to destroy the experiment instance */
SIMULA_EXPORT_C_EXP void fcnExpDel(void *obj,const char* name=NULL)
{
	if ( obj )
	{
		delete obj;
	}
}

/* Function to get information of the model*/
SIMULA_EXPORT_C_EXP void getInformationExperiment(char *experimentName, bool& isDebug)
{
	sprintf(experimentName,"%s","exp1");
	isDebug = false;
}

#ifndef SIMULA_USE_DECK_SYMBOLS
	extern "C" void* getDeck(int t) {return 0;};
	extern "C" void delDeck(void *obj) {return;}
#endif //SIMULA_USE_DECK_SYMBOLS


