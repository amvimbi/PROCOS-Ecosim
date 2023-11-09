//==================================================================
// Code generated automatically
// Description: Experiment class header 
//==================================================================
#ifndef CO2__FP_CritProp_compute_default_exp1_H
#define CO2__FP_CritProp_compute_default_exp1_H
#include "CO2.+f+p_+crit+prop_compute.default.h"

class CO2__FP_CritProp_compute_default_exp1: public CO2__FP_CritProp_compute_default
{
public:
	CO2__FP_CritProp_compute_default_exp1(const char* mgr= NULL,const char* installDir=NULL,bool dmode=false);
public:
	virtual void runExperiment();
	virtual void addExptVariables();
	virtual void initDefaultsLibraryGlobalsInExp();
	virtual void initDefaultsExp();
	virtual void evalBoundsExp(double TIME);
	virtual void initDelaysExpt();
	virtual void initInternalModels();
	virtual bool gcs(unsigned int* size, unsigned int* chunkSize, const char*** chunkedStr);

	static CO2__FP_CritProp_compute_default_exp1* s_current;
private:


// Experiment variables

// experiment functions

};
#endif
