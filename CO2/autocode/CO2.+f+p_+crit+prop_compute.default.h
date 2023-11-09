//==================================================================
// Code generated automatically
// Description: Partition class header
//==================================================================
#ifndef CO2__FP_CritProp_compute_default_H
#define CO2__FP_CritProp_compute_default_H

#include <INTEG_simula.h>

class CO2__FP_CritProp_compute_default: public INTEG_simula
{
public:
	// methods to cope with the continuous part
	CO2__FP_CritProp_compute_default(const char *mgr=NULL,const char *dirInstall=NULL,bool dmode=false);
	virtual void copyBack( double dyn[], double ldr[] );
	virtual void initBlocks( double dyn[], double ldr[], double *_time );
	virtual void fres( double *_time, double dyn[], double der[], double res[] );

	// methods to cope with the discrete part
	virtual void checkAsserts( double *_time);
	virtual void constraints( double *_time, double ev[],double dyn[],double ldr[]);
	virtual void evalWhen(double *_time, bool w[],bool cont[] );
	virtual void executeWhen(double *_time,  int index );
	virtual void evalZones( double *_time,  int branchZone[],bool cont[] );
	virtual INTEG_simula::t_initEvent* initEvents(int& nEvents,int& nWhen,int& nZones,int& nConstraints, const char**& whenTxt,const char**& zoneTxt,int *&zoneTxtIndex);
	virtual void initDelays();
	virtual void initInternalModels();

	// Pointer use for numerical wrapper
	static CO2__FP_CritProp_compute_default* s_current;
public: 
	//EL functions declaration
	double CO2__CRYO_PF_CritProp_CORR(const int & chem,const int & fprop,int * ier);
	double CRYOLIB__HEPAK_He_X_vs_JN1_JN2(const int & JN1,const int & JN2,const double & X1,const double & X2);
	double MATH__max(const double & x,const double & y);
	double CRYOLIB__HEPAK_He_prop_vs_ph(const double & p,const double & h,const int & n_prop);
	double CRYOLIB__CRYO_PF_prop_vs_ph(const int & chem,const double & p_bar,const double & h,const int & fprop,int * ier,int * jx,int * jy);

	virtual int sparseGetEquiv( int nvars,int v[] );

	 bool gst(unsigned int* size, unsigned int* chunkSize, const char*** chunkedStr);

	 bool gsi(unsigned int* size, unsigned int* chunkSize, const char*** chunkedStr);

	virtual bool gcs(unsigned int* size, unsigned int* chunkSize, const char*** chunkedStr);
};
#endif
