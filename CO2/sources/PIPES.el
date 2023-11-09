/*------------------------------------------------------------------------------
 LIBRARY: CO2
 FILE: PIPES
 AUTHOR: Viren Bhanot
 COMPANY: CERN
 DESCRIPTION: Discretized streams and pipe models
 CREATION DATE: 20/05/2017
------------------------------------------------------------------------------*/
USE CRYOLIB
USE MATH
USE PORTS_LIB
USE THERMAL

/*
 * INDEX
 * 1. REFRIGERANT STREAMS
 *	1. DiscretizedVolume
 *	2. CylindricalStream
 *	3. FluidStream_ph
 *	4. FluidStream_phr
 *	6. FluidStream_ru
 *	7. AnnularStream
 *	8. AnnularStream_ph
 *	9. AnnularStream_ru
 * 2. PIPES
 *	1. Pipe_ph
 *	2. Pipe_phr
 *	3. Pipe_ru
 */
 
--------------------------------------------------------------------------------
// REFRIGERANT STREAMS
--------------------------------------------------------------------------------
ABSTRACT COMPONENT DiscretizedVolume (
	INTEGER n = 3 "Number of control volumes"
)
"Base class for discretized fluid streams"
// Does not specify geometry
// Does not assign geometry dependent port variables (eg P_aux)
// Handles refrigerant flow only (no wall calculations)
PORTS
	IN fluid f_in
	OUT fluid f_out
DATA
	REAL L = 1 UNITS u_m "Length (Height for vertical components)"
	REAL z_in = 0 UNITS u_m "Elevation of inlet port wrt user-defined base"
	REAL dz = 0 UNITS u_m "Elevation change, aka z_in - z_out"
	REAL m0 = 0.005 UNITS u_kg_s "Nominal (design) mass flow rate"
	REAL dP0 = 0.01 UNITS u_bar "Nominal (design) pressure drop"
	REAL dP_small = 0.0001 UNITS u_bar "dP below which low mass flow rate smoothing applied"
	ENUM InitialConditions init  = Ph "Specify initial condition variables"
	REAL P0_i = 20 UNITS u_bar "Initial inlet pressure"
	REAL h0_i = 250000 UNITS u_J_kg "Initial inlet specific enthalpy"	
	REAL T0_i = 255 UNITS u_K "Initial inlet temperature"
	REAL x0_i = 0.5 "Initial inlet vapour quality"
	REAL P0_o = 20 UNITS u_bar "Initial outlet pressure"
	REAL h0_o = 250000 UNITS u_J_kg "Initial outlet specific enthalpy"	
	REAL T0_o = 255 UNITS u_K "Initial outlet temperature"
	REAL x0_o = 0.5 "Initial outlet vapour quality"
DECLS
	REAL cp[n] UNITS u_J_kgK "Specific heat capacity"
	REAL delT[n] UNITS u_K "Superheat/Subcooling (positive is superheat)"
	REAL dP UNITS u_bar "Total pressure drop in stream"
	REAL drho_dh[n]
	REAL drho_dP[n]
	REAL gamma[n] "Void fraction"
	REAL h[n] UNITS u_J_kg "Specific enthalpy"
	REAL k[n] UNITS u_W_mK "Thermal conductivity"
	REAL M UNITS u_kg "Refrigerant mass"
	REAL m[n+1] UNITS u_kg_s "Mass flow rate"
	REAL mh[n+1] UNITS u_W "Enthalpy flow rate"
	REAL mu[n] UNITS u_Pas "Dynamic viscosity"
	REAL P[n] UNITS u_bar "Pressure"
	ENUM Phase phase[n] "Fluid phase"
	REAL P_Pa[n] UNITS u_Pa "Pressure, in Pa"
	REAL Pr[n] "Prandtl number"
	REAL rho[n] UNITS u_kg_m3 "Density"
	REAL sigma[n] UNITS u_N_m "Surface tension"
	REAL T[n] UNITS u_K "Temperature"
	REAL Tsat[n] UNITS u_K "Saturation temperature"
	REAL u[n] UNITS u_J "Specific internal energy"
	REAL vel[n] UNITS u_m_s "Fluid velocity"
	REAL Vol UNITS u_dm3 "Internal volume in Litre"
	REAL vsound[n] UNITS u_m_s "Mach number"
	REAL x[n] "Vapour quality"
	DISCR REAL z[n+1] UNITS u_m "Elevation of the control volume interfaces"
	INTEGER ier,ipx,ipy // error codes
	DISCR REAL delz
	REAL dP_static[n] UNITS u_bar
OBJECTS
	SaturationProperties sat[n]
	RefGeometryRecord geo[n]
TOPOLOGY
	PATH f_in TO f_out
INIT
	f_out.is_C = TRUE
	linspace(P0_i,P0_o,n,P,1)
	linspace(z_in,z_in-dz,n+1,z)
	delz = dz/n
	assignInitDiscretized(n,init,f_in.fluid,P0_i,P0_o,h0_i,h0_o,T0_i,T0_o,x0_i,x0_o,P,h,rho,u,x)
DISCRETE
	ASSERT (abs(dz)<=L) FATAL "Elevation change cannot exceed tube length (abs(dz) must be < L)"
CONTINUOUS
	EXPAND_BLOCK (i IN 1,n)
		vel[i] = 0.5*(m[i]+ m[i+1])/(rho[i]*geo[i].Ap) // central difference
		P[i] = P_Pa[i]*1E-5
	END EXPAND_BLOCK
	dP = P[1] - P[n]
	M = SUM(i IN 1,n; geo[i].V*rho[i])
	Vol = SUM(i IN 1,n; geo[i].V*1000)
<mi> 	m[1] = f_in.m
<mo> 	m[n+1] = f_out.m
	//EXPAND (i IN 2,n) m[i] := m0 * regRoot2(P[i-1] - P[i] + (0.5*(rho[i-1]+rho[i]) * 9.81 * (delz))/1e5, dP_small) / sqrt(dP0/n)
	EXPAND (i IN 2,n) m[i] := m0 * regRoot2(P[i-1] - P[i] + (rho[i-1]*9.81*delz)/1e5, dP_small) / sqrt(dP0/n)
	//EXPAND (i IN 2,n) m[i] := m0  * regRoot2(P[i-1]-P[i],dP_small)/sqrt(dP0/n) + sqrt((0.5*(rho[i-1]+rho[i])*9.81*(delz))/1e5)
	
	EXPAND (i IN 1,n) dP_static[i] := rho[i]*9.81*(z[i]-z[i+1])/1e5
	EXPAND (i IN 1,n) delT[i] = T[i] - Tsat[i]

	f_in.fluid = f_out.fluid
	f_in.n_fluid = f_out.n_fluid
	f_in.Gcrit = CRYO_Gcrit_fun(f_in.fluid,P[1]*1e5,rho[1],T[1],0,f_in.v,vsound[1],ier,ipy)
<hb> 	f_in.hb = h[1]
<Pi> 	f_in.P = P[1] //- 0.5*rho[1]*9.81*(z[1]-z[2])/1e5 // add static pressure if z2>z1
	f_in.rho = rho[1]
	f_out.Gcrit = CRYO_Gcrit_fun(f_in.fluid,P[n]*1e5,rho[n],T[n],0,f_out.v,vsound[n],ier,ipy)
<hf> 	f_out.hf = h[n]
<Po>	f_out.P = P[n] //+ 0.5*rho[n]*9.81*(z[n]-z[n+1])/1e5
	f_out.rho = rho[n]
END COMPONENT



--------------------------------------------------------------------------------
// CYLINDRICAL STREAMS
--------------------------------------------------------------------------------
ABSTRACT COMPONENT CylindricalStream IS_A DiscretizedVolume (
	BOOLEAN ConstantRefHTC = TRUE "If FALSE, use correlations",
	BOOLEAN DynamicHTC = TRUE "If TRUE, apply smoothing between correlation HTCs"
)
"Adds thermal ports (one per control volume) and htc calculations to DiscretizedVolume"
// does not add governing equations
PORTS
	OUT thermal(n=n) tp_out // n nodes
DATA
	REAL D = 0.01 UNITS u_m "Diameter"
	ENUM SinglePhaseHTC choice_1p = UseDittusBoelter
	ENUM TwoPhaseHTC choice_2p = UseChen
	REAL alpha_l = 500 UNITS u_W_m2K "Liquid HTC, for ConstantHTC==TRUE"
	REAL alpha_tp = 1000 UNITS u_W_m2K "2Phase HTC, for ConstantHTC==TRUE"
	REAL alpha_g = 400 UNITS u_W_m2K "Vapour HTC, for ConstantHTC==TRUE"
	REAL m_steady = 0.05 UNITS u_kg_s "Steady-state mass flow rate for which HTC calculated"
	REAL tauAlpha = 3 UNITS u_s "Time constant for heat transfer coefficient calculations"
DECLS
	DISCR REAL Pcrit UNITS u_bar
	DISCR REAL MW "Molar mass"
	REAL alphaHat[n] UNITS u_W_m2K "HTC calculated from correlations"
	REAL alpha[n] UNITS u_W_m2K "Low pass filtered HTC"
	REAL Q[n] UNITS u_W "Heat transfer rate"
	REAL q[n] UNITS u_W_m2 "Heat flux"
	REAL Qr_total UNITS u_W "Overall heat transfer"
	REAL m_avg[n] UNITS u_kg_s "Averaged mass flow rate for HTC calculations"
	REAL Pred[n] UNITS u_bar "Reduced pressure"
INIT
	FOR (i IN 1,n)
		geo[i].EndSurfaces = 0
		geo[i].setCylindricalGeometry(D,L/n)
		IF (x[i]<=0) THEN
			alpha[i] = alpha_l
			alphaHat[i] = alpha_l
		ELSEIF(x[i]>=1) THEN
			alpha[i] = alpha_g
			alphaHat[i] = alpha_g
		ELSE
			alpha[i] = alpha_tp
			alphaHat[i] = alpha_tp
		END IF
		P_Pa[i] = P[i]*1e5
	END FOR
	Pcrit = 73.77 // CRYO_PF_CritProp(f_in.fluid,fprop_pressure,ier)
	MW = getMolWt(f_in.fluid)
CONTINUOUS
	EXPAND_BLOCK (i IN 1,n)
		m_avg[i] = 0.5*(m[i] + m[i+1])
		Pred[i] = P[i]/Pcrit
		q[i] = abs(alphaHat[i]*(T[i] - tp_out.Tk[i]))
	END EXPAND_BLOCK
<alpha>	EXPAND_BLOCK (i IN 1,n)
		IF (ConstantRefHTC) INSERT
			alphaHat[i] = max(0,getSmoothHTC(alpha_g,alpha_l,alpha_tp,m_avg[i],m_steady,x[i]))
		ELSE
			alphaHat[i] = getCorrelationHTC(f_in.fluid,choice_1p,choice_2p,geo[i].Ap,cp[i],D,k[i],m_avg[i],mu[i],MW,P[i],Pred[i],q[i],sigma[i],T[i],tp_out.Tk[i],x[i],sat[i])
		END IF
		IF (DynamicHTC) INSERT
			alpha[i]' = 1/tauAlpha * (alphaHat[i] - alpha[i])
		ELSE
			alpha[i] = alphaHat[i]
		END IF
	END EXPAND_BLOCK
<Q>	EXPAND_BLOCK (i IN 1,n)
		Q[i] = alpha[i]*geo[i].As*(T[i]-tp_out.Tk[i]) // see note in THERMAL file on sign convention
		tp_out.q[i] = Q[i]
	END EXPAND_BLOCK
	Qr_total = SUM(i IN 1,n; Q[i])

	f_in.A = geo[1].Ap
	f_in.P_aux = P[1] + geo[1].V*rho[1]*vsound[1]*1E-5/geo[1].As
	f_in.I = 0.5*geo[1].V/geo[1].As**2
	f_in.v = f_in.m/(rho[1]*geo[1].Ap)
	f_out.A = geo[n].Ap
	f_out.P_aux = P[n] + geo[n].V*rho[n]*vsound[n]*1E-5/geo[n].As
	f_out.I = 0.5*geo[n].V/geo[n].As**2
	f_out.v = f_out.m/(rho[n]*geo[n].Ap)
END COMPONENT



COMPONENT FluidStream_Ph IS_A CylindricalStream (
	ENUM VoidFractionModel vfModel = Homogeneous "Void fraction model for two-phase flow"
)
"Ph state variables. No tube wall, only refrigerant stream"
// Also adds gravity pressure drop and two-phase flow models
DECLS
	REAL h_corr[n] UNITS u_J_kg "Correction term for two-phase port enthalpy"
	REAL h_flow[n] UNITS u_J_kg "Flow-weighted enthalpy for two-phase flow"
CONTINUOUS
	// Two-phase flow slip-ratio based correction
	EXPAND_BLOCK (i IN 1,n)
		h_corr[i] = getCorrectionEnthalpy(vfModel,h[i],rho[i],sat[i])
		h_flow[i] = h[i] + h_corr[i]
	END EXPAND_BLOCK

<mhi>	mh[1] = m[1]*donor_cell(m[1],f_in.hf,h_flow[1])
<mho>	mh[n+1] = m[n+1]*donor_cell(m[n+1],h_flow[n],f_out.hb)
<mh>	EXPAND(i IN 2,n) mh[i] = m[i]*donor_cell(m[i],h_flow[i-1],h_flow[i])
<state> EXPAND (i IN 1,n) getStatePh(f_in.fluid,P[i],h[i],cp[i],drho_dP[i],drho_dh[i],gamma[i],k[i],mu[i],phase[i],Pr[i],rho[i],sigma[i],T[i],Tsat[i],u[i],vsound[i],x[i],sat[i],ier,ipx,ipy)
	
<gov>	EXPAND_BLOCK (i IN 1,n)
		geo[i].V*(drho_dP[i]*P_Pa[i]'+drho_dh[i]*h[i]') = m[i] - m[i+1]
		geo[i].V*((h[i]*drho_dP[i]-1)*P_Pa[i]' + (h[i]*drho_dh[i]+rho[i])*h[i]') = mh[i] - mh[i+1] - Q[i] // - 9.81*(z[i]-z[i+1])*rho[i]*vel[i]*geo[i].V/L
	END EXPAND_BLOCK
	
<:hb> 	f_in.hb = h_flow[1]
<:hf> 	f_out.hf = h_flow[n]
END COMPONENT



COMPONENT FluidStream_Phr IS_A FluidStream_Ph
"Phr state variables. No tube wall, only refrigerant stream"
DECLS
	REAL dh_dP[n]
	REAL dh_drho[n]
	PRIVATE REAL tmp[n] // dummy density
CONTINUOUS
<:gov>	EXPAND_BLOCK (i IN 1,n)
		dh_drho[i] = 1/drho_dh[i]
		dh_dP[i] = -drho_dP[i]/drho_dh[i]
		geo[i].V*rho[i]' = m[i] - m[i+1]
		geo[i].V*((rho[i]*dh_dP[i]-1)*P_Pa[i]' + (dh_drho[i]+h[i])*rho[i]') - mh[i] + mh[i+1] + Q[i] = 0 // TODO: add gravity term
		h[i]' = dh_dP[i]*P_Pa[i]' + dh_drho[i]*rho[i]'
	END EXPAND_BLOCK
<:state> EXPAND (i IN 1,n) getStatePh(f_in.fluid,P[i],h[i],cp[i],drho_dP[i],drho_dh[i],gamma[i],k[i],mu[i],phase[i],Pr[i],tmp[i],sigma[i],T[i],Tsat[i],u[i],vsound[i],x[i],sat[i],ier,ipx,ipy)
END COMPONENT



COMPONENT FluidStream_ru IS_A FluidStream_Ph
"rho-u state variables. No thermal wall"
CONTINUOUS
<:gov>	EXPAND_BLOCK (i IN 1,n)
		geo[i].V*rho[i]' = m[i] - m[i+1]
		geo[i].V*(rho[i]*u[i]'+u[i]*rho[i]') = mh[i] - mh[i+1] - Q[i] // TODO: add gravity term
	END EXPAND_BLOCK
<:state> EXPAND (i IN 1,n) getStateRU(f_in.fluid,rho[i],u[i],cp[i],drho_dP[i],drho_dh[i],gamma[i],h[i],k[i],mu[i],P[i],phase[i],Pr[i],sigma[i],T[i],Tsat[i],vsound[i],x[i],sat[i],ier,ipx,ipy)
END COMPONENT



// dP/dt CONSTANT
COMPONENT FluidStream_dpdt IS_A FluidStream_Ph
"Uniform pressure derivative cylindrical stream"
// REFERENCE: Qiao and Laughman, 2018, Comparison of  Approximate Momentum Equations in Dynamic Models of Vapor Compression Systems
DECLS
	REAL dP_dt "Pressure derivative"
INIT
	dP_dt = 0
CONTINUOUS
	dP_dt = P[1]'
<:gov>	EXPAND_BLOCK (i IN 1,n)
		geo[i].V*(drho_dP[i]*dP_dt+drho_dh[i]*h[i]') = m[i] - m[i+1]
		geo[i].V*((h[i]*drho_dP[i]-1)*dP_dt + (h[i]*drho_dh[i]+rho[i])*h[i]') = mh[i] - mh[i+1] // TODO: add gravity and heat source terms
	END EXPAND_BLOCK
END COMPONENT



--------------------------------------------------------------------------------
// ANNULAR STREAMS
--------------------------------------------------------------------------------
ABSTRACT COMPONENT AnnularStream IS_A DiscretizedVolume (
	BOOLEAN ConstantRefHTC = TRUE "False if using correlations",
	BOOLEAN DynamicHTC = TRUE "If TRUE, applies a low-pass filter to HTC"
)
"Connects to thermal port on both its inner and outer surfaces"
PORTS
	IN thermal(n=n) tp_in "Annulus inner wall"
	OUT thermal(n=n) tp_out "Annulus outer wall"
DATA
	REAL Do = 0.018 UNITS u_m "Outer diameter"
	REAL Di = 0.012 UNITS u_m "Inner diameter"
	ENUM SinglePhaseHTC choice_1p = UseDittusBoelter
	ENUM TwoPhaseHTC choice_2p = UseChen
	REAL tauAlpha = 3 UNITS u_s "Time constant for low-pass HTC filter"
	REAL alpha_l = 1000 UNITS u_W_m2K "Liquid HTC"
	REAL alpha_tp = 1000 UNITS u_W_m2K "Two Phase HTC"
	REAL alpha_g = 1000 UNITS u_W_m2K "Vapour HTC"
	REAL m_steady = 0.05 UNITS u_kg_s "Steady state mass flow rate for which HTCs calculated"
DECLS
	DISCR REAL Pcrit UNITS u_bar
	DISCR REAL MW "Molar Weight"
	REAL alpha_i[n], alpha_o[n] UNITS u_W_m2K
	REAL alphaHat_i[n], alphaHat_o[n] UNITS u_W_m2K
	REAL Q_i[n] UNITS u_W "Heat transfer at inner surface"
	REAL Q_o[n] UNITS u_W "Heat transfer at outer surface"
	REAL q_i[n] UNITS u_W_m2 "Heat flux on inner surface"
	REAL q_o[n] UNITS u_W_m2 "Heat flux on outer surface"
	REAL Qi_total, Qo_total, Qr_total UNITS u_W
	PRIVATE REAL m_avg[n] UNITS u_kg_s
	PRIVATE REAL Pred[n] UNITS u_bar "Reduced pressure"
	DISCR REAL tmp1
INIT
	FOR (i IN 1,n)
		geo[i].EndSurfaces = 0
		geo[i].setAnnularGeometry(Do,(Do-Di)*0.5,L/n)
		IF (x[i]<=0) THEN
			alpha_i[i] = alpha_g
			alpha_o[i] = alpha_g
			alphaHat_i[i] = alpha_g
			alphaHat_o[i] = alpha_g
		ELSEIF (x[i]>=1) THEN
			alpha_i[i] = alpha_l
			alpha_o[i] = alpha_l
			alphaHat_i[i] = alpha_l
			alphaHat_o[i] = alpha_l
		ELSE
			alpha_i[i] = alpha_tp
			alpha_o[i] = alpha_tp
			alphaHat_i[i] = alpha_tp
			alphaHat_o[i] = alpha_tp
		END IF
		P_Pa[i] = P[i]*1e5
	END FOR
	Pcrit = 73.77 //CRYO_PF_CritProp(f_in.fluid,fprop_pressure,ier)
	MW = getMolWt(f_in.fluid)
	tmp1 = geo[1].As
DISCRETE
	ASSERT (Do > Di) FATAL "Inner diameter must be smaller than outer diameter"
CONTINUOUS
	EXPAND_BLOCK (i IN 1,n)
		m_avg[i] = 0.5*(abs(m[i]) + abs(m[i+1]))
		Pred[i] = P[i]/Pcrit
		q_i[i] = abs(alphaHat_i[i]*(tp_in.Tk[i]-T[i]))
		q_o[i] = abs(alphaHat_o[i]*(T[i] - tp_out.Tk[i]))
	END EXPAND_BLOCK
<alpha>	EXPAND_BLOCK (i IN 1,n)
		IF (ConstantRefHTC) INSERT
			alphaHat_i[i] = getSmoothHTC(alpha_g,alpha_l,alpha_tp,m_avg[i],m_steady,x[i])
			alphaHat_o[i] = alphaHat_i[i]
		ELSE
			alphaHat_i[i] = getCorrelationHTC(f_in.fluid,choice_1p,choice_2p,geo[i].Ap,cp[i],Do-Di,k[i],m_avg[i],mu[i],MW,P[i],Pred[i],q_i[i],sigma[i],T[i],tp_in.Tk[i],x[i],sat[i])
			alphaHat_o[i] = getCorrelationHTC(f_in.fluid,choice_1p,choice_2p,geo[i].Ap,cp[i],Do-Di,k[i],m_avg[i],mu[i],MW,P[i],Pred[i],q_o[i],sigma[i],T[i],tp_out.Tk[i],x[i],sat[i])
		END IF
		IF (DynamicHTC) INSERT
			alpha_i[i]' = 1/tauAlpha * (alphaHat_i[i] - alpha_i[i])
			alpha_o[i]' = 1/tauAlpha * (alphaHat_o[i] - alpha_o[i])
		ELSE
			alpha_i[i] = alphaHat_i[i]
			alpha_o[i] = alphaHat_o[i]
		END IF
	END EXPAND_BLOCK
<Q>	EXPAND_BLOCK (i IN 1,n)
		Q_i[i] = alpha_i[i]*geo[i].As_i*(tp_in.Tk[i]-T[i]) //positive when heat flowing in
		Q_o[i] = alpha_o[i]*geo[i].As_o*(T[i]-tp_out.Tk[i]) // positive when heat flowing out
		tp_in.q[i] = Q_i[i]
		tp_out.q[i] = Q_o[i]
	END EXPAND_BLOCK
	Qi_total = SUM(i IN 1,n; Q_i[i])
	Qo_total = SUM(i IN 1,n; Q_o[i])
	Qr_total = Qi_total + Qo_total

 	f_in.A = geo[1].Ap
	f_in.P_aux = P[1] + geo[1].V*rho[1]*vsound[1]*1E-5/geo[1].As
	f_in.I = 0.5*geo[1].V/geo[1].As**2
	f_in.v = f_in.m/(rho[1]*geo[1].Ap)
 	f_out.A = geo[n].Ap
	f_out.P_aux = P[n] + geo[n].V*rho[n]*vsound[n]*1E-5/geo[n].As
	f_out.I = 0.5*geo[n].V/geo[n].As**2
	f_out.v = f_out.m/(rho[n]*geo[n].Ap)
END COMPONENT



COMPONENT AnnularStream_Ph IS_A AnnularStream (
	ENUM VoidFractionModel vfModel = Smith "Slip-ratio based void fraction model"
)
"Annular stream with pressure-enthalpy state variables and void fraction two-phase models"
DECLS
	REAL h_corr[n] UNITS u_J_kg "Difference between flow-weighted and density-weighted enthalpy"
	REAL h_flow[n] UNITS u_J_kg "Flow-weighted enthalpy"
CONTINUOUS
	EXPAND_BLOCK (i IN 1,n)
		h_corr[i] = getCorrectionEnthalpy(vfModel,h[i],rho[i],sat[i])
		h_flow[i] = h[i] + h_corr[i]
	END EXPAND_BLOCK
	mh[1] = m[1]*donor_cell(m[1],f_in.hf,h_flow[1])
	mh[n+1] = m[n+1]*donor_cell(m[n+1],h_flow[n],f_out.hb)
	EXPAND(i IN 2,n) mh[i] = m[i]*donor_cell(m[i],h_flow[i-1],h_flow[i])

<gov>	EXPAND_BLOCK (i IN 1,n)
		geo[i].V*(drho_dP[i]*P_Pa[i]'+drho_dh[i]*h[i]') = m[i] - m[i+1]
		geo[i].V*((h[i]*drho_dP[i]-1)*P_Pa[i]' + (h[i]*drho_dh[i]+rho[i])*h[i]') = mh[i] - mh[i+1] + Q_i[i] - Q_o[i]
	END EXPAND_BLOCK
<state>	EXPAND (i IN 1,n) getStatePh(f_in.fluid,P[i],h[i],cp[i],drho_dP[i],drho_dh[i],\
			gamma[i],k[i],mu[i],phase[i],Pr[i],rho[i],sigma[i],
			T[i],Tsat[i],u[i],vsound[i],x[i],sat[i],ier,ipx,ipy)
<:hb>	f_in.hb = h_flow[1]
<:hf>	f_out.hf = h_flow[n]
END COMPONENT



COMPONENT AnnularStream_ru IS_A AnnularStream_Ph
"rho,u state variables"
CONTINUOUS
<:gov>	EXPAND_BLOCK (i IN 1,n)
		geo[i].V*rho[i]' = m[i] - m[i+1]
		geo[i].V*(rho[i]*u[i]'+u[i]*rho[i]') = mh[i]-mh[i+1] + Q_i[i]-Q_o[i]
	END EXPAND_BLOCK
<:state> EXPAND (i IN 1,n) getStateRU(f_in.fluid,rho[i],u[i],cp[i],drho_dP[i],drho_dh[i],\
		gamma[i],h[i],k[i],mu[i],P[i],phase[i],Pr[i],sigma[i],
		T[i],Tsat[i],vsound[i],x[i],sat[i],ier,ipx,ipy)
END COMPONENT



--------------------------------------------------------------------------------
// PIPES
--------------------------------------------------------------------------------
ABSTRACT COMPONENT Pipe (
	INTEGER n = 3 "Number of control volumes",
	BOOLEAN ConstantAirHTC = TRUE "If FALSE, uses ChurchillChu natural convective correlation",
	BOOLEAN ConstantRefHTC = TRUE "False if using correlations",
	BOOLEAN DynamicRefHTC = TRUE "If TRUE, applies a low-pass filter to ref HTC"
)
"Only handles DATA variables. TOPOLOGY left to child components"
// This abstract component does not have any physical significance.
// Merely a coding convenience so that I don't have to duplicate code.
PORTS
	IN fluid f_in
	OUT fluid f_out
	IN analog_signal(n=1) T_amb
DATA
	ENUM Material mat = SS_304
	REAL Do = 0.01 UNITS u_m "Outer diameter"
	REAL Thw = 0.001 UNITS u_m "Wall thickness"
	REAL L = 1 UNITS u_m "Length"
	REAL z_in = 0 UNITS u_m "Elevation of the inlet port wrt some user-defined base"
	REAL dz = 0 UNITS u_m "Elevation change between inlet and outlet port, z_in - z_out"

	ENUM SinglePhaseHTC choice_1p = UseDittusBoelter
	ENUM TwoPhaseHTC choice_2p = UseChen
	REAL tau = 3 UNITS u_s "Time constant for Dynamic HTC"
	REAL alpha_l = 1000 UNITS u_W_m2K "Liquid HTC, for ConstantHTC==TRUE"
	REAL alpha_tp = 1000 UNITS u_W_m2K "Two Phase HTC, for ConstantHTC==TRUE"
	REAL alpha_g = 1000 UNITS u_W_m2K "Vapour HTC, for ConstantHTC==TRUE"
	REAL m_steady = 0.05 UNITS u_kg_s "Steady state refrigerant mass flow rate at which constant HTCs calculated"
	REAL m0 = 0.005 UNITS u_kg_s "Nominal (design) mass flow rate"
	REAL dP0 = 0.01 UNITS u_bar "Nominal (design) pressure drop"
	REAL alpha_a0 = 50 UNITS u_W_m2K "Air side constant heat transfer coefficient"

	ENUM InitialConditions init = Ph "Specify initial condition variables"
	REAL P0_i = 20 UNITS u_bar "Initial inlet pressure"
	REAL h0_i = 250000 UNITS u_J_kg "Initial inlet specific enthalpy"	
	REAL x0_i = 0.5 "Initial inlet vapour quality"
	REAL T0_i = 255 UNITS u_K "Initial inlet temperature"
	REAL P0_o = 20 UNITS u_bar "Initial outlet pressure"
	REAL h0_o = 250000 UNITS u_J_kg "Initial outlet specific enthalpy"	
	REAL x0_o = 0.5 "Initial outlet vapour quality"
	REAL T0_o = 255 UNITS u_K "Initial outlet temperature"
	REAL Tw0_i = 245 UNITS u_K "Initial wall temperature - first segment"
	REAL Tw0_o = 245 UNITS u_K "Initial wall temperature - last segment"
TOPOLOGY
	AnnularWall_naturalConvection (
		n=n,
		ConstantAirHTC=ConstantAirHTC,
		accountForRadiation=FALSE
	) Wall (
		mat = mat,
		Do = Do,
		Th = Thw,
		L = L,
		alpha_a0 = alpha_a0,
		Tw0_i = Tw0_i,
		Tw0_o = Tw0_o
	)	
	CONNECT T_amb TO Wall.T_amb
END COMPONENT



COMPONENT Pipe_Ph IS_A CO2.Pipe (
	ENUM VoidFractionModel vfModel = Homogeneous "Two-phase slip-ratio based void-fraction model"
)
"Insulated pipe. P,h state variables"
TOPOLOGY
	FluidStream_Ph (
		n=n,
		ConstantRefHTC=ConstantRefHTC,
		DynamicHTC=DynamicRefHTC,
		vfModel = vfModel
	) Stream (
		D = Do-2*Thw,
		L = L,
		dz = dz,
		z_in = z_in,
		init = init,
		P0_i = P0_i,
		h0_i = h0_i,
		x0_i = x0_i,
		T0_i = T0_i,
		P0_o = P0_o,
		h0_o = h0_o,
		T0_o = T0_o,
		x0_o = x0_o,
		choice_1p = choice_1p,
		choice_2p = choice_2p,
		alpha_l = alpha_l,
		alpha_tp = alpha_tp,
		alpha_g = alpha_g,
		m_steady = m_steady,
		m0 = m0,
		dP0 = dP0
	)
	CONNECT f_in TO Stream
	CONNECT Stream TO f_out
	CONNECT Stream.tp_out TO Wall.tp_in
END COMPONENT



COMPONENT Pipe_Phr IS_A CO2.Pipe
"Insulated pipe. P,h,rho state variables"
TOPOLOGY
	FluidStream_Phr (
		n = n,
		ConstantRefHTC = ConstantRefHTC
	) Stream (
		D = Do-2*Thw,
		L = L,
		dz = dz,
		z_in = z_in,
		init = init,
		P0_i = P0_i,
		h0_i = h0_i,
		x0_i = x0_i,
		T0_i = T0_i,
		P0_o = P0_o,
		h0_o = h0_o,
		T0_o = T0_o,
		x0_o = x0_o,
		choice_1p = choice_1p,
		choice_2p = choice_2p,
		alpha_l = alpha_l,
		alpha_tp = alpha_tp,
		alpha_g = alpha_g,
		m_steady = m_steady,
		m0 = m0,
		dP0 = dP0
	)
	CONNECT f_in TO Stream
	CONNECT Stream TO f_out
	CONNECT Stream.tp_out TO Wall.tp_in
END COMPONENT




COMPONENT Pipe_ru IS_A CO2.Pipe
"1D Pipe. rho,u state variables."
TOPOLOGY
	FluidStream_ru (
		n = n,
		ConstantRefHTC = ConstantRefHTC
	) Stream (
		D = Do-2*Thw,
		L = L,
		dz = dz,
		z_in = z_in,
		init = init,
		P0_i = P0_i,
		h0_i = h0_i,
		x0_i = x0_i,
		T0_i = T0_i,
		P0_o = P0_o,
		h0_o = h0_o,
		T0_o = T0_o,
		x0_o = x0_o,
		choice_1p = choice_1p,
		choice_2p = choice_2p,
		alpha_l = alpha_l,
		alpha_tp = alpha_tp,
		alpha_g = alpha_g,
		m_steady = m_steady,
		m0 = m0, 
		dP0 = dP0
	)
	CONNECT f_in TO Stream
	CONNECT Stream TO f_out
	CONNECT Stream.tp_out TO Wall.tp_in
END COMPONENT