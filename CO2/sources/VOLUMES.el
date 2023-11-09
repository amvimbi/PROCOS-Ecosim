/*-------------------------------------------------------------------------------
 LIBRARY: CO2
 FILE: VOLUMES
 AUTHOR: Viren Bhanot
 COMPANY: CERN
 DESCRIPTION: Lumped capacitance models
 CREATION DATE: 10/10/2016
-------------------------------------------------------------------------------*/
USE CRYOLIB
USE MATH
USE PORTS_LIB
USE THERMAL



/*
 * INDEX
 * 
 * 1. BASE COMPONENTS
 *	1. LumpedAdiabaticVolume
 *	2. LumpedVolume_heatTransfer
 * 2. ADIABATIC LUMPED VOLUMES
 *	1. AdiabaticVol_Ph
 *	2. AdiabaticVol_Ph2P
 *	3. AdiabaticVol_Phr
 *	4. AdiabaticVol_ru
 * 3. LUMPED VOLUMES WITH HEAT TRANSFER
 *	1. CV_Ph
 *	2. CV_Ph2P
 *	3. CV_Phr
 *	4. CV_ru
 */



--------------------------------------------------------------------------------
// BASE COMPONENTS
--------------------------------------------------------------------------------
ABSTRACT COMPONENT LumpedAdiabaticVolume
"Base class for lumped adiabatic cylindrical control volumes"
// Handles geometry, initial conditions and port parameters
// Does not handle governing equations or heat transfer calculations
PORTS
	IN fluid f_in
	OUT fluid f_out
DATA
	REAL D = 0.01 UNITS u_m RANGE Eps,Inf "Diameter"
	REAL L = 1 UNITS u_m RANGE Eps,Inf "Length"
	REAL z_in = 0 UNITS u_m "Elevation of inlet port relative to user-defined base"
	REAL dz = 0 UNITS u_m "Elevation change between inlet and outlet, z_in - z_out, negative for upward pipe"
	ENUM InitialConditions init = Ph "Initial condition specification"
	REAL P0 = 20 UNITS u_bar "Initial pressure"
	REAL h0 = 250000 UNITS u_J_kg "Initial specific enthalpy"
	REAL T0 = 255 UNITS u_K "Initial temperature"
	REAL x0 = 0.5 UNITS no_units "Initial vapour quality"
DECLS
	DISCR REAL Vol UNITS u_dm3 "Internal volume (Liters)"
	REAL cp UNITS u_J_kgK "Specific heat capacity"
	REAL drho_dh "Partial derivative of density wrt specific enthalpy"
	REAL drho_dP "Partial derivative of density wrt pressure"
	REAL gamma "Void fraction"
	REAL G UNITS u_kg_sm2 "Mass flux"
	REAL h UNITS u_J_kg "Specific enthalpy"
	REAL k UNITS u_W_mK "Thermal conductivity"
	REAL M UNITS u_kg "Refrigerant mass"
	REAL mu UNITS u_Pas "Dynamic viscosity"
	REAL P UNITS u_bar "Pressure"
	REAL Pr "Prandtl"
	REAL rho UNITS u_kg_m3 "Density"
	REAL sigma UNITS u_N_m "Surface tension"
	REAL T UNITS u_K "Temperature"
	REAL Tsat UNITS u_K "Saturation temperature"
	REAL delT UNITS u_K "Superheating/Subcooling (positive is superheat)"
	REAL u UNITS u_J_kg "Specific internal energy"
	REAL vsound UNITS u_m_s "Mach number"
	REAL vel UNITS u_m_s "Fluid velocity"
	REAL x "Vapour quality"
	ENUM Phase phase
	PRIVATE INTEGER ier,ipx,ipy // error codes
OBJECTS
	RefGeometryRecord geo
	SaturationProperties sat
TOPOLOGY
	PATH f_in TO f_out
INIT
	geo.EndSurfaces = 0 // Volume is open on both ends
	geo.setCylindricalGeometry(D,L)
	Vol = geo.V*1000 // Liters
	f_out.is_C = TRUE // Outlet port is *Capacitive* (must connect to resistive component at outlet)
	assignInitLumped(init,f_in.fluid,P0,h0,T0,x0,P,h,rho,u)
DISCRETE
	ASSERT (abs(dz)<=L) FATAL "Elevation change (dz) cannot be greater than length of volume (L)"
CONTINUOUS
	vel = 0.5*(f_in.m + f_out.m) / (rho*geo.Ap)
	G = 0.5*(f_in.m + f_out.m) / (geo.Ap)
<M> 	M = rho*geo.V
	delT = T - Tsat
	// PORT ASSIGNMENTS	
	f_in.fluid = f_out.fluid
	f_in.n_fluid = f_out.n_fluid
	f_in.A = geo.Ap
	f_in.Gcrit = CRYO_Gcrit_fun(f_in.fluid,P*1E5,rho,T,0,f_in.v,vsound,ier,ipy)
<hb> 	f_in.hb = h - 0.5*9.8*dz
	f_in.I = 0.5*geo.V/geo.As**2
<Pi> 	f_in.P = P - 0.5*rho*9.81*dz/1e5 
<Pia> 	f_in.P_aux = P + geo.V*rho*vsound*1e-5/geo.As
	f_in.rho = rho
	f_in.v = f_in.m/(rho*geo.Ap)
	f_out.A = geo.Ap
	f_out.Gcrit = CRYO_Gcrit_fun(f_in.fluid,P*1e5,rho,T,0,f_out.v,vsound,ier,ipy)
<hf> 	f_out.hf = h + 0.5*9.8*dz 
	f_out.I = 0.5*geo.V/geo.As**2
<Po> 	f_out.P = P + 0.5*rho*9.81*dz/1e5 
<Pao> 	f_out.P_aux = P + geo.V*rho*vsound*1e-5/geo.As
	f_out.rho = rho
	f_out.v = f_out.m/(rho*geo.Ap)	
END COMPONENT



ABSTRACT COMPONENT LumpedVolume_heatTransfer IS_A LumpedAdiabaticVolume (
	BOOLEAN ConstantHTC = TRUE "False if using correlations",
	BOOLEAN DynamicHTC = TRUE "If TRUE, applies a low-pass filter to HTC calculations - recommended with FALSE ConstantHTC"
)
"Adds heat transfer calculations to LumpedAdiabaticVolume"
// still no governing equations
// alphaHat is the correlation value, alpha is the low-pass filtered final value
PORTS
	OUT thermal(n=1) tp_out
DATA
	ENUM SinglePhaseHTC choice_1p = UseDittusBoelter
	ENUM TwoPhaseHTC choice_2p = UseChen
	REAL tauAlpha = 3 UNITS u_s "Time constant for dynamic HTC calculations"
	REAL alpha_l = 750 UNITS u_W_m2K "Liquid HTC, for ConstantHTC==TRUE"
	REAL alpha_tp = 1000 UNITS u_W_m2K "Two Phase HTC, for ConstantHTC==TRUE"
	REAL alpha_g = 500 UNITS u_W_m2K "Vapour HTC, for ConstantHTC==TRUE"
	REAL m_steady = 0.05 UNITS u_kg_s "Steady state mass flow rate for which HTCs calculated"
DECLS
	REAL Q UNITS u_W "Heat transfer rate"
	REAL q UNITS u_W_m2 "Heat flux"
	REAL alpha UNITS u_W_m2K "Low-pass filtered HTC value"
	REAL alphaHat UNITS u_W_m2K "HTC calculated from correlation"
	REAL m_avg UNITS u_kg_s "Average mass flow rate"
	REAL MW "Molecular weight"
	REAL P_red "Reduced pressure"
	DISCR REAL Pcrit UNITS u_bar
INIT
	Pcrit = 73.77 //CRYO_PF_CritProp(f_in.fluid,fprop_pressure,ier)

	IF (x0<=0) THEN
		alpha = alpha_l
		alphaHat = alpha_l
	ELSEIF (x0>=1) THEN
		alpha = alpha_g
		alphaHat = alpha_l
	ELSE
		alpha = alpha_tp
		alphaHat = alpha_tp
	END IF
CONTINUOUS
	m_avg = (f_in.m + f_out.m)/2
	P_red = P/(Pcrit*1E-5)
	q = abs(alphaHat*(T-tp_out.Tk[1]))
	MW = getMolWt (f_in.fluid)
	IF (ConstantHTC) INSERT
		alphaHat = max(1,getSmoothHTC(alpha_g,alpha_l,alpha_tp,m_avg,m_steady,x))
	ELSE
		alphaHat = getCorrelationHTC(f_in.fluid,choice_1p,choice_2p,geo.Ap,cp,D,k,m_avg,mu,MW,P,P_red,q,sigma,T,tp_out.Tk[1],x,sat)
	END IF
	IF (DynamicHTC) INSERT
		alpha' = 1/tauAlpha * (alphaHat - alpha)
	ELSE
		alpha = alphaHat
	END IF
<Q> 	Q = alpha*geo.As*(T - tp_out.Tk[1])
	tp_out.q[1] = Q
END COMPONENT






--------------------------------------------------------------------------------
// ADIABATIC LUMPED VOLUMES
--------------------------------------------------------------------------------
COMPONENT AdiabVol_Ph IS_A LumpedAdiabaticVolume (
	ENUM VoidFractionModel voidFractionModel = Homogeneous
)
"Lumped adiabatic control volume with P,h state variables"
// Includes void fraction model
DECLS
	REAL h_corr UNITS u_J_kg "Corrected enthalpy"
	REAL h_flow UNITS u_J_kg "Flow-weighted enthalpy"
	REAL hin UNITS u_J_kg "Inlet enthalpy"
	REAL hout UNITS u_J_kg "Outlet enthalpy"
	PRIVATE REAL P_Pa UNITS u_Pa
INIT
	P_Pa = P0*1E5 // Solver works better when P in Pa because h is in J/kg
CONTINUOUS
	P = P_Pa*1E-5
	
	h_corr = getCorrectionEnthalpy(voidFractionModel,h,rho,sat)
	h_flow = h + h_corr
	hin = donor_cell(f_in.m,f_in.hf,h_flow)
	hout = donor_cell(f_out.m,h_flow,f_out.hb)

<mas>	geo.V*(drho_dP*P_Pa'+drho_dh*h') = f_in.m - f_out.m
<ene> 	geo.V*((h*drho_dP-1)*P_Pa' + (h*drho_dh+rho)*h') - 9.81*(dz-z_in)*rho*vel/L = f_in.m*hin - f_out.m*hout
<st> 	getStatePh(f_in.fluid,P,h,cp,drho_dP,drho_dh,gamma,k,mu,phase,Pr,rho,sigma,T,Tsat,u,vsound,x,sat,ier,ipx,ipy)

<:hb> 	f_in.hb = h_flow - 0.5*9.8*dz // Lumped implementation of gravity equation
<:hf> 	f_out.hf = h_flow + 0.5*9.8*dz
END COMPONENT



COMPONENT AdiabVol_Phr IS_A AdiabVol_Ph
"Lumped adiabatic volume with P,h,rho state variables"
DECLS
	REAL tmp // dummy density
	REAL dh_drho
	REAL dh_dP
CONTINUOUS
	dh_drho = 1/drho_dh
	dh_dP = -(drho_dP)/(drho_dh)
<:mas> 	geo.V*rho' = f_in.m - f_out.m
<:ene> 	geo.V*((rho*dh_dP - 1)*P_Pa' + (dh_drho + h)*rho') -  9.81*(dz-z_in)*rho*vel/L = f_in.m*hin - f_out.m*hout
	h' = dh_dP*P_Pa' + dh_drho*rho'
<:st> 	getStatePh(f_in.fluid,P,h,cp,drho_dP,drho_dh,gamma,k,mu,phase,Pr,tmp,sigma,T,Tsat,u,vsound,x,sat,ier,ipx,ipy)
END COMPONENT



COMPONENT AdiabVol_ru IS_A LumpedAdiabaticVolume
"Lumped adiabatic volume with rho,u state variables"
CONTINUOUS
<mas> 	geo.V*rho' = f_in.m - f_out.m
<ene> 	geo.V*(rho*u'+u*rho') - 9.81*(dz-z_in)*rho*vel/L = f_in.mh - f_out.mh
	getStateRU(f_in.fluid,rho,u,cp,drho_dP,drho_dh,gamma,h,k,mu,P,phase,Pr,sigma,T,Tsat,vsound,x,sat,ier,ipx,ipy)
END COMPONENT






--------------------------------------------------------------------------------
// LUMPED VOLUMES WITH THERMAL PORT
--------------------------------------------------------------------------------
COMPONENT CV_Ph IS_A LumpedVolume_heatTransfer (
	ENUM VoidFractionModel voidFractionModel = Homogeneous
)
"Lumped volume with Ph state variables and thermal port heat transfer"
// Model derived from Laughman et al. 'A Comparison of Transient Heat Pump 
// Cycle Models Using Alternative Flow Descriptions', STBE, 2015". 

DECLS
	REAL h_corr UNITS u_J_kg "Corrected enthalpy"
	REAL h_flow UNITS u_J_kg "Flow-weighted enthalpy"
	REAL hin UNITS u_J_kg "Inlet enthalpy"
	REAL hout UNITS u_J_kg "Outlet enthalpy"
	REAL P_Pa UNITS u_Pa
	REAL As UNITS u_m2 "Surface area"
INIT
	P_Pa = P0*1E5
DISCRETE
	ASSERT (abs(dz-z_in)<=L) FATAL "The elevation of the tube cannot be higher than the tube length itself"
CONTINUOUS
	As = geo.As
	P = P_Pa*1E-5
	h_corr = getCorrectionEnthalpy(voidFractionModel,h,rho,sat)	
	h_flow = h + h_corr
	hin = donor_cell(f_in.m,f_in.hf,h_flow)
	hout = donor_cell(f_out.m,h_flow,f_out.hb)

<mas> 	geo.V*(drho_dP*P_Pa'+drho_dh*h') = f_in.m - f_out.m
<ene> 	geo.V*((h*drho_dP-1)*P_Pa' + (h*drho_dh+rho)*h') - 9.81*(dz-z_in)*rho*vel/L = f_in.m*hin - f_out.m*hout - Q 
<st> 	getStatePh(f_in.fluid,P,h,cp,drho_dP,drho_dh,gamma,k,mu,phase,Pr,rho,sigma,T,Tsat,u,vsound,x,sat,ier,ipx,ipy)

<:hb> 	f_in.hb = h_flow - 0.5*9.8*dz
<:hf> 	f_out.hf = h_flow + 0.5*9.8*dz
END COMPONENT



COMPONENT CV_Phr IS_A CV_Ph
"Lumped volume with Phr state variables and thermal port heat transfer"
DECLS
	REAL tmp
	REAL dh_drho
	REAL dh_dP
CONTINUOUS
	dh_drho = 1/drho_dh
	dh_dP = -(drho_dP)/(drho_dh)
<:mas> 	geo.V*rho' = f_in.m - f_out.m
<:ene> 	geo.V*((rho*dh_dP - 1)*P_Pa' + (dh_drho + h)*rho') -  9.81*(dz-z_in)*rho*vel/L = f_in.m*hin - f_out.m*hout - Q
	h' = dh_dP*P_Pa' + dh_drho*rho'
<:st> 	getStatePh(f_in.fluid,P,h,cp,drho_dP,drho_dh,gamma,k,mu,phase,Pr,tmp,sigma,T,Tsat,u,vsound,x,sat,ier,ipx,ipy)
END COMPONENT



COMPONENT CV_ru IS_A LumpedVolume_heatTransfer
"Lumped volume with rho,u state variables and thermal port heat transfer"
CONTINUOUS
<mas> 	geo.V*rho' = f_in.m - f_out.m
<ene> 	geo.V*(rho*u'+u*rho') - 9.81*(dz-z_in)*rho*vel/L  = f_in.mh - f_out.mh - Q
	getStateRU(f_in.fluid,rho,u,cp,drho_dP,drho_dh,gamma,h,k,mu,P,phase,Pr,sigma,T,Tsat,vsound,x,sat,ier,ipx,ipy)
END COMPONENT