/*-----------------------------------------------------------------------------------------
 LIBRARY: CO2
 FILE: FUNCTIONS
 AUTHOR: Viren Bhanot
 COMPANY: CERN
 DESCRIPTION: Mathematical functions for thermofluid calculations
 CREATION DATE: 10/10/2016
-----------------------------------------------------------------------------------------*/
USE MATH
USE CRYOLIB
USE PORTS_LIB

//"C++" FUNCTION a,REAL rp_hPx(IN REAL P, IN REAL x) IN "ecosimproRefprop.lib"
//"C++" FUNCTION REAL rp_rhoPx(IN REAL P, IN REAL x) IN "ecosimproRefprop.lib"
//"C++" FUNCTION REAL rp_all(IN REAL P_P OUT REAL hl, OUT REAL hv, OUT REAL rhol, OUT REAL rhov, OUT REAL d1, OUT REAL d2, OUT REAL drhol_dP, OUT REAL drhol_dh, OUT REAL drhov_dP, OUT REAL drhov_dh) IN "ecosimproRefprop.lib"

/*
 * INDEX
 * 
 * 1. Initial Conditions Assignment
 * 	1. assignInitLumped
 * 	2. assignInitDiscretized
 * 	3. assignInit2D
 * 2. Modelica Functions
 * 	1. semiLinear
 * 	2. regRoot
 * 	3. evaluatePoly3_derivativeAtZero
 * 	4. regRoot2_utility
 * 	5. regRoot2
 * 	6. spliceFunction
 * 3. linspace
 * 4. getCorrectionEnthalpy
 */



--------------------------------------------------------------------------------
// GLOBALS
--------------------------------------------------------------------------------
ENUM InitialConditionsAll = {
	P_h, // Pressure-Enthalpy
	P_level, // Pressure-Liquid level
	P_rho, // Pressure-Density
	P_sub, // Pressure-Subcooling
	P_sup, // Pressure-Superheating
	P_T, // Pressure-Temperature
	P_x, // Pressure-Vapour quality
	T_h, // Temperature-Enthalpy
	T_level, // Temperature-Liquid level
	T_rho, // Temperature-Density
	T_x // Temperature-Vapour quality
}
ENUM InitialConditions = {Ph,PT,Px}
SET_OF(InitialConditionsAll) InitialConditionsAccu = {P_h,P_level,P_T,P_x}
SET_OF(InitialConditionsAll) InitialConditions2FM = {P_level,P_x,T_level,T_x}
ENUM VoidFractionModel = {
	Homogeneous,
	Baroczy, // Equation doesn't match Baroczy data for void fraction above 0.9. Perhaps annular flow
	Levy,
	LockhartMartinelli,
	Smith,
	Thom, // Steam-water
	Turner, // Turbulent-turbulent flow, poor predictions 
	Zivi // Hated by Andrea, Butterworth says good for pressure drop/heat transfer prediction and for condensation
}



FUNCTION NO_TYPE linspace(REAL inlet, REAL outlet, INTEGER n, OUT REAL y[n], INTEGER type = 1)
"Linearly space an array between inlet and outlet values (derived from matlab)"
/* 
 * Linear distributions are used for obtaining cohesive initial distribution
 * of thermodynamic parameters. Several options have been implemented.
 * For values of inlet=1, outlet=3, n=3, outputs are:
 * 	- type 1 -> 1, 2, 3
 * 	- type 2 -> 1, 1.667, 2.333 // final ghost cell uses outlet=3 value
 * 	- type 3 -> 1.667, 2.333, 3 // initial ghost cell uses inlet=1 value
 */
BODY
	IF (n==1) THEN 
		y[1] = 0.5*(inlet + outlet)
		RETURN
	END IF
	FOR (i IN 1,n)
		IF (type == 2) THEN y[i] = inlet + (i-1)/n*(outlet-inlet)
		ELSEIF (type == 3) THEN y[i] = inlet + i/n*(outlet-inlet)
		ELSE y[i] = inlet + (i-1)*(outlet-inlet)/(n-1)
		END IF
	END FOR
END FUNCTION



--------------------------------------------------------------------------------
// INITIAL CONDITIONS
--------------------------------------------------------------------------------
FUNCTION NO_TYPE assignInitLumped (
	IN ENUM InitialConditions init,
	IN ENUM ChemName fluid,
	IN REAL P0,
	IN REAL h0,
	IN REAL T0,
	IN REAL x0,
	OUT REAL P,
	OUT REAL h,
	OUT REAL rho,
	OUT REAL u
) 
"Assigns initial conditions for lumped volume components"
// Note that all of P,h,rho and u are assigned, because some components might
// use rho,u while others might use P,h etc.
DECLS
	INTEGER ier,ipx,ipy
BODY
	P = P0
	IF (init==Ph) THEN
		h = h0
		rho = CRYO_PF_prop_vs_ph(fluid,P0,h0,fprop_density,ier,ipx,ipy)
		u = CRYO_PF_prop_vs_ph(fluid,P0,h0,fprop_energy,ier,ipx,ipy)
	ELSEIF (init==PT) THEN
		h = CRYO_PF_prop_vs_pT(fluid,P0,T0,fprop_enthalpy,ier,ipx,ipy)
		rho = CRYO_PF_prop_vs_pT(fluid,P0,T0,fprop_density,ier,ipx,ipy)
		u = CRYO_PF_prop_vs_pT(fluid,P0,T0,fprop_energy,ier,ipx,ipy)
	ELSE
		h = CRYO_PF_prop_vs_Px(fluid,P0,x0,fprop_enthalpy,ier,ipx,ipy)
		rho = CRYO_PF_prop_vs_Px(fluid,P0,x0,fprop_density,ier,ipx,ipy)
		u = CRYO_PF_prop_vs_Px(fluid,P0,x0,fprop_energy,ier,ipx,ipy)
	END IF
END FUNCTION



FUNCTION NO_TYPE assignInitAccu (
	IN ENUM InitialConditionsAccu init,
	IN ENUM ChemName fluid,
	IN REAL P0,
	IN REAL h0,
	IN REAL T0,
	IN REAL x0,
	IN REAL Level0,
	OUT REAL P,
	OUT REAL h,
	OUT REAL rho,
	OUT REAL u
)
"Assign initial conditions for homogeneous accumulators"
DECLS
	REAL lvlFraction "Fractional liquid level"
	REAL gamma "Void fraction"
	REAL rhog, rhof UNITS u_kg_m3 "Saturated densities"
	REAL xx
	INTEGER ier,ipx,ipy
BODY
	P = P0
	lvlFraction = min(1,max(0,Level0/100))
	IF (init==P_h) THEN
		h = h0
		rho = CRYO_PF_prop_vs_ph(fluid,P0,h0,fprop_density,ier,ipx,ipy)
		u = CRYO_PF_prop_vs_ph(fluid,P0,h0,fprop_energy,ier,ipx,ipy)
	ELSEIF (init==P_T) THEN
		h = CRYO_PF_prop_vs_pT(fluid,P0,T0,fprop_enthalpy,ier,ipx,ipy)
		rho = CRYO_PF_prop_vs_pT(fluid,P0,T0,fprop_density,ier,ipx,ipy)
		u = CRYO_PF_prop_vs_pT(fluid,P0,T0,fprop_energy,ier,ipx,ipy)
	ELSEIF (init==P_x) THEN
		h = CRYO_PF_prop_vs_Px(fluid,P0,x0,fprop_enthalpy,ier,ipx,ipy)
		rho = CRYO_PF_prop_vs_Px(fluid,P0,x0,fprop_density,ier,ipx,ipy)
		u = CRYO_PF_prop_vs_Px(fluid,P0,x0,fprop_energy,ier,ipx,ipy)
	ELSEIF (init==P_level) THEN
		rhof = CRYO_PF_prop_vs_Px(fluid,P0,0,fprop_density,ier,ipx,ipy)
		rhog = CRYO_PF_prop_vs_Px(fluid,P0,1,fprop_density,ier,ipx,ipy)
		gamma = (1 - lvlFraction)
		rho = rhof*(1 - gamma) + gamma*rhog
		xx = gamma*rhog/rho
		h = CRYO_PF_prop_vs_Px(fluid,P0,xx,fprop_enthalpy,ier,ipx,ipy)
		u = CRYO_PF_prop_vs_Px(fluid,P0,xx,fprop_energy,ier,ipx,ipy)
	END IF
END FUNCTION



FUNCTION NO_TYPE assignInitAccu2FM (
	IN ENUM InitialConditions2FM init,
	IN ENUM ChemName fluid,
	IN REAL P0 "Initial pressure in bar",
	IN REAL Tsat0 "Initial saturation temperature in K",
	IN REAL Level0 "Initial liquid level as percentage",
	IN REAL x0 "Initial vapour quality",
	OUT REAL P "Pressure in bar",
	OUT REAL hf "Liquid enthalpy",
	OUT REAL hg "Vapour enthalpy",
	OUT REAL y "Liquid level as fraction of 1"
)
"Assigns initial conditions for two-fluid accumulator models"
DECLS
	REAL rhof
	REAL rhog
	INTEGER ier,jx,jy // error codes
BODY
	IF (init==P_level) THEN
		P = P0
		hf = rp_hPx(P0*1e5,0)
		hg = rp_hPx(P0*1e5,1)
		y = max(0,min(1,Level0/100))
	ELSEIF (init==P_x) THEN
		P = P0
		rhof = rp_rhoPx(P0*1e5,0)
		rhog = rp_rhoPx(P0*1e5,1)
		y = 1 / (1 + (1-x0)/x0*(rhog/rhof))
		hf = rp_hPx(P0*1e5,0)
		hg = rp_hPx(P0*1e5,1)
	ELSEIF (init==T_level) THEN
		P = CRYO_PF_prop_vs_Tx(fluid,Tsat0,x0,fprop_pressure,ier,jx,jy)
		hf = rp_hPx(P0*1e5,0)
		hg = rp_hPx(P0*1e5,1)
		y = max(0,min(1,Level0/100))
	ELSEIF (init==T_x) THEN
		P = CRYO_PF_prop_vs_Tx(fluid,Tsat0,x0,fprop_pressure,ier,jx,jy)
		rhof = rp_rhoPx(P0*1e5,0)
		rhog = rp_rhoPx(P0*1e5,1)
		y = 1 / (1 + (1-x0)/x0*(rhog/rhof))
		hf = rp_hPx(P0*1e5,0)
		hg = rp_hPx(P0*1e5,1)
	END IF
END FUNCTION



FUNCTION NO_TYPE assignInitDiscretized (
	IN INTEGER n "Number of control volumes",
	IN ENUM InitialConditions init,
	IN ENUM ChemName fluid,
	IN REAL P0_i,
	IN REAL P0_o,
	IN REAL h0_i,
	IN REAL h0_o,
	IN REAL T0_i,
	IN REAL T0_o,
	IN REAL x0_i,
	IN REAL x0_o,
	OUT REAL P[n],
	OUT REAL h[n],
	OUT REAL rho[n],
	OUT REAL u[n],
	OUT REAL x[n],
	IN INTEGER type = 1 "linspace type argument"
)
"Assigns initial conditions for a 1D discretized volume"
DECLS
	REAL T[n]
	REAL Tsat[n]
	INTEGER ier,ipx,ipy
BODY
	linspace(P0_i,P0_o,n,P,type) // see linspace below for explanation
	IF (init==Ph) THEN
		linspace(h0_i,h0_o,n,h,type)
		FOR (i IN 1,n)
			rho[i] = CRYO_PF_prop_vs_ph(fluid,P[i],h[i],fprop_density,ier,ipx,ipy)
			u[i] = CRYO_PF_prop_vs_ph(fluid,P[i],h[i],fprop_energy,ier,ipx,ipy)
			x[i] = CRYO_PF_prop_vs_ph(fluid,P[i],h[i],fprop_quality,ier,ipx,ipy)
		END FOR
	ELSEIF (init==PT) THEN
		linspace(T0_i,T0_o,n,T,type)
		FOR (i IN 1,n)
			Tsat[i] = CRYO_PF_Tsat_vs_p(fluid,P[i],ier,ipy)
			IF (T[i]<=Tsat[i]) THEN
				x[i] = 0
			ELSEIF(T[i]>Tsat[i]) THEN
				x[i] = 1
			END IF
			h[i] = CRYO_PF_prop_vs_pT(fluid,P[i],T[i],fprop_enthalpy,ier,ipx,ipy)
			rho[i] = CRYO_PF_prop_vs_pT(fluid,P[i],T[i],fprop_density,ier,ipx,ipy)
			u[i] = CRYO_PF_prop_vs_pT(fluid,P[i],T[i],fprop_energy,ier,ipx,ipy)
			
		END FOR
	ELSE
		linspace(x0_i,x0_o,n,x,type)
		FOR (i IN 1,n)
			h[i] = CRYO_PF_prop_vs_Px(fluid,P[i],x[i],fprop_enthalpy,ier,ipx,ipy)
			rho[i] = CRYO_PF_prop_vs_Px(fluid,P[i],x[i],fprop_density,ier,ipx,ipy)
			u[i] = CRYO_PF_prop_vs_Px(fluid,P[i],x[i],fprop_energy,ier,ipx,ipy)
		END FOR
	END IF
END FUNCTION



FUNCTION NO_TYPE assignInit2D (
	IN INTEGER nBank,
	IN INTEGER nTube,
	IN ENUM InitialConditions init,
	IN ENUM ChemName fluid,
	IN REAL P0_i,
	IN REAL P0_o,
	IN REAL h0_i,
	IN REAL h0_o,
	IN REAL T0_i,
	IN REAL T0_o,
	IN REAL x0_i,
	IN REAL x0_o,
	OUT REAL P[nBank,nTube],
	OUT REAL h[nBank,nTube]
)
"Assign initial conditions for a 2D discretized volume"
DECLS
	REAL T[nBank,nTube]
	REAL x[nBank,nTube]
	INTEGER ier,ipx,ipy
BODY
	FOR (i IN 1,nBank)
		FOR (j IN 1,nTube)
			P[i,j] = P0_i + ((i-1)*nTube + j - 1) / (nBank*nTube) * (P0_o - P0_i)
		END FOR
	END FOR
	FOR (i IN 1,nBank)
		FOR (j IN 1,nTube)
			IF (init==Ph) THEN
				h[i,j] = h0_i + ((i-1)*nTube + j) / (nBank*nTube) * (h0_o - h0_i)
			ELSEIF (init==PT) THEN
				T[i,j] = T0_i + ((i-1)*nTube + j) / (nBank*nTube) * (T0_o - T0_i)
				h[i,j] = CRYO_PF_prop_vs_pT(fluid,P[i,j],T[i,j],fprop_enthalpy,ier,ipx,ipy)
			ELSE
				x[i,j] = x0_i + ((i-1)*nTube + j) / (nBank*nTube) * (x0_o - x0_i)
				h[i,j] = CRYO_PF_prop_vs_Px(fluid,P[i,j],x[i,j],fprop_enthalpy,ier,ipx,ipy)
			END IF
		END FOR
	END FOR		
END FUNCTION



--------------------------------------------------------------------------------
// MODELICA FUNCTIONS
--------------------------------------------------------------------------------
FUNCTION REAL semiLinear(REAL x, REAL pos, REAL neg)
"Function for handling reverse flows. Returns x*pos for x>=0, and x*neg for x<0"
DECLS
	REAL y
BODY
	IF(x>=0) THEN y = pos*x 
	ELSE y = neg*x
	END IF
	RETURN y
END FUNCTION



FUNCTION REAL regRoot(REAL x, REAL delta=0.01)
"Function for calculating square root, with finite slope around x=0"
BODY
	RETURN x/(x*x+delta*delta)**0.25
END FUNCTION



FUNCTION REAL evaluatePoly3_derivativeAtZero(REAL x,REAL x1,REAL y1,REAL y1d,REAL y0d)
"Fritsch-Carlson Order 3 Polynomial (derived from modelica). Basically, it increases monotonically"
DECLS
	REAL a1, a2, a3, xx
BODY	
	a1 = x1*y0d
	a2 = 3*y1-x1*y1d-2*a1
	a3 = y1-a2-a1
	xx = x/x1
	RETURN xx*(a1+xx*(a2+xx*a3))
END FUNCTION



FUNCTION REAL regRoot2_utility(REAL x,REAL x1,REAL k1,REAL k2,BOOLEAN use_yd0,REAL yd0)
"Used in regRoot2 function"
DECLS
	REAL y0d, y1, y2, y1d, y2d, x2, w, w1, w2, ret
BODY
	x2 = -x1*(k2/k1)
	IF (x<=x2) THEN
		ret = -sqrt(k2*abs(x))
	ELSE
		y1 = sqrt(k1*x1)
		y2 =-sqrt(k2*abs(x2))
		y1d = sqrt(k1/x1)/2
		y2d = sqrt(2/abs(x2))/2
		IF (use_yd0) THEN
			y0d = yd0
		ELSE
			w = x2/x1
			y0d = ((3*y2-x2*y2d)/w - (3*y1-x1*y1d)*w)
		END IF
		w1 = sqrt(8.75*k1/x1)
		w2 = sqrt(8.75*k2/abs(x2))
		y0d = min(y0d,0.9*min(w1,w2))
		IF (x>=0) THEN
			ret = y1*evaluatePoly3_derivativeAtZero(x/x1,1,1,y1d*x1/y1,y0d*x1/y1)
		ELSE
			ret = y1*evaluatePoly3_derivativeAtZero(x/x1,x2/x1,y2/y1,y2d*x1/y1,y0d*x1/y1)
		END IF
	END IF
	RETURN ret
END FUNCTION



FUNCTION REAL regRoot2 (REAL x, REAL x_small=0.01, REAL k1=1, REAL k2=1, BOOLEAN use_yd0=FALSE, REAL yd0=1)
"Function for calculating square root using Fritsch Carlson monotonically increasing 3rd order polynomial around x=0 region"
DECLS
	REAL ret
BODY
	IF (x>=x_small) THEN
		ret = sqrt(k1*x)
	ELSEIF (x<=-x_small) THEN
		ret = -sqrt(k2*abs(x))
	ELSEIF (k1>=k2) THEN
		ret = regRoot2_utility(x,x_small,k1,k2,use_yd0,yd0)
	ELSE
		ret = -regRoot2_utility(-x,x_small,k2,k1,use_yd0,yd0)
	END IF
	RETURN ret
END FUNCTION



FUNCTION REAL regPowGen(REAL x, REAL x_small=0.01, REAL a=0.5, REAL b=1)
"Way better than regRoot2 at smoothing low flow regions. a=0.5/b=1 and a=0.5/b=3 are recommended combos"
BODY
	RETURN x**b * (x**2 + x_small**2)**((a-b)/2)
END FUNCTION


FUNCTION REAL spliceFunction (IN REAL pos, IN REAL neg, IN REAL x, IN REAL delx)
DECLS
	REAL sX
	REAL sX1
	REAL w
BODY
	sX1 = x/delx
	sX = sX1 * asin(1)
	IF (sX1 <=-0.999999999) THEN
	    w = 0
	ELSEIF (sX1 >= 0.999999999) THEN
	    w = 1
	ELSE
	    w = (tanh(tan(sX)) + 1)/2
	END IF
	RETURN w*pos + (1-w)*neg
END FUNCTION



FUNCTION REAL getCorrectionEnthalpy (
	IN ENUM VoidFractionModel choice,
	IN REAL h,
	IN REAL rho,
	IN SaturationProperties sat
) 
"Returns a correction term for staggered-grid thermal cell boundary enthalpy, to handle void fraction"
/*
 * Returns a correction term that is added to the in-situ enthalpy (cell enthalpy).
 * Bauer (1999) has a method where the _mass flow rate_ is corrected for by calculating
 * a phase velocity. This leads to problems when using correlations that don't 
 * calculate slip ratio explicitly. Instead, Hongtao Qiao's (Qiao, 2014, Phd thesis) method prefers 
 * to correct for the enthalpy directly. See his 2015 Flash Tank paper, and his paper with Laughman.
 *
 * Note that the symbols for flow quality and static quality are inverted in the FTVI 
 * paper/thesis vs his papers at MERL with Laughman. The correlations themselves are based on 
 * Buterworth's 1975 paper.
 */
DECLS
	REAL c,q,r,s // correlation constants
	REAL gamma // void fraction
	REAL hfg
	REAL S
	REAL x "Static quality"
	REAL xHat "Flow weighted quality"
	REAL h_corr "Enthalpy correction term"
BODY
	IF (sat.rho_g == sat.rho_l) THEN
		h_corr = 0
	ELSE
		hfg = sat.h_g - sat.h_l
		gamma = max(0,min(1,(rho-sat.rho_l)/(sat.rho_g-sat.rho_l)))
		IF (sat.h_g==sat.h_l) THEN
			xHat = 0
		ELSE
			xHat = max(0,min(1,(h-sat.h_l)/(sat.h_g-sat.h_l)))
		END IF
		IF (choice == Baroczy) THEN
			c = 1
			q = 0.74
			r = 0.65
			s = 0.13
		ELSEIF (choice == LockhartMartinelli) THEN
			c = 0.28
			q = 0.64
			r = 0.36
			s = 0.07
		ELSEIF (choice == Thom) THEN
		// Good for steam-water mixtures
			c = 1
			q = 1
			r = 0.89
			s = 0.18
		ELSEIF (choice == Turner) THEN
			c = 1
			q = 0.72
			r = 0.40
			s = 0.08
		ELSEIF (choice == Smith) THEN
			c = 0.79
			q = 0.78
			r = 0.58
			s = 0
		ELSEIF (choice == Zivi) THEN
		// Good for condensation (NB: Andrea hates this model)
			c = 1
			q = 1
			r = 0.67
			s = 0
		ELSE // Homogeneous assumption
			c = 1
			q = 1
			r = 1
			s = 0
		END IF
		x = max(0,min(1,gamma**(1/q) / (gamma**(1/q) + (1/c*(sat.rho_l/sat.rho_g)**r * (sat.mu_g/sat.mu_l)**s)**(1/q) * (1-gamma)**(1/q))))
		h_corr = (x-xHat)*hfg
	END IF
	RETURN h_corr
END FUNCTION



--------------------------------------------------------------------------------
// PLC RELATED FUNCTIONS
--------------------------------------------------------------------------------
// LOGICAL OPERATIONS
FUNCTION BOOLEAN LogicalAND (
	IN BOOLEAN bool1,
	IN BOOLEAN bool2
)
DECLS
	BOOLEAN y
BODY
	IF (bool1 AND bool2) THEN
		y = TRUE
	ELSE
		y = FALSE
	END IF
	RETURN y
END FUNCTION



FUNCTION BOOLEAN LogicalNOT (
	IN BOOLEAN bool
)
DECLS
	BOOLEAN y
BODY
	IF (bool) THEN
		y = FALSE
	ELSE
		y = TRUE
	END IF
	RETURN y
END FUNCTION



FUNCTION REAL TransitionLogic (
	IN BOOLEAN Tr0,
	IN BOOLEAN Tr1,
	IN BOOLEAN Tr2,
	IN BOOLEAN Tr3
)
"Logic for transitions between stepper positions"
// We only go up through the stepper, and not down.
// Instead, if a STOP condition is activated, the
// PCO should control for it.
DECLS
	REAL y
BODY
	IF (Tr3) THEN
		y = 3
		RETURN y
	ELSEIF (Tr2) THEN
		y = 2
		RETURN y
	ELSEIF (Tr1) THEN
		y = 1
		RETURN y
	ELSE 
		y = 0
	END IF	
	RETURN y
END FUNCTION



FUNCTION REAL TransitionLogicDemo (
	IN BOOLEAN Tr0,
	IN BOOLEAN Tr1,
	IN BOOLEAN Tr2,
	IN BOOLEAN Tr3,
	IN BOOLEAN Tr4,
	IN BOOLEAN Tr5,
	IN BOOLEAN Tr6,
	IN BOOLEAN Tr7
)
"Logic for transitions for DEMO"
// Same as previous TransitionLogic component
// but just has more steps.
DECLS
	REAL y
BODY
	IF (Tr7) THEN
		y = 7
	ELSEIF (Tr6) THEN
		y = 6
	ELSEIF (Tr5) THEN
		y = 5
		RETURN y
	ELSEIF (Tr4) THEN
		y = 4
		RETURN y
	ELSEIF (Tr3) THEN
		y = 3
		RETURN y
	ELSEIF (Tr2) THEN
		y = 2
		RETURN y
	ELSEIF (Tr1) THEN
		y = 1
		RETURN y
	ELSE
		y = 0
	END IF	
	RETURN y
END FUNCTION



/*
FUNCTION valveFlowRate (
	IN ChemName fluid,
	IN REAL Cv,
	IN REAL P1,
	IN REAL h1,
	IN REAL P2,
	IN REAL h2
)
"Calculate flow rate through valve, with "
DECLS
	REAL dP
	REAL Gl
	REAL Gv
	CONST REAL N1 = 14.42
	REAL Pin
	REAL Pout
	REAL hin
	REAL rhoin
	REAL Tin
	REAL xin
	REAL vDotl
	REAL vDotv
	REAL vDotv_crit
BODY
	dP = P1 - P2
	Pin = donor_cell(dP,P1,P2)
	Pout = donor_cell(dP,P2,P1)
	hin = donor_cell(dP,h1,h2)
	rhoin = CRYO_PF_prop_vs_ph(fluid,Pin,hin,fprop_density,ier,xx,yy)
	Tin = CRYO_PF_prop_vs_ph(fluid,Pin,hin,fprop_temperature,ier,xx,yy)
	xin = CRYO_PF_prop_vs_ph(fluid,Pin,hin,fprop_quality,ier,xx,yy)
	
	rhoRef = CRYO_PF_prop_vs_pT(f_in.fluid,1.01325,293.15,fprop_density,ier,xx,yy)
	Gl = rhoRef/998 // against water
	Gv = rhoRef/1.205 // against air
	
	vDotv = N2*Cv*Pin * (1-2*abs(dP)/(3*Pin)) * regRoot2(dP/(Pin*G*Tin),1e-4) // low dP case
	vDotv_crit= sign(dP) * 0.471 * N2*Cv*Pin * regRoot2(1/(G*Tin),1e-4) // high dP case (choked flow)
	vDotl = N1*Cv*regRoot2(dP/G,0.01) // liquid case
	
	m = ZONE (x>=1) vDot * 1e-3/60 * rho * (Tin/273.15) * (1/Pin) // convert slpm to kg/s
	OTHERS vDot * 1e-3/60 * rho // convert lpm to kg/s

END FUNCTION
*/