/*-----------------------------------------------------------------------------------------
 LIBRARY: CO2
 FILE: FLUIDPROPS
 AUTHOR: Viren Bhanot
 COMPANY: CERN
 CREATION DATE: 14/06/2017
-----------------------------------------------------------------------------------------*/

/*
 * Extends from Cryolib state_vs_Ph and state_vs_ru methods.
 * I use a different method for calculating partial derivatives in two-phase region
 */

USE MATH
USE THERMAL
USE THERMO_TABLE_INTERP
USE CRYOLIB



CONST REAL CO2_CRITENTH = 332245.6524734457 UNITS u_J_kg "CO2 critical enthalpy" // Cryolib is returning error for hcrit in v2.0.6
CONST REAL CO2_CRITPR = 73.77 UNITS u_bar "CO2 critical pressure" // Cryolib is returning garbage for Pcrit in v2.0.6
CONST REAL CO2_CRITTEMP = 304.128 UNITS u_K "CO2 critical temperature" // Cryolib is returning garbage for Pcrit in v2.0.6


FUNCTION NO_TYPE getStateRU (
	IN ENUM ChemName chem "Chemical constituents of the Working fluid",
	IN REAL rho UNITS u_kg_m3 "Partial Density of Chemical Const. 1",
	IN REAL u "Specific internal energy",
	OUT REAL cp "Specific heat at constant pressure",
	OUT REAL drho_dP "drho/dp at constant h",
	OUT REAL drho_dh "drho/dh at constant P",
	OUT REAL gamma "Void fraction",
	OUT REAL h "Enthalpy",
	OUT REAL k "Thermal conductivity of single phase",
	OUT REAL mu "Viscosity of single phase",
	OUT REAL p_bar "Pressure of Chemical Const. 1",
	OUT ENUM Phase phase "Phase",
	OUT REAL Pr "Prandtl number",
	OUT REAL sigma"Surface tension",
	OUT REAL T "Temperature",  
	OUT REAL Tsat "Saturation temperature",
	OUT REAL vsound "Speed of sound",
	OUT REAL x "Quality",
	OUT SaturationProperties sat "Class containing saturation properties",
	OUT INTEGER ier "Error flag",
	OUT INTEGER jx "Pointer to last table interval",
	OUT INTEGER jy "Pointer to last table interval"
)
DECLS
	REAL beta // isobaric thermal expansion coefficient
	REAL beta_g
	REAL beta_l
	REAL dhl_dP,dhg_dP
	REAL dP_dT
	REAL dT_dP
	REAL dv_dh,dv_dP
	REAL dvl_dP,dvg_dP
	REAL dxHat_dP
	REAL kappa // isothermal compressibility coefficient
	REAL kappa_g 
	REAL kappa_l
	REAL P
	REAL Pcrit
	REAL prop[3,20]
	REAL prop_crit[20]
	REAL prop_tril[20]
	REAL Tcrit
	REAL v_l,v_g
	REAL vprop
	REAL xHat
	INTEGER ig // Real fluid number
BODY
	ig = setofPos(ChemName, chem)
	IF(itab[ig] == 0) THEN
		IF(path_prop =="") THEN
			path_prop = expandFilePath(Pathprops)
			FOR (i IN ChemName)   // loop to transfer defined CEA labels to FORTRAN
				read_label(CEA_GAS_Name[i]) 
			END FOR
			// reading of CEA coefficients of all declared ChemName
			read_cea(setofSize(ChemName), path_prop, MW_ch, ier)
		END IF
		WRITE("\n   ******    Reading properties file %s.dat ******\n",ChemFileName[chem])
		table_read(" ",ig, path_prop, ChemFileName[chem], ig, ier, fml, UsrFluiddH[chem])
		ng_tot = ng_tot + 1
		itab[ig] = ng_tot
	END IF
	ASSERT (ier < 100) FATAL "StateRU: problems reading thermo properties file"
	crit_cond(ig, prop_crit, ier)
	trip_l_cond(ig, prop_tril, ier)
	Pcrit = prop_crit[1]
	Tcrit = prop_crit[2]

	// Get thermodynamic properties:
	thermo_prop(ig, 3, rho, 7, u, 111, prop, x, ier, INTERP_ORDER, jx, jy) // note: returns x too
	
	P = prop[1,1]
	T = prop[1,2]
	// Phase detection
	IF (P < Pcrit) THEN
		IF (x <= 0) THEN phase = liquid
		ELSEIF (x >= 1) THEN phase = vapor
		ELSE phase = two_phase
		END IF
		Tsat = prop[2,2]
		sat.rho_l = prop[2,3]
		sat.h_l = prop[2,5]
		sat.cp_l = prop[2,11]
		sat.mu_l = prop[2,20]
		sat.k_l = prop[2,19]
		beta_l = prop[2,14]
		kappa_l = prop[2,16]
		sat.rho_g = prop[3,3]
		sat.h_g = prop[3,5]
		sat.cp_g = prop[3,11]
		sat.mu_g = prop[3,20]
		sat.k_g = prop[3,19]
		beta_g = prop[3,14]
		kappa_g = prop[3,16]  
	ELSE
		IF (T >= Tcrit) THEN
			phase = supercritical
			x = 1.
		ELSE
			phase = liquid
			x = 0.
		END IF
		Tsat = 0 // above critical point
		sat.rho_l = prop_crit[3]
		sat.h_l = prop_crit[5]
		sat.cp_l = prop_crit[11]
		sat.mu_l = prop_crit[20]
		sat.k_l = prop_crit[19]
		beta_l = prop_crit[14]
		kappa_l = prop_crit[16]
		sat.rho_g = sat.rho_l
		sat.h_g = sat.h_l
		sat.cp_g = sat.cp_l
		sat.mu_g = sat.mu_l
		sat.k_g = sat.k_l
		beta_g = beta_l
		kappa_g = kappa_l
	END IF
	//  Property calculation as function of void fraction
	IF (x <= 0) THEN
		gamma = 0.
		sigma = prop[2,18]
		cp = prop[1,11]    
		mu = prop[1,20]
		k = prop[1,19]
		beta = prop[1,14]
		kappa = prop[1,16]
		drho_dh = -beta * rho / cp
		drho_dP = (rho * kappa + beta / cp * (1. - beta *T))
		vsound = prop[1,17]
	ELSEIF (x >= 1) THEN
		gamma = 1.
		sigma = prop[2,18]
		cp = prop[1,11]    
		mu = prop[1,20]
		k = prop[1,19]
		beta = prop[1,14]
		kappa = prop[1,16]
		drho_dh = -beta * rho / cp
		drho_dP = (rho * kappa + beta / cp * (1. - beta *T))
		vsound = prop[1,17]
	ELSE
		gamma = x / sat.rho_g / (x/sat.rho_g + (1.-x)/sat.rho_l)
		sigma = prop[2,18]
		// Approximated computation of vsound as equivalent two phase mixture
		vsound = max(0.01, 1/ssqrt(rho*(gamma/sat.rho_g/prop[3,17]**2 + (1-gamma)/sat.rho_l/prop[2,17]**2)))
		// Transport properties for two phase case
		cp = x*sat.cp_g + (1-x)*sat.cp_l
		k = x*sat.k_g + (1-x)*sat.k_l
		mu = x*sat.mu_g + (1-x)*sat.mu_l				
		// Tummescheit/Bauer/Thorade etc. method for density partial derivatives
		xHat = (h - sat.h_l) / (sat.h_g - sat.h_l)
		v_g = 1/sat.rho_g
		v_l = 1/sat.rho_l
		dT_dP = T * (v_g - v_l) / (sat.h_g - sat.h_l)
		dhl_dP = v_l*(1-beta_l*T) + sat.cp_l*dT_dP
		dhg_dP = v_g*(1-beta_g*T) + sat.cp_g*dT_dP
		dxHat_dP = (xHat*dhg_dP + (1-xHat)*dhl_dP) / (sat.h_l - sat.h_g)
		dv_dh = (v_g - v_l) / (sat.h_g - sat.h_l)
		dvl_dP = beta_l*v_l*dT_dP - kappa_l*v_l
		dvg_dP = beta_g*v_g*dT_dP - kappa_g*v_g
		dv_dP = dvl_dP + dxHat_dP*(v_g-v_l) + xHat*(dvg_dP - dvl_dP)
		drho_dh = -rho**2 * dv_dh
		drho_dP = -rho**2 * dv_dP
	END IF
	p_bar = P/1e5  
	h = u + P/rho
	Pr = cp * mu/k
	RETURN
END FUNCTION



FUNCTION NO_TYPE getStatePh (
	IN ENUM ChemName chem "Working chemical constituent",
	IN REAL p_bar "Pressure",
	IN REAL h "Enthalpy",	
	OUT REAL cp "Specific heat at constant pressure",
	OUT REAL drho_dP "drho/dp at constant h",
	OUT REAL drho_dh "drho/dh at constant P",
	OUT REAL gamma "Void fraction",
	OUT REAL k "Thermal conductivity of single phase",
	OUT REAL mu "Viscosity of single phase",
	OUT ENUM Phase phase  "Phase",
	OUT REAL Pr "Prandtl number",
	OUT REAL rho "Density",
	OUT REAL sigma "Surface tension",
	OUT REAL T "Temperature",
	OUT REAL Tsat "Saturation temperature",
	OUT REAL u "Specific internal energy",
	OUT REAL vsound "Speed of sound",
	OUT REAL x "Quality",
	OUT SaturationProperties sat "Saturation properties class",
	OUT INTEGER ier "Error flag",
	OUT INTEGER ipx	"last index for x interpolation",
	OUT INTEGER ipy "last index for y interpolation"
)
DECLS
	REAL beta // isobaric thermal expansion coefficient
	REAL beta_g
	REAL beta_l
	REAL dhl_dP,dhg_dP
	REAL dP_dT
	REAL dT_dP
	REAL dv_dh,dv_dP
	REAL dvl_dP,dvg_dP
	REAL dxHat_dP
	REAL kappa // isothermal compressibility coefficient
	REAL kappa_g
	REAL kappa_l
	REAL P UNITS u_Pa "Pressure"
	REAL Pcrit UNITS u_Pa "critical pressure"
	REAL prop[3,20]
	REAL prop_crit[20] "Fluid properties at critical point"
	REAL prop_tril[20] "Fluid properties at triple liquid point"
	REAL Tcrit "critical Temperature (K)"
	REAL v_l,v_g
	REAL vprop
	REAL xHat
	INTEGER ig "Real fluid number"
	INTEGER jx,jy
BODY
	P = p_bar*1e5	
	ig = setofPos(ChemName, chem) // specifies index of fluid in list of ChemName
	IF(itab[ig] == 0) THEN
		IF(path_prop =="") THEN
			path_prop = expandFilePath(Pathprops)
			FOR (i IN ChemName) // loop to transfer defined CEA labels to FORTRAN
				read_label(CEA_GAS_Name[i]) 
			END FOR
			// reading of CEA coefficients of all declared ChemName
			read_cea(setofSize(ChemName), path_prop, MW_ch, ier)
		END IF
		WRITE("\nReading properties file %s.dat\n",ChemFileName[chem])
		table_read(" ",ig, path_prop, ChemFileName[chem], ig, ier, fml, UsrFluiddH[chem])
		ng_tot = ng_tot + 1
		itab[ig] = ng_tot
	END IF
	ASSERT (ier < 100) FATAL "Problems reading thermo properties file"
	crit_cond(ig, prop_crit, ier)
	trip_l_cond(ig, prop_tril, ier)
	Pcrit = prop_crit[1]
	Tcrit = prop_crit[2]
	thermo_prop(ig, 1, P, 5, h, 111, prop, x, ier, INTERP_ORDER, ipx, ipy)
	T = prop[1,2]
	rho = prop[1,3]
	u = prop[1,7]
	// Saturation properties
	IF (P < Pcrit) THEN // Subcritical
		IF (x <= 0) THEN phase = liquid
		ELSEIF (x >= 1) THEN phase = vapor
		ELSE phase = two_phase
		END IF 
		Tsat = prop[2,2]
		sat.rho_l = prop[2,3]
		sat.h_l = prop[2,5]
		sat.cp_l = prop[2,11]
		sat.mu_l = prop[2,20]
		sat.k_l = prop[2,19]
		beta_l = prop[2,14]
		kappa_l = prop[2,16]
		sat.rho_g = prop[3,3]
		sat.h_g = prop[3,5]
		sat.cp_g = prop[3,11]
		sat.mu_g = prop[3,20]
		sat.k_g = prop[3,19]
		beta_g = prop[3,14]
		kappa_g = prop[3,16]  
	ELSE
		IF (T >= Tcrit) THEN 
			phase = supercritical
			x = 1.
		ELSE
			phase = liquid
			x = 0.
		END IF
		Tsat = 0
		sat.rho_l = prop_crit[3]
		sat.h_l = prop_crit[5]
		sat.cp_l = prop_crit[11]
		sat.mu_l = prop_crit[20]
		sat.k_l = prop_crit[19]
		beta_l = prop_crit[14]
		kappa_l = prop_crit[16]	
		// As soon as P>Pcrit, liquid and vapour saturation properties are equal to each other:
		sat.rho_g = sat.rho_l
		sat.h_g = sat.h_l
		sat.cp_g = sat.cp_l
		sat.mu_g = sat.mu_l
		sat.k_g = sat.k_l
		beta_g = beta_l
		kappa_g = kappa_l
	END IF
	// Other properties
	IF (x <= 0) THEN
		gamma = 0.
		sigma = prop[2,18]
		cp = prop[1,11] 
		mu = prop[1,20]
		k = prop[1,19]
		beta = prop[1, 14]
		kappa = prop[1, 16]
		drho_dh = -beta * rho / cp
		drho_dP = (rho*kappa + beta/cp*(1.-beta*T))
		vsound = prop[1,17]
	ELSEIF (x >= 1) THEN
		gamma = 1.
		sigma = 0
		cp = prop[1,11]  
		mu = prop[1,20]
		k = prop[1,19]
		beta = prop[1, 14]
		kappa = prop[1, 16]
		drho_dh = -beta * rho / cp
		drho_dP = (rho*kappa + beta/cp*(1-beta*T))
		vsound = prop[1,17]
	ELSE
		// Two-Phase
		gamma = x / sat.rho_g / (x/sat.rho_g + (1.-x)/sat.rho_l)
		sigma = prop[2,18]
		THERMO_TABLE_INTERP.thermo_prop_vs_xy(ig,3,rho,2,T,17,vsound,ier,P,x,INTERP_ORDER,jx,jy)
		cp = x*sat.cp_g + (1-x)*sat.cp_l
		k = x*sat.k_g + (1-x)*sat.k_l
		mu = x*sat.mu_g + (1-x)*sat.mu_l		
		// Method described in Thorade/Sadaat EAS 2013:
		xHat = (h-sat.h_l)/(sat.h_g-sat.h_l)
		v_g = 1/sat.rho_g
		v_l = 1/sat.rho_l
		dT_dP = T*(v_g - v_l)/(sat.h_g - sat.h_l)
		dhl_dP = v_l*(1-beta_l*T) + sat.cp_l*dT_dP
		dhg_dP = v_g*(1-beta_g*T) + sat.cp_g*dT_dP
		dxHat_dP = (xHat*dhg_dP + (1-xHat)*dhl_dP) / (sat.h_l - sat.h_g)
		dv_dh = (v_g - v_l)/(sat.h_g - sat.h_l)
		dvl_dP = beta_l*v_l*dT_dP - kappa_l*v_l
		dvg_dP = beta_g*v_g*dT_dP - kappa_g*v_g
		dv_dP = dvl_dP + dxHat_dP*(v_g-v_l) + xHat*(dvg_dP - dvl_dP)
		drho_dh = -rho**2 * dv_dh
		drho_dP = -rho**2 * dv_dP
	END IF
	Pr = cp * mu/k
	RETURN
END FUNCTION



FUNCTION NO_TYPE getThermodynamicState_ph (
	IN ENUM ChemName chem,
	IN REAL p_bar,
	IN REAL h,
	OUT REAL rho,
	OUT REAL s,
	OUT REAL T,
	OUT REAL u,
	OUT REAL x,
	OUT ENUM Phase phase,
	OUT INTEGER iex,
	OUT INTEGER ipx,
	OUT INTEGER ipy
)
"Smaller than getStatePh. Used in ThermodynamicState class. Note that this returns the entropy also"
DECLS
	REAL P UNITS u_Pa "Pressure"
	REAL Pcrit UNITS u_Pa "critical pressure"
	REAL prop[3,20]
	REAL prop_crit[20] "Fluid properties at critical point"
	REAL prop_tril[20] "Fluid properties at triple liquid point"
	REAL Tcrit "critical Temperature (K)"
	REAL h_l, h_g
	INTEGER ig "Real fluid number"
	INTEGER jx,jy
BODY
	P = p_bar*1e5	
	ig = setofPos(ChemName, chem)
	IF(itab[ig] == 0) THEN
		IF(path_prop =="") THEN
			path_prop = expandFilePath(Pathprops)
			FOR (i IN ChemName) // loop to transfer defined CEA labels to FORTRAN
				read_label(CEA_GAS_Name[i]) 
			END FOR
			// reading of CEA coefficients of all declared ChemName
			read_cea(setofSize(ChemName), path_prop, MW_ch, ier)
		END IF
		table_read(" ",ig, path_prop, ChemFileName[chem], ig, ier, fml, UsrFluiddH[chem])
		ng_tot = ng_tot + 1
		itab[ig] = ng_tot
	END IF
	ASSERT (ier < 100) FATAL "Problems reading thermo properties file"
	crit_cond(ig, prop_crit, ier)
	trip_l_cond(ig, prop_tril, ier)
	Pcrit = prop_crit[1]
	Tcrit = prop_crit[2]
	thermo_prop(ig, 1, P, 5, h, 111, prop, x, ier, INTERP_ORDER, ipx, ipy)
	T = prop[1,2]
	rho = prop[1,3]
	u = prop[1,7]
	s = prop[1,6]
	// Saturation properties
	IF (P < Pcrit) THEN // Subcritical
		h_l = prop[2,5]
		h_g = prop[3,5]
		x = (h-h_l)/(h_g-h_l)
		IF (x <= 0) THEN phase = liquid
		ELSEIF (x >= 1) THEN phase = vapor
		ELSE phase = two_phase
		END IF
	ELSE
		IF (T >= Tcrit) THEN 
			phase = supercritical
		ELSE
			phase = liquid
		END IF
		h_l = prop_crit[5]
		h_g = h_l
		x = 1
	END IF
	RETURN
END FUNCTION



FUNCTION NO_TYPE getThermodynamicState_ru (
	IN ENUM ChemName chem,
	IN REAL rho,
	IN REAL u,
	OUT REAL P,
	OUT REAL h,
	OUT REAL s,
	OUT REAL T,
	OUT REAL x,
	OUT ENUM Phase phase,
	OUT INTEGER iex,
	OUT INTEGER ipx,
	OUT INTEGER ipy
)
"Smaller than getStateRU. Used in ThermodynamicState class. Note that this returns the entropy also"
DECLS
	REAL Pcrit UNITS u_Pa "critical pressure"
	REAL prop[3,20]
	REAL prop_crit[20] "Fluid properties at critical point"
	REAL prop_tril[20] "Fluid properties at triple liquid point"
	REAL Tcrit "critical Temperature (K)"
	REAL h_l, h_g
	INTEGER ig "Real fluid number"
	INTEGER jx,jy
BODY
	ig = setofPos(ChemName, chem)
	IF(itab[ig] == 0) THEN
		IF(path_prop =="") THEN
			path_prop = expandFilePath(Pathprops)
			FOR (i IN ChemName) // loop to transfer defined CEA labels to FORTRAN
				read_label(CEA_GAS_Name[i]) 
			END FOR
			// reading of CEA coefficients of all declared ChemName
			read_cea(setofSize(ChemName), path_prop, MW_ch, ier)
		END IF
		table_read(" ",ig, path_prop, ChemFileName[chem], ig, ier, fml, UsrFluiddH[chem])
		ng_tot = ng_tot + 1
		itab[ig] = ng_tot
	END IF
	ASSERT (ier < 100) FATAL "Problems reading thermo properties file"
	crit_cond(ig, prop_crit, ier)
	trip_l_cond(ig, prop_tril, ier)
	Pcrit = prop_crit[1]
	Tcrit = prop_crit[2]
	thermo_prop(ig, 3, rho, 7, u, 111, prop, x, ier, INTERP_ORDER, jx, jy)
	P = prop[1,1]
	T = prop[1,2]
	s = prop[1,6]
	// Saturation properties
	IF (P < Pcrit) THEN // Subcritical
		h_l = prop[2,5]
		h_g = prop[3,5]
		x = (h-h_l)/(h_g-h_l)
		IF (x <= 0) THEN phase = liquid
		ELSEIF (x >= 1) THEN phase = vapor
		ELSE phase = two_phase
		END IF
	ELSE
		IF (T >= Tcrit) THEN 
			phase = supercritical
		ELSE
			phase = liquid
		END IF
		h_l = prop_crit[5]
		h_g = h_l
		x = 1
	END IF
	RETURN
END FUNCTION



FUNCTION REAL getMolWt (IN ENUM ChemName chem)
DECLS
	REAL M
	INTEGER ig
	INTEGER ier
BODY
	ig = setofPos(ChemName, chem)
	IF(itab[ig] == 0) THEN
		IF(path_prop =="") THEN
			path_prop = expandFilePath(Pathprops)
			FOR (i IN ChemName) // loop to transfer defined CEA labels to FORTRAN
				read_label(CEA_GAS_Name[i]) 
			END FOR
			// reading of CEA coefficients of all declared ChemName
			read_cea(setofSize(ChemName), path_prop, MW_ch, ier)
		END IF
	END IF
	IF (chem == R410A) THEN
		M = 72.58541424
	ELSEIF (chem == UsrDefined2) THEN
		M = 97.6 // R404A
	ELSE
		M = MW_ch[chem]
	END IF
	RETURN M
END FUNCTION



FUNCTION NO_TYPE getPartialDers (
	IN ENUM ChemName chem "Working chemical constituent",
	IN REAL p_bar "Pressure",
	IN REAL h "Enthalpy",
	OUT REAL drho_dP "drho/dp at constant h",
	OUT REAL drho_dh "drho/dh at constant P",
	OUT REAL x,
	OUT INTEGER ier,
	OUT INTEGER ipx,
	OUT INTEGER ipy
)
"Get partial derivatives of density wrt P and h"
DECLS
	REAL beta // isobaric thermal expansion coefficient
	REAL beta_g
	REAL beta_l
	REAL cp
	REAL dhl_dP,dhg_dP
	REAL dP_dT
	REAL dT_dP
	REAL dv_dh,dv_dP
	REAL dvl_dP,dvg_dP
	REAL dxHat_dP
	REAL kappa // isothermal compressibility coefficient
	REAL kappa_g
	REAL kappa_l
	REAL P UNITS u_Pa "Pressure"
	REAL Pcrit UNITS u_Pa "critical pressure"
	REAL prop[3,20]
	REAL prop_crit[20] "Fluid properties at critical point"
	REAL prop_tril[20] "Fluid properties at triple liquid point"
	REAL rho
	REAL T
	REAL Tcrit "critical Temperature (K)"
	REAL u
	REAL v_l,v_g
	REAL vprop
	REAL xHat
	INTEGER ig "Real fluid number"
OBJECTS
	SaturationProperties sat
BODY
	P = p_bar*1e5	
	// Set path
	ig = setofPos(ChemName, chem) // specifies index of fluid in list of ChemName
	IF(itab[ig] == 0) THEN
		IF(path_prop =="") THEN
			path_prop = expandFilePath(Pathprops)
			FOR (i IN ChemName) // loop to transfer defined CEA labels to FORTRAN
				read_label(CEA_GAS_Name[i]) 
			END FOR
			// reading of CEA coefficients of all declared ChemName
			read_cea(setofSize(ChemName), path_prop, MW_ch, ier)
		END IF
		WRITE("\nReading properties file %s.dat\n",ChemFileName[chem])
		table_read(" ",ig, path_prop, ChemFileName[chem], ig, ier, fml, UsrFluiddH[chem])
		ng_tot = ng_tot + 1
		itab[ig] = ng_tot
	END IF
	ASSERT (ier < 100) FATAL "Problems reading thermo properties file"
	
	// Critical/Triple properties
	crit_cond(ig, prop_crit, ier)
	trip_l_cond(ig, prop_tril, ier)
	Pcrit = prop_crit[1]
	Tcrit = prop_crit[2]
	
	// Get property matrix
	thermo_prop(ig, 1, P, 5, h, 111, prop, x, ier, 3, ipx, ipy)
	
	T = prop[1,2]
	rho = prop[1,3]
	u = prop[1,7]
	
	// Saturation properties
	IF (P < Pcrit) THEN
		// Subcritical
		sat.rho_l = prop[2,3]
		sat.h_l = prop[2,5]
		sat.cp_l = prop[2,11]

		sat.rho_g = prop[3,3]
		sat.h_g = prop[3,5]
		sat.cp_g = prop[3,11]
	ELSE
		// Supercritical
		IF (T >= Tcrit) THEN
			x = 1.
		ELSE
			x = 0.
		END IF
		sat.rho_l = prop_crit[3]
		sat.h_l = prop_crit[5]
		sat.cp_l = prop_crit[11]
		sat.rho_g = sat.rho_l
		sat.h_g = sat.h_l
		sat.cp_g = sat.cp_l
	END IF
	
	IF (x <= 0) THEN
	// Liquid
		cp = prop[1,11]
		beta = prop[1,14]
		kappa = prop[1,16]
		drho_dh = -beta * rho / cp
		drho_dP = (rho*kappa + beta/cp*(1.-beta*T))
	ELSEIF (x >= 1) THEN
	// Vapour
		cp = prop[1,11]
		beta = prop[1,14]
		kappa = prop[1,16]
		drho_dh = -beta * rho / cp
		drho_dP = (rho*kappa + beta/cp*(1-beta*T))
	ELSE
	// Two-Phase
		xHat = (h-sat.h_l)/(sat.h_g-sat.h_l)
		v_g = 1/sat.rho_g
		v_l = 1/sat.rho_l
		dT_dP = T*(v_g - v_l)/(sat.h_g - sat.h_l)
		dhl_dP = v_l*(1-beta_l*T) + sat.cp_l*dT_dP
		dhg_dP = v_g*(1-beta_g*T) + sat.cp_g*dT_dP
		dxHat_dP = (xHat*dhg_dP + (1-xHat)*dhl_dP) / (sat.h_l - sat.h_g)
		dv_dh = (v_g - v_l)/(sat.h_g - sat.h_l)
		dvl_dP = beta_l*v_l*dT_dP - kappa_l*v_l
		dvg_dP = beta_g*v_g*dT_dP - kappa_g*v_g
		dv_dP = dvl_dP + dxHat_dP*(v_g-v_l) + xHat*(dvg_dP - dvl_dP)
		drho_dh = -rho**2 * dv_dh
		drho_dP = -rho**2 * dv_dP
	END IF
	RETURN
END FUNCTION