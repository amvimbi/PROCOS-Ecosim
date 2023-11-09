/*-----------------------------------------------------------------------------------------
 LIBRARY: CO2
 FILE: CORRELATIONS
 AUTHOR: Viren Bhanot
 COMPANY: CERN
 CREATION DATE: 24/01/2017
 
 DESCRIPTION:
 Heat transfer coefficient correlations for single phase, two phase (boiling and 
 condensation) and air. Smoothing functions are available to smooth transitions between 
 single phase and two phase. Pressure drop correlations are not included in this library
 although they are available in Cryolib if necessary.
-----------------------------------------------------------------------------------------*/
USE MATH
USE CRYOLIB
USE THERMO_TABLE_INTERP



/*
 * INDEX:
 * 1. SINGLE PHASE CONVECTIVE
 * 	1. Gnielinski
 * 	2. DittusBoelter
 * 	3. Colburn
 * 2. TWO PHASE CONVECTIVE
 * 	1. Chen (simplified function call compared to Cryolib)
 * 	2. GungorWinterton (1986 and 1987)
 * 	3. JungRadermacher
 * 	4. Kandlikar
 * 	5. Dobson (Condensing)
 * 3. POOL BOILING
 * 	1. Cooper
 * 4. AIR HEAT TRANSFER
 * 	1. Kim-Webb
 * 	2. Chang-Wang
 * 	3. Churchill-Chu
 * 	4. McAdams
 * 	5. Churchill-Bernstein
 * 5. OVERALL HTC CALCULATIONS
 * 	1. getSmoothHTC
 	2. getSinglePhaseHTC
 * 	3. getTwoPhaseHTC
 * 	4. getCorrelationHTC
 */



--------------------------------------------------------------------------------
// SINGLE-PHASE CONVECTIVE HEAT TRANSFER
--------------------------------------------------------------------------------
ENUM SinglePhaseHTC = {UseGnielinski,UseDittusBoelter,UseColburn}



FUNCTION REAL getSinglePhaseHTC (
	IN ENUM SinglePhaseHTC choice,
	IN REAL A UNITS u_m2 "Cross sectional area",
	IN REAL cp UNITS u_J_kgK "Fluid specific heat capacity",
	IN REAL D UNITS u_m "Hydraulic diameter",
	IN REAL k UNITS u_W_mK "Fluid thermal conductivity",
	IN REAL mu UNITS u_Pas "Fluid dynamic viscosity",
	IN REAL m UNITS u_kg_s "Mass flow rate"
)
"Single-phase HTC calculations"
DECLS
 	REAL alpha
	REAL f
	REAL G "Mass flux"
	REAL Nu "Nusselt"
	REAL Pr "Prandtl"
	REAL Re "Reynolds"
BODY
	G = max(1E-6,abs(m)) / A
	Re = G * D / mu
	Pr = mu * cp / k
	IF (choice==UseDittusBoelter) THEN
		Nu = 0.023 * Re**0.8 * Pr**0.4
	ELSEIF (choice==UseGnielinski) THEN
		f = (1.58 * log(max(1e4,Re)) - 3.28)**-2
		Nu = 0.5 * f * (max(1E4,Re)-1000) * Pr / (1 + 12.7*(0.5*f)**0.5 * (Pr**0.667-1))
	ELSEIF (choice==UseColburn) THEN
		Nu = 0.023 * Pr**0.33 * Re**0.8
	END IF
	alpha = Nu * k / D
	RETURN max(1,alpha)
END FUNCTION



--------------------------------------------------------------------------------
// BOILING CORRELATIONS
--------------------------------------------------------------------------------
ENUM TwoPhaseHTC = {UseChen,UseDenglerAdams,UseGungorWinterton1986,UseGungorWinterton1987,UseJungRadermacher,UseKandlikar,UseLiuWinterton} // boiling correlations



FUNCTION REAL Chen (
	IN REAL A "Cross-sectional area",
	IN REAL D "Hydraulic diameter",
	IN REAL m "Mass flow rate",
	IN REAL q "Heat flux",
	IN REAL sigma "SUrface tension",
	IN REAL Tsat "Saturation temperature",
	IN REAL Tw "Wall temperature",
	IN REAL x "Vapour quality",
	OUT SaturationProperties sat
)
"Original Chen boiling correlation"
/*
 * Chen, J.C., 1966. Correlation for boiling heat transfer to saturated ﬂuids 
 * in convective ﬂow. Ind. Eng. Chem. Process Des. Dev. 5(3), 322-329.
 * Also refer to Chen and Fang "A note on the Chen correlation ofsaturated flow 
 * boiling heat transfer" 2014
 * FLUIDS: Mainly Water. Methanol, cyclohexane and pentane
 * OPERATION RANGE: 0.01<x<0.71, 54<G<4070, 6.3<q<2397.5E3
 *
 * NOTE: For now, I've used max(0,(Tw-Tsat))**0.99 in alpha_nuc otherwise the 
 * correlation crashes if we suddenly transition to condensation. This 
 * implementation is physically nonsense, but needed for numerical robustness)
 */
DECLS
	REAL alpha "HTC"
	REAL alpha_con "Convective boiling HTC"
	REAL alpha_nuc "Nucleate boiling HTC"
	REAL alpha_l "Liquid HTC"
	REAL G "Mass flux"
	REAL hfg "Latent heat"
	REAL B, F, S "Two-phase multipliers"
	REAL Re, Re_l "Reynolds"
	REAL vfg
	REAL xm "Range-limited vapour quality, constrained between 0 and 1"
	REAL Xtt_inv "Reciprocal Martinelli parameter"
BODY
	xm = max(0, min(0.99,x))
	G = abs(m)/A
	// Single-phase forced convective contribution
	Xtt_inv = (xm/(1-xm))**0.9 * (sat.rho_l/sat.rho_g)**0.5 * (sat.mu_g/sat.mu_l)**0.1
	Re_l = (1-xm)*G*D/sat.mu_l
	F = (1 + (2.35*(Xtt_inv + 0.213)**0.736)**20)**(1/20)
	alpha_l = getSinglePhaseHTC(UseDittusBoelter,A,sat.cp_l,D,sat.k_l,sat.mu_l,m)
	alpha_con = F * alpha_l
	// Nucleate boiling contribution
	hfg = max(sat.h_g - sat.h_l, 1E-6)
	vfg = 1/sat.rho_g - 1/sat.rho_l
	Re = Re_l * F**1.25
	S = 1 / (1 + 2.53E-6 * Re**1.17)
	B = 0.00122 * (sat.k_l**0.79 * sat.cp_l**0.45 * sat.rho_l**0.49 * S) / (max(sigma, 1e-4)**0.5 * sat.mu_l**0.29 * sat.rho_g**0.24)
	alpha_nuc = B * max(0,(Tw-Tsat))**0.99 * hfg**0.51 / (max(1e-4,vfg)*Tsat)**0.75
	alpha = alpha_con + alpha_nuc
	//NOTE: CRYOLIB fn call -> alpha = htc_tube_boiling_Chen (m,x,A,D,20,Tsat,Tsat,sat.h_l,sat.h_g,sat.k_l,sat.rho_l,sat.rho_g,sat.mu_l,sat.mu_g,sat.cp_l,sigma,Tw)
	RETURN alpha
END FUNCTION



FUNCTION REAL GungorWinterton (
	IN REAL A UNITS u_m2 "Cross-sectional area",
	IN REAL D UNITS u_m "Hydraulic diameter",
	IN REAL m UNITS u_kg_s "Mass flow rate",
	IN REAL MW UNITS u_kg "Molecular mass",
	IN REAL P_red "Reduced pressure",
	IN REAL q UNITS u_W_m2 "Heat flux",
	IN REAL x "Vapour quality",
	OUT SaturationProperties sat,
	IN BOOLEAN Year86=FALSE
)
"Includes both 1986 and 1987 correlations"
DECLS
	REAL hfg "Latent heat"
	REAL Xtt "Martinelli parameter"
	REAL G "Mass flux"
	REAL Bo "Boiling number"
	REAL Pr_l "Liquid Prandtl"
	REAL Re_l "Liquid Reynolds"
	REAL alpha_l "Liquid HTC"
	REAL Fr "Froude number"
	REAL E, E2, S, S2 // Multipliers
	REAL alpha_pool "Pool boiling HTC"
	REAL xm "Numerically corrected vapour quality"
	REAL alpha "HTC"
BODY
	hfg = max(sat.h_g - sat.h_l, 1E-6)
	xm = max(0.001, min(0.999, x)) // correlation fails if x is exactly zero or one
	Xtt = ((1-xm)/xm)**0.9 * (sat.rho_g/sat.rho_l)**0.5 * (sat.mu_l/sat.mu_g)**0.1
	G = max(1E-6,m) / A
	Bo = abs(q) / (G*hfg)
	Pr_l = sat.mu_l * sat.cp_l / sat.k_l
	Re_l = (1-xm)*G*D / sat.mu_l
	alpha_l = 0.023 * Re_l**0.8 * Pr_l**0.4 * sat.k_l / D
	IF (Year86==TRUE) THEN
		Fr = G ** 2 / (sat.rho_l**2 * D * 9.81)
		E = 1 + 24000 * Bo**1.16 + 1.37 * (1/Xtt)**0.86
		S = 1 / (1 + 1.15E-6 * E**2 * Re_l**1.17)
		IF (Fr < 0.05) THEN
			E2 = Fr ** (0.1 - 2.0 * Fr)
			S2 = Fr ** 0.5
			E = E * E2
			S = S * S2
		END IF
		alpha_pool = 55 * P_red**0.12 * (MW)**(-0.5) * abs(q)**0.67 * (-log10(P_red))**(-0.55)
		alpha = E*alpha_l + S*alpha_pool
	ELSE
		E = 1 + 3000 * Bo**0.86 + 1.12 * (xm/(1-xm))**0.75 * (sat.rho_l/sat.rho_g)**0.41
		alpha = E * alpha_l
	END IF
	RETURN alpha
END FUNCTION



FUNCTION REAL LiuWinterton (
	IN REAL A UNITS u_m2 "Cross-sectional area",
	IN REAL D UNITS u_m "Hydraulic diameter",
	IN REAL m UNITS u_kg_s "Mass flow rate",
	IN REAL MW UNITS u_kg "Molecular mass",
	IN REAL P_red "Reduced pressure",
	IN REAL q UNITS u_W_m2 "Heat flux",
	IN REAL x "Vapour quality",
	OUT SaturationProperties sat,
	IN BOOLEAN Year86=FALSE
)
/*
 * Liu and Winterton, "A general correlation for saturated and subcooled flow boiling
 * in tubes and annuli, based on a nucleate pool boiling equation",
 * Int. J Heat Mass Transfer, 1991.
 */
DECLS
	REAL alpha UNITS u_W_m2K "Overall HTC"
	REAL alpha_l UNITS u_W_m2K "Liquid HTC"
	REAL alpha_pool UNITS u_W_m2K "Nucleate boiling HTC"
	REAL F "Convective enhancement factor"
	REAL Fr_l "Liquid Froude"
	REAL G UNITS u_kg_sm2 "Mass flux"
	REAL Pr_l "Liquid Prandtl"
	REAL Re_l "Liquid Reynolds"
	REAL S "Supression factor"
	REAL xm "Numerically corrected vapour quality"
BODY
	xm = max(0.01,min(x,0.99))
	G = max(1E-6,m) / A
	Pr_l = sat.mu_l * sat.cp_l / sat.k_l
	Re_l = (1-xm)*G*D / sat.mu_l
	F = (1.0 + xm * Pr_l * (sat.rho_l/sat.rho_g - 1.0)) ** 0.35
	S = (1.0 + 0.055*F**0.1 * Re_l**0.16) ** -1.0
	Fr_l = G ** 2 / (sat.rho_l**2 * D * 9.81)
	IF (Fr_l<0.05) THEN
		F = F * Fr_l**(0.1 - 2.0*Fr_l)
		S = S * sqrt(Fr_l)
	END IF
	alpha_l = 0.023 * Re_l**0.8 * Pr_l**0.4 * sat.k_l/D
	alpha_pool = 55 * P_red**0.12 * (MW)**(-0.5) * abs(q)**0.67 * (-log10(P_red))**(-0.55)
	alpha = sqrt((F*alpha_l)**2 + (S*alpha_pool)**2)
	RETURN alpha
END FUNCTION



FUNCTION REAL JungRadermacher (
	IN REAL A UNITS u_m2 "Cross-sectional area",
	IN REAL D UNITS u_m "Hydraulic diameter",
	IN REAL m UNITS u_kg_s "Mass flow rate",
	IN REAL q UNITS u_W_m2 "Heat flux",
	IN REAL sigma "Surface tension",
	IN REAL T UNITS u_K "Fluid temperature",
	IN REAL x "Vapour quality",
	OUT SaturationProperties sat
)
"Designed primarily for refrigerant mixtures"
DECLS 
	REAL hfg
	REAL Xtt
	REAL G
	REAL Bo
	REAL Pr_l
	REAL Re_l
	REAL alpha_l
	REAL N
	REAL Fp
	REAL bd
	REAL beta
	REAL alpha_sa
	REAL xm
BODY
	hfg = sat.h_g - sat.h_l
	beta = 35
	xm = max(0.01,min(0.99,x))
	Xtt = ((1-xm)/xm)**0.9 * (sat.rho_g/sat.rho_l)**0.5 * (sat.mu_l/sat.mu_g)**0.1
	G = max(1e-6,m) / A
	Pr_l = sat.mu_l * sat.cp_l / sat.k_l
	Re_l = G * D / sat.mu_l
	alpha_l = 0.023 * Re_l**0.8 * Pr_l**0.4 * sat.k_l / D
	Fp = 2.37 * (0.29 + 1/Xtt)**0.85
	bd = 0.0146 * beta * (2 * sigma / (9.81 * (sat.rho_l - sat.rho_g)))**0.5
	alpha_sa = 207 * sat.k_l/bd * (abs(q)*bd / (sat.k_l*T))**0.745 * (sat.rho_g/sat.rho_l)**0.581 * Pr_l**0.533
	Bo = abs(q) / (G*hfg)
	N = 4048 * Xtt**1.22 * Bo**1.13
	RETURN N*alpha_sa + Fp*alpha_l
END FUNCTION



FUNCTION REAL Kandlikar (
	IN REAL A UNITS u_m2 "Cross-sectional area",
	IN REAL D UNITS u_m "Hydraulic diameter",
	IN REAL m UNITS u_kg_s "Mass flow rate",
	IN REAL q UNITS u_W_m2 "Heat flux",
	IN REAL x "Vapour quality",
	OUT SaturationProperties sat
)
DECLS
	REAL alpha_l "Liquid HTC"
	REAL alpha_nuc "Nucleate boiling HTC"
	REAL alpha_con "Forced convection HTC"
	REAL alpha "Overall HTC"
	REAL Bo "Boiling number"
	REAL C5_con, C5_nuc
	REAL Co "Convection number"
	REAL Ffl "Fluid dependent parameter"
	REAL Fr_l "Liquid Froude"
	REAL G "Mass flux"
	REAL hfg "Latent heat"
	REAL Pr_l "Liquid Prandtl"
	REAL Re_l "Liquid Reynolds"
	REAL xm
	CONST REAL C1_con = 1.1360
	CONST REAL C2_con = -0.9
	CONST REAL C3_con = 667.2
	CONST REAL C4_con = 0.7
	CONST REAL C1_nuc = 0.6683
	CONST REAL C2_nuc = -0.2
	CONST REAL C3_nuc = 1058.0
	CONST REAL C4_nuc = 0.7
BODY
	hfg = sat.h_g - sat.h_l
	xm = max(0.01, min(0.99, x))
	G = max(1e-6,m) / A
	Re_l = (1-xm) * G * D / sat.mu_l
	Pr_l = sat.mu_l * sat.cp_l / sat.k_l
	Co = ((1-xm)/xm)**0.8 * (sat.rho_g/sat.rho_l)**0.5
	Bo = q / (G*hfg)
	Fr_l = G**2 / (9.81*D*sat.rho_l**2)	
	IF (Fr_l < 0.04) THEN
		C5_con = 0.3
		C5_nuc = 0.2 // horizontal tubes
	ELSE
		C5_con = 0
		C5_nuc = 0 // vertical tubes
	END IF
	Ffl = 1 // fluid dependent parameter
	alpha_l = 0.023 * Re_l**0.8 * Pr_l**0.4 * sat.k_l / D
	alpha_con = alpha_l * (C1_con * Co**C2_con * (25*Fr_l)**C5_con + C3_con * Bo**C4_con * Ffl)
	alpha_nuc = alpha_l * (C1_nuc * Co**C2_nuc * (25*Fr_l)**C5_nuc + C3_nuc * Bo**C4_nuc * Ffl)
	alpha = alpha_nuc
	RETURN alpha
END FUNCTION



FUNCTION REAL DenglerAdams (
	IN REAL A UNITS u_m2 "Cross sectional area",
	IN REAL D UNITS u_m "Hydraulic diameter",
	IN REAL m UNITS u_kg_s "Mass flow rate",
	IN REAL x "Vapour quality",
	OUT SaturationProperties sat
)
"Two-phase HTC with single phase correction multiplier"
DECLS
	REAL xm "Numerically corrected vapour quality"
	REAL Xtt "Martinelli papraameter"
	REAL alpha_sp "Liquid HTC"
	REAL alpha "Two-Phase HTC"
	REAL a = 3.5 "Curve fit parameter"
	REAL b = 0.5 "Curve fit parameter"
BODY
	alpha_sp = getSinglePhaseHTC(UseDittusBoelter,A,sat.cp_l,D,sat.k_l,sat.mu_l,m)
	xm = max(0.01, min(0.99, x))
	Xtt = ((1 - xm) / xm)**0.9 * (sat.rho_g/sat.rho_l)**0.5 * (sat.mu_l/sat.mu_g)**0.1
	alpha = alpha_sp * a * (1/Xtt)**b
	RETURN alpha
END FUNCTION



--------------------------------------------------------------------------------
// CONDENSATION CORRELATIONS
--------------------------------------------------------------------------------
FUNCTION REAL Dobson (
	IN REAL A UNITS u_m2 "Cross-sectional area",
	IN REAL D UNITS u_m "Hydraulic diameter",
	IN REAL m UNITS u_kg_s "Mass flow rate",
	IN REAL q UNITS u_W_m2 "Heat flux",
	IN REAL x "Vapour quality",
	OUT SaturationProperties sat
)
"Dobson condensation HTC"
DECLS
	REAL alpha
	REAL G
	REAL Nu
	REAL Pr_l
	REAL Re_l
	REAL xm "Numerically corrected vapour quality"
	REAL Xtt
BODY
	xm = max(0.01,min(0.99,x))
	Xtt = ((1-xm)/xm)**0.9 * (sat.rho_l/sat.rho_g)**0.5 * (sat.mu_l /sat.mu_g)**0.1
	G = max(1e-6, m) / A
	Pr_l = sat.mu_l * sat.cp_l / sat.k_l
	Re_l = G * D * (1-xm) / sat.mu_l
	Nu = max(3.66, 0.023 * Re_l**0.8 * Pr_l**0.4 * (1 + 2.22 / Xtt**0.89))
	alpha = sat.k_l * Nu / D
	RETURN alpha
END FUNCTION



--------------------------------------------------------------------------------
// OVERALL TWO-PHASE HTC
--------------------------------------------------------------------------------
FUNCTION REAL getTwoPhaseHTC (
	IN ENUM TwoPhaseHTC choice,
	IN REAL A UNITS u_m2 "Cross-sectional area",
	IN REAL D UNITS u_m "Hydraulic diameter",
	IN REAL m UNITS u_kg_s "Mass flow rate",
	IN REAL MW UNITS u_kg "Fluid MOlecular mass",
	IN REAL P_red "Reduced pressure",
	IN REAL q "Heat flux",
	IN REAL sigma "Surface tension",
	IN REAL T UNITS u_K "Temperature",
	IN REAL Tw UNITS u_K "Wall temperature",
	IN REAL x "Vapour quality",
	OUT SaturationProperties sat,
	REAL Tb = 0.1 UNITS u_K "ΔT threshold for boiling HTC", // optional input
	REAL Tc = -0.1 UNITS u_K "ΔT threshold for condensation" // optional input
)
"Calculate HTC in two-phase region. Accounts for both boiling and condensation"
// slightly hacky implementation that prevents having to calculate BOTH boiling and condensation at each time step.
DECLS
	REAL alpha_boil UNITS u_W_m2K "Boiling HTC"
	REAL alpha_cond UNITS u_W_m2K "Condensation HTC"
	REAL alpha UNITS u_W_m2K "Overall HTC"
	REAL delT UNITS u_K
BODY
	delT = Tw - T
	
	// BOILING
	IF delT>Tc THEN	
		IF choice==UseChen THEN
			alpha_boil = Chen(A,D,m,q,sigma,T,Tw,x,sat)
		ELSEIF choice==UseDenglerAdams THEN
			alpha_boil = DenglerAdams(A,D,m,x,sat)
		ELSEIF choice==UseGungorWinterton1986 THEN
			alpha_boil = GungorWinterton(A,D,m,MW,P_red,q,x,sat,TRUE)
		ELSEIF choice==UseGungorWinterton1987 THEN
			alpha_boil = GungorWinterton(A,D,m,MW,P_red,q,x,sat,FALSE)
		ELSEIF choice==UseJungRadermacher THEN
			alpha_boil = JungRadermacher(A,D,m,q,sigma,T,x,sat)
		ELSEIF choice==UseKandlikar THEN
			alpha_boil = Kandlikar(A,D,m,q,x,sat)
		ELSEIF choice==UseLiuWinterton THEN
			alpha_boil = LiuWinterton(A,D,m,MW,P_red,q,x,sat)
		END IF
	ELSE
		alpha_boil = 0
	END IF
	
	// CONDENSATION
	IF delT<Tb THEN
		alpha_cond = Dobson(A,D,m,q,x,sat)
	ELSE
		alpha_cond = 0
	END IF
	
	// COMBINED
	IF delT>Tb THEN
		alpha = alpha_boil
	ELSEIF delT<Tc THEN
		alpha = alpha_cond
	ELSE
		alpha = alpha_boil*CRYOLIB.Smooth((delT-Tc)/(Tb-Tc)) + alpha_cond*(1 - CRYOLIB.Smooth((delT-Tc)/(Tb-Tc)))
	END IF
	RETURN alpha
END FUNCTION



--------------------------------------------------------------------------------
// POOL BOILING
--------------------------------------------------------------------------------
FUNCTION REAL Cooper (IN REAL delT, IN REAL M, IN REAL P_red, IN REAL q_flux)
"Pool boiling heat transfer coefficient"
BODY
	RETURN q_flux**0.67 * 55 * P_red**0.12 * M**-0.5 * (-log10(P_red))**-0.55
END FUNCTION



--------------------------------------------------------------------------------
// OVERALL HTC CALCULATIONS
--------------------------------------------------------------------------------
FUNCTION REAL getSmoothHTC (
	REAL alpha_g UNITS u_W_m2K "Vapour HTC",
	REAL alpha_l UNITS u_W_m2K "Liquid HTC",
	REAL alpha_tp UNITS u_W_m2K "Two-Phase HTC",
	REAL m_actual UNITS u_kg_s "Actual mass flow rate",
	REAL m_steady UNITS u_kg_s "Design mass flow rate",
	REAL x "Vapour quality",
	REAL x_low=0.1,
	REAL x_high=0.9,
	BOOLEAN mDotAdjustment = TRUE
)
"Spline interpolate HTC between fluid phases and optionally correct for mass flow rate"
// mDotSmoothing is only for Constant HTC conditions. When HTC calculated using
// correlations, no need to adjust for mass flow rate
DECLS
	REAL alpha
BODY
	/*
	IF (x<=0) THEN
		alpha = alpha_l
	ELSEIF (x>=1) THEN
		alpha = alpha_g
	ELSEIF (x>0 AND x<=x_low) THEN
		alpha = (-2*(alpha_tp-alpha_l)/x_low**3 * x**3) + 3/x_low**2*(alpha_tp-alpha_l)*x**2 + alpha_l
	ELSEIF (x>=x_high AND x<1) THEN
		alpha = 2*(alpha_g-alpha_tp)/(x_high-1)**3 * x**3 - 3*(x_high+1)*(alpha_g - alpha_tp) / (x_high-1)**3 * x**2 + ((-6 * (alpha_g - alpha_tp) / (x_high - 1) ** 3) + 6 * (x_high + 1) * (alpha_g - alpha_tp) / (x_high - 1) ** 3) * x + alpha_g - 2 * (alpha_g - alpha_tp) / (x_high - 1) ** 3 + 3 * (x_high + 1) * (alpha_g - alpha_tp) / (x_high - 1) ** 3 - ((-6 * (alpha_g - alpha_tp) / (x_high - 1) ** 3) + 6 * (x_high + 1) * (alpha_g - alpha_tp) / (x_high - 1) ** 3)
	ELSE
		alpha = alpha_tp
	END IF
	*/
	IF x<=0.5*(x_low+x_high) THEN
		alpha = spliceFunction(alpha_tp,alpha_l,x-x_low/2,x_low/2)
	ELSE
		alpha = spliceFunction(alpha_g,alpha_tp,x-(x_high+(1-x_high)/2),(1-x_high)/2)
	END IF
	IF (mDotAdjustment) THEN
		alpha = abs(m_actual/m_steady)**0.8 * alpha
	END IF
	RETURN alpha
END FUNCTION



FUNCTION REAL getCorrelationHTC (
	IN ENUM ChemName fluid "Refrigerant",
	IN ENUM SinglePhaseHTC choice_1p,
	IN ENUM TwoPhaseHTC choice_2p,
	IN REAL A UNITS u_m2 "Cross sectional area",
	IN REAL cp UNITS u_J_kgK "Fluid specific heat capacity", //single phase
	IN REAL D UNITS u_m "Hydraulic diameter",
	IN REAL k UNITS u_W_mK "Fluid thermal conductivity", //single phase
	IN REAL m UNITS u_kg_s "Mass flow rate",
	IN REAL mu UNITS u_Pas "Dynamic viscosity", //single phase
	IN REAL MW UNITS u_kg "Fluid Molecular weight",
	IN REAL P UNITS u_bar "Fluid pressure",
	IN REAL P_red "Reduced pressure",
	IN REAL q UNITS u_W_m2 "Heat flux",
	IN REAL sigma UNITS u_N_m "Surface tension",
	IN REAL T UNITS u_K "Fluid temperature",
	IN REAL Tw UNITS u_K "Wall temperature",
	IN REAL x "vapour quality",
	OUT SaturationProperties sat,
	IN REAL x_low = 0.1 "Transition region between liquid and two phase HTC", // optional input
	IN REAL x_high = 0.9 "Transition region between two-phase and vapour HTC" // optional input
)
"One stop shop for getting refrigerant-side HTC in every region"
// Note that the HTCs of all three regions are calculated at each time
// step, so simulations might be slowed down significantly"
DECLS
	REAL alpha_sp UNITS u_W_m2K
	REAL alpha_tp UNITS u_W_m2K
	REAL alpha UNITS u_W_m2K
BODY
	alpha_sp = getSinglePhaseHTC(choice_1p,A,cp,D,k,mu,abs(m))
	alpha_tp = getTwoPhaseHTC(choice_2p,A,D,abs(m),MW,P_red,q,sigma,T,Tw,x,sat)
	alpha = max(1,getSmoothHTC(alpha_sp,alpha_sp,alpha_tp,0,0,x,x_low,x_high,FALSE))
	RETURN alpha
END FUNCTION



--------------------------------------------------------------------------------
// AIR HEAT TRANSFER CORRELATIONS
--------------------------------------------------------------------------------
ENUM AirHTC = {UseKimWebb, UseChangWang, UseConstant} // fin-tube HX



FUNCTION REAL KimWebb (
	REAL Re,
	REAL Pr,
	REAL Gc,
	REAL CpAir,
	IN AirFlowPassageRecord geo
)
"For round tube plate fin heat exchangers. Only valid for < 3 banks of tubes"
DECLS
	REAL j1
	REAL j
	REAL St
	REAL alpha,alpha1
BODY
	ASSERT (geo.nBank<=3) FATAL "Kim-Webb correlation only applicable for up to 3 rows of tubes"
	j = 0.163 * Re**-0.369 * (geo.Pt/geo.Pl)**0.106 * (geo.Fp/geo.Dc)**0.0138 * (geo.Pt/geo.Dc)**0.13
	j = j * 1.043 * (Re**-0.14 * (geo.Pt/geo.Pl)**-0.564 * (geo.Fp/geo.Dc)**-0.123 * (geo.Pt/geo.Dc)**1.17)**(3-geo.nBank)
	St = j/Pr**(2/3)
	alpha1 = St*Gc*CpAir
	IF (geo.staggered==FALSE) THEN
		alpha = 0.7*alpha1  // Inline Tube HTC Factor
	ELSE 
		alpha = alpha1
	END IF
	RETURN alpha
END FUNCTION



FUNCTION REAL ChangWang (
	IN REAL Re,
	IN REAL Pr,
	IN REAL Gc,
	IN REAL cp_air,
	IN AirFlowPassageRecord geo
)
"For round tube, plate fin heat exchangers"
DECLS
	REAL j
	REAL Dh
	REAL P1, P2, P3, P4, P5, P6
	REAL St // Stanton
BODY
	Dh = 4*geo.A_min*geo.HX_l/geo.Ao
	IF (geo.nBank==1) THEN
		P1 = 1.9 - 0.23*log(Re)
		P2 = -0.236 + 0.126*log(Re)
		j = 0.108 * Re**-0.29 * (geo.Pt/geo.Pl)**P1 * (geo.Fp/geo.Dc)**-1.084 * (geo.Fp/Dh)**-0.786 * (geo.Fp/geo.Pt)**P2
	ELSE
		P3 = -0.361 - 0.042*geo.nBank/log(Re) + 0.158*log(geo.nBank*(geo.Fp/geo.Dc)**0.41)
	   	P4 = -1.224 - 0.076 * (geo.Pl/Dh)**1.42 / log(Re)
	   	P5 = -0.083 + 0.058*geo.nBank/log(Re)
	   	P6 = -5.735+1.2*log(Re/geo.nBank)
		j = 0.086*Re**P3 * geo.nBank**P4 * (geo.Fp/geo.Dc)**P5 * (geo.Fp/Dh)**P6 * (geo.Fp/geo.Pt)**-0.93
		St = j/Pr**(2/3)
	END IF
	St = j/Pr**(2/3)
	RETURN St*Gc*cp_air
END FUNCTION



FUNCTION REAL ChurchillChu (
	IN REAL Ra,
	IN REAL Pr,
	IN REAL k_air,
	IN REAL L "Characteristic length wrt gravity direction"
)
"Free convection on a vertical plate (both laminar and turbulent)"
DECLS
	REAL Nu
BODY
	ASSERT (Ra < 1E12) WARNING "ChurchillChu correlation is valid only up to Ra < 1E12"
	IF (Ra<1E9) THEN
		Nu = (0.680 + (0.670*Ra**(1/4)) / (1 + (0.492/Pr)**(9/16))**(4/9)) ** 2
	ELSE
		Nu = (0.825 + (0.387*Ra**(1/6)) / (1 + (0.492/Pr)**(9/16))**(8/27)) ** 2
	END IF
	RETURN k_air * Nu / L
END FUNCTION



FUNCTION REAL McAdams (
	IN REAL Ra_D,
	IN REAL Ra_L,
	IN REAL k_air,
	IN REAL L
)
"McAdams correlation for natural convection on hot horizontal circular plates"
DECLS
	REAL Nu
BODY
	ASSERT (Ra_L > 1E5 AND Ra_L < 3E10) WARNING "McAdams correlation is only valid within 1E5 < Ra < 3E10"
	IF (Ra_L < 2E7) THEN
		Nu = 0.14*Ra_D**0.5
	ELSE
		Nu = 0.54*Ra_D**0.24
	END IF
	RETURN k_air * Nu / L
END FUNCTION


FUNCTION REAL ChurchillBernstein (
	IN REAL Re,
	IN REAL Pr,
	IN REAL k_air,
	IN REAL L
)
"Churchill Bernstein correlation for forced convection on cylinder in crossflow"
DECLS
	REAL Nu
BODY
	ASSERT (Re < 1E4) WARNING "Churchill Bernstein is only valid until Re < 1E4"
	Nu = 0.3 + (0.62 * Re**0.5 * Pr**(1/3)) / (1 + (0.4/Pr)**(2/3))**0.25
	RETURN k_air * Nu / L
END FUNCTION