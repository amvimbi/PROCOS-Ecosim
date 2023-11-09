/*-----------------------------------------------------------------------------------------
 LIBRARY: CO2
 FILE: TWOFLUID
 AUTHOR: Viren Bhanot
 COMPANY: CERN
 DESCRIPTION: Two fluid models for more detailed Accumulator simulations
 CREATION DATE: 05/04/2023
-----------------------------------------------------------------------------------------*/
USE CRYOLIB
USE MATH
USE PORTS_LIB
USE THERMAL
USE THERMO_TABLE_INTERP


/* 
 * CONTENTS:
 * 1. Tank_TFM: Base model for TFM tanks (no inlet ports)
 * 
 */


// External Refprop Functions
"C++" FUNCTION REAL rp_rhoPh(IN REAL a, IN REAL b) IN "ecosimproRefprop.lib"
"C++" FUNCTION NO_TYPE rp_partialDers(IN REAL a, IN REAL b, OUT REAL c, OUT REAL d) IN "ecosimproRefprop.lib"
"C++" FUNCTION REAL rp_rhoPx(IN REAL a, IN REAL b) IN "ecosimproRefprop.lib"
"C++" FUNCTION REAL rp_hPx(IN REAL a, IN REAL b) IN "ecosimproRefprop.lib"
"C++" FUNCTION NO_TYPE rp_all(IN REAL P, IN REAL hhf, IN REAL hhg, OUT REAL rhof, OUT REAL rhog, OUT REAL hfs, OUT REAL hgs, OUT REAL ddphf, OUT REAL ddhpf, OUT REAL ddphg, OUT REAL ddhpg) IN "ecosimproRefprop.lib"

--------------------------------------------------------------------------------
// TWO FLUID MODELS
--------------------------------------------------------------------------------

/*
 * NOTE ON FLUID PROPERTIES:
 * All two-fluid models will only use Refprop for calculating fluid properties
 * This is because the properties change drastically around the saturation 
 * conditions and two fluid models *always* hover near saturation conditions. 
 * In particular, the partial derivative calculations can be off if using the 
 * look-up table.
 */
 
/*
 * NOTE ON OUTLET PORT VOID FRACTION
 * It is useful to think of the outlet port as it's own pipe. Based on the liquid
 * level in the accumulator vs the height of the pipe, there is a certain submersion
 * level of the outlet pipe. This leads to a certain void fraction inthe outlet
 * pipe. There is also a flow vapour quality, which is the ratio of the flow rate 
 * of vapour to the total flow rate. Using these two parameters, we can calculate
 * how much flow the vapour phase in the accu is losing, vs how much the liquid 
 * phase is losing.
*/


ABSTRACT COMPONENT TwoFluidTank (
	INTEGER nports_out = 1,
	BOOLEAN AccountforGravity = FALSE "If TRUE, adds static pressure at fluid ports",
	REAL D_out = 0.05 UNITS u_m "Outlet port diameter",
	REAL z_out = 0.02 UNITS u_m "Elevation of bottom of outlet port, 0 if underneath the vessel"
)
"Shell of the tank with only geometry and materials parameters"
PORTS 
	OUT fluid f_out[nports_out]
DATA
	// Geometry
	REAL D = 1 UNITS u_m "Accumulator inner diameter"
	REAL L = 2 UNITS u_m "Accumulator height"
	
	// Design parameters
	REAL tau_l = 1 UNITS u_s "Time constant for evaporation"
	REAL tau_v = 1 UNITS u_s "Time constant for condensation"
	// Initial Conditions
	ENUM InitialConditions2FM init = P_level "Specify initial condition variables"
	REAL P0 = 20 UNITS u_bar "Initial pressure"
	REAL h0 = 250000 UNITS u_J_kg "Initial enthalpy"
	REAL Tsat0 = 253.644825 UNITS u_K "Initial saturation temperature"
	REAL x0 = 0.1 UNITS no_units "Initial vapour quality"
	REAL Level0 = 20 UNITS u_pct "Initial liquid level in %"
DECLS
	REAL dP_out[nports_out] UNITS u_Pa "Outlet pressure change due to gravitational effects"
	DISCR REAL V UNITS u_m3 "Volume in m^3"
	DISCR REAL Vol UNITS u_dm3 "Volume in Liters"
	REAL z_l UNITS u_m "Elevation of liquid in the port with respect to z_out"
	REAL Ap, Apv, Apl UNITS u_m2 "Cross-sectional areas in port; total, vapor area, liquid area"
	DISCR REAL Pcrit UNITS u_bar "Critical pressure"
	//REAL rhocrit UNITS u_kg_m3 "Critical density"
	REAL drhol_dP,drhol_dh,drhov_dP,drhov_dh "Partial density derivatives"
	REAL gamma "Void fraction"
	REAL gamma_out "Void fraaction in outlet port"
	REAL h UNITS u_J_kg "Averaged enthalpy"
	REAL hl, hv UNITS u_J_kg "Liquid and vapour phase enthalpies"
	REAL hls, hvs "Liquid and vapour saturation enthalpies"
	REAL level "Percentage of liquid column inside the vessel"
	REAL level_frac "Fractional liquid level in the vessel"
	REAL Ml, Mv UNITS u_kg "Liquid and vapour phase masses"
	REAL M UNITS u_kg "Total refrigerant charge"
	REAL m_evap, m_cond UNITS u_kg_s "Evaporation and condensation flow rates"
	REAL mout_l, mout_v UNITS u_kg_s "Liquid and vapour portions of the outlet mass flow rate"
	REAL hout_l, hout_v UNITS u_kg_s "Liquid and vapour portions of the outlet specific enthalpy"
	REAL mhout_l, mhout_v UNITS u_W "Liquid and vapour portions of the outlet enthalpy flow rate"
	REAL P UNITS u_bar "Pressure"
	REAL P_Pa UNITS u_Pa "Pressure in Pascal"
	REAL rhol, rhov UNITS u_kg_m3 "Liquid and vapour densities"
	REAL rho UNITS u_kg_m3 "Density of liquid/vapor mix"
	REAL Tl, Tv UNITS u_K "Temperatures of liquid and vapour phases"
	REAL Tsat UNITS u_C
	REAL Vl, Vv UNITS u_m3 "Volume occupied by phases"
	REAL x "Vapour quality"
	REAL x_out "Vapour quality in outlet port depicted as its own pipe"
	//REAL vsound_l UNITS u_m_s "Velocity of sound in the liquid phase"
	REAL d1,d2 // dummy variables
	INTEGER ier,jx,jy "Error codes"
OBJECTS
	RefGeometryRecord geo
INIT
	geo.EndSurfaces = 2
	geo.setCylindricalGeometry(D,L)
	
	V = geo.V
	Vol = 1000*V
	//Initialization  
	Vl = (Level0/100)*L* MATH.PI*(D/2)**2
	P_Pa = P0*1e5
	assignInitAccu2FM(init,f_out[1].fluid,P0,Tsat0,Level0,x0,P,hl,hv,level_frac)
	Pcrit = CRYO_PF_CritProp_CORR(f_out[1].fluid,fprop_pressure,ier)                                
	//rhocrit = CRYO_PF_CritProp_CORR(f_out[1].fluid,fprop_density,ier) 
	f_out[1].is_C = TRUE // sets Accumulator as capacitive component
DISCRETE
	ASSERT (z_out >= 0 AND z_out <= L - D_out) FATAL "Outlet pipe height is greater than vessel height"	
CONTINUOUS

	P := P_Pa * 1e-5
	rp_all(P_Pa,hl,hv,rhol,rhov,d1,d2,drhol_dP,drhol_dh,drhov_dP,drhov_dh) // dummy vars d1/d2 used to prevent crash in supercrticial case
	
	Tsat = CRYO_PF_prop_vs_Px(f_out[1].fluid,P,0,fprop_temperature,ier,jx,jy)-273.15
	
	hls := rp_hPx(P_Pa,0)
	hvs := rp_hPx(P_Pa,1)
	
	Tl := CRYO_PF_prop_vs_ph(f_out[1].fluid,P,hl,fprop_temperature,ier,jx,jy)
	Tv := CRYO_PF_prop_vs_ph(f_out[1].fluid,P,hv,fprop_temperature,ier,jx,jy)
	
	Vv := geo.V - Vl
	Ml := rhol*Vl
	Mv := rhov*Vv
	M := Ml + Mv
	x := ZONE (P<Pcrit) Mv/M OTHERS 1
	rho := M/V
	h := x*hv + (1-x)*hl
	level_frac := Vl/V // fractional liquid level
	gamma := 1-level_frac // void fraction
	level := level_frac*100
	
	Ap := MATH.PI*(D_out/2)**2
	
	// Phase change flow rates
	m_evap := ZONE(hl<hls) 0 OTHERS tau_l*rhol*Vl*(hl-hls)/(hvs-hls)
	m_cond := ZONE(hv>hvs) 0 OTHERS tau_v*rhov*Vv*(hvs-hv)/(hvs-hls)
	
	// Port balance equations
	IF (z_out==0) INSERT
		z_l := D_out // Port submerged
	ELSE
		z_l := level - z_out
	END IF
	characterize_partially_submerged_circular_surface(z_l,D_out,Apl,Apv,gamma_out)
	x_out := gamma_out*rhov/((1-gamma_out)*rhol + gamma_out*rhov)
	
	mout_l := (1-x_out)*f_out[1].m
	mout_v := x_out*f_out[1].m 
	
	hout_l := ZONE(x_out>0 AND x_out<1 AND P<Pcrit) CRYO_PF_prop_vs_Px(f_out[1].fluid,f_out[1].P,0,fprop_enthalpy,ier,jx,jy) OTHERS f_out[1].hb
	hout_v := ZONE(x_out>0 AND x_out<1 AND P<Pcrit) CRYO_PF_prop_vs_Px(f_out[1].fluid,f_out[1].P,1,fprop_enthalpy,ier,jx,jy) OTHERS f_out[1].hb
	
	mhout_l := semiLinear(f_out[1].m,0,(1-x_out)*(hout_l-hl))      
	mhout_v := semiLinear(f_out[1].m,0,x_out*(hout_v-hv))	
	
	// Governing equations
<mal>	Vl*(drhol_dP*P_Pa' + drhol_dh*hl') + rhol*Vl' = m_cond - m_evap - mout_l									 	        
<mav>	Vv*(drhov_dP*P_Pa' + drhov_dh*hv') - rhov*Vl' = m_evap - m_cond - mout_v											
<enl>	Vl*(rhol*hl' - P_Pa') = m_cond*(hls-hl) - m_evap*(hvs-hl) - mhout_l
<env>	Vv*(rhov*hv' - P_Pa') = m_evap*(hvs-hv) - m_cond*(hls-hv) - mhout_v

	// Port equations
	dP_out[1] = (rhov*(L - max(z_out,level)) + rhol*(max(z_out,level) - z_out))*9.806
	f_out[1].P = ZONE(AccountforGravity) P + dP_out[1]/1e5 OTHERS P
	f_out[1].hf = x_out*hv+(1-x_out)*hl
	f_out[1].rho = gamma_out*rhov + (1-gamma_out)*rhol
	f_out[1].A = Ap
	
	f_out[1].Gcrit = 1 
	f_out[1].P_aux = P
	f_out[1].I = 0.5*geo.V/geo.As**2
	f_out[1].v = (1-gamma_out)*(1-x_out)*f_out[1].m/(rhol*Apl) + gamma_out*x_out*f_out[1].m/(rhov*Apv)
		
	//HOW TO DEAL WITH THE REST OF THE PORT PARAMETERS?	
	/*IF (z_out==0) INSERT
	// (A) Port always submerged
		vsound_l = CRYO_PF_prop_vs_ph(f_out[1].fluid,f_out[1].P,hl,fprop_vsound,ier,jx,jy) //  hl???  +  f_out[1].hf or f_out[1].hb
		f_out[1].Gcrit = 1 -- CRYO_Gcrit_fun(f_out[1].fluid,f_out.P,rhol,Tl,0,f_out[1].v,vsound_l,ier,jy)
		f_out[1].P_aux = P -- P + geo.V*rhol*vsound_l*1E-5/geo.As
		f_out[1].I = 0.5*geo.V/geo.As**2
		f_out[1].v = f_out[1].m/(rhol*f_out[1].A) 
	ELSE	
	// (B) Port possibly partially submerged
		f_out[1].Gcrit = 1 --CRYO_Gcrit_fun(f_out[1].fluid,P*1E5,rho,T,0,f_out[1].v,vsound,ier,jy)                                        // *** rho ?
		f_out[1].P_aux = P --P + geo.V*rho*vsound*1E-5/geo.As
		f_out[1].I = 0.5*geo.V/geo.As**2
		f_out[1].v = (1-gamma_out)*(1-x_out)*f_out[1].m/(rhol*Apl) + gamma_out*x_out*f_out[1].m/(rhov*Apv)  // Mixture velocity = (1-gamma_out)* velocity_liquid + gamma_out* velocity_vapor	
	END IF*/
END COMPONENT



COMPONENT TwoFluidBase_Prescribed IS_A TwoFluidTank	
PORTS 		 
	IN analog_signal(n=1) K_heat
	IN analog_signal(n=1) K_cool
DATA
	REAL alphaL_ref = 10000 UNITS u_W_m2K "Liquid to shell heat transfer coefficient"
	REAL alphaV_ref = 5000 UNITS u_W_m2K "Liquid to shell heat transfer coefficient"
	REAL alphaLV = 10 UNITS u_W_m2K "Heat transfer coefficient between liquid and vapour phases"
	REAL Q_heat_max = 45000 UNITS u_W "Maximum heating capacity"
	REAL Q_cool_max = 45000 UNITS u_W "Maximum cooling capacity"
	REAL tauAlphaL = 3 UNITS u_s "Time constant for low-pass HTC calculations for Liquid"
	REAL tauAlphaV = 3 UNITS u_s "Time constant for low-pass HTC calculations for Vapor"
DECLS
	CLOSE nports_out = 1
	REAL Pred UNITS u_bar "Reduced pressure"
	REAL Alv UNITS u_m2 "Contact area between liquid and vapor"
	REAL Q_coolPow UNITS u_W
	REAL Q_heatPow UNITS u_W
	REAL Qvl UNITS u_W "Heat transfer from vapor to liquid"
CONTINUOUS
	Pred := P/(Pcrit*1e-5)
	
	//Compute area of contact between liquid and vapor
	Alv := MATH.PI*(D/2)**2
				
	// Heat Transfer terms		
	Q_heatPow := K_heat.signal[1]/100*Q_heat_max
	Q_coolPow := K_cool.signal[1]/100*Q_cool_max
	
	Qvl := alphaLV * Alv * (Tv - Tl) // Heat transfer from vapor to liquid
	
<:enl>	Vl*(rhol*hl' - P_Pa') = m_cond*(hls-hl) - m_evap*(hvs-hl) - mhout_l + Q_heatPow + Qvl
<:env>	Vv*(rhov*hv' - P_Pa') = m_evap*(hvs-hv) - m_cond*(hls-hv) - mhout_v - Q_coolPow - Qvl
END COMPONENT



COMPONENT TwoFluidAccumulator_Prescribed (
	BOOLEAN AccountforGravity = FALSE,
	REAL D_out = 0.05 UNITS u_m "Out Port diameter",
	REAL z_out = 0.02 UNITS u_m "Out Port bottom elevation"
)
PORTS
	IN analog_signal(n=1) K_heat
	IN analog_signal(n=1) K_cool
	OUT fluid f_out
DATA 
	// Wall
	ENUM PipeMat mat = SS_304 "Shell material"
	REAL Do = 0.6 UNITS u_m "Shell outer diameter"
	REAL Th = 0.01 UNITS u_m "Shell thickness"
	REAL L = 1.8 UNITS u_m "Shell height"
	// Design parameters
	REAL tau_l = 1 UNITS u_s "Time constant for evaporation"
	REAL tau_v = 1 UNITS u_s "Time constant for condensation"
	REAL alphaL_ref = 10000 UNITS u_W_m2K "Liquid to shell heat transfer coefficient"
	REAL alphaV_ref = 5000 UNITS u_W_m2K "Liquid to shell heat transfer coefficient"
	REAL alphaLV = 10 UNITS u_W_m2K "Heat transfer coefficient between liquid and vapour phases"
	REAL Q_heat_max = 45000 UNITS u_W "Maximum heating capacity"
	REAL Q_cool_max = 45000 UNITS u_W "Maximum cooling capacity"
	REAL tauAlphaL = 3 UNITS u_s "Time constant for low-pass HTC calculations for Liquid"
	REAL tauAlphaV = 3 UNITS u_s "Time constant for low-pass HTC calculations for Vapor"
	// Initial Conditions
	ENUM InitialConditions2FM init = P_level "Specify initial condition variables"
	REAL P0 = 20 UNITS u_bar "Initial pressure"
	REAL h0 = 250000 UNITS u_J_kg "Initial enthalpy"
	REAL Tsat0 = 253.644825 UNITS u_K "Initial saturation temperature"
	REAL x0 = 0.1 UNITS no_units "Initial vapour quality"
	REAL Level0 = 20 UNITS u_pct "Initial liquid level in %"
	REAL Tw0 = 243.15 UNITS u_K "Initial wall temperature"
TOPOLOGY
	LumpedShell_insulated wall (
		mat = mat,
		Th = Th,
		Do = Do,	        
		L = L,
		Tw0 = Tw0
	)
	TwoFluidBase_Prescribed (AccountforGravity = AccountforGravity, D_out = D_out,	z_out = z_out) ref (
		D = Do - 2*Th, 
		L = L,
		tau_l = tau_l,
		tau_v = tau_v,
		init = init,
		P0 = P0,
		h0 = h0,
		Tsat0 = Tsat0,
		x0 = x0,
		Level0 = Level0,
		alphaL_ref = alphaL_ref,
		alphaV_ref = alphaV_ref,
		alphaLV = alphaLV,
		Q_heat_max = Q_heat_max,
		Q_cool_max = Q_cool_max,
		tauAlphaL = tauAlphaL,
		tauAlphaV = tauAlphaV
	)
	CONNECT K_heat TO ref.K_heat
	CONNECT K_cool TO ref.K_cool
	CONNECT ref.f_out[1] TO f_out
END COMPONENT














--*******************************************************************************
// From TEST_CODE.el
--*******************************************************************************
/*


--------------------------------------------------------------------------------
// ACCUMULATOR MODELS
--------------------------------------------------------------------------------
COMPONENT t_Accumulator_TFM (
	BOOLEAN refprop = TRUE "TRUE: use REFPROP equations. Might be slower",
	INTEGER choice = 1 "1: FullForm, 2: AlternativeFull, 3: Incompressible"
)
"Two-fluid accumulator model with outlet port"
PORTS
	IN analog_signal(n=1) Heat
	IN analog_signal(n=1) Cool
	OUT fluid f_out
	OUT analog_signal(n=1) Pressure "Pressure (bar)"
	OUT analog_signal(n=1) Level "Level (fractional)"
DATA
	REAL D = 0.1 UNITS u_m "Diameter"
	REAL L = 1.27324 UNITS u_m "Height"
	REAL Qheat_max = 10000 UNITS u_W "Maximum heating power"
	REAL Qcool_max = 45000 UNITS u_W "Maximum cooling power"
	REAL P0 = 20 UNITS u_bar "Initial pressure"
	REAL y0 = 0.3148 "Initial level"
DECLS
	CONST REAL alphafg = 0 UNITS u_W_m2K "Interphase heat transfer"
	CONST REAL eps = 1e-4 UNITS u_m3 "Minimum allowed volume for any phase"
	CONST REAL xf0 = 0.001
	CONST REAL xg0 = 0.999
	CONST REAL zHeat = 0.01
	CONST REAL zCool = 0.99
	DISCR REAL Pcrit UNITS u_bar "Critical pressure"
	DISCR REAL V UNITS u_m3 "Volume"
	DISCR REAL Ac UNITS u_m2 "Cross section area"
	REAL drhof_dP, drhog_dP
	REAL drhof_dh, drhog_dh
	REAL P UNITS u_bar "Pressure"
	REAL P_Pa UNITS u_Pa "Pascal pressure"
	REAL hf, hg
	REAL hfs, hgs, hfs_dummy, hgs_dummy
	REAL rhof, rhog UNITS u_kg_m3
	REAL M, Mg, Mf UNITS u_kg
	REAL Vf, Vg UNITS u_m3
	REAL m_evap, m_cond
	REAL xf, xg, x
	REAL Qh, Qc
	REAL Qhf, Qhg, Qcf, Qcg
	REAL x_out
	REAL mhf_out, mhg_out
	REAL y
	REAL level
	REAL rho
	REAL gamma
	REAL h
	REAL sf, sg
	REAL Qfg
	REAL Tf,Tg
	REAL hl_out, hv_out
	INTEGER ier,jx,jy "Error codes"
INIT
	Ac = MATH.PI * (D/2)**2
	V = Ac*L
	P_Pa = P0*1E5
	Pcrit = CRYO_PF_CritProp_CORR(f_out.fluid,fprop_pressure,ier) //CRYO_PF_CritProp(f_out.fluid,fprop_pressure,ier)*1e-5  
	Vf = V*y0
	Vg = V-Vf
	y = y0
	IF (refprop) THEN
		hf = rp_hPx(P0*1e5,0)
		hg = rp_hPx(P0*1e5,1)
	ELSE
		hf = CRYO_PF_prop_vs_Px(f_out.fluid,P0,0,fprop_enthalpy,ier,jx,jy)
		hg = CRYO_PF_prop_vs_Px(f_out.fluid,P0,1,fprop_enthalpy,ier,jx,jy)
	END IF
CONTINUOUS
	IF (refprop) INSERT
		rp_all(P_Pa,hf,hg,rhof,rhog,hfs_dummy,hgs_dummy,drhof_dP,drhof_dh,drhog_dP,drhog_dh)
		hfs := ZONE(P<Pcrit) hfs_dummy OTHERS CO2_CRITENTH
		hgs := ZONE(P<Pcrit) hgs_dummy OTHERS CO2_CRITENTH
		xf := max(0,min(1,(hf-hfs)/(hgs-hfs)))
		xg := max(0,min(1,(hg-hfs)/(hgs-hfs)))
	ELSE
		getPartialDers(f_out.fluid,P,hf,drhof_dP,drhof_dh,xf,ier,jx,jy)
		getPartialDers(f_out.fluid,P,hg,drhog_dP,drhog_dh,xg,ier,jx,jy)
		hfs = ZONE (P<Pcrit) CRYO_PF_prop_vs_Px(f_out.fluid,P,0,fprop_enthalpy,ier,jx,jy) OTHERS CO2_CRITENTH
		hgs = ZONE (P<Pcrit) CRYO_PF_prop_vs_Px(f_out.fluid,P,1,fprop_enthalpy,ier,jx,jy) OTHERS CO2_CRITENTH
		rhof = CRYO_PF_prop_vs_ph(f_out.fluid,P,hf,fprop_density,ier,jx,jy)
		rhog = CRYO_PF_prop_vs_ph(f_out.fluid,P,hg,fprop_density,ier,jx,jy)
	END IF
	
	P = P_Pa * 1e-5
	
	Tf := CRYO_PF_prop_vs_ph(f_out.fluid,P,hf,fprop_temperature,ier,jx,jy)
	Tg := CRYO_PF_prop_vs_ph(f_out.fluid,P,hg,fprop_temperature,ier,jx,jy)
	sf := CRYO_PF_prop_vs_ph(f_out.fluid,P,hf,fprop_entropy,ier,jx,jy)
	sg := CRYO_PF_prop_vs_ph(f_out.fluid,P,hg,fprop_entropy,ier,jx,jy)

	Vf := V*y
	Vg := V - Vf
	level := Vf/V
	gamma = 1-level
	
	Mg := rhog*Vg
	Mf := rhof*Vf
	M := Mf + Mg
	x := Mg/M
	
	h := x*hg + (1-x)*hf
	rho:= (1-gamma)*rhof + gamma*rhog

	m_evap := ZONE(Vf<eps OR P>Pcrit) 0 OTHERS max(0,rhof*Vf*(xf-xf0))
	m_cond := ZONE(Vg<eps OR P>Pcrit) 0 OTHERS max(0,rhog*Vg*(xg0-xg))

	Qh := Heat.signal[1]/100 * Qheat_max
	Qc := Cool.signal[1]/100 * Qcool_max
	
	Qhf = ZONE(level>=zHeat) Qh OTHERS 0
	Qhg = ZONE(level>=zHeat) 0 OTHERS Qh
	Qcf = ZONE(level<=zCool) 0 OTHERS Qc
	Qcg = ZONE(level<=zCool) Qc OTHERS 0
	
	Qfg := alphafg*Ac*(Tg-Tf)
		
	x_out := ZONE(f_out.m>=0) 0 OTHERS max(0,min(1,CRYO_PF_prop_vs_ph(f_out.fluid,f_out.P,f_out.hb,fprop_quality,ier,jx,jy)))
	hl_out := ZONE(x_out<1 AND x_out>0 AND P<Pcrit) CRYO_PF_prop_vs_Px(f_out.fluid,f_out.P,0,fprop_enthalpy,ier,jx,jy)
		OTHERS f_out.hb
	hv_out := ZONE(x_out>0 AND x_out<1 AND P<Pcrit) CRYO_PF_prop_vs_Px(f_out.fluid,f_out.P,1,fprop_enthalpy,ier,jx,jy)
		OTHERS f_out.hb

<mf>	Vf*drhof_dP*P_Pa' + rhof*V*y' + Vf*drhof_dh*hf' = m_cond - m_evap - f_out.m*(1-x_out)
<mg>	Vg*drhog_dP*P_Pa' - rhog*V*y' + Vg*drhog_dh*hg' = m_evap - m_cond - f_out.m*x_out

<en>	IF (choice==1) INSERT
		// Full Form
		mhf_out = ZONE(f_out.m>=0) f_out.m*hf OTHERS f_out.m*(1-x_out) * hl_out
		mhg_out = ZONE(f_out.m>=0) 0 OTHERS f_out.m*x_out * hv_out
		(rhof*hf-P_Pa)*V*y' + Vf*((hf*drhof_dP-1)*P_Pa' + (hf*drhof_dh+rhof)*hf') = m_cond*hfs - m_evap*hgs + Qhf - Qcf - mhf_out + Vf*P_Pa' + Qfg
		(P_Pa-rhog*hg)*V*y' + Vg*((hg*drhog_dP-1)*P_Pa' + (hg*drhog_dh+rhog)*hg') = m_evap*hgs - m_cond*hfs + Qhg - Qcg - mhg_out + Vg*P_Pa' - Qfg
	ELSEIF (choice==2) INSERT
		// Alternative Full Form
		mhf_out = ZONE(f_out.m>=0) f_out.m*(hf-(hf-P_Pa/rhof)) OTHERS f_out.m*(1-x_out)*(hl_out-(hf-P_Pa/rhof))
		mhg_out = ZONE(f_out.m>=0) 0 OTHERS f_out.m*x_out*(hv_out-(hg-P_Pa/rhog))
		Vf*(P_Pa/rhof*drhof_dP-1)*P_Pa' + Vf*(P_Pa/rhof*drhof_dh+rhof)*hf' = m_cond*(hfs-(hf-P_Pa/rhof)) - m_evap*(hgs-(hf-P_Pa/rhof)) - mhf_out + Qhf - Qcf + Vf*P_Pa' + Qfg
		Vg*(P_Pa/rhog*drhog_dP-1)*P_Pa' + Vg*(P_Pa/rhog*drhog_dh+rhog)*hg' = m_evap*(hgs-(hg-P_Pa/rhog)) - m_cond*(hfs-(hg-P_Pa/rhog)) - mhg_out + Qhg - Qcg + Vg*P_Pa' - Qfg
	ELSE
		// Incompressible Form
		mhf_out = ZONE(f_out.m>=0) f_out.m*(hf-hf) OTHERS f_out.m*(1-x_out)*(hl_out-hf)
		mhg_out = ZONE(f_out.m>=0) 0 OTHERS f_out.m*x_out*(hv_out-hg)
		Vf*(rhof*hf' - P_Pa') = m_cond*(hfs-hf) - m_evap*(hgs-hf) - mhf_out + Qhf - Qcf + Vf*P_Pa' + Qfg
		Vg*(rhog*hg' - P_Pa') = m_evap*(hgs-hg) - m_cond*(hfs-hg) - mhg_out + Qhg - Qcg + Vg*P_Pa' - Qfg
	END IF
	
	Level.signal[1] = level*100
	Pressure.signal[1] = P
	f_out.P = P
	f_out.P_aux = P
	f_out.hf = hf
	f_out.rho = rhof
	f_out.v = 1
	f_out.A = 1
	f_out.Gcrit = 1
END COMPONENT



COMPONENT t_AccumulatorTFM_v IS_A t_Accumulator_TFM
PORTS
	IN fluid f_in "Spray inlet"
CONTINUOUS
<:mf>	Vf*drhof_dP*P_Pa' + rhof*V*y' + Vf*drhof_dh*hf' = m_cond - m_evap - f_out.m*(1-x_out) + f_in.m
<:mg>	Vg*drhog_dP*P_Pa' - rhog*V*y' + Vg*drhog_dh*hg' = m_evap - m_cond - f_out.m*x_out

<:en>	IF (choice==1) INSERT
		// Full Form
		mhf_out = ZONE(f_out.m>=0) f_out.m*hf OTHERS f_out.m*(1-x_out) * hl_out
		mhg_out = ZONE(f_out.m>=0) 0 OTHERS f_out.m*x_out * hv_out
		(rhof*hf-P_Pa)*V*y' + Vf*((hf*drhof_dP-1)*P_Pa' + (hf*drhof_dh+rhof)*hf') = m_cond*hfs - m_evap*hgs + Qh - mhf_out + Vf*P_Pa' + Qfg + f_in.m*(hfs-hf)
		(P_Pa-rhog*hg)*V*y' + Vg*((hg*drhog_dP-1)*P_Pa' + (hg*drhog_dh+rhog)*hg') = m_evap*hgs - m_cond*hfs - Qc - mhg_out + Vg*P_Pa' - Qfg - f_in.m*(hfs-f_in.hf)
	ELSEIF (choice==2) INSERT
		// Alternative Full Form Derivation
		mhf_out = ZONE(f_out.m>=0) f_out.m*(hf-(hf-P_Pa/rhof)) OTHERS f_out.m*(1-x_out)*(hl_out-(hf-P_Pa/rhof))
		mhg_out = ZONE(f_out.m>=0) 0 OTHERS f_out.m*x_out*(hv_out-(hg-P_Pa/rhog))
		Vf*(P_Pa/rhof*drhof_dP-1)*P_Pa' + Vf*(P_Pa/rhof*drhof_dh+rhof)*hf' = m_cond*(hfs-(hf-P_Pa/rhof)) - m_evap*(hgs-(hf-P_Pa/rhof)) - mhf_out + Qh + Vf*P_Pa' + Qfg
		Vg*(P_Pa/rhog*drhog_dP-1)*P_Pa' + Vg*(P_Pa/rhog*drhog_dh+rhog)*hg' = m_evap*(hgs-(hg-P_Pa/rhog)) - m_cond*(hfs-(hg-P_Pa/rhog)) - mhg_out - Qc + Vg*P_Pa' - Qfg
	ELSE
		// Incompressible Form
		mhf_out = ZONE(f_out.m>=0) f_out.m*(hf-hf) OTHERS f_out.m*(1-x_out)*(hl_out-hf)
		mhg_out = ZONE(f_out.m>=0) 0 OTHERS f_out.m*x_out*(hv_out-hg)
		Vf*(rhof*hf' - P_Pa') = m_cond*(hfs-hf) - m_evap*(hgs-hf) - mhf_out + Qh + Vf*P_Pa' + Qfg
		Vg*(rhog*hg' - P_Pa') = m_evap*(hgs-hg) - m_cond*(hfs-hg) - mhg_out - Qc + Vg*P_Pa' - Qfg
	END IF
	
	f_in.fluid = f_out.fluid
	f_in.n_fluid = f_out.n_fluid
	f_in.P = P
	f_in.P_aux = P
	f_in.hb = hg
	f_in.rho = rhog
	f_in.v = 1
	f_in.A = 1
	f_in.Gcrit = 1
END COMPONENT



COMPONENT test_AccumulatorTFM (
	BOOLEAN useRefprop = FALSE "If TRUE, expect slower performance"
)
"Two-fluid accumulator model with outlet port"
PORTS
	OUT fluid f_out
	IN analog_signal(n=1) Heat
	IN analog_signal(n=1) Cool
	OUT analog_signal(n=1) Level "Liquid level in %"
	OUT analog_signal(n=1) Pressure "Pressure in Pa"
DATA

	REAL D = 0.5 UNITS u_m "Diameter"
	REAL L = 1.96 UNITS u_m "Height"
	REAL z_out = 0.1 UNITS u_m "Elevation of bottom of outlet port from Accumulator bottom"
	REAL D_out = 0.01 UNITS u_m "Diameter of outlet port"
	REAL z_heat = 0.1 UNITS u_m "Elevation of heater from Accumulator bottom"
	REAL z_cool = 1.8 UNITS u_m "Elevation of cooling coil from Accumulator bottom"

	REAL tau_f = 1 UNITS u_Hz "Time constant for evaporation flow rate, leave default if unsure"
	REAL tau_g = 1 UNITS u_Hz "Time constant for condensation flow rate"
	REAL Qheat_max = 100e3 UNITS u_W "Maximum heating power"
	REAL Qcool_max = 100e3 UNITS u_W "Maximum cooling power"

	REAL P0 = 20 UNITS u_bar "Initial pressure"
	REAL y0 = 0.1 "Initial level as fraction of 1"
DECLS
	CONST REAL xf0 = 0.01
	CONST REAL xg0 = 0.999
	CONST REAL eps = 1e-4 UNITS u_m3 "Minimum volume for either phase"

	DISCR REAL V UNITS u_m3 "Volume"
	DISCR REAL Vol UNITS u_dm3 "Volume in Litres"
	DISCR REAL Ac UNITS u_m2 "Cross-sectional area"
	DISCR REAL Pcrit UNITS u_bar "Critical pressure"
	DISCR REAL rhocrit UNITS u_kg_m3 "Critical density"
	
	REAL drhof_dP, drhog_dP, drhof_dh, drhog_dh "Partial density derivatives"
	REAL gamma "Void fraction"
	REAL h UNITS u_J_kg "Averaged enthalpy"
	REAL hf, hg "Liquid and vapour phase enthalpies"
	REAL hfs, hgs "Saturation liquid and vapour enthalpies"
	REAL lvl UNITS no_units "Fractional liquid level, from 0 to 1"
	REAL level_abs UNITS u_m "Absolute liquid level in meter"
	REAL M UNITS u_kg "Refrigerant charge"
	REAL Mf, Mg UNITS u_kg "Mass of the two phases"
	REAL m_evap, m_cond "Evaporation and condensation flow rates"
	REAL mout_f, mout_g "Portions of inlet flow rate going to liquid and vapour phases" 
	REAL mhout_f, mhout_g "Portion of inlet enthalpy flow rate going to liquid and vapour phases"
	REAL P UNITS u_bar "Pressure"
	REAL P_Pa UNITS u_Pa "Pressure in Pascal"
	REAL Qh, Qc UNITS u_W "Overall heating and cooling powers"
	REAL Qfh, Qgh, Qfc, Qgc UNITS u_W "Portion of heating/cooling powers being applied to liquid/vapour phases"
	REAL rho UNITS u_kg_m3 "Averaged density"
	REAL rhof, rhog UNITS u_kg_m3 "Densities"
	REAL Tf, Tg UNITS u_K "Fluid temperatures"
	REAL Tsat UNITS u_K "Saturation temperature"
	REAL Vf, Vg UNITS u_m3 "Volume occupied by phases"
	REAL x "Vapour quality"
	REAL xf, xg "Local vapour quality within each phase"
	REAL dum1,dum2,dum3,dum4 // Dummy variables
	INTEGER ier,jx,jy "Error codes"

	REAL H, U
	REAL uf, ug UNITS u_J_kg
INIT
	Ac = MATH.PI * (D/2)**2
	V = Ac*L
	Vol = V*1000

	P_Pa = P0*1e5
	Vf = V*y0
	Vg = V-Vf
	hf = rp_hPx(P0*1e5,0)
	hg = rp_hPx(P0*1e5,1)

	Pcrit = CRYO_PF_CritProp(f_out.fluid,fprop_pressure,ier)*1e-5 // same between refprop and lookup table
	rhocrit = CRYO_PF_CritProp(f_out.fluid,fprop_density,ier)
	f_out.is_C = TRUE
CONTINUOUS
	P = P_Pa * 1e-5
	V = Vf + Vg
	
	lvl = spliceFunction(Vf/V,0,Pcrit-P,0.01)
	gamma = 1-lvl
	level_abs = lvl*L
	
	Mf := rhof*Vf
	Mg := rhog*Vg
	M := Mf + Mg
	x := Mg/M
	
	Tsat := CRYO_PF_prop_vs_Px(f_out.fluid,P,0,fprop_temperature,ier,jx,jy) // spits out non-garbage at P>Pcrit
	rho = lvl*rhof + gamma*rhog
	h = (Mf*hf + Mg*hg)/M
	
	SEQUENTIAL
	IF (P<Pcrit) THEN
		rp_all(P_Pa,hf,hg,rhof,rhog,dum1,dum2,drhof_dP,drhof_dh,drhog_dP,drhog_dh) // save some time and calculate all together
	ELSE
		rp_partialDers(P_Pa,hf,drhof_dP,drhof_dh)
		rp_partialDers(P_Pa,hg,drhog_dP,drhog_dh)
		rhof = rp_rhoPh(P_Pa,hf)
		rhog = rp_rhoPh(P_Pa,hg)
	END IF
	END SEQUENTIAL
	hfs := spliceFunction(rp_hPx(P_Pa,0), CO2_CRITENTH, Pcrit-P, 0.1)
	hgs := spliceFunction(rp_hPx(P_Pa,1), CO2_CRITENTH, Pcrit-P, 0.1)
	
	Tf:= CRYO_PF_prop_vs_ph(f_out.fluid,P,hf,fprop_temperature,ier,jx,jy)
	Tg:= CRYO_PF_prop_vs_ph(f_out.fluid,P,hg,fprop_temperature,ier,jx,jy)
	
	div_safe(hf-hfs, hgs-hfs, dum3) // avoids division by zero errors when hgs==hfs
	div_safe(hg-hfs, hgs-hfs, dum4)
	xf := ZONE(P<Pcrit) max(0,min(1,dum3)) OTHERS 0
	xg := ZONE(P<Pcrit) max(0,min(1,dum4)) OTHERS 1
	
	m_evap := ZONE(Vf<eps OR P>Pcrit) 0 OTHERS max(0,tau_f*rhof*Vf*(xf-xf0))
	m_cond := ZONE(Vg<eps OR P>Pcrit) 0 OTHERS max(0,tau_g*rhog*Vg*(xg0-xg))
	
	Qh := Heat.signal[1]/100 * Qheat_max
	Qc := Cool.signal[1]/100 * Qcool_max
	
	Qfh := spliceFunction(Qh, 0, level_abs-z_heat, 0.001*L)
	Qgh := spliceFunction(0, Qh, level_abs-z_heat, 0.001*L)
	Qfc := spliceFunction(Qc, 0, level_abs-z_cool, 0.001*L)
	Qgc := spliceFunction(0, Qc, level_abs-z_cool, 0.001*L)

	// (assume any two-phase inlet flow is well-mixed and not stratified)
	mout_f = spliceFunction(f_out.m, 0, level_abs-(z_out+D_out/2), D_out/2)
	mout_f + mout_g = f_out.m
	mhout_f = donor_cell(mout_f, mout_f*f_out.hf, mout_f*f_out.hb)
	mhout_g = donor_cell(mout_g, mout_g*f_out.hf, mout_g*f_out.hb)
	
<maf>	Vf*(drhof_dP*P_Pa' + drhof_dh*hf') + rhof*Vf' = m_cond - m_evap - mout_f
<mag>	Vg*(drhog_dP*P_Pa' + drhog_dh*hg') - rhog*Vf' = m_evap - m_cond - mout_g
	
<enf>	(rhof*hf-P_Pa)*Vf' + Vf*((hf*drhof_dP-1)*P_Pa' + (hf*drhof_dh+rhof)*hf') + P_Pa*Vf' = m_cond*hfs - m_evap*hgs + Qfh - Qfc - f_out.mh
<eng>	(P_Pa-rhog*hg)*Vf' + Vg*((hg*drhog_dP-1)*P_Pa' + (hg*drhog_dh+rhog)*hg') - P_Pa*Vf' = m_evap*hgs - m_cond*hfs + Qgh - Qgc
	
	uf := CRYO_PF_prop_vs_ph(f_out.fluid,P,hf,fprop_energy,ier,jx,jy)
	ug := CRYO_PF_prop_vs_ph(f_out.fluid,P,hg,fprop_energy,ier,jx,jy)
	H = hf*Mf + hg*Mg
	U = uf*Mf + ug*Mg
	
	Pressure.signal[1] = P
	Level.signal[1] = lvl*100
	f_out.P = P
	f_out.P_aux = P
	f_out.hf = spliceFunction(hf, hg, level_abs-(z_out+D_out/2), D_out/2)
	f_out.rho = rhof
	f_out.v = 1.0
	f_out.A = MATH.PI * (D_out/2.0)**2.0
	f_out.Gcrit = 1.0
END COMPONENT


COMPONENT test_AccumulatorTFM_pini
PORTS
	OUT fluid f_out
	IN analog_signal(n=1) Heat
	IN analog_signal(n=1) Cool
	OUT analog_signal(n=1) Level "Liquid level in %"
	OUT analog_signal(n=1) Pressure "Pressure in Pa"
DATA
	REAL D = 0.5 UNITS u_m "Diameter"
	REAL L = 1.96 UNITS u_m "Height"
	REAL D_out = 0.005 UNITS u_m "Outlet port diameter"
	REAL z_out = 1.7 UNITS u_m "Outlet port height"
	REAL Qheat_max = 100e3 UNITS u_W "Maximum heating power"
	REAL Qcool_max = 100e3 UNITS u_W "Maximum cooling power"
	REAL tau_f = 1 UNITS u_s "Time constant for evaporation"
	REAL tau_g = 1 UNITS u_s "Time constant for condensation"
	REAL alpha_fg = 10 UNITS u_W_m2K "Heat transfer coefficient between liquid and vapour phases"
	REAL P0 = 20 UNITS u_bar "Initial pressure"
	REAL y0 = 0.1 "Initial level as fraction of 1"
DECLS
	DISCR REAL V UNITS u_m3 "Volume"
	DISCR REAL Vol UNITS u_dm3 "Volume in Litres"
	DISCR REAL Ac UNITS u_m2 "Cross-sectional area"
	DISCR REAL Pcrit UNITS u_bar "Critical pressure"
	DISCR REAL rhocrit UNITS u_kg_m3 "Critical density"
	REAL drhof_dP,drhof_dh,drhog_dP,drhog_dh
	REAL h UNITS u_J_kg "Averaged enthalpy"
	REAL hf, hg UNITS u_J_kg "Liquid and vapour phase enthalpies"
	REAL hfs, hgs "Saturation liquid and vapour enthalpies"
	REAL lvl, level
	REAL Mf, Mg UNITS u_kg "Mass of the two phases"
	REAL M UNITS u_kg "Total mass in vessel"
	REAL m_evap, m_cond "Evaporation and condensation flow rates"
	REAL mh_out UNITS u_W "Outlet port enthalpy flow rate"
	REAL Tsat UNITS u_K
	REAL P UNITS u_bar "Pressure"
	REAL P_Pa UNITS u_Pa "Pressure in Pascal"
	REAL Qh, Qc UNITS u_W "Overall heating and cooling powers"
	REAL Qfg UNITS u_W "Heat transfer between liquid and vapour phases"
	REAL rhof, rhog UNITS u_kg_m3 "Densities"
	REAL Tf, Tg UNITS u_K "Temperatures of liquid and vapour phases"
	REAL Vf, Vg UNITS u_m3 "Volume occupied by phases"
	REAL x "Vapour quality"
	REAL d1,d2
	INTEGER ier,jx,jy "Error codes"
	
	REAL delta = 1000 UNITS u_J_kg "Minimum distance from the saturated line" // prevent phase disappearing
	
INIT
	Ac = MATH.PI * (D/2)**2
	V = Ac*L
	Vol = V*1000

	P_Pa = P0*1e5
	Vf = V*y0
	Vg = V-Vf
	hf = rp_hPx(P0*1e5,0)
	hg = rp_hPx(P0*1e5,1)
	rhof = rp_rhoPx(P0*1e5,0)
	rhog = rp_rhoPx(P0*1e5,1)
	Mf = rhof*Vf
	Mg = rhog*Vg

	Pcrit = 73.77
	rhocrit = CRYO_PF_CritProp(f_out.fluid,fprop_density,ier)
	f_out.is_C = TRUE
CONTINUOUS
	lvl = Vf/V // fractional
	level = lvl*L // in meters
	P := P_Pa * 1e-5
	hfs := rp_hPx(P_Pa,0)
	hgs := rp_hPx(P_Pa,1)
<st>	rp_all(P_Pa,hf,hg,rhof,rhog,d1,d2,drhof_dP,drhof_dh,drhog_dP,drhog_dh)	
	Tf = CRYO_PF_prop_vs_ph(f_out.fluid,P,hf,fprop_temperature,ier,jx,jy)
	Tg = CRYO_PF_prop_vs_ph(f_out.fluid,P,hg,fprop_temperature,ier,jx,jy)
	Tsat = CRYO_PF_prop_vs_Px(f_out.fluid,P,0,fprop_temperature,ier,jx,jy)
	
	x := Mg/M
	h := x*hg + (1-x)*hf
	M = Mf + Mg

	m_evap := ZONE(hf<hfs+delta) 0 OTHERS tau_f*rhof*Vf*(hf-hfs)/(hgs-hfs) // Keep a minimal distance from the edge
	m_cond := ZONE(hg>hgs-delta) 0 OTHERS tau_g*rhog*Vg*(hgs-hg)/(hgs-hfs)
	
	Qh := Heat.signal[1]/100 * Qheat_max
	Qc := Cool.signal[1]/100 * Qcool_max
	
	mh_out = semiLinear(f_out.m,0,f_out.hb-hf)
	
	Qfg = alpha_fg * Ac * (Tf-Tg) // +ive => liquid hotter => heat transfer from liquid to vapour
	
<maf>	Vf*(drhof_dP*P_Pa' + drhof_dh*hf') + rhof*Vf' = m_cond - m_evap - f_out.m
<mag>	Vg*(drhog_dP*P_Pa' + drhog_dh*hg') - rhog*Vf' = m_evap - m_cond
<enf>	Vf*(rhof*hf' - P_Pa') = m_cond*(hfs-hf) - m_evap*(hgs-hf) - mh_out + Qh - Qfg
<eng>	Vg*(rhog*hg' - P_Pa') = m_evap*(hgs-hg) - m_cond*(hfs-hg) - Qc + Qfg
	
	Vf = Mf/rhof
	Vg = Mg/rhog
	V = Vf + Vg
	
	Pressure.signal[1] = P
	Level.signal[1] = lvl*100
	f_out.P = P
	f_out.P_aux = P
	f_out.hf = spliceFunction(hf,hg,level-(z_out+0.5*D_out),D_out/2.0)
	f_out.rho = rhof
	f_out.I = 1
	f_out.A = MATH.PI * (0.05/2)**2
	f_out.Gcrit = 1
	f_out.v = f_out.m/(rhof*f_out.A)	
END COMPONENT



COMPONENT test_AccumulatorTFM_spray IS_A test_AccumulatorTFM_pini
PORTS
	IN fluid f_in
DATA
	REAL z_in = 1.8 UNITS u_m "Elevation of spray port"
	REAL D_in = 0.01 UNITS u_m "Diameter of spray port"
CONTINUOUS
<:maf>	Vf*(drhof_dP*P_Pa' + drhof_dh*hf') + rhof*Vf' = m_cond - m_evap + f_in.m - f_out.m
<:mag>	Vg*(drhog_dP*P_Pa' + drhog_dh*hg') - rhog*Vf' = m_evap - m_cond
<:enf>	rhof*Vf*hf' - Vf*P_Pa' = m_cond*(hfs-hf) - m_evap*(hgs-hf) - mh_out + Qh - Qfg + f_in.m*(hfs-hf)
<:eng>	rhog*Vg*hg' - Vg*P_Pa' = m_evap*(hgs-hg) - m_cond*(hfs-hg) - Qc + Qfg - f_in.m*(hfs-f_in.hf)
	
	f_in.fluid = f_out.fluid
	f_in.n_fluid = f_out.n_fluid
	f_in.P = P
	f_in.P_aux = P
	f_in.hb = hg
	f_in.rho = rhog
	f_in.v = 1
	f_in.A = 1
	f_in.Gcrit = 1
END COMPONENT


COMPONENT test_AccumulatorTFM_flowthrough IS_A test_AccumulatorTFM_pini
PORTS
	IN fluid f_in
	OUT analog_signal(n=1) SaturationTemp
DATA
	REAL z_in = 1.8 UNITS u_m "Elevation of spray port"
	REAL D_in = 0.01 UNITS u_m "Diameter of spray port"
CONTINUOUS
<:maf>	Vf*(drhof_dP*P_Pa' + drhof_dh*hf') + rhof*Vf' = m_cond - m_evap + f_in.m - f_out.m
<:mag>	Vg*(drhog_dP*P_Pa' + drhog_dh*hg') - rhog*Vf' = m_evap - m_cond
<:enf>	rhof*Vf*hf' - Vf*P_Pa' = m_cond*(hfs-hf) - m_evap*(hgs-hf) - mh_out + Qh - Qfg + f_in.m*(f_in.hf-hf)
<:eng>	rhog*Vg*hg' - Vg*P_Pa' = m_evap*(hgs-hg) - m_cond*(hfs-hg) - Qc + Qfg - f_in.m*(hfs-f_in.hf)
	
	SaturationTemp.signal[1] = CRYO_PF_prop_vs_Px(f_out.fluid,P,0,fprop_temperature,ier,jx,jy) - 273.15
	
	f_in.fluid = f_out.fluid
	f_in.n_fluid = f_out.n_fluid
	f_in.P = P
	f_in.P_aux = P
	f_in.hb = hg
	f_in.rho = rhog
	f_in.v = 1
	f_in.A = 1
	f_in.Gcrit = 1
END COMPONENT



COMPONENT test_AccumulatorTFM_mix IS_A test_AccumulatorTFM_spray
"Spray cooling model where the spray mixes completely with the vapour phase"
// unlike previous model, here the spray is assumed to evaporate completely while it mixes with the vapour phase
// therefore, the mass balance accounting for the spray is done in the vapour phase.
DECLS
	REAL mh_in UNITS u_W "Inlet port enthalpy flow rate"
CONTINUOUS
	mh_in = semiLinear(f_in.m,f_in.m*(f_in.hf-hg),0)
	
<:maf>	Vf*(drhof_dP*P_Pa' + drhof_dh*hf') + rhof*Vf' = m_cond - m_evap - f_out.m
<:mag>	Vg*(drhog_dP*P_Pa' + drhog_dh*hg') - rhog*Vf' = m_evap - m_cond + f_in.m
<:enf>	rhof*Vf*hf' - Vf*P_Pa' = m_cond*(hfs-hf) - m_evap*(hgs-hf) + Qh - Qfg - mh_out
<:eng>	rhog*Vg*hg' - Vg*P_Pa' = m_evap*(hgs-hg) - m_cond*(hfs-hg) - Qc + Qfg + mh_in
END COMPONENT		

*/

