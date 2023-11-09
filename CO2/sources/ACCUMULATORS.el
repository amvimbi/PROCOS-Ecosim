/*------------------------------------------------------------------------------
 LIBRARY: CO2
 FILE: ACCUMULATOR
 AUTHOR: Viren Bhanot
 COMPANY: CERN
 DESCRIPTION: Accumulator Models
 CREATION DATE: 18/10/2016
------------------------------------------------------------------------------*/
USE CRYOLIB
USE MATH
USE PORTS_LIB
USE THERMAL


/*
 * NOTE 1: INHERITANCE/COMPOSITION LIMITATIONS
 * Complication involved in trying to use composition instead of inheritance.
 * In EcosimPro, it is not possible to do the following:
 *
 *	COMPONENT A (INTEGER n = 1)
 *	DATA
 * 		REAL x[n] = 2
 *	END COMPONENT
 *	
 *	COMPONENT B IS_A A
 * 	DECLS
 *		CLOSE n = 2
 *	END COMPONENT
 *
 *	COMPONENT C
 *	TOPOLOGY
 *		B b (x[1] = 5) // doesn't work
 *	END COMPONENT
 *
 * In practice, this means that to assign values to D_in[nports_in] etc of the 
 * parent component (AbsTank), I have to manually declare D_in[1] in the child
 * component (Accumulator_Receiver). Even then, unless the "[1]" matches the
 * "nports_in", I will get a compile error.
 *
 * NOTE 2: INIT BLOCK LIMITATION
 * Cannot set f_out[i].is_C in INIT block
 * Because one cannot have a variable index (i) in the INIT block.
 * Thus, everytime the nports_out number is changed, we should 
 * manually assign the is_C of the relevant port. (The nports_in number will be 
 * handled as usual by downstream component)
 */



--------------------------------------------------------------------------------
// BASE COMPONENTS
--------------------------------------------------------------------------------
ABSTRACT COMPONENT Tank (
	INTEGER nports_out = 1 "Number of outlet ports",
	BOOLEAN AccountforGravity = FALSE "If TRUE, adds static pressure at fluid ports",
	BOOLEAN ConstantRefHTC = TRUE "If FALSE, uses Cooper pool-boiling correlation",
	BOOLEAN DynamicRefHTC = FALSE "If TRUE, adds low-pass fliter to HTC calculations"
)
"Base refrigerant-side component for phase separators with only outlet ports"
// No geometry assignment at this stage
// No heat transfer calculations or thermal wall
// No inlet ports
PORTS
	OUT fluid f_out[nports_out]
	OUT thermal(n=1) tp_out
DATA
	REAL D = 0.1 UNITS u_m "Diameter"
	REAL L = 1 UNITS u_m "Height"
	REAL z_bottom = 0 UNITS u_m "Elevation of bottom of tank wrt user-defined coordinate system"
	REAL D_out[nports_out] = 0.001 UNITS u_m "Outlet port diameters"
	REAL z_out[nports_out] = 0.02 UNITS u_m "Outlet port heights"
	REAL alpha_ref = 10000 UNITS u_W_m2K "Refrigerant-to-shell HTC"
	REAL tauAlpha = 3 UNITS u_s "Time constant for HTC correlation low-pass filter"
	ENUM InitialConditionsAccu init = P_h "Specify initial conditions"
	REAL P0 = 20 UNITS u_bar "Initial refrigerant pressure"
	REAL h0 = 250000 UNITS u_J_kg "Initial refrigerant enthalpy"
	REAL T0 = 255 UNITS u_K "Initial refrigerant temperature"
	REAL x0 = 0.1 "Initial vapour quality"
	REAL Level0 = 20 UNITS u_pct "Initial liquid level (0-100%)"
DECLS
	REAL alpha_o UNITS u_W_m2K "Low-pass filtered outer wall refrigerant-side HTC"
	REAL alphaHat_o UNITS u_W_m2K "Outer wall refrigerant HTC from correlation"
	REAL cp UNITS u_J_kgK "Specific heat capacity"
	REAL delT UNITS u_K "Superheating/Subcooling (positive is superheat)"
	REAL dP_out[nports_out] UNITS u_Pa "Outlet pressure change due to gravitational effects"
	REAL drho_dh "Partial density derivative wrt enthalpy"
	REAL drho_dP "Partial density derivative wrt pressure"
	REAL gamma "Void fraction"
	REAL h UNITS u_J_kg "Specific enthalpy"
	REAL level UNITS u_m "Height of liquid level above tank bottom"
	REAL level_frac "Fraction of vessel filled"
	REAL lvl UNITS u_m "Parameter used to ensure numerical robustness in level calculations"
	REAL k UNITS u_W_mK "Thermal conductivity"
	REAL M UNITS u_kg "Refrigerant mass"
	REAL m_total_out UNITS u_kg_s "Total outlet mass flow rate"
	REAL mh_total_out UNITS u_W "Total outlet enthalpy flow rate"
	REAL mu UNITS u_Pas "Dynamic viscosity"
	DISCR REAL MW "Molecular weight"
	REAL P UNITS u_bar "Pressure"
	REAL P_Pa UNITS u_Pa "Pressure in Pascal"
	DISCR REAL Pcrit "Critical pressure"
	ENUM Phase phase "Fluid phase"
	REAL Pred "Reduced pressure"
	REAL Pr "Prandtl number"
	REAL Q_o UNITS u_W "Heat transfer at the outer surface"
	REAL q_o UNITS u_W_m2 "Heat flux at outer surface"
	REAL rho UNITS u_kg_m3 "Density"
	REAL sigma UNITS u_N_m "Surface tension"
	REAL T UNITS u_K "Temperature"
	REAL Tsat UNITS u_K "Saturation temperature"
	REAL u UNITS u_J_kg "Internal energy"
	REAL Vol UNITS u_dm3 "Internal volume (L)"
	REAL vsound UNITS u_m_s "Mach number"
	REAL x "Vapour quality"
	REAL z_level UNITS u_m "Elevation of the liquid level wrt user-defined coordinate system"
	REAL z_top UNITS u_m "Elevation of the top of the tank wrt user-defined coordinate system"
	INTEGER ier,jx,jy // error codes
OBJECTS
	RefGeometryRecord geo
	SaturationProperties sat
INIT
	alpha_o = alpha_ref
	alphaHat_o = alpha_ref
	assignInitAccu(init,f_out[1].fluid,P0,h0,T0,x0,Level0,P,h,rho,u)
	FOR (i IN 1,nports_out)
		dP_out[i] = 0
	END FOR
	f_out[1].is_C = TRUE
	geo.EndSurfaces = 2 // closed at both ends
	MW = getMolWt(f_out[1].fluid)
	Pcrit = 73.77 // CRYO_PF_CritProp(f_out[1].fluid,fprop_pressure,ier)
	P_Pa = P0*1E5
DISCRETE
	EXPAND (i IN 1,nports_out) ASSERT (z_out[i] <= L - D_out[i]/2) FATAL "Outlet pipe height is greater than vessel height"
CONTINUOUS
<Vol>	Vol = geo.V*1000
<M> 	M = rho*geo.V // modified later when cooling spiral takes up some of the volume
	//Level
	div_safe((rho-sat.rho_g)*L, sat.rho_l-sat.rho_g, lvl) // calculate 'lvl' for numerical robustness. See div_safe in MATHS library.
	delT = T - Tsat
	level = ZONE (P>=Pcrit) L OTHERS min(L,max(0,lvl))
	z_top = z_bottom + L
	z_level = z_bottom + level
	level_frac = level/L
	// Pressures
	P_Pa = P * 1e5
	Pred = P/(Pcrit*1e-5)
	// Heat Transfer
	IF (ConstantRefHTC) INSERT
		alphaHat_o = alpha_ref
	ELSE
		alphaHat_o = Cooper(abs(T-tp_out.Tk[1]),MW,Pred,q_o)
	END IF
	IF (DynamicRefHTC) INSERT
		alpha_o' = 1/tauAlpha * (alphaHat_o - alpha_o)
	ELSE
		alpha_o = alphaHat_o
	END IF
	Q_o = alpha_o*geo.As_o*(T - tp_out.Tk[1])
	tp_out.q[1] = Q_o
	q_o = abs(alpha_o*(T-tp_out.Tk[1]))
	// Enthalpy flow rates
<mo>	m_total_out = SUM(i IN 1,nports_out; f_out[i].m)
<mho>	mh_total_out = SUM(i IN 1,nports_out;f_out[i].mh)
	// Governing equations
<mas>	geo.V*(drho_dP*P_Pa'+drho_dh*h') = - m_total_out
<ene>	geo.V*((h*drho_dP-1)*P_Pa' + (h*drho_dh+rho)*h') = - mh_total_out - Q_o
	// Thermodynamic state
<st>	getStatePh(f_out[1].fluid,P,h,cp,drho_dP,drho_dh,gamma,k,mu,phase,Pr,rho,sigma,T,Tsat,u,vsound,x,sat,ier,jx,jy)	
	// Port assignments
<hf>	EXPAND (i IN 1,nports_out) f_out[i].hf = ZONE (phase==two_phase) spliceFunction(sat.h_l,sat.h_g,level/(z_out[i]+0.5*D_out[i])-1,D_out[i]/(2.0*(z_out[i]+D_out[i]/2.0))) OTHERS h	
	IF (nports_out>0) INSERT
		EXPAND_BLOCK (i IN 1,nports_out)
			dP_out[i] = (sat.rho_g*(z_top - max(z_out[i],z_level)) + sat.rho_l*(max(z_out[i],z_level) - z_out[i]))*9.806
			f_out[i].A = MATH.PI*(D_out[i]/2)**2
			f_out[i].Gcrit = CRYO_Gcrit_fun(f_out[i].fluid,P*1E5,rho,T,0,f_out[i].v,vsound,ier,jy)
			f_out[i].rho = CRYO_PF_prop_vs_ph(f_out[1].fluid, f_out[i].P,f_out[i].h,fprop_density,ier,jx,jy)
			f_out[i].P = ZONE(AccountforGravity) P + dP_out[i]/1e5 OTHERS P
			f_out[i].P_aux = P + geo.V*rho*vsound*1E-5/geo.As
			f_out[i].I = 0.5*geo.V/geo.As**2
			f_out[i].v = f_out[i].m/(rho*geo.Ap)
		END EXPAND_BLOCK
	END IF
END COMPONENT



COMPONENT Tank_addInlet IS_A Tank (
	INTEGER nports_in = 1
)
"Add inlet ports to Tank model"
// still no geometry assigned
PORTS
	IN fluid f_in[nports_in]
DATA
	REAL D_in[nports_in] = 0.001 UNITS u_m "Inlet port diameters"
	REAL z_in[nports_in] = 0.02 UNITS u_m "Inlet port heights"
DECLS
	REAL dP_in[nports_in] UNITS u_bar "Inlet pressure change due to gravity"
	REAL m_total_in UNITS u_kg_s "Total inlet mass flow rate"
	REAL mh_total_in UNITS u_W "Total inlet enthalpy flow rate"
DISCRETE
	EXPAND_BLOCK(i IN 1,nports_in)
		ASSERT (f_in[i].is_C == FALSE) FATAL "The inlet port(s) must be connected to Resistive components"
		ASSERT (z_in[i] <= L - D_in[i]/2) FATAL "Inlet pipe height cannot be greater than vessel height"
	END EXPAND_BLOCK
CONTINUOUS
	m_total_in = SUM(i IN 1,nports_in; f_in[i].m)
	mh_total_in = SUM(i IN 1,nports_in; f_in[i].mh)
<:mas>	geo.V*(drho_dP*P_Pa'+ drho_dh*h') = m_total_in - m_total_out
<:ene>	geo.V*((h*drho_dP-1)*P_Pa' + (h*drho_dh+rho)*h') = mh_total_in - mh_total_out - Q_o

	EXPAND (i IN 2,nports_out) f_out[i].fluid = f_out[1].fluid
	EXPAND (i IN 2,nports_out) f_out[i].n_fluid = f_out[1].n_fluid
<hb>	EXPAND (i IN 1,nports_in) f_in[i].hb = ZONE (phase==two_phase) spliceFunction(sat.h_l,sat.h_g,level/(z_in[i]+0.5*D_in[i])-1,D_in[i]/(2.0*(z_in[i]+D_in[i]/2.0))) OTHERS h
	EXPAND (i IN 1,nports_in) f_in[i].fluid = f_out[1].fluid
	EXPAND (i IN 1,nports_in) f_in[i].n_fluid = f_out[1].n_fluid
	IF (nports_in>0) INSERT
		EXPAND_BLOCK (i IN 1,nports_in)
			dP_in[i] = (sat.rho_g*(z_top - max(z_in[i],z_level)) + sat.rho_l*(max(z_in[i],z_level) - z_in[i]))*9.806
			f_in[i].A = MATH.PI*(D_in[i]/2)**2
			f_in[i].Gcrit = CRYO_Gcrit_fun(f_in[i].fluid,P*1E5,rho,T,0,f_in[i].v,vsound,ier,jy)
			f_in[i].rho = CRYO_PF_prop_vs_ph(f_in[i].fluid, f_in[i].P,f_in[i].h,fprop_density,ier,jx,jy)
			f_in[i].P = ZONE(AccountforGravity) P + dP_in[i] OTHERS P
			f_in[i].P_aux = P + geo.V*rho*vsound*1E-5/geo.As
			f_in[i].I = 0.5*geo.V/geo.As**2
			f_in[i].v = f_in[i].m/(rho*geo.Ap)
		END EXPAND_BLOCK
	END IF
END COMPONENT






--------------------------------------------------------------------------------
// HEAT PUMP COMPONENTS
--------------------------------------------------------------------------------
COMPONENT AccumulatorBase IS_A Tank_addInlet
"Refrigerant side base component for j-type accumulators"
DECLS
	CLOSE nports_in = 1
	CLOSE nports_out = 1
INIT
	geo.setCylindricalGeometry(D,L)
END COMPONENT



COMPONENT FlashTankBase IS_A Tank_addInlet
"Refrigerant-side base component for Flash tanks"
DECLS
	CLOSE nports_in = 1
	CLOSE nports_out = 2
INIT
	geo.setCylindricalGeometry(D,L)
END COMPONENT



COMPONENT SuctionChamberVolume IS_A Tank_addInlet
"Adds a thermal port on the inside to account for Motor heat transfer." 
// Also assigns annular geometry to tank
PORTS
	IN thermal(n=1) tp_in
DATA
	REAL Di = 0.05 UNITS u_m2 "Inner diameter"
DECLS
	CLOSE nports_in = 1
	CLOSE nports_out = 1
	CLOSE AccountforGravity = FALSE
	
	REAL alphaHat_i UNITS u_W_m2K
	REAL alpha_i UNITS u_W_m2K
	REAL q_i UNITS u_W_m2
	REAL Q_i UNITS u_W
	REAL x_in
	REAL d1, d2 // dummy variables
INIT
	geo.setAnnularGeometry(D,0.5*(D-Di),L)
	alphaHat_i = alpha_ref
	alpha_i = alpha_ref
CONTINUOUS
	x_in = CRYO_PF_prop_vs_ph(f_in[1].fluid,f_in[1].P,f_in[1].h,fprop_quality,ier,jx,jy)
	IF (ConstantRefHTC) INSERT
		alphaHat_i = alpha_ref
	ELSE
		alphaHat_i = Cooper(abs(T-tp_in.Tk[1]),MW,Pred,q_i)
	END IF
	IF (DynamicRefHTC) INSERT
		alpha_i' = 1/tauAlpha * (alphaHat_i - alpha_i)
	ELSE
		alpha_i = alphaHat_i
	END IF
	Q_i = alpha_i*geo.As_i*(tp_in.Tk[1] - T)
	tp_in.q[1] = Q_i
	q_i = abs(alpha_i*(tp_in.Tk[1]-T))
	
<:ene>	geo.V*((h*drho_dP-1)*P_Pa' + (h*drho_dh+rho)*h') = mh_total_in - mh_total_out + Q_i - Q_o
	d1 := geo.V
	d2 := mh_total_in - mh_total_out + Q_i - Q_o
END COMPONENT



COMPONENT Accumulator_Receiver (
	BOOLEAN ConstantRefHTC = TRUE "If FALSE, uses Cooper pool boiling correlation",
	BOOLEAN DynamicRefHTC = FALSE "If TRUE, applies low-pass filter to HTC calculations",
	BOOLEAN ConstantAirHTC = TRUE "If FALSE, uses Churchill-Chu natural convective correlation"
)
"Suction-line accumulator accumulator. Can be used as liquid-line receiver also."
PORTS
	IN fluid f_in
	OUT fluid f_out
	IN analog_signal(n=1) T_amb "Ambient temperature"
DATA
	ENUM Material mat = SS_304 "Shell material"
	REAL D = 0.1329 UNITS u_m "Shell outre diameter"
	REAL Th = 0.005 UNITS u_m "Shell thickness"
	REAL L = 0.3429  UNITS u_m "Shell Height"
	REAL D_in[1] = 0.01 UNITS u_m "Inlet pipe diameter" //see note at the top of this file
	REAL z_in[1] = 0.28 UNITS u_m "Inlet pipe height"
	REAL D_out[1] = 0.01 UNITS u_m "Outlet pipe diameter"
	REAL z_out[1] = 0.28 UNITS u_m "Outlet pipe height"

	REAL alpha_a = 25 UNITS u_W_m2K "Shell to air heat transfer coefficient"
	REAL alpha_ref = 200 UNITS u_W_m2K "Refrigerant to shell heat transfer coefficient"
	REAL tauAlpha = 3 UNITS u_s "Time constant for HTC correlation low-pass filter"

	ENUM InitialConditionsAccu init = P_h "Initial condition specification"
	REAL P0 = 20 UNITS u_bar "Initial pressure"
	REAL h0 = 250000 UNITS u_J_kg "Initial specific enthalpy"
	REAL T0 = 255 UNITS u_K "Initial temperature"
	REAL x0 = 0.5 "Initial vapour quality"
	REAL Level0 = 20 UNITS u_pct "Initial liquid level in %"
	REAL Tw0 = 300 UNITS u_K "Shell initial temperature"
TOPOLOGY
	AccumulatorBase (
		ConstantRefHTC = ConstantRefHTC,
		DynamicRefHTC = DynamicRefHTC
	) Accumulator (
		D = D-2*Th,
		L = L,
		D_in = D_in,
		z_in = z_in,
		D_out = D_out,
		z_out = z_out,
		alpha_ref = alpha_ref,
		tauAlpha = tauAlpha,
		init = init,
		P0 = P0,
		h0 = h0,
		T0 = T0,
		x0 = x0,
		Level0 = Level0
	)
	LumpedShell(ConstantAirHTC=ConstantAirHTC) Wall (
		mat = mat,
		Do = D,
		Th = Th,
		L = L,
		Tw0 = Tw0,
		alpha_a0 = alpha_a
	)
	CONNECT T_amb TO Wall.T_amb
	CONNECT f_in TO Accumulator.f_in[1]
	CONNECT Accumulator.f_out[1] TO f_out
	CONNECT Accumulator.tp_out TO Wall.tp_in
END COMPONENT



COMPONENT FlashTank (
	BOOLEAN ConstantAirHTC = TRUE "If FALSE, uses ChurchillChu natural convection correlation",
	BOOLEAN ConstantRefHTC = TRUE "If FALSE, uses Cooper pool boiling correlation",
	BOOLEAN DynamicRefHTC = FALSE "If TRUE, applies low-pass filter to HTC calculations",
	// these need to be here rather than in DATA to avoid syntax errors:
	REAL D_vap = 0.01 UNITS u_m "Vapour pipe diameter",
	REAL H_vap = 0.28 UNITS u_m "Vapour pipe height",	
	REAL D_liq = 0.01 UNITS u_m "Liquid pipe diameter",
	REAL H_liq = 0.0001 UNITS u_m "Liquid pipe height"
)
"Flash tank component for phase separation"
PORTS
	IN fluid f_in
	OUT fluid f_vap
	OUT fluid f_liq
	IN analog_signal(n=1) T_amb
DATA
	ENUM Material mat = SS_304 "Shell material"
	REAL D = 0.1329 UNITS u_m "Shell outer diameter"
	REAL Th = 0.005 UNITS u_m "Shell thickness"
	REAL L = 0.3429  UNITS u_m "Length (or Height for vertical components)"
	REAL D_in[1] = 0.01 UNITS u_m "Inlet pipe diameter"
	REAL z_in[1] = 0.28 UNITS u_m "Inlet pipe height"
	REAL alpha_a = 25 UNITS u_W_m2K "Shell to air heat transfer coefficient"
	REAL alpha_ref = 100 UNITS u_W_m2K "Shell to refrigerant heat transfer coefficient"
	REAL tauAlpha = 3 UNITS u_s "Time constant for HTC correlation low-pass filter"
	ENUM InitialConditionsAccu init = P_h "Initial condition specification"
	REAL P0 = 20 UNITS u_bar "Initial pressure"
	REAL h0 = 250000 UNITS u_J_kg "Initial specific enthalpy"
	REAL T0 = 255 UNITS u_K "Initial temperature"
	REAL x0 = 0.5 "Initial vapour quality"
	REAL Level0 = 20 UNITS u_pct "Initial liquid level in %"
	REAL Tw0 = 300 UNITS u_K "Shell initial temperature"
TOPOLOGY
	FlashTankBase (
		ConstantRefHTC = ConstantRefHTC,
		DynamicRefHTC = DynamicRefHTC
	) FlashTank (
		D = D-2*Th,
		L = L,
		D_in = D_in,
		z_in = z_in,
		D_out = {D_vap,D_liq},
		z_out = {H_vap,H_liq},
		alpha_ref = alpha_ref,
		tauAlpha = tauAlpha,
		init = init,
		P0 = P0,
		h0 = h0,
		T0 = T0,
		x0 = x0,
		Level0 = Level0
	)
	LumpedShell (ConstantAirHTC=ConstantAirHTC) Wall (
		Do = D,
		Th = Th,
		L = L,
		alpha_a0 = alpha_a
	)
	CONNECT T_amb TO Wall.T_amb
	CONNECT f_in TO FlashTank.f_in[1]
	CONNECT FlashTank.f_out[1] TO f_vap
	CONNECT FlashTank.f_out[2] TO f_liq
	CONNECT FlashTank.tp_out TO Wall.tp_in
END COMPONENT



COMPONENT FlashTank_insulated (
	BOOLEAN ConstantRefHTC = TRUE "If FALSE, uses Cooper pool boiling correlation",
	BOOLEAN DynamicRefHTC = FALSE "If TRUE, uses a low-pass filter for HTC calculation",
	REAL D_vap = 0.01 UNITS u_m "Vapour pipe diameter",
	REAL H_vap = 0.28 UNITS u_m "Vapour pipe height",	
	REAL D_liq = 0.01 UNITS u_m "Liquid pipe diameter",
	REAL H_liq = 0.0001 UNITS u_m "Liquid pipe height"
)
"Flash tank component for phase separation"
PORTS
	IN fluid f_in
	OUT fluid f_vap
	OUT fluid f_liq
DATA
	ENUM Material mat = SS_304 "Shell material"
	REAL D = 0.1329 UNITS u_m "Shell outer diameter"
	REAL Th = 0.005 UNITS u_m "Shell thickness"
	REAL L = 0.3429  UNITS u_m "Length (or Height for vertical components)"
	REAL D_in[1] = 0.01 UNITS u_m "Inlet pipe diameter"
	REAL z_in[1] = 0.28 UNITS u_m "Inlet pipe height"
	REAL alpha_ref = 100 UNITS u_W_m2K "Shell to refrigerant heat transfer coefficient"
	REAL tauAlpha = 3 UNITS u_s "Time constant for low-pass HTC calculations"
	ENUM InitialConditionsAccu init = P_h "Initial condition specification"
	REAL P0 = 20 UNITS u_bar "Initial pressure"
	REAL h0 = 250000 UNITS u_J_kg "Initial specific enthalpy"
	REAL T0 = 255 UNITS u_K "Initial temperature"
	REAL x0 = 0.5 "Initial vapour quality"
	REAL Level0 = 20 UNITS u_pct "Initial liquid level in %"
	REAL Tw0 = 300 UNITS u_K "Shell initial temperature"
TOPOLOGY
	FlashTankBase (ConstantRefHTC = ConstantRefHTC, DynamicRefHTC = DynamicRefHTC) FlashTank (
		D = D-2*Th,
		L = L,
		D_in = D_in,
		z_in = z_in,
		D_out = {D_vap,D_liq},
		z_out = {H_vap,H_liq},
		alpha_ref = alpha_ref,
		tauAlpha = tauAlpha,
		init = init,
		P0 = P0,
		h0 = h0,
		T0 = T0,
		x0 = x0,
		Level0 = Level0
	)
	LumpedShell_insulated Wall (
		Do = D,
		Th = Th,
		L = L
	)
	CONNECT f_in TO FlashTank.f_in[1]
	CONNECT FlashTank.f_out[1] TO f_vap
	CONNECT FlashTank.f_out[2] TO f_liq
	CONNECT FlashTank.tp_out TO Wall.tp_in
END COMPONENT






--------------------------------------------------------------------------------
// CO2 COMPONENTS
--------------------------------------------------------------------------------
/* BASE COMPONENTS
 * Components with 'Base' in their name make up the fluid control volume
 * of the accumulator models. These are later used in creating composite
 * accumulator models with walls, heaters, chiller coils etc.
 */
 
ABSTRACT COMPONENT CO2AccumulatorBase IS_A Tank
"Refrigerant side base component for CO2 accumulators"
// Only specifies geometry and ports
// Does not specify heater and cooling coil components
DECLS
	CLOSE nports_out = 1
INIT
	geo.setCylindricalGeometry(D,L)
END COMPONENT



COMPONENT CO2AccumulatorBase_prescribed IS_A CO2AccumulatorBase
"Prescribed cooling and heating input values. Gravity neglected and two-phase fluid assumed homogeneous."
PORTS
	IN analog_signal(n=1) K_heat
	IN analog_signal(n=1) K_cool
DATA
	REAL Q_cool_max = 1000 UNITS u_W "Maximum cooling capacity"
	REAL Q_heat_max = 1000 UNITS u_W "Maximum heating capacity"
DECLS
	REAL Q_coolPow UNITS u_W
	REAL Q_heatPow UNITS u_W
	REAL hliq, hvap
CONTINUOUS
	Q_heatPow = K_heat.signal[1]/100*Q_heat_max
	Q_coolPow = K_cool.signal[1]/100*Q_cool_max
<:ene>	geo.V*((h*drho_dP-1)*P_Pa' + (h*drho_dh+rho)*h') = Q_heatPow - Q_coolPow - mh_total_out - Q_o
	hliq = sat.h_l
	hvap = sat.h_g
END COMPONENT



COMPONENT CO2AccumulatorBase_heater IS_A CO2AccumulatorBase
"Adds cartridge heater, maintains prescribed cooling"
// by prescribed cooling, I mean max cooling power is user-specified
PORTS
	IN thermal(n=1) tp_in_eh
	IN analog_signal(n=1) K_cool
DATA
	REAL D_heater = 0.01 UNITS u_m "Heater diameter"
	REAL L_heater = 0.1 UNITS u_m "Heater length"
	REAL alpha_heater_l = 1000 UNITS u_W_m2K "Liquid Heat transfer coefficient from heater to refrigerant"
	REAL alpha_heater_tp = 1000 UNITS u_W_m2K "Two-phase HTC heater to refrigerant" // TODO: should be splined straight from liquid to vapour. No such thing as two-phase
	REAL alpha_heater_v = 100 UNITS u_W_m2K "Vapour HTC heater to refrigerant"
	REAL Q_cool_max = 1000 UNITS u_W "Maximum cooling capacity"
DECLS
	REAL Q_eh UNITS u_W "Heat transfer from heater"
	REAL Q_coolPow UNITS u_W "Cooling capacity"
	DISCR REAL As_heater UNITS u_m2 "Heater outer surface area"
	REAL alpha_heater UNITS u_W_m2K "Actual HTC heater to refrigerant"
INIT
	As_heater = 2*MATH.PI*(D_heater/2)*L_heater + MATH.PI*(D_heater/2)**2
CONTINUOUS
	alpha_heater = getSmoothHTC(alpha_heater_v,alpha_heater_l,alpha_heater_tp,0,0,x,0.1,0.9,FALSE)
	Q_eh = alpha_heater * As_heater * (tp_in_eh.Tk[1] - T)
	tp_in_eh.q[1] = Q_eh
	Q_coolPow = K_cool.signal[1]/100*Q_cool_max
<:ene>	geo.V*((h*drho_dP-1)*P_Pa' + (h*drho_dh+rho)*h') = - mh_total_out - Q_o + Q_eh - Q_coolPow
END COMPONENT



COMPONENT CO2AccumulatorBase_externChiller IS_A AccumulatorBase
"CO2 Accumulator with connection for cold liquid injection at the top"
// Adds a cartridge heater to the Accumulator_Receiver component
PORTS
	IN thermal(n=1) tp_in_eh "Cartridge electrical heater thermal port"
DATA
	REAL D_heater = 0.1 UNITS u_m "Heater diameter"
	REAL L_heater = 1 UNITS u_m "Heater length"
	REAL alpha_heater = 1000 UNITS u_W_m2K "Heat transfer coefficient with cartridge heater"
DECLS
	REAL Q_eh UNITS u_W "Heat flow from electrical heater"
	DISCR REAL As_heater UNITS u_m2 "Heater surface area"
INIT
	As_heater = 2 * MATH.PI * (D_heater/2) * L_heater
CONTINUOUS
	Q_eh = alpha_heater * As_heater * (tp_in_eh.Tk[1] - T) // Positive if heat flowing in from heater
	tp_in_eh.q[1] = Q_eh
<:ene> 	geo.V*((h*drho_dP-1)*P_Pa' + (h*drho_dh+rho)*h') = mh_total_in - mh_total_out + Q_eh
END COMPONENT



COMPONENT CO2AccumulatorBase_flowThrough IS_A Tank_addInlet
"Base model for flow-through accumulator"
// adds an inlet port that supplies return fluid from the transfer line outlet
PORTS
	IN thermal(n=1) tp_in_eh "Cartridge heater thermal port"
	IN thermal(n=1) tp_in_ch "Chiller coil thermal port"
DATA
	REAL D_heater = 0.1 UNITS u_m "Heater diameter"
	REAL L_heater = 1 UNITS u_m "Heater length"
	REAL D_coil = 0.1 UNITS u_m "Chiller coil inner diameter"
	REAL L_coil = 1 UNITS u_m "Chiller coil length"
	
	REAL alpha_heater_v = 1000 UNITS u_W_m2K "Heater-to-refrigerant vapour HTC"
	REAL alpha_heater_tp = 10000 UNITS u_W_m2K "Heater-to-refrigerant two-phase HTC"
	REAL alpha_heater_l = 3500 UNITS u_W_m2K "Heater-to-refrigerant liquid HTC"
	
	REAL alpha_coil = 1000 UNITS u_W_m2K "Cooling-spiral-to-refrigerant heat transfer coefficient"
DECLS
	CLOSE nports_out = 1
	CLOSE nports_in = 1
	REAL alpha_heater UNITS u_W_m2K "Heater-to-refrigerant HTC"
	REAL Q_ch UNITS u_W "Chiller-to-rerigerant heat transfer"
	REAL Q_eh UNITS u_W "Heater-to-refrigerant heat transfer"
	DISCR REAL As_heater UNITS u_m2 "Heater surface area"
	DISCR REAL As_coil UNITS u_m2 "Chiller coil surface area"
	DISCR REAL V_coil UNITS u_m3 "Chiller coil volume"
INIT
	geo.setAccumulatorGeometry(D,L)
	As_heater = 2 * MATH.PI * (D_heater/2) * L_heater
	As_coil = 2 * MATH.PI * (D_coil/2) * L_coil
	V_coil = MATH.PI * (D_coil/2)**2 * L_coil
CONTINUOUS
	Q_ch = alpha_coil * As_coil * (tp_in_ch.Tk[1] - T)
	tp_in_ch.q[1] = Q_ch
	alpha_heater = getSmoothHTC(alpha_heater_v,alpha_heater_l,alpha_heater_tp,0,0,x,0.1,0.9,FALSE)
	Q_eh = alpha_heater * As_heater * (tp_in_eh.Tk[1] - T)
	tp_in_eh.q[1] = Q_eh
<:ene>	geo.V * ((h*drho_dP-1)*P_Pa' + (h*drho_dh+rho)*h') = mh_total_in - mh_total_out + Q_eh + Q_ch - Q_o
END COMPONENT



COMPONENT CO2AccumulatorRefrigerantVolume (
	BOOLEAN AccountForGravity = FALSE "If TRUE, adds static pressure at fluid ports"
)
"Complete CO2Accumulator model with thermal ports for chiller coil and heater heat transfer"
// Unlike CO2AccumluatorBase_externChiller, this model does not include calculations for the 
// spiral within it. Instead, the spiral is assumed to be handled elsewhere, while the model
// just deals with the port interface exchanging heat transfer.
PORTS
	OUT fluid f_out
	IN thermal(n=1) tp_in_ch "Chiller thermal port"
	IN thermal(n=1) tp_in_eh "Cartridge electrical heater thermal port"
	OUT thermal(n=1) tp_out "Heat transfer with thermal shell"
DATA
	REAL D = 0.1 UNITS u_m "Diameter"
	REAL L = 1 UNITS u_m "Height"
	REAL z_bottom = 0 UNITS u_m "Elevation of bottom of tank wrt user-defined coordinate system" 
	REAL D_out = 0.01 UNITS u_m "Outlet port diameter"
	REAL z_out = 0.01 UNITS u_m "Outlet port height"
	REAL D_coil = 0.009525 UNITS u_m "Chiller coil outer diameter"
	REAL L_coil = 11 UNITS u_m "Chiller coil length"
	REAL D_heater = 0.1 UNITS u_m "Heater outer diameter"
	REAL L_heater = 1 UNITS u_m "Heater length"

	REAL alpha_ref = 200 UNITS u_W_m2K "Shell-to-refrigerant heat transfer coefficient"
	REAL alpha_coil = 1000 UNITS u_W_m2K "Chiller-coil-to-refrigerant heat transfer coefficient"
	REAL alpha_heater = 1000 UNITS u_W_m2K "Heater-to-refrigerant heat transfer coefficient"

	ENUM InitialConditionsAccu init = P_h "Specify initial conditions"
	REAL P0 = 20 UNITS u_bar "Initial refrigerant pressure"
	REAL h0 = 250000 UNITS u_J_kg "Initial refrigerant enthalpy"
	REAL T0 = 255 UNITS u_K "Initial refrigerant temperature"
	REAL x0 = 0.1 "Initial vapour quality"
	REAL Level0 = 20 UNITS u_pct "Initial liquid level in %"
DECLS
	DISCR REAL As_coil UNITS u_m2 "Outer surface area of coil"
	DISCR REAL As_heater UNITS u_m2 "Heater outer surface area"
	PRIVATE REAL cp UNITS u_J_kgK
	REAL dP_out UNITS u_Pa "Outlet pressure change due to gravitational effects"
	REAL drho_dh "Partial density derivative wrt enthalpy"
	REAL drho_dP "Partial density derivative wrt pressure"
	REAL gamma "Void fraction"
	REAL h UNITS u_J_kg "Specific enthalpy"
	PRIVATE REAL k UNITS u_J_kgK "Thermal conductivity"
	REAL lvl "Dummy level variable for safe division"
	REAL level UNITS u_m "Height of liquid level above tank bottom"
	REAL level_frac "Liquid level as fraction of 1"
	REAL M UNITS u_kg "Refrigerant charge"
	PRIVATE REAL mu UNITS u_Pas
	REAL P UNITS u_bar "Pressure"
	DISCR REAL Pcrit UNITS u_bar
	REAL P_Pa UNITS u_Pa "Pressure in Pascal"
	REAL Q_o UNITS u_W "Heat transfer to shell"
	REAL Q_ch UNITS u_W "Heat flow from chiller"
	REAL Q_eh UNITS u_W "Heat flow from electrical heater"
	REAL rho UNITS u_kg_m3 "Density"
	REAL sigma UNITS u_N_m "Surface tension"
	REAL T UNITS u_K "Temperature"
	REAL Tsat UNITS u_K "Saturation temperature"
	REAL delT UNITS u_K "Superheating/Subcooling (positive is superheat)"
	REAL u UNITS u_J_kg "Internal energy"
	DISCR REAL V_coil UNITS u_m3 "Volume taken up by cooling spiral"
	REAL vsound UNITS u_m_s "Mach number"
	REAL x "Vapour quality"
	REAL z_level UNITS u_m "Elevation of the liquid level wrt user-defined coordinate system"
	REAL z_top UNITS u_m "Elevation of the top of the tank wrt user-defined coordinate system"
	REAL Vol "Internal volume"
	ENUM Phase phase
	PRIVATE REAL Pr "Prandtl number"
	PRIVATE INTEGER ier,jx,jy // error codes
OBJECTS
	RefGeometryRecord geo
	SaturationProperties sat
INIT
	geo.setCylindricalGeometry(D,L)
	geo.EndSurfaces = 2
	assignInitAccu(init,f_out.fluid,P0,h0,T0,x0,Level0,P,h,rho,u)
	Pcrit = 73.77 // HardCoded // CRYO_PF_CritProp (f_out.fluid,fprop_pressure,ier)
	P_Pa = P0*1E5
	f_out.is_C = TRUE	
	As_coil = 2 * MATH.PI * (D_coil/2) * L_coil
	V_coil = MATH.PI * (D_coil/2)**2 * L_coil
	As_heater = 2 * MATH.PI * (D_heater/2) * L_heater
CONTINUOUS
	Vol = geo.V*1000 // Liters - useful in debugging
	delT = T - Tsat
	M = rho * (geo.V-V_coil)
	div_safe((rho-sat.rho_g)*L , sat.rho_l-sat.rho_g , lvl) // prevents division by zero errors
	level = ZONE (P>=Pcrit) L OTHERS min(L,max(0,lvl))
	z_top = z_bottom + L
	z_level = z_bottom + level
	level_frac = level/L
	f_out.hf = ZONE (phase==two_phase) spliceFunction(sat.h_l,sat.h_g,level/(z_out+0.5*D_out)-1,D_out/(2.0*(z_out+D_out/2.0)))
			OTHERS h
	Q_o = alpha_ref * geo.As_o * (T-tp_out.Tk[1])
	Q_ch = alpha_coil * As_coil * (tp_in_ch.Tk[1] - T)
	Q_eh = alpha_heater * As_heater * (tp_in_eh.Tk[1] - T) // positive if heat flowing in from heater
	tp_out.q[1] = Q_o
	tp_in_ch.q[1] = Q_ch
	tp_in_eh.q[1] = Q_eh

	P_Pa = P*1E5
	geo.V*(drho_dP*P_Pa'+drho_dh*h') = -f_out.m
	geo.V*((h*drho_dP-1)*P_Pa' + (h*drho_dh+rho)*h') = -f_out.mh - Q_o + Q_eh + Q_ch
	getStatePh(f_out.fluid,P,h,cp,drho_dP,drho_dh,gamma,k,mu,phase,Pr,rho,sigma,T,Tsat,u,vsound,x,sat,ier,jx,jy)

	dP_out = (sat.rho_g*(z_top - max(z_out,z_level)) + sat.rho_l*(max(z_out,z_level) - z_out))*9.806
	f_out.A = MATH.PI*(D_out/2)**2
	f_out.Gcrit = CRYO_Gcrit_fun(f_out.fluid,P*1e5,rho,T,0,f_out.v,vsound,ier,jy)
	f_out.rho = rho
	f_out.P = ZONE (AccountForGravity) P + dP_out/1e5 OTHERS P
	f_out.P_aux = P + geo.V*rho*vsound*1E-5/geo.As
	f_out.I = 0.5*geo.V/geo.As**2
	f_out.v = f_out.m/(rho*geo.Ap)
END COMPONENT



--------------------------------------------------------------------------------
// COMPLETE COMPONENTS
--------------------------------------------------------------------------------
ABSTRACT COMPONENT CO2Accumulator (
	BOOLEAN ConstantRefHTC = TRUE "If FALSE, uses Cooper pool boiling correlation",
	BOOLEAN DynamicRefHTC = FALSE "If TRUE, applies a low-pass filter to HTC calculations",
	BOOLEAN AccountforGravity = TRUE "If TRUE, applies static head to pressure calcs",
	BOOLEAN TsatInCelsius = TRUE "If FALSE SaturationTemp port outputs Tsat in Kelvin"
)
"Includes all the data and ports that EVERY CO2Accumulator will need"
// Makes the code for the individual component models smaller
PORTS
	OUT fluid f_out
	OUT analog_signal(n=1) SaturationTemp
	OUT analog_signal(n=1) Pressure
	OUT analog_signal(n=1) Level "Liquid level inside vessel in % of volume"
DATA
	ENUM PipeMat mat = SS_304 "Shell material"
	REAL D = 1 UNITS u_m "Shell diameter"
	REAL Th = 0.01 UNITS u_m "Shell thickness"
	REAL L = 1.5 UNITS u_m "Shell height"
	REAL alpha_ref = 1000 UNITS u_W_m2K "Refrigerant to shell heat transfer coefficient"
	REAL tauAlpha = 3 UNITS u_s "Time constant for low-pass HTC calculations"

	ENUM InitialConditionsAccu init = P_h
	REAL P0 = 20 UNITS u_bar "Initial pressure"
	REAL h0 = 250000 UNITS u_J_kg "Initial enthalpy"
	REAL T0 = 255 UNITS u_K "Initial temperature"
	REAL x0 = 0.1 UNITS no_units "Initial vapour quality"
	REAL Level0 = 20 UNITS u_pct "Initial liquid level in %"
	REAL Tw0 = 243.15 UNITS u_K "Initial wall temperature"
TOPOLOGY
	LumpedShell_insulated wall (
		mat = mat,
		Do = D,
		Th = Th,
		L = L,
		Tw0 = Tw0
	)
END COMPONENT



COMPONENT CO2Accumulator_prescribed IS_A CO2Accumulator (
	REAL D_out = 0.001 UNITS u_m "Outlet pipe diameter",
	REAL z_out = 0.005 UNITS u_m "Outlet pipe elevation above bottom"
)
"CO2 Accumulator with heating and cooling powers prescribed"
PORTS
	IN analog_signal(n=1) K_heat "Percentage of max heating power"
	IN analog_signal(n=1) K_cool "Percentage of max cooling power"
DATA
	REAL Q_cool_max = 1000 UNITS u_W "Maximum cooling capacity"
	REAL Q_heat_max = 1000 UNITS u_W "Maximum heating capacity"
DECLS
	INTEGER ier,jx
TOPOLOGY
	CO2AccumulatorBase_prescribed (AccountforGravity = AccountforGravity, ConstantRefHTC = ConstantRefHTC, DynamicRefHTC = DynamicRefHTC) ref (
		D = D-2*Th,
		L = L,
		D_out = {D_out},
		z_out = {z_out},
		alpha_ref = alpha_ref,
		tauAlpha = tauAlpha,
		Q_cool_max = Q_cool_max,
		Q_heat_max = Q_heat_max,
		init = init,
		P0 = P0,
		h0 = h0,
		T0 = T0,
		x0 = x0,
		Level0 = Level0
	)
	CONNECT ref.f_out[1] TO f_out
	CONNECT K_heat TO ref.K_heat
	CONNECT K_cool TO ref.K_cool
	CONNECT ref.tp_out TO wall.tp_in
CONTINUOUS
	SaturationTemp.signal[1] = ZONE (TsatInCelsius) CRYO_PF_Tsat_vs_p(f_out.fluid,ref.P,ier,jx)-273.15 OTHERS CRYO_PF_Tsat_vs_p(f_out.fluid,ref.P,ier,jx)
	Pressure.signal[1] = ref.P
	Level.signal[1] = ref.level_frac*100
END COMPONENT


COMPONENT CO2Accumulator_prescribed_new (
	REAL D_out = 0.01 UNITS u_m "Outlet port diameter",
	REAL z_out = 0.01 UNITS u_m "Outlet port height",
	BOOLEAN TsatInCelsius = TRUE "If FALSE, outputs saturation temperature in Kelvin"
)
"CO2 Accumulator with cartridge heater component and cooling spiral"
PORTS
	OUT fluid f_out
	IN analog_signal(n=1) K_heat
	IN analog_signal(n=1) K_cool
	OUT analog_signal(n=1) SaturationTemp
	OUT analog_signal(n=1) Pressure
	OUT analog_signal(n=1) Level
DATA
	ENUM PipeMat mat = SS_304 "Shell material"
	REAL D = 0.16828 UNITS u_m "Shell outer diameter"
	REAL Th = 0.007 UNITS u_m "Shell thickness"
	REAL L = 1.5 UNITS u_m "Shell height"
	REAL Q_heat_max = 1000 UNITS u_W "Maximum heating power"
	REAL Q_cool_max = 1000 UNITS u_W "Maximum cooling power"
	REAL alpha_ref = 1000 UNITS u_W_m2K "HTC between refrigerant to shell"
	ENUM InitialConditionsAccu init = P_h
	REAL P0 = 20 UNITS u_bar "Initial pressure"
	REAL h0 = 250000 UNITS u_J_kg "Initial enthalpy"
	REAL T0 = 255 UNITS u_K "Initial temperature"
	REAL x0 = 0.1 UNITS no_units "Initial vapour quality"
	REAL Level0 = 20 UNITS u_pct "Initial liquid level in %"
	REAL Tw0 = 293.15 UNITS u_K "Initial accumulator wall temperature"	
DECLS
	PRIVATE INTEGER ier,jx,jy
TOPOLOGY
	CO2AccumulatorBase_prescribed ref (
		D = D-2*Th,
		L = L,
		D_out = D_out,
		z_out = z_out,
		alpha_ref = alpha_ref,
		Q_heat_max = Q_heat_max,
		Q_cool_max = Q_cool_max,
		init = init,
		P0 = P0,
		h0 = h0,
		T0 = T0,
		x0 = x0,
		Level0 = Level0
	)
	LumpedShell_insulated wall (
		mat = mat,
		Do = D,
		Th = Th,
		L = L,
		Tw0 = Tw0
	)
	CONNECT ref.f_out[1] TO f_out
	CONNECT ref.tp_out TO wall.tp_in
	CONNECT K_cool TO ref.K_cool
	CONNECT K_heat TO ref.K_heat
CONTINUOUS
	SaturationTemp.signal[1] = ZONE (TsatInCelsius) ref.T -273.15 OTHERS ref.T
	Pressure.signal[1] = ref.P
	Level.signal[1] = ref.level_frac*100
END COMPONENT



COMPONENT CO2Accumulator_heater IS_A CO2Accumulator (
	REAL D_out = 0.001 UNITS u_m "Outlet pipe diameter",
	REAL z_out = 0.005 UNITS u_m "Outlet pipe elevation above bottom"
)
"CO2 Accumulator with cartridge heater component"
PORTS
	IN analog_signal(n=1) K_heat "Heater power control signal (in %)"
	IN analog_signal(n=1) K_cool "Cooling power control signal (in %)"
	OUT analog_signal(n=1) T_heater "Heater temperature"
DATA
	ENUM PipeMat heaterMat = SS_304 "Heater material"
	REAL D_heater = 0.01 UNITS u_m "Heater diameter"
	REAL L_heater = 0.25 UNITS u_m "Heater length"
	REAL alpha_heater_l = 1000 UNITS u_W_m2K "Liquid HTC heater to refrigerant"
	REAL alpha_heater_tp = 1000 UNITS u_W_m2K "Two-phase HTC heater to refrigerant"
	REAL alpha_heater_v = 100 UNITS u_W_m2K "Vapour HTC heater to refrigerant"
	REAL Q_cool_max = 1000 UNITS u_W "Maximum cooling power"
	REAL Q_heat_max = 1000 UNITS u_W "Maximum heater power"
	REAL Tw0_heater = 293.15 UNITS u_K "Initial heater temperature"
DECLS
	INTEGER ier, jx // error codes
TOPOLOGY
	CO2AccumulatorBase_heater (ConstantRefHTC=ConstantRefHTC, DynamicRefHTC=DynamicRefHTC, AccountforGravity=AccountforGravity) ref (
		D = D-2*Th,
		L = L,
		D_out = {D_out},
		z_out = {z_out},
		alpha_ref = alpha_ref,
		alpha_heater_l = alpha_heater_l,
		alpha_heater_tp = alpha_heater_tp,
		alpha_heater_v = alpha_heater_v,
		tauAlpha = tauAlpha,
		Q_cool_max = Q_cool_max,
		init = init,
		P0 = P0,
		h0 = h0,
		T0 = T0,
		x0 = x0,
		Level0 = Level0
	)
	CartridgeHeaterCtrl (n=1) heater (
		mat = heaterMat,
		D = D_heater,
		L = L_heater,
		Q_max = Q_heat_max,
		Tw0_i = Tw0_heater,
		Tw0_o = Tw0_heater
	)
	CONNECT K_heat TO heater.K
	CONNECT K_cool TO ref.K_cool
	CONNECT heater.tp_out TO ref.tp_in_eh
	CONNECT ref.tp_out TO wall.tp_in
	CONNECT ref.f_out[1] TO f_out
CONTINUOUS
	SaturationTemp.signal[1] = ZONE (TsatInCelsius) CRYO_PF_Tsat_vs_p(f_out.fluid,ref.P,ier,jx)-273.15 OTHERS CRYO_PF_Tsat_vs_p(f_out.fluid,ref.P,ier,jx)
	Pressure.signal[1] = ref.P
	Level.signal[1] = ref.level_frac*100
	T_heater.signal[1] = heater.T_avg - 273.15 // Output in °C
END COMPONENT



COMPONENT CO2Accumulator_externChiller IS_A CO2Accumulator
PORTS
	IN fluid f_in
	IN analog_signal(n=1) K_heat "Heater power in %"
DATA
	ENUM PipeMat heaterMat = SS_304 "Heater material"
	REAL D_heater = 0.01 UNITS u_m "Heater diameter"
	REAL L_heater = 0.25 UNITS u_m "Heater length"
	REAL D_in[1] = 0.01 UNITS u_m "Inlet pipe diameter" // See note at the top of this file
	REAL z_in[1] = 0.28 UNITS u_m "Inlet pipe height"
	REAL D_out[1] = 0.01 UNITS u_m "Outlet pipe diameter"
	REAL z_out[1] = 0.01 UNITS u_m "Outlet pipe height"
	REAL alpha_heater = 100 UNITS u_W_m2K "Heater to refrigerant HTC"
	REAL Q_heat_max = 1000 UNITS u_W "Maximum heater power"
TOPOLOGY
	CO2AccumulatorBase_externChiller (ConstantRefHTC=TRUE) Accumulator (
		D = D-2*Th,
		L = L,
		D_in = D_in,
		z_in = z_in,
		D_out = D_out,
		z_out = z_out,
		init = init,
		alpha_ref = alpha_ref,
		alpha_heater = alpha_heater,
		P0 = P0,
		h0 = h0,
		T0 = T0,
		x0 = x0,
		Level0 = Level0
	)
	CartridgeHeaterCtrl (n=1) heater (
		mat = heaterMat,
		D = D_heater,
		L = L_heater,
		Q_max = Q_heat_max,
		Tw0_i = Tw0,
		Tw0_o = Tw0
	)
	CONNECT f_in TO Accumulator.f_in[1]
	CONNECT Accumulator.f_out[1] TO f_out
	CONNECT K_heat TO heater.K
	CONNECT heater.tp_out TO Accumulator.tp_in_eh
	CONNECT Accumulator.tp_out TO wall.tp_in
CONTINUOUS
	Pressure.signal[1] = Accumulator.P
	IF (TsatInCelsius) INSERT
		SaturationTemp.signal[1] = Accumulator.Tsat - 273.15
	ELSE
		SaturationTemp.signal[1] = Accumulator.Tsat
	END IF
	Level.signal[1] = Accumulator.level_frac*100
END COMPONENT



COMPONENT CO2Accumulator_flowThrough IS_A CO2Accumulator
PORTS
	IN fluid f_in "Flow through port"
	IN fluid ch_in "Chiller coil inlet port"
	OUT fluid ch_out "Chiller coil outlet port"
	IN analog_signal(n=1) K_heat "Heater power in %"
	OUT analog_signal(n=1) T_heater "Heater temperature"
DATA
	ENUM PipeMat heaterMat = SS_304 "Heater material"
	ENUM PipeMat coilMat = SS_304 "Chiller coil material"
	REAL D_in[1] = 0.01 UNITS u_m "Inlet pipe diameter" // See note at the top of this file
	REAL z_in[1] = 0.28 UNITS u_m "Inlet pipe height"
	REAL D_out[1] = 0.01 UNITS u_m "Outlet pipe diameter"
	REAL z_out[1] = 0.01 UNITS u_m "Outlet pipe height"
	REAL D_heater = 0.01 UNITS u_m "Heater diameter"
	REAL L_heater = 0.25 UNITS u_m "Heater length"
	REAL D_coil = 0.01 UNITS u_m "Chiller coil outer diameter"
	REAL Th_coil = 0.001 UNITS u_m "Chiller coil wall thickness"
	REAL L_coil = 1 UNITS u_m "Chiller coil spiral total length"
	REAL Q_heat_max = 1000 UNITS u_W "Maximum heater power"
	REAL alpha_heater_v = 500 UNITS u_W_m2K "Heater-to-refrigerant vapour heat transfer coefficient"
	REAL alpha_heater_tp = 15000 UNITS u_W_m2K "Heater-to-refrigerant two-phase heat transfer coefficient"
	REAL alpha_heater_l = 5000 UNITS u_W_m2K "Heater-to-refrigerant liquid heat transfer coefficient"
	
	REAL alpha_coil = 1000 UNITS u_W_m2K "Coil-to-refrigerant heat transfer coefficient"
	
	REAL alpha_l_coil = 2500 UNITS u_W_m2K "Liquid HTC for chiller fluid in coil"
	REAL alpha_tp_coil = 10000 UNITS u_W_m2K "Two-phase HTC for chiller fluid in coil"
	REAL alpha_g_coil = 1500 UNITS u_W_m2K "Vapour HTC for chiller fluid in coil"
	REAL m_steady_coil = 0.2 UNITS u_kg_s "Flow rate for whichi coil fluid HTCs calculated"
	
	ENUM InitialConditions init_coil = Ph
	REAL P0_coil = 10 UNITS u_bar "Initial chiller spiral pressure"
	REAL h0_coil = 250000 UNITS u_J_kg "Initial chiller spiral enthalpy"
	REAL T0_coil = 293.15 UNITS u_K "Initial chiller spiral temperature"
	REAL x0_coil = 0.5 "Initial chiller spiral vapour quality"
	
	REAL Tw0_coil = 293.15 UNITS u_K  "Initial chiller spiral wall temperature"
	REAL Tw0_heater = 293.15 UNITS u_K "Initial heater wall temperature"
TOPOLOGY
	CO2AccumulatorBase_flowThrough (ConstantRefHTC=TRUE) Accumulator (
		D = D-2*Th,
		L = L,
		D_in = D_in,
		z_in = z_in,
		D_out = D_out,
		z_out = z_out,
		D_heater = D_heater,
		L_heater = L_heater,
		D_coil = D_coil-2*Th_coil,
		L_coil = L_coil,
		init = init,
		alpha_ref = alpha_ref,
		alpha_heater_v = alpha_heater_v,
		alpha_heater_tp = alpha_heater_tp,
		alpha_heater_l = alpha_heater_l,
		alpha_coil = alpha_coil,
		P0 = P0,
		h0 = h0,
		T0 = T0,
		x0 = x0,
		Level0 = Level0
	)
	CartridgeHeaterCtrl (n=1) heater (
		mat = heaterMat,
		D = D_heater,
		L = L_heater,
		Q_max = Q_heat_max,
		Tw0_i = Tw0_heater,
		Tw0_o = Tw0_heater
	)
	CV_Ph (ConstantHTC=TRUE, voidFractionModel=Smith) coolingSpiral (
		D = D_coil-2*Th_coil,
		L = L_coil,
		alpha_l = alpha_l_coil,
		alpha_tp = alpha_tp_coil,
		alpha_g = alpha_g_coil,
		m_steady = m_steady_coil,
		init = init_coil,
		P0 = P0_coil,
		h0 = h0_coil,
		T0 = T0_coil,
		x0 = x0_coil
	)
	AnnularWall_twoPort (n=1,Countercurrent = FALSE) coolingSpiralWall (
		mat = coilMat,
		Do = D_coil,
		Th = Th_coil,
		L = L_coil,
		Tw0_i = Tw0_coil,
		Tw0_o = Tw0_coil
	)	
	CONNECT f_in TO Accumulator.f_in[1]
	CONNECT Accumulator.f_out[1] TO f_out
	CONNECT ch_in TO coolingSpiral.f_in
	CONNECT coolingSpiral.f_out TO ch_out
	CONNECT coolingSpiral.tp_out TO coolingSpiralWall.tp_in
	CONNECT coolingSpiralWall.tp_out TO Accumulator.tp_in_ch
	CONNECT K_heat TO heater.K
	CONNECT heater.tp_out TO Accumulator.tp_in_eh
	CONNECT Accumulator.tp_out TO wall.tp_in
CONTINUOUS
	Pressure.signal[1] = Accumulator.P
	IF (TsatInCelsius) INSERT
		SaturationTemp.signal[1] = Accumulator.Tsat - 273.15
		T_heater.signal[1] = heater.T_avg - 273.15
	ELSE
		SaturationTemp.signal[1] = Accumulator.Tsat
		T_heater.signal[1] = heater.T_avg
	END IF
	Level.signal[1] = Accumulator.level_frac*100
END COMPONENT



COMPONENT CO2Accumulator_complete (
	REAL D_out = 0.01 UNITS u_m "Outlet port diameter",
	REAL z_out = 0.01 UNITS u_m "Outlet port height",
	BOOLEAN TempsInCelsius = TRUE "If FALSE, outputs heater and saturation temperatures in Kelvin"
)
"CO2 Accumulator with cartridge heater component and cooling spiral"
PORTS
	OUT fluid f_out
	IN fluid ch_in
	OUT fluid ch_out
	IN analog_signal(n=1) K "Heater power control signal"
	OUT analog_signal(n=1) SaturationTemp
	OUT analog_signal(n=1) Pressure
	OUT analog_signal(n=1) Level
	OUT analog_signal(n=1) T_heater
	OUT analog_signal(n=1) Q_ch
DATA
	ENUM PipeMat mat = SS_304 "Shell material"
	ENUM PipeMat coilMat = Copper "Spiral coil material"
	ENUM PipeMat heaterMat = SS_304 "Heater material"
	REAL D = 0.16828 UNITS u_m "Shell outer diameter"
	REAL Th = 0.007 UNITS u_m "Shell thickness"
	REAL L = 1.5 UNITS u_m "Shell height"
	REAL D_coil = 0.009525 UNITS u_m "Chiller coil outer diameter"
	REAL Th_coil = 0.000889 UNITS u_m "Chiller coil wall thickness"
	REAL L_coil = 11 UNITS u_m "Chiller coil length"
	REAL D_heater = 0.01 UNITS u_m "Heater diameter"
	REAL L_heater = 0.25 UNITS u_m "Heater length"

	REAL Q_max_heat = 1000 UNITS u_W "Maximum power of cartridge heater"
	REAL alpha_ref = 1000 UNITS u_W_m2K "HTC between refrigerant to shell"
	REAL alpha_coil = 1000 UNITS u_W_m2K "HTC for CO2 to spiral heat transfer"
	REAL alpha_heater = 1000 UNITS u_W_m2K "HTC for CO2 to heater heat transfer"
	REAL alpha_l_coil = 1000 UNITS u_W_m2K "Chiller fluid liquid HTC"
	REAL alpha_tp_coil = 1000 UNITS u_W_m2K "Chiller fluid two-phase HTC"
	REAL alpha_g_coil = 1000 UNITS u_W_m2K "Chiller vapour HTC"
	REAL m_steady_coil = 0.02 UNITS u_kg_s "Mass flow raten at which chiller HTC calculated"	

	ENUM InitialConditionsAccu init = P_h
	REAL P0 = 20 UNITS u_bar "Initial pressure"
	REAL h0 = 250000 UNITS u_J_kg "Initial enthalpy"
	REAL T0 = 255 UNITS u_K "Initial temperature"
	REAL x0 = 0.1 UNITS no_units "Initial vapour quality"
	REAL Level0 = 20 UNITS u_pct "Initial liquid level in %"
	
	ENUM InitialConditions init_coil = Ph
	REAL P0_coil = 10 UNITS u_bar "Initial chiller spiral pressure"
	REAL h0_coil = 250000 UNITS u_J_kg "Initial chiller spiral enthalpy"
	REAL T0_coil = 293.15 UNITS u_K "Initial chiller spiral temperature"
	REAL x0_coil = 0.5 "Initial chiller spiral vapour quality"
	
	REAL Tw0_wall = 293.15 UNITS u_K "Initial accumulator wall temperature"
	REAL Tw0_coil = 293.15 UNITS u_K  "Initial chiller spiral wall temperature"
	REAL Tw0_heater = 293.15 UNITS u_K "Initial heater wall temperature"
	
DECLS
	PRIVATE INTEGER ier,jx,jy
TOPOLOGY
	CO2AccumulatorRefrigerantVolume ref (
		D = D-2*Th,
		L = L,
		D_coil = D_coil,
		L_coil = L_coil,
		D_out = D_out,
		z_out = z_out,
		
		alpha_ref = alpha_ref,
		alpha_coil = alpha_coil,
		alpha_heater = alpha_heater,
		
		init = init,
		P0 = P0,
		h0 = h0,
		T0 = T0,
		x0 = x0,
		Level0 = Level0
	)
	CV_Ph (ConstantHTC=TRUE, voidFractionModel=Smith) coolingSpiral (
		D = D_coil-2*Th_coil,
		L = L_coil,
		alpha_l = alpha_l_coil,
		alpha_tp = alpha_tp_coil,
		alpha_g = alpha_g_coil,
		m_steady = m_steady_coil,
		init = init_coil,
		P0 = P0_coil,
		h0 = h0_coil,
		T0 = T0_coil,
		x0 = x0_coil
	)
	AnnularWall_twoPort (n=1,Countercurrent = FALSE) coolingSpiralWall (
		mat = coilMat,
		Do = D_coil,
		Th = Th_coil,
		L = L_coil,
		Tw0_i = Tw0_coil,
		Tw0_o = Tw0_coil
	)	
	LumpedShell_insulated wall (
		mat = mat,
		Do = D,
		Th = Th,
		L = L,
		Tw0 = Tw0_wall
	)
	CartridgeHeaterCtrl (n=1) heater (
		mat = heaterMat,
		D = D_heater,
		L = L_heater,
		Q_max = Q_max_heat,
		Tw0_i = Tw0_heater,
		Tw0_o = Tw0_heater
	)
	CONNECT ref.f_out TO f_out
	CONNECT ch_in TO coolingSpiral.f_in
	CONNECT ch_out TO coolingSpiral.f_out
	CONNECT coolingSpiral.tp_out TO coolingSpiralWall.tp_in
	CONNECT coolingSpiralWall.tp_out TO ref.tp_in_ch
	CONNECT heater.tp_out TO ref.tp_in_eh
	CONNECT K TO heater.K
	CONNECT ref.tp_out TO wall.tp_in
CONTINUOUS
	SaturationTemp.signal[1] = ZONE (TempsInCelsius) CRYO_PF_Tsat_vs_p(f_out.fluid,ref.P,ier,jx)-273.15 OTHERS CRYO_PF_Tsat_vs_p(f_out.fluid,ref.P,ier,jx)
	Pressure.signal[1] = ref.P
	Level.signal[1] = ref.level_frac*100
	T_heater.signal[1] = ZONE(TempsInCelsius) heater.Tw[1] - 273.15 OTHERS heater.Tw[1]
	Q_ch.signal[1] = -coolingSpiral.Q
END COMPONENT



--------------------------------------------------------------------------------
// DAMPER MODEL
--------------------------------------------------------------------------------
COMPONENT DamperBase IS_A CO2AccumulatorBase
"Base component for Damper. Assumes no ambient heat transfer"
PORTS
	IN thermal(n=1) tp_in_eh "Heater thermal port"
DATA
	REAL D_heater = 0.1 UNITS u_m "Heater outer diameter"
	REAL L_heater = 1 UNITS u_m "Heater length"
	REAL alpha_heater = 1000 UNITS u_W_m2K "Heat transfer coefficient with cartridge heater"
DECLS
	REAL Q_eh UNITS u_W "Heat flow from electrical heater"
	DISCR REAL As_heater UNITS u_m2 "Heater outer surface area"
INIT
	As_heater = 2 * MATH.PI * (D_heater/2) * L_heater
CONTINUOUS
	Q_eh = alpha_heater * As_heater * (tp_in_eh.Tk[1] - T) // positive if heat flow FROM heater
	tp_in_eh.q[1] = Q_eh
<:mas>	geo.V*rho' = - f_out[1].m
<:ene>	geo.V*(u*rho'+ rho*u') = - f_out[1].mh + Q_eh
<:st>	getStateRU(f_out[1].fluid,rho,u,cp,drho_dP,drho_dh,gamma,h,k,mu,P,phase,Pr,\
			sigma,T,Tsat,vsound,x,sat,ier,jx,jy)
END COMPONENT



COMPONENT Damper
PORTS
	OUT fluid f_out
	OUT analog_signal(n=1) T_heater "Heater temperature"
	IN analog_signal(n=1) DO "DigitalOutput object for heater power"
DATA
	ENUM PipeMat mat = SS_304 "Shell material"
	ENUM PipeMat heaterMat = SS_304 "Heater material"
	REAL D = 0.10 UNITS u_m "Outer diameter"
	REAL Thw = 0.001 UNITS u_m "Shell thickness"
	REAL L = 0.15 UNITS u_m "Height"
	REAL D_heater = 0.01 UNITS u_m "Heater diameter"
	REAL L_heater = 0.25 UNITS u_m "Heater length"

	REAL alpha_heater = 1000 UNITS u_W_m2K "HTC for CO2 to heater heat transfer"
	REAL Q_max = 1000 UNITS u_W "Maximum heating power"
	ENUM InitialConditionsAccu init = P_h
	REAL P0 = 20 UNITS u_bar "Initial pressure"
	REAL h0 = 250000 UNITS u_J_kg "Initial enthalpy"
	REAL T0 = 255 UNITS u_K "Initial temperature"
	REAL x0 = 0.1 UNITS no_units "Initial vapour quality"
	REAL Level0 = 20 UNITS u_pct "Initial liquid level in %"
	REAL Tw0 = 293.15 UNITS u_K "Initial heater and wall temperature"
DECLS
	PRIVATE INTEGER ier,jx,jy
TOPOLOGY
	DamperBase ref (
		D = D-2*Thw,
		L = L,
		alpha_heater = alpha_heater,
		init = init,
		P0 = P0,
		h0 = h0,
		T0 = T0,
		x0 = x0,
		Level0 = Level0
	)
	CartridgeHeaterCtrl (n=1) heater (
		mat = heaterMat,
		D = D_heater,
		L = L_heater,
		Q_max = Q_max,
		Tw0_i = Tw0,
		Tw0_o = Tw0
	)
	LumpedShell_insulated wall (
		mat = mat,
		Do = D,
		Th = Thw,
		L = L,
		Tw0 = Tw0
	)
	CONNECT ref.f_out[1] TO f_out
	CONNECT heater.tp_out TO ref.tp_in_eh
	CONNECT ref.tp_out TO wall.tp_in
	CONNECT DO TO heater.K
CONTINUOUS
	T_heater.signal[1] = heater.T_avg - 273.15 // hardcoded in °C
END COMPONENT



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


ABSTRACT COMPONENT TankTwoFluid
"Base component for two-fluid accumulator"
PORTS
	OUT fluid f_out
DATA
	// Geometry
	REAL D = 0.5 UNITS u_m "Diameter"
	REAL L = 1.96 UNITS u_m "Height"
	REAL D_out = 0.01 UNITS u_m "Outlet port diameter"
	REAL z_out = 0.01 UNITS u_m "Elevation of bottom of outlet port, 0 if underneath the vessel"
	// Design parameters
	REAL tau_l = 1 UNITS u_s "Time constant for evaporation"
	REAL tau_v = 1 UNITS u_s "Time constant for condensation"
	// Initial conditions
	ENUM InitialConditions2FM init = P_level "Specify initial condition variables"
	REAL P0 = 20 UNITS u_bar "Initial pressure"
	REAL Tsat0 = 20 UNITS u_C "Initial saturation temperature (°C)"
	REAL Level0 = 20 UNITS u_pct "Initial liquid level (%)"
	REAL x0 = 0.1 "Initial vapour quality"
DECLS
	DISCR REAL V UNITS u_m3 "Volume"
	DISCR REAL Vol UNITS u_dm3 "Volume in Litres"
	DISCR REAL Ac UNITS u_m2 "Cross-sectional area"
	DISCR REAL Pcrit UNITS u_bar "Critical pressure"
	DISCR REAL rhocrit UNITS u_kg_m3 "Critical density"
	REAL drhol_dP,drhol_dh,drhov_dP,drhov_dh "Partial density derivatives"
	REAL gamma "Void fraction"
	REAL h UNITS u_J_kg "Averaged enthalpy"
	REAL hl, hv UNITS u_J_kg "Liquid and vapour phase enthalpies"
	REAL hls, hvs "Liquid and vapour saturation enthalpies"
	REAL level UNITS u_m "Height of liquid column inside the vessel"
	REAL level_frac "Fractional liquid level in the vessel"
	REAL Ml, Mv UNITS u_kg "Liquid and vapour phase masses"
	REAL M UNITS u_kg "Total refrigerant charge"
	REAL m_evap, m_cond UNITS u_kg_s "Evaporation and condensation flow rates"
	REAL mout_l, mout_v UNITS u_kg_s "Liquid and vapour portions of the outlet flow rate"
	REAL mhout_l, mhout_v UNITS u_J_kg
	REAL P UNITS u_bar "Pressure"
	REAL P_Pa UNITS u_Pa "Pressure in Pascal"
	REAL rhol, rhov UNITS u_kg_m3 "Liquid and vapour densities"
	REAL Tl, Tv UNITS u_K "Temperatures of liquid and vapour phases"
	REAL Vl, Vv UNITS u_m3 "Volume occupied by phases"
	REAL x "Vapour quality"
	REAL d1,d2 // dummy variables
	INTEGER ier,jx,jy "Error codes"
OBJECTS
	RefGeometryRecord geo
INIT
	geo.setAccumulatorGeometry(D,L) // TODO: account for ellpisoid-cap surface area
	Vol = geo.V*1000
	assignInitAccu2FM(init,f_out.fluid,P0,Tsat0+273.15,Level0,x0,P,hl,hv,level_frac)
	Pcrit = 73.77 //CRYO_PF_CritProp(f_out.fluid,fprop_pressure,ier)*1e-5
	rhocrit = CRYO_PF_CritProp(f_out.fluid,fprop_density,ier)
	f_out.is_C = TRUE // sets Accumulator as capacitive component
CONTINUOUS
	P := P_Pa * 1e-5
	rp_all(P_Pa,hl,hv,rhol,rhov,d1,d2,drhol_dP,drhol_dh,drhov_dP,drhov_dh) // dummy vars d1/d2 used to prevent crash in supercrticial case
	hls := rp_hPx(P_Pa,0)
	hvs := rp_hPx(P_Pa,1)
	Tl := CRYO_PF_prop_vs_ph(f_out.fluid,P,hl,fprop_temperature,ier,jx,jy)
	Tv := CRYO_PF_prop_vs_ph(f_out.fluid,P,hv,fprop_temperature,ier,jx,jy)	
	Vv := geo.V - Vl
	Ml := rhol*Vl
	Mv := rhov*Vv
	M := Ml + Mv
	x := ZONE (P<Pcrit) Mv/M OTHERS 1
	h := x*hv + (1-x)*hl
	level_frac := Vl/V // fractional liquid level
	gamma := 1-level_frac
	level := level_frac*L
	// Phase change flow rates
	m_evap := ZONE(hl<hls) 0 OTHERS tau_l*rhol*Vl*(hl-hls)/(hvs-hls)
	m_cond := ZONE(hv>hvs) 0 OTHERS tau_v*rhov*Vv*(hvs-hv)/(hvs-hls)
	// Port balance equations
	mout_l = ZONE(level>z_out) f_out.m OTHERS 0
	mout_v = ZONE(level>z_out) 0 OTHERS f_out.m
	mhout_l = semiLinear(f_out.m,0,f_out.hb-hl)
	mhout_v = 0 // assume always submerged TODO:improve
	// Governing equations
<mal>	Vl*(drhol_dP*P_Pa' + drhol_dh*hl') + rhol*Vl' = m_cond - m_evap
<mav>	Vv*(drhov_dP*P_Pa' + drhov_dh*hv') - rhov*Vl' = m_evap - m_cond
<enl>	Vl*(rhol*hl' - P_Pa') = m_cond*(hls-hl) - m_evap*(hvs-hl) - mhout_l
<env>	Vv*(rhov*hv' - P_Pa') = m_evap*(hvs-hv) - m_cond*(hls-hv) - mhout_v
	// Port equations	
	f_out.P = P
	f_out.P_aux = P
	f_out.hf = hl
	f_out.rho = rhol
	f_out.v = 1
	f_out.A = 1
	f_out.Gcrit = 1
END COMPONENT



ABSTRACT COMPONENT AccumulatorTwoFluid
PORTS
	OUT fluid f_out
	OUT analog_signal(n=1) Pressure
	OUT analog_signal(n=1) SaturationTemp
	OUT analog_signal(n=1) Level
DATA
	ENUM PipeMat mat = SS_304 "Shell material"
	REAL D = 1 UNITS u_m "Shell diameter"
	REAL Th = 0.01 UNITS u_m "Shell thickness"
	REAL L = 1.5 UNITS u_m "Shell height"
	REAL D_out = 0.01 UNITS u_m "Outlet port diameter"
	REAL z_out = 0.01 UNITS u_m "Elevation of bottom of outlet port, 0 if underneath the vessel"
	// Design parameters
	REAL tau_l = 1 UNITS u_s "Time constant for evaporation"
	REAL tau_v = 1 UNITS u_s "Time constant for condensation"
	REAL alpha_ref = 1000 UNITS u_W_m2K "Refrigerant to shell heat transfer coefficient"
	REAL tauAlpha = 3 UNITS u_s "Time constant for low-pass HTC calculations"
	// Initial Conditions
	ENUM InitialConditions2FM init = P_level
	REAL P0 = 20 UNITS u_bar "Initial pressure"
	REAL h0 = 250000 UNITS u_J_kg "Initial enthalpy"
	REAL x0 = 0.1 UNITS no_units "Initial vapour quality"
	REAL Level0 = 20 UNITS u_pct "Initial liquid level in %"
	REAL Tw0 = 243.15 UNITS u_K "Initial wall temperature"
TOPOLOGY
	LumpedShell_insulated wall (
		mat = mat,
		Do = D,
		Th = Th,
		L = L,
		Tw0 = Tw0
	)
END COMPONENT



COMPONENT AccumulatorTwoFluid_prescribed IS_A AccumulatorTwoFluid
PORTS 
	IN analog_signal(n=1) K_heat
	IN analog_signal(n=1) K_cool
DATA
	REAL Q_cool_max = 1000 UNITS u_W "Maximum cooling capacity"
	REAL Q_heat_max = 1000 UNITS u_W "Maximum heating capacity"
TOPOLOGY
	
END COMPONENT