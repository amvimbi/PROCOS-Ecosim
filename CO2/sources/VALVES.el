/*-----------------------------------------------------------------------------------------
 LIBRARY: CO2
 FILE: VALVES
 AUTHOR: Viren Bhanot
 COMPANY: CERN
 DESCRIPTION: Valve Models
 CREATION DATE: 22/08/2016
-----------------------------------------------------------------------------------------*/
USE CRYOLIB
USE CONTROL
USE MATH
USE PLC
USE THERMAL
USE PORTS_LIB



/* 
 * INDEX
 * 1. dP
 * 2. GenericOrifice
 * 3. ControlValve
 * 4. TXV
 * 5. BackPressureRegulator
 * 6. CheckValve
 * 7. ReversingValve
 */



ENUM FlowEquation = {Swagelok, DischargeCoefficient, Quadratic}
ENUM dPSmoothing = {dP_regRoot2, dP_regPowGen1, dP_regPowGen3}


FUNCTION NO_TYPE getSwagelokFlow (
	IN ENUM ChemName fluid,
	IN REAL Cv,
	IN REAL P1 UNITS u_bar,
	IN REAL h1 UNITS u_J_kg,
	IN REAL P2 UNITS u_bar,
	IN REAL h2 UNITS u_J_kg,
	OUT REAL mDot UNITS u_kg_s,
	OUT REAL vDot "Volumetric flow rate (L/min)"
)
"Returns mass and volumetric flow rates through Swagelok valves"
// smooths the flow rate through the two-phase region by weighting with vapour quality
// includes choked flow calculation
DECLS
	REAL dP UNITS u_bar
	REAL Gl, Gv "Specific gravity"
	REAL hin UNITS u_J_kg
	REAL mDotl, mDotv
	REAL Pin UNITS u_bar
	REAL Pout UNITS u_bar
	REAL rhoin UNITS u_kg_m3
	REAL rhoRef UNITS u_kg_m3
	REAL Tin UNITS u_K "Inlet fluid temperature"
	REAL x "Inlet vapour quality"
	REAL vDotl, vDotv
	CONST REAL N1 = 14.42
	CONST REAL N2 = 6950
	INTEGER a,b,c // error codes
BODY
	dP = P1 - P2
	Pin = donor_cell(dP,P1,P2)
	hin = donor_cell(dP,h1,h2)
	Pout = donor_cell(dP,P2,P1)
	rhoin = CRYO_PF_prop_vs_ph(fluid,Pin,hin,fprop_density,a,b,c)
	Tin = CRYO_PF_prop_vs_ph(fluid,Pin,hin,fprop_temperature,a,b,c)
	x = max(0,min(CRYO_PF_prop_vs_ph(fluid,Pin,hin,fprop_quality,a,b,c),1))
	rhoRef = CRYO_PF_prop_vs_pT(fluid,1.01325,293.15,fprop_density,a,b,c)
	Gl = rhoin/998
	Gv = rhoRef/1.205
	// Weighted average of vapour and liquid flow rates
	// Makes a smooth transition between single and two-phase
	IF(Pin < CO2_CRITPR) THEN
		// Vapour flow
		IF (Pout>=0.5*Pin) THEN
			vDotv = N2*Cv*Pin * (1-2*abs(dP)/(3*Pin)) * regRoot2(dP,1e-4) / sqrt(Pin*Gv*Tin) // low dP
		ELSE
			vDotv = sign(dP) * 0.471 * N2*Cv*Pin * regRoot2(1/(Gv*Tin),1e-4) // high dP (choked flow)
		END IF
		mDotv =  vDotv * 1e-3/60 * rhoin * (Tin/273.15) * (1/Pin)
		// Liquid flow
		vDotl = N1*Cv*regRoot2(dP/Gl,1e-3)
		mDotl = vDotl * 1e-3/60 * rhoin
		// Overall flow
		vDot = x*vDotv + (1-x)*vDotl
		mDot = x*mDotv + (1-x)*mDotl
	ELSE
		// Liquid flow
		vDot = N1*Cv*regRoot2(dP/Gl,1e-3)
		mDot = vDot * 1e-3/60 * rhoin
	END IF
END FUNCTION



COMPONENT dP IS_A CRYOLIB.AbsJunction (
	ENUM dPSmoothing choice = dP_regPowGen3 "Smoothing method. See note in VALVES.el for explanation"
)
"Static pressure drop component. Used often between two capacitive components to ensure cohesive boundary conditions"
DATA
	REAL m0 = 0.01 UNITS u_kg_s "Nominal mass flow rate"
	REAL dP0 = 0.01 UNITS u_bar "Nominal pressure drop"
	REAL dP_small = 0.001 UNITS u_bar "Threshold for low-mass-flow-rate smoothing"
CONTINUOUS
<m>	IF (choice==dP_regRoot2) INSERT
		m := m0 * regRoot2(f_in.P-f_out.P,dP_small)/sqrt(dP0)
	ELSEIF (choice==dP_regPowGen1) INSERT
		m := m0 * regPowGen(f_in.P-f_out.P,dP_small)/sqrt(dP0)
	ELSE
		m := m0 * regPowGen(f_in.P-f_out.P,dP_small,0.5,3)/sqrt(dP0)
	END IF
END COMPONENT



COMPONENT GenericOrifice IS_A AbsJunction (
	BOOLEAN SwagelokEquation = FALSE "If FALSE uses Discharge Coefficient equation"
)
"Short tube orifice (fixed diameter opening)"
/* Verification
 * To verify this model, Swagelok's own CV plots can be used. Three cases are studied
 * 1. Liquid flow
 * 2. Gas flow with low dP
 * 3. Gas flow with high dP
 */
DATA
	REAL D = 0.002 UNITS u_m "Orifice Diameter"
	REAL Cv = 0. UNITS no_units "Flow Coefficient, or Discharge coefficient if not Swagelok"
DECLS
	REAL A UNITS u_m2 "Opening area"
	REAL dP UNITS u_bar "Pressure drop across valve"
	REAL Gl, Gv "Specific gravity"
	REAL hin UNITS u_J_kg "Inlet specific enthalpy"
	REAL Pin UNITS u_bar "Inlet pressure"
	REAL Pout UNITS u_bar "Outlet pressures"
	REAL rho UNITS u_kg_m3 "Inlet fluid density"
	REAL rhoRef UNITS u_kg_m3 "Reference density at 1 atm and 20°C"
	REAL Tin UNITS u_K "Inlet temperature"
	REAL vDot "Volumetric flow rate (L/min)"
	REAL x "Vapour quality"
	CONST REAL N1 = 14.42 "Constant for Swagelok liquid CV for Bar and L/min"
	CONST REAL N2 = 6950 "Constant for Swagelok gas CV for Bar and L/min"
	INTEGER ier,xx,yy "Error codes"
CONTINUOUS
	dP := f_in.P - f_out.P
	Pin = donor_cell(dP,f_in.P,f_out.P)
	Pout = donor_cell(dP,f_out.P,f_in.P)
	hin = donor_cell(dP,f_in.hf,f_in.hb)
	rho = donor_cell(dP,f_in.rho,f_out.rho)
	Tin = CRYO_PF_prop_vs_ph(f_in.fluid,Pin,hin,fprop_temperature,ier,xx,yy)
	x = CRYO_PF_prop_vs_ph(f_in.fluid,Pin,hin,fprop_quality,ier,xx,yy)
	rhoRef = CRYO_PF_prop_vs_pT(f_in.fluid,1.01325,293.15,fprop_density,ier,xx,yy)
	Gl = rho/998
	Gv = rhoRef/1.205
<m> 	IF (SwagelokEquation) INSERT
		getSwagelokFlow(f_in.fluid,Cv,f_in.P,f_in.h,f_out.P,f_out.h,m,vDot)
	ELSE
		m = Cv * D**2 * regRoot2(dP,0.01) * sqrt(rho*1e5) // discharge coefficient equation
		vDot = m / rho * (60/1000) // Liter per minute
	END IF
	A = MATH.PI*D**2/4
END COMPONENT

COMPONENT ControlValve IS_A GenericOrifice (
	BOOLEAN delayOpening = FALSE "If TRUE, adds a time constant for changes in opening"
)
"Generic valve with controllable orifice opening"
PORTS
	IN analog_signal(n=1) Pos "Position request (%)"
	OUT analog_signal(n=1) feedback "Feedback position (%)"
DATA
	REAL opening0 = 100 UNITS u_pct "Initial opening"
	REAL tau = 3 UNITS u_s "Time constant for opening changes"
DECLS
	REAL mDot_dummy, vDot_dummy
	REAL opening UNITS u_pct "Valve opening"
INIT
	opening = opening0
CONTINUOUS
	feedback.signal[1] = opening
	IF (delayOpening) INSERT
		opening' = (Pos.signal[1] - opening)/tau
	ELSE
		opening = Pos.signal[1]
	END IF
<:m> 	IF (SwagelokEquation) INSERT
		getSwagelokFlow(f_in.fluid,Cv,f_in.P,f_in.h,f_out.P,f_out.h,mDot_dummy,vDot_dummy)
		m = opening/100 * mDot_dummy
		vDot = opening/100 * vDot_dummy
	ELSE
		m = opening/100 * Cv * D**2 * regRoot2(dP,0.01) * sqrt(rho*1e5)
		vDot = m / rho * (60/1000)
	END IF
END COMPONENT


COMPONENT TXV IS_A AbsJunction
"Thermostatic expansion valve. Includes delay in response of the sensor bulb"
PORTS
	IN analog_signal(n=1) Psuc
	IN analog_signal(n=1) Tsuc
DATA
	REAL D_orif = 3.86E-3 UNITS u_m "Orifice diameter"
	REAL D_rod = 2.69E-3 UNITS u_m "Rod (Stem) diameter"
	REAL D_pin = 4.74E-3 UNITS u_m "Pin (head) diameter"
	REAL H_Pin = 1.25E-3 UNITS u_m "Pin height"
	REAL Amin = 1E-6 UNITS u_m2 "Minimum opening area"
	REAL Cv = 0.28 UNITS no_units "Flow Coefficient"
	REAL tau = 90 UNITS u_s	"Bulb time constant"
	REAL offset = 1.5E5 UNITS u_Pa "Initial spring force"
	REAL Gfac = 1.5E8 UNITS no_units "G-factor"
	REAL Tb0 = 299 UNITS u_K "Initial bulb temperature"
DECLS
	REAL Amax
	REAL Ahor
	REAL Aeff
	REAL dP
	REAL h
	REAL kv
	REAL opening
	REAL P
	REAL Pb
	REAL rho
	REAL T	
	REAL Tb
	REAL y "Actual valve stem displacement (cannot be negative)"
	REAL yr "Calculated displacement of the valve stem - higher yr means valve is more open"
	REAL deltaP
	DISCR REAL y0
	INTEGER ier, jx, jy
INIT 
	Tb = Tb0
	y0 = (D_orif-D_rod)/(D_pin-D_rod)*H_Pin		
CONTINUOUS
	dP = f_in.P - f_out.P
	P = IF(dP>=0) f_in.P ELSE f_out.P
	h = IF(dP>=0) f_in.h ELSE f_out.h
	rho = IF(dP>=0) f_in.rho ELSE f_out.rho
	T = CRYO_PF_prop_vs_ph(f_in.fluid,P,h,fprop_temperature,ier,jx,jy)
	Pb = CRYO_PF_psat_vs_T(f_in.fluid,Tb,ier,jy)
	deltaP = Pb - Psuc.signal[1]
	yr = (Pb*1E5 - Psuc.signal[1]*1E5 - offset) / Gfac // force balance on the diaphragm	
	y = max(0,yr)
	Amax = MATH.PI/4*(D_orif**2 - D_rod**2) // maximum possible opening area
	Ahor = MATH.PI/4*D_orif**2 - MATH.PI*(D_rod/2 + max(0,(y0-y)) * (D_pin-D_rod)/(2*H_Pin))**2 // partially-closed area calculation
	Aeff = max(Amin,min(Amax,Ahor))
	kv = Cv*Aeff
	m = kv * regRoot2(dP,0.1) * sqrt(rho*1E5)
	Tb' = (Tsuc.signal[1] - Tb)/tau
	opening = Aeff/Amax
END COMPONENT



ENUM ChokeCalculation = {IEC,Dhar}

COMPONENT DetailedTXV IS_A AbsJunction (
	ENUM ChokeCalculation method = IEC "IEC or Dhar"
)
"Detailed Thermostatic expansion valve. Includes choked flow calculations"
PORTS
	IN analog_signal(n=1) Psuc
	IN analog_signal(n=1) Tsuc
DATA
	REAL D_orif = 3.86E-3 UNITS u_m "Orifice diameter"
	REAL D_rod = 2.69E-3 UNITS u_m "Rod (Stem) diameter"
	REAL D_pin = 4.74E-3 UNITS u_m "Pin (head) diameter"
	REAL H_Pin = 1.25E-3 UNITS u_m "Pin height"
	REAL Amin = 1E-6 UNITS u_m2 "Minimum opening area"
	REAL Cv = 0.28 UNITS no_units "Flow Coefficient"
	REAL tau = 90 UNITS u_s	"Bulb time constant"
	REAL offset = 1.5E5 UNITS u_Pa "Initial spring force"
	REAL Gfac = 1.5E8 UNITS no_units "G-factor"
	REAL Tb0 = 299 UNITS u_K "Initial bulb temperature"
DECLS
	REAL Amax
	REAL Ahor
	REAL Aeff
	REAL dP
	REAL h
	REAL kv
	REAL opening
	REAL P
	REAL Pb
	REAL rho
	REAL T	
	REAL Tb
	REAL y "Actual valve stem displacement (cannot be negative)"
	REAL yr "Calculated displacement of the valve stem - higher yr means valve is more open"
	REAL deltaP
	DISCR REAL y0
	INTEGER ier, jx, jy
	REAL Pup, Pdown
	CONST REAL xt = 0.87
	CONST REAL Gamma = 1.3
	CONST REAL N6 = 2.73
	REAL Y
	REAL FGamma
	REAL chokelimit
	REAL K
	REAL FlowCoeff
	REAL x
INIT 
	Tb = Tb0
	y0 = (D_orif-D_rod)/(D_pin-D_rod)*H_Pin
CONTINUOUS
	// OPENING AREA CALCULATIONS
	h = IF(dP>=0) f_in.h ELSE f_out.h
	T = CRYO_PF_prop_vs_ph(f_in.fluid,P,h,fprop_temperature,ier,jx,jy)
	Pb = CRYO_PF_psat_vs_T(f_in.fluid,Tb,ier,jy)
	deltaP = Pb - Psuc.signal[1]
	yr = (Pb*1E5 - Psuc.signal[1]*1E5 - offset) / Gfac // force balance on the diaphragm	
	y = max(0,yr)
	Amax = MATH.PI/4*(D_orif**2 - D_rod**2) // maximum possible opening area
	Ahor = MATH.PI/4*D_orif**2 - MATH.PI*(D_rod/2 + max(0,(y0-y)) * (D_pin-D_rod)/(2*H_Pin))**2 // partially-closed area calculation
	Aeff = max(Amin,min(Amax,Ahor))
	FlowCoeff = Cv*Aeff
	Tb' = (Tsuc.signal[1] - Tb)/tau
	opening = Aeff/Amax
	// FLOW RATE CALCULATIONS
	dP = f_in.P - f_out.P
	Pup = IF (dP>=0) f_in.P*1e5 ELSE f_out.P*1e5
	Pdown = IF (dP>=0) f_out.P*1e5 ELSE f_in.P*1e5
	rho = IF (dP>=0) f_in.rho ELSE f_out.rho
	K = CRYO_PF_prop_vs_ph(f_in.fluid,Pup*1e-5,CRYO_PF_prop_vs_Px(f_in.fluid,Pup*1e-5,1,fprop_enthalpy,ier,jx,jy)+10e3,fprop_cp,ier,jx,jy)/1000
	P = IF (dP>=0) f_in.P ELSE f_out.P
	IF (method==IEC) INSERT
		x = dP*1e5/Pup
		chokelimit = FGamma*xt
		m = IF (abs(x)<chokelimit) 1/3600*FlowCoeff*N6*Y*regRoot2(x*Pup/1e3*rho,0.001) ELSE 1/3600*FlowCoeff*2/3*N6*regRoot2(sign(dP*1e5)*FGamma*xt*Pup/1e3*rho, 0.001)
	ELSE
		x = Pdown/Pup
		chokelimit = (2/(K+1))**(K/(K-1))
		m = ZONE(abs(x)>chokelimit) FlowCoeff*regRoot2(sign(dP)*(1-x)*Pup*rho, 0.001) \
			OTHERS FlowCoeff*regRoot2(sign(dP)*(1-chokelimit)*Pup*rho, 0.001)
	END IF
	FGamma = Gamma/1.40
	Y = 1 - sign(dP)*x/(3*FGamma*xt) // Don't use abs() chattering
END COMPONENT



COMPONENT BackPressureRegulator IS_A AbsJunction
"Test model based on the TXV diaphragm model"
DATA
	REAL D_orif = 4E-3 UNITS u_m "Orifice diameter"
	REAL D_rod = 2E-3 UNITS u_m "Rod (Stem) diameter"
	REAL D_pin = 4E-3 UNITS u_m "Pin (head) diameter"
	REAL H_Pin = 1E-3 UNITS u_m "Pin height"
	REAL Amin = 1E-8 UNITS u_m2 "Minimum opening area"
	REAL Cv = 0.20 UNITS no_units "Flow Coefficient"
	REAL offset = 5E5 UNITS u_Pa "Initial spring force"
	REAL Gfac = 2E8 UNITS no_units "G-factor"
DECLS
	REAL Amax
	REAL Ahor
	REAL Aeff
	REAL dP
	REAL kv
	REAL opening,per_opening
	REAL y "Actual valve stem displacement (cannot be negative)"
	REAL yr "Calculated displacement of the valve stem - higher yr means valve is more open"
	DISCR REAL y0
	INTEGER ier, jx, jy
INIT 
	y0 = (D_orif-D_rod)/(D_pin-D_rod)*H_Pin
CONTINUOUS
	dP = f_in.P - f_out.P
	yr = (f_in.P*1e5 - f_out.P*1e5 - offset) / Gfac
	y = max(0,yr)
	Amax = MATH.PI/4 * (D_orif**2 - D_rod**2)
	Ahor = MATH.PI/4*D_orif**2 - MATH.PI*(D_rod/2 + max(0,(y0-y)) * (D_pin-D_rod)/(2*H_Pin))**2
	Aeff = max(Amin,min(Amax,Ahor))
	kv = Cv*Aeff
	m = kv*regRoot2(dP,0.01) * sqrt(f_in.rho*1e5)
	opening = Aeff/Amax
	per_opening = opening*100
END COMPONENT



COMPONENT CheckValve IS_A dP (
	ENUM FlowEquation eqnChoice = Swagelok
)
"On-Off Valve. Default flow direction is In to Out."
PORTS
	IN bool_signal(n=1) K "on/off control signal"
	OUT analog_signal(n=1) K_feedback "Position feedback"
DATA
	REAL Cv = 1 "Flow coefficient"
	REAL Cd = 0.5 "Discharge coefficient"
	REAL D = 0.002 UNITS u_m "Orifice diameter"
	REAL crack = 0 UNITS u_bar "Cracking pressure"
	REAL tau = 1 UNITS u_s "Time constant of the valve"
	REAL phi_min = 0.01 "Minimum opening - 0 to 1"
	REAL phi_init = 1.0 "Initial opening - 0 to 1"
DECLS
	REAL dP UNITS u_bar
	REAL Pin UNITS u_bar
	REAL Pout UNITS u_bar
	REAL hin
	REAL G "Specific gravity"
	REAL phi "Actual opening"
	REAL phi_target "Targeted opening"
	REAL rho
	REAL rhoRef "Density of refrigerant at reference conditions"
	REAL Tin
	REAL x "Vapour quality"
	REAL vDot
	REAL mDot_dummy, vDot_dummy
	CONST REAL N1 = 14.42 "Swagelok CV constant"
	CONST REAL N2 = 6950
	REAL feedback "Actual valve position"
	INTEGER ier,jx,jy
INIT
	phi = phi_init
DISCRETE
	ASSERT(phi_min>=0 AND phi_min<=1) FATAL "Minimum opening must be between 0 and 1"
	ASSERT(phi_init>=0 AND phi_init<=1) FATAL "Initial opening must be between 0 and 1"
CONTINUOUS
	feedback = phi*100
	phi_target := ZONE (K.signal[1]) 1.0 OTHERS phi_min // K!=0 is On mode, during which valve is fully open
	phi' := (phi_target - phi) / tau
	dP := f_in.P - f_out.P
	rho := donor_cell(dP,f_in.rho,f_out.rho)
	Pin := donor_cell(dP,f_in.P,f_out.P)
	hin := donor_cell(dP,f_in.h,f_out.h)
	Tin := CRYO_PF_prop_vs_ph(f_in.fluid,Pin,hin,fprop_temperature,ier,jx,jy)
	Pout := donor_cell(dP,f_out.P,f_in.P)
	x := CRYO_PF_prop_vs_ph(f_in.fluid,Pin,hin,fprop_quality,ier,jx,jy)
	rhoRef := CRYO_PF_prop_vs_pT(f_in.fluid,1.01325,293.15,fprop_density,ier,jx,jy)
	G := rho/998 // OTHERS rhoRef/1.205
 	IF (eqnChoice==Swagelok) INSERT
		getSwagelokFlow(f_in.fluid,Cv,f_in.P,f_in.h,f_out.P,f_out.h,mDot_dummy,vDot_dummy)
		m = phi * mDot_dummy
		vDot = phi * vDot_dummy
	ELSEIF (eqnChoice==DischargeCoefficient) INSERT
		m = phi * Cv*D**2 * regRoot2(dP,dP_small) * sqrt(rho*1e5)
		vDot = m/rho
		mDot_dummy = m
		vDot_dummy = vDot
	ELSE
		m = phi * m0/sqrt(dP0) * regRoot2(max(0,dP-crack),dP_small)
		vDot = m/rho
		mDot_dummy = m
		vDot_dummy = vDot
	END IF
	K_feedback.signal[1] = feedback
END COMPONENT

COMPONENT FloatingBallValve IS_A dP
"Cryo-Floating-Ball Valcve Default flow direction is In to Out."
PORTS
	IN bool_signal(n=1) K "on/off control signal"
	IN bool_signal (n=1) Dir "Orientation WRT input output ports. True => Leak hole faces the inlet, False means otherwise."
	OUT analog_signal(n=1) K_feedback "Position feedback"
DATA
	REAL Cv = 1 "Flow coefficient"
	REAL Cvl = 0.0005 "Flow coefficient"
	REAL tau = 1 UNITS u_s "Time constant of the valve"
	REAL phi_min = 0.01 "Minimum opening - 0 to 1"
	REAL phi_init = 1.0 "Initial opening - 0 to 1"
	REAL dP_leak = 3 UNITS u_bar "dP above which leak"
DECLS
	REAL dP UNITS u_bar
	REAL Cv_act
	REAL phi_act
	REAL Pin UNITS u_bar
	REAL Pout UNITS u_bar
	REAL hin
	REAL G "Specific gravity"
	REAL phi "Actual opening"
	REAL phi_l "Leak opening"
	REAL rho
	REAL rhoRef "Density of refrigerant at reference conditions"
	REAL Tin
	REAL x "Vapour quality"
	REAL vDot
	REAL mDot_dummy, vDot_dummy
	CONST REAL N1 = 14.42 "Swagelok CV constant"
	CONST REAL N2 = 6950
	REAL feedback "Actual valve position"
	INTEGER ier,jx,jy 
	
INIT
	
	
DISCRETE
	ASSERT(phi_min>=0 AND phi_min<=1) FATAL "Minimum opening must be between 0 and 1"
	ASSERT(phi_init>=0 AND phi_init<=1) FATAL "Initial opening must be between 0 and 1"
	
	
CONTINUOUS
	feedback = phi*100 //Express feedback as a percentage
	phi := 1//ZONE (K.signal[1]) 1.0 OTHERS 0 // K!=0 is On mode, during which valve is fully open
	//phi' := (phi - phi) / tau //
	phi_l := 0.05//ZONE (NOT K.signal[1]) 0.0005 OTHERS 0
	dP := f_in.P - f_out.P
	rho := donor_cell(dP,f_in.rho,f_out.rho)
	Pin := donor_cell(dP,f_in.P,f_out.P)
	hin := donor_cell(dP,f_in.h,f_out.h)
	Tin := CRYO_PF_prop_vs_ph(f_in.fluid,Pin,hin,fprop_temperature,ier,jx,jy)
	Pout := donor_cell(dP,f_out.P,f_in.P)
	x := CRYO_PF_prop_vs_ph(f_in.fluid,Pin,hin,fprop_quality,ier,jx,jy)
	rhoRef := CRYO_PF_prop_vs_pT(f_in.fluid,1.01325,293.15,fprop_density,ier,jx,jy)
	G := rho/998 // OTHERS rhoRef/1.205
	Cv_act = Cv
	/*Cv_act := ZONE (K.signal[1]) Cv\
				 ZONE (NOT K.signal[1] AND dP>=0 AND NOT Dir.signal[1]) Cvl\
				 ZONE (NOT K.signal[1] AND dP<=0 AND Dir.signal[1]) Cvl\
				 OTHERS 0*/
	
 	phi_act := ZONE (NOT K.signal[1] AND dP>=dP_leak AND NOT Dir.signal[1]) phi_l\
				 ZONE (NOT K.signal[1] AND dP<=-dP_leak AND Dir.signal[1]) phi_l\
				 OTHERS phi
	
	getSwagelokFlow(f_in.fluid,Cv_act,f_in.P,f_in.h,f_out.P,f_out.h,mDot_dummy,vDot_dummy)
<:m>	m = phi_act * mDot_dummy
	vDot = phi_act * vDot_dummy
	
	/*
	EXPAND_BLOCK(K.signal[1] AND LStatus==1)
				getSwagelokFlow(f_in.fluid,Cvl,f_in.P,f_in.h,f_out.P,f_out.h,mDot_dummy,vDot_dummy)
				m = phi_l * mDot_dummy
				vDot = phi_l * vDot_dummy
	END EXPAND_BLOCK
	
	EXPAND_BLOCK((NOT k) AND (LStatus == 0))
				m = 0.00005
				vDot = m/rho
				vDot_dummy = m
				vDot_dummy = vDot
	END EXPAND_BLOCK
	*/
	
	/*
	ELSEIF (K.signal[1] == FALSE AND LeakStatus == 1.) INSERT
			getSwagelokFlow(f_in.fluid,Cvl,f_in.P,f_in.h,f_out.P,f_out.h,mDot_dummy,vDot_dummy)
			m = phi * mDot_dummy
			vDot = phi * vDot_dummy
		ELSE 
			m = phi * Cvl*D**2 * regRoot2(dP,dP_small) * sqrt(rho*1e5)
			vDot = m/rho
			mDot_dummy = m
			vDot_dummy = vDot
		END IF
	ELSEIF (K.signal[1] == FALSE AND LeakStatus == 0.) INSERT
		m = 0.0001
	*/
	K_feedback.signal[1] = feedback	
END COMPONENT

COMPONENT ReversingValve
"Reverse flow between two input and two output streams"
PORTS
	IN fluid f_in_high // from compressor discharge
	OUT fluid f_out_high // to outdoor unit
	IN fluid f_in_low // from indoor unit
	OUT fluid f_out_low // to compressor suction
	IN bool_signal(n=1) K "1 for normal operation, 0 for reversed position"
DATA
	ENUM InitialConditions init = Ph
	REAL P0_l = 20 UNITS u_bar
	REAL h0_l = 150000 UNITS u_J_kg
	REAL T0_l = 255 UNITS u_K
	REAL x0_l = 0.5
	REAL P0_h = 20 UNITS u_bar
	REAL h0_h = 150000 UNITS u_J_kg
	REAL T0_h = 255 UNITS u_K
	REAL x0_h = 0.5
	REAL m0 = 0.01 UNITS u_kg_s "Nominal (design) mass flow rate"
	REAL dP0 = 0.01 UNITS u_bar "Nominal (design) pressure drop"
	REAL dP_small = 0.0001 UNITS u_bar "dP below which to apply low mass flow rate smoothing"
	REAL phi_init = 1.0 "Initial opening of valves that ARE open"
	REAL phi_min = 0.01 "Minimum valve opening"
	REAL tau = 1.0 UNITS u_s "Time constant for open/close toggling"
TOPOLOGY
	AdiabVol_Ph HighSideCV(D=0.01,L=0.01,init=init,P0=P0_h,h0=h0_h,T0=T0_h,x0=x0_h)
	AdiabVol_Ph LowSideCV(D=0.01,L=0.01,init=init,P0=P0_l,h0=h0_l,T0=T0_l,x0=x0_l)
	CheckValve (eqnChoice=Quadratic) HighSideValve(m0=m0,dP0=dP0,dP_small=dP_small,phi_init=1,phi_min=phi_min,tau=tau)
	CheckValve (eqnChoice=Quadratic) LowSideValve(m0=m0,dP0=dP0,dP_small=dP_small,phi_init=1,phi_min=phi_min,tau=tau)
	CheckValve (eqnChoice=Quadratic) LowSideValve1(m0=m0,dP0=dP0,dP_small=dP_small,phi_init=phi_min,phi_min=phi_min,tau=tau)
	CheckValve (eqnChoice=Quadratic) LowSideValve2(m0=m0,dP0=dP0,dP_small=dP_small,phi_init=phi_min,phi_min=phi_min,tau=tau)
	Gate_NOT reverse

	CONNECT HighSideCV.f_in TO f_in_high
	CONNECT HighSideCV.f_out TO HighSideValve.f_in
	CONNECT HighSideCV.f_out TO LowSideValve1.f_in
	CONNECT LowSideCV.f_in TO f_in_low
	CONNECT LowSideCV.f_out TO LowSideValve.f_in	
	CONNECT LowSideValve1.f_out TO LowSideValve.f_in
	CONNECT LowSideValve.f_out TO f_out_low
	CONNECT HighSideValve.f_out TO f_out_high
	CONNECT LowSideValve2.f_out TO f_out_low
	CONNECT HighSideValve.f_out TO LowSideValve2.f_in

	CONNECT K TO HighSideValve.K
	CONNECT K TO LowSideValve.K
	CONNECT K TO reverse.s_in
	CONNECT reverse.s_out TO LowSideValve1.K
	CONNECT reverse.s_out TO LowSideValve2.K	
INIT
	f_out_high.is_C = FALSE
	f_out_low.is_C = FALSE
END COMPONENT



COMPONENT ThreeWayValve (
	ENUM FlowEquation eqnChoice = Swagelok
)
PORTS
	IN fluid f_in
	OUT fluid f_out1 
	OUT fluid f_out2
	IN bool_signal(n=1) K "If TRUE, flow to f_out1"
DATA
	// Design parameters
	REAL Cv = 90
	REAL Cd = 15
	REAL D = 0.002 UNITS u_m "Orifice Diameter"
	REAL m0 = 0.01 UNITS u_kg_s "Nominal flow rate for Quadratic relation"
	REAL dP0 = 0.01 UNITS u_bar "Nominal pressure drop for Quadratic relation"
	REAL dP_small = 0.001 UNITS u_bar "Region at which Low Mass Flow Rate smoothing to apply"
	REAL phi_init = 1.0 "Initial position between 0 and 1. 0 means flow to f_out2"
	REAL phi_min = 0.01
	REAL tau = 1.0 UNITS u_s "Time constant for open/close toggling"
	ENUM InitialConditions init = Ph
	REAL P0 = 20 UNITS u_bar
	REAL h0 = 150000 UNITS u_J_kg
	REAL T0 = 255 UNITS u_K
	REAL x0 = 0.5
TOPOLOGY
	AdiabVol_Ph CV(D=0.04,L=0.05,init=init,P0=P0,h0=h0,T0=T0,x0=x0) // hardcoded geometry
	CheckValve (eqnChoice=eqnChoice) Valve1(Cv=Cv,Cd=Cd,D=D,m0=m0,dP0=dP0,dP_small=dP_small,phi_init=1,phi_min=phi_min,tau=tau)
	CheckValve (eqnChoice=eqnChoice) Valve2(Cv=Cv,Cd=Cd,D=D,m0=m0,dP0=dP0,dP_small=dP_small,phi_init=1,phi_min=phi_min,tau=tau)
	Gate_NOT reverse

	CONNECT f_in TO CV.f_in
	CONNECT CV.f_out TO Valve1.f_in
	CONNECT CV.f_out TO Valve2.f_in
	CONNECT Valve1.f_out TO f_out1
	CONNECT Valve2.f_out TO f_out2
	CONNECT K TO Valve1.K
	CONNECT K TO reverse.s_in
	CONNECT reverse.s_out TO Valve2.K
INIT
	f_out1.is_C = FALSE
	f_out2.is_C = FALSE
END COMPONENT



COMPONENT NonReturnValve IS_A AbsJunction 
"One way valve / No reverse flow"
DATA
	REAL Cv = 0.1 UNITS no_units "Flow Coefficient, or Discharge coefficient if not Swagelok"
	REAL Pcrack = 1 UNITS u_bar
DECLS
	REAL dP UNITS u_bar "Pressure drop across valve"
	REAL hin UNITS u_J_kg "Inlet specific enthalpy"
	REAL Pin UNITS u_bar "Inlet pressure"
	REAL Pout UNITS u_bar "Outlet pressures"
	REAL vDot "Volumetric flow rate (L/min)"
	REAL m_dummy
	REAL x "Vapour quality"
	INTEGER ier,xx,yy "Error codes"
CONTINUOUS
	dP := f_in.P - f_out.P
	Pin = donor_cell(dP,f_in.P,f_out.P)
	Pout = donor_cell(dP,f_out.P,f_in.P)
	hin = donor_cell(dP,f_in.hf,f_in.hb)
	x = CRYO_PF_prop_vs_ph(f_in.fluid,Pin,hin,fprop_quality,ier,xx,yy)
	getSwagelokFlow(f_in.fluid,Cv,f_in.P,f_in.h,f_out.P+Pcrack,f_out.h,m_dummy,vDot)
	m = ZONE (dP-Pcrack>=0) m_dummy OTHERS 0
END COMPONENT

