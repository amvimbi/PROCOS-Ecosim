/*-----------------------------------------------------------------------------------------
 LIBRARY: CO2
 FILE: PUMPS
 AUTHOR: Viren Bhanot
 COMPANY: CERN
 DESCRIPTION: Pumps and Compressors
 CREATION DATE: 01/11/2016
-----------------------------------------------------------------------------------------*/
USE CRYOLIB
USE MATH
USE PORTS_LIB
USE THERMAL



--------------------------------------------------------------------------------
// COMPRESSOR MODELS
--------------------------------------------------------------------------------
ABSTRACT COMPONENT AbsCompressor
"Base class for compressor and Pump models"
// NOTE: Cryolib also has AbsCompressor
PORTS
	IN fluid f_in
	OUT fluid f_out
DECLS
	REAL m UNITS u_kg_s "Refrigerant mass flow rate"
	REAL rho UNITS u_kg_m3 "Inlet density"
	REAL s UNITS u_J_K "Outlet entropy"
	REAL hs UNITS u_J_kg "Isentropic outlet enthalpy"
	PRIVATE INTEGER ier,ipx,ipy "Error codes"
TOPOLOGY
	PATH f_in TO f_out
INIT
	f_out.is_C = FALSE
CONTINUOUS
	rho = CRYO_PF_prop_vs_ph(f_in.fluid,f_in.P,f_in.hf,fprop_density,ier,ipx,ipy)
	s = CRYO_PF_prop_vs_ph(f_in.fluid,f_in.P,f_in.hf,fprop_entropy,ier,ipx,ipy)
	hs = CRYO_PF_prop_vs_ps(f_in.fluid,f_out.P,s,fprop_enthalpy,ier,ipx,ipy)

	f_in.fluid = f_out.fluid
	f_in.n_fluid = f_out.n_fluid
	f_in.mf = f_out.mf
	f_in.mb = f_out.mb
	f_in.mb = 0 // backwards flow impossible
	f_in.mf = IF (m<0) 0 ELSE m // backwards flow impossible
<hb>	f_in.hb = f_out.hb
END COMPONENT



COMPONENT Compressor IS_A CO2.AbsCompressor
"Simple compressor model, assumes constant efficiencies"
// RPM itself is an input parameter rather than a percentage
PORTS
	IN analog_signal(n=1) RPM "Rotational Speed"
	OUT analog_signal(n=1) Q_loss
DATA
	REAL disp = 27.42E-6 UNITS u_m3 "Displacement volume"
	REAL etaVol = 0.94 "Volumetric efficiency"
	REAL etaIse = 0.69 "Isentropic efficiency"
	REAL etaMot = 0.90 "Motor efficiency"
DECLS
	REAL power UNITS u_W "Power consumption"
	REAL motorLoss UNITS u_W "Motor power losses, dissipated as heat"
	REAL T UNITS u_K
CONTINUOUS
	T = CRYO_PF_prop_vs_ph(f_in.fluid,f_in.P,f_in.hf,fprop_temperature,ier,ipx,ipy)
<hf> 	f_out.hf = f_in.hf + (hs - f_in.hf) / etaIse
<m>	m = max(1e-4,(RPM.signal[1]/60)*disp*etaVol*rho)
<W>	power = m * (f_out.hf - f_in.hf) / etaMot
<loss>	motorLoss = (1 - etaMot) * power
	Q_loss.signal[1] = motorLoss
END COMPONENT



COMPONENT Compressor_vfd IS_A CO2.AbsCompressor
"Compressor model with Variable Frequency Drive"
// 'frequencies' here refer to rotational speed of crankshaft
// The data required for this model is available in Bitzer's SSPG8 software
PORTS
	IN bool_signal(n=1) start "Start signal"
	IN analog_signal(n=1) speed "Compressor speed signal in % (0-100)"
	OUT analog_signal(n=1) Q_loss "Power loss in the form of heat"
DATA
	INTEGER nCyl = 4 "No. of cylinders"
	REAL bore = 40 UNITS u_mm "Bore diameter"
	REAL stroke = 39.3 UNITS u_mm "Stroke length"
	REAL MaxFrequency = 82 UNITS u_Hz "PLC imposed maximum frequency allowed"
	TABLE_1D FrequencyRange = {{0,100},{25,87}} "Frequencies at 0 and 100% input"
	TABLE_1D MotorRPM = {{50,60},{1450,1750}} "Actual rotational speeds of motor at given frequencies"
	REAL etaVol = 0.95 "Volumetric efficiency"
	REAL etaIse = 0.75 "Isentropic efficiency"
	REAL etaMot = 0.90 "Motor efficiency"
DECLS
	DISCR REAL disp UNITS u_m3 "Displacement volume per revolution of the crank"
	REAL RPM UNITS u_rpm "Revolution speed of compressor crank"
	REAL nu UNITS u_Hz "Frequency at the inlet signal"
	REAL power UNITS u_W "Power consumption"
	REAL motorLoss UNITS u_W "Motor power losses, dissipated as heat"
	REAL T UNITS u_K "Refrigerant temperature"
INIT
	disp = nCyl * MATH.PI * (bore*1E-3/2)**2 * (stroke*1E-3) // displacement volume per revolution
CONTINUOUS
	nu = min(MaxFrequency,linearInterp1D(FrequencyRange,speed.signal[1])) // cannot be more than user specified maximum
	RPM = linearInterp1D(MotorRPM,nu)
	T = CRYO_PF_prop_vs_ph(f_in.fluid,f_in.P,f_in.hf,fprop_temperature,ier,ipx,ipy)
<m>	m = ZONE(start.signal[1]) speed.signal[1]/100 * (RPM/60) * disp * etaVol * rho OTHERS 0
<hf>	f_out.hf = f_in.hf + (hs - f_in.hf) / etaIse
<W>	power = m * (f_out.hf - f_in.hf) / etaMot
<loss>	motorLoss = (1 - etaMot) * power
	Q_loss.signal[1] = motorLoss
END COMPONENT



COMPONENT CompressorAHRI IS_A CO2.AbsCompressor (
	BOOLEAN SuperheatCorrection = FALSE "If true, uses the Dabiri-Rice method to correct for superheat"
)
"AHRI 10-coefficient polynomial compressor model"
DATA
	REAL cm[10] "Coefficients for mass flow rate"
	REAL cW[10] "Coefficients for power consumption"
	REAL Tsup = 5 "Superheat for which map calculated, either °C or °F"
DECLS
	REAL power UNITS u_W "Compressor power consumption"
	REAL Td UNITS u_K "Saturated discharge dewpoint temperature"
	REAL Tin UNITS u_K "Actual suction temperature, using map superheat Tsup"
	REAL Ts UNITS u_K "Saturated suction dewpoint temperature"
	REAL vMap UNITS u_m3_kg "Specific volume of suction fluid"
	REAL vActual UNITS u_m3_kg "Actual specific volume at saturation temperature"
	
CONTINUOUS
	Ts = CRYO_PF_prop_vs_Px(f_in.fluid,f_in.P,1,fprop_temperature,ier,ipx,ipy)
	Td = CRYO_PF_prop_vs_Px(f_in.fluid,f_out.P,1,fprop_temperature,ier,ipx,ipy)
	Tin = Ts + Tsup
<m>	m = cm[1] + cm[2]*Ts + cm[3]*Td + cm[4]*Ts**2 + cm[5]*Ts*Td + cm[6]*Td**2 + cm[7]*Ts**3 + cm[8]*Ts**2*Td + cm[9]*Ts*Td**2 + cm[10]*Td**3
<W>	power = cW[1] + cW[2]*Ts + cW[3]*Td + cW[4]*Ts**2 + cW[5]*Ts*Td + cW[6]*Td**2 + cW[7]*Ts**3 + cW[8]*Ts**2*Td + cW[9]*Ts*Td**2 + cW[10]*Td**3
END COMPONENT



COMPONENT ScrollCompressor (
	BOOLEAN ConstantRefHTC = TRUE,
	BOOLEAN ConstantAirHTC = TRUE
)
"Inlet must connect to resistive component. Outlet to capacitive"
PORTS
	IN fluid f_in
	OUT fluid f_out
	IN analog_signal(n=1) T_amb
	IN analog_signal(n=1) RPM "Rotational Speed"
	OUT analog_signal(n=1) Q_suc
	OUT analog_signal(n=1) Q_dis
DATA
	ENUM Material mat = SS_304 "Compressor shell material"
	REAL Do_suc = 0.1396 UNITS u_m "Suction chamber shell outer diameter"
	REAL Th_suc = 0.005 UNITS u_m "Suction chamber shell thickness"
	REAL L_suc = 0.5 UNITS u_m "Suction chamber height"
	REAL D_motor = 0.08 UNITS u_m "Motor diameter - motor is inside suction chamber, and as tall as the suction chamber itself"
	REAL Do_dis = 0.1396 UNITS u_m "Discharge chamber shell outer diameter"
	REAL Th_dis = 0.005 UNITS u_m "Discharge chamber shell thickness"
	REAL L_dis = 0.06 UNITS u_m "Discharge chamber height"
	REAL D_i[1] = 0.01 UNITS u_m "SUction chamber inlet port height"
	REAL H_i[1] = 0.4 UNITS u_m "Suction chamber inlet port height"
	REAL D_o[1] = 0.01 UNITS u_m "Suction chamber outlet port height"
	REAL H_o[1] = 0.4 UNITS u_m "Suction chamber outlet port height"

	REAL disp = 27.42E-6 UNITS u_m3 "Displacement volume"
	REAL etaVol = 0.94 "Volumetric efficiency"
	REAL etaIse = 0.69 "Isentropic efficiency"
	REAL etaMot = 0.90 "Motor efficiency"

	REAL alpha_a_suc = 10 UNITS u_W_m2K "Suction chamber natural-convective HTC"
	REAL alpha_a_dis = 10 UNITS u_W_m2K "Discharge chamber natural-convective HTC"

	REAL alpha_l_dis = 150 UNITS u_W_m2K "Discharge chamber - Liquid heat transfer coefficient"
	REAL alpha_tp_dis = 200 UNITS u_W_m2K "Discharge chamber - Two phase heat transfer coefficient"
	REAL alpha_g_dis = 150 UNITS u_W_m2K "Discharge chamber - Vapour heat transfer coefficient"
	REAL alpha_suc = 200 UNITS u_W_m2K "Suction chamber heat transfer coefficient"
	REAL m_steady = 0.05 UNITS u_kg_s "Steady state refrigerant mass flow rate at which HTCs were calculated"
	REAL m0 = 0.01 UNITS u_kg_s "Nominal mass flow rate"
	REAL dP0 = 0.02 UNITS u_bar "Nominal pressure drop"

	ENUM InitialConditionsAccu init_suc = P_h "Initial conditions specification for suction chamber"
	REAL P0_suc = 17 UNITS u_bar "Initial suction pressure"
	REAL h0_suc = 300000 UNITS u_J_kg "Initial suction enthalpy"
	REAL T0_suc = 300 UNITS u_K "Initial suction temperature"
	ENUM InitialConditions init_dis = Ph "Initial conditions specification for discharge chamber"
	REAL P0_dis = 17 UNITS u_bar "Initial discharge pressure"
	REAL h0_dis = 300000 UNITS u_J_kg "Initial discharge enthalpy"
	REAL T0_dis = 300 UNITS u_K "Initial discharge temperature"
	REAL Tw0_suc = 300 UNITS u_K "Initial suction chamber shell temperature"
	REAL Tw0_dis = 300 UNITS u_K "Initial discharge chamber shell temperature"
DECLS
	DISCR REAL Mw UNITS u_kg "Metal weight"
	CONST REAL D_ii = 0.01
TOPOLOGY
	CO2.Compressor Compressor (
		disp = disp,
		etaVol = etaVol,
		etaIse = etaIse,
		etaMot = etaMot
	)
	SuctionChamberVolume SuctionChamber (
		D = Do_suc-2*Th_suc,
		Di = D_motor,
		L = L_suc,
		D_in = D_i,
		z_in = H_i,
		D_out = D_o,
		z_out = H_o,
		alpha_ref = alpha_suc,
		init = init_suc,
		P0 = P0_suc,
		h0 = h0_suc,
		T0 = T0_suc
	)
	LumpedShell (ConstantAirHTC=ConstantAirHTC,accountForRadiation=FALSE) SuctionChamberWall (
		mat = mat,
		Do = Do_suc,
		Th = Th_suc,
		L = L_suc,
		alpha_a0 = alpha_a_suc,
		Tw0 = Tw0_suc
	)
	CartridgeHeater(n=1) Motor (
		D = D_motor,
		mat = mat,
		L = L_suc,
		Tw0_i = Tw0_suc,
		Tw0_o = Tw0_suc		
	)	
	CV_Ph DischargeChamber (
		D = Do_dis-2*Th_dis,
		L = L_dis,
		alpha_l = alpha_l_dis,
		alpha_tp = alpha_tp_dis,
		alpha_g = alpha_g_dis,
		m_steady = m_steady,
		init = init_dis,
		P0 = P0_dis,
		h0 = h0_dis,
		T0 = T0_dis
	)
	LumpedShell (ConstantAirHTC=ConstantAirHTC,accountForRadiation=TRUE) DischargeChamberWall (
		mat = mat,
		Do = Do_dis,
		Th = Th_dis,
		L = L_dis,
		alpha_a0 = alpha_a_dis,
		Tw0 = Tw0_dis
	)
	dP dp (m0=m0, dP0=dP0) // small pressure drop element at compressor outlet
	
	CONNECT f_in TO SuctionChamber.f_in[1]
	CONNECT SuctionChamber.f_out[1] TO Compressor
	CONNECT Compressor TO DischargeChamber
	CONNECT DischargeChamber TO dp
	CONNECT dp TO f_out
	CONNECT SuctionChamber.tp_out TO SuctionChamberWall.tp_in
	CONNECT Motor.tp_out TO SuctionChamber.tp_in
	CONNECT DischargeChamber.tp_out TO DischargeChamberWall.tp_in
	CONNECT Compressor.Q_loss TO Motor.Q_in
	CONNECT RPM TO Compressor.RPM
	CONNECT T_amb TO SuctionChamberWall.T_amb
	CONNECT T_amb TO DischargeChamberWall.T_amb
INIT
	Mw = SuctionChamberWall.Mw + DischargeChamberWall.Mw + Motor.Mw
DISCRETE
	ASSERT (init_dis!=Px) FATAL "Px initial conditions not supported for Compressor model"
	ASSERT (init_suc!=P_x AND init_suc!=P_level) FATAL "Two phase initial conditions not supported for Compressor model"
CONTINUOUS	
	Q_suc.signal[1] = SuctionChamberWall.Qa_total
	Q_dis.signal[1] = DischargeChamberWall.Qa_total
END COMPONENT



COMPONENT EconCompressorBase
"Base class for compressor with vapour-injection port"
PORTS
	IN fluid f_suc
	IN fluid f_econ
	OUT fluid f_dis
	IN analog_signal(n=1) RPM
	OUT analog_signal(n=1) Q_loss
DATA
	REAL disp = 29.5E-6 UNITS u_m3 "Displacement volume"
	REAL RPM_nom = 3500 UNITS u_rpm "Nominal (design) revolution speed"
	REAL etaMech = 0.90 "Mechanical efficiency"
	REAL etaMot = 0.90 "Motor efficiency"
	REAL y = 1.35 "Polytropic index (=Cp/Cv)"
	REAL z = 1.21 "First stage volume ratio"
	REAL mc[2] = {-0.1396,1.9894} "Mass flow curve fit coefficients"
	REAL pc[4] = {1.9,0.61,1.9,0.61} "Power curve fit coefficients"
DECLS
	REAL P_suc UNITS u_bar
	REAL P_econ UNITS u_bar
	REAL P_dis UNITS u_bar
	REAL P_int UNITS u_bar
	REAL h_suc UNITS u_J_kg
	REAL h_mix UNITS u_J_kg
	REAL h_int UNITS u_J_kg
	REAL h_dis UNITS u_J_kg
	REAL m_suc UNITS u_kg_s
	REAL m_econ UNITS u_kg_s
	REAL m_dis UNITS u_kg_s
	REAL s_suc
	REAL rho_suc UNITS u_kg_m3
	REAL rho_econ UNITS u_kg_m3
	REAL rho_dis UNITS u_kg_m3
	REAL rho_int UNITS u_kg_m3
	REAL rho_mix UNITS u_kg_m3
	REAL vDot_mix UNITS u_m3_s
	REAL vDot_suc UNITS u_m3_s
	REAL power UNITS u_W
	REAL motorLoss UNITS u_W
	
	CONST REAL Cv_econ = 1.5e-6
	CONST REAL Cv_econ_backflow = 5e-6
	REAL m_econ1, m_econ2
	REAL mh_suc, mh_econ // dummy
	INTEGER ier,iex,iey // error codes
	INTEGER tmp1,tmp2,tmp3 // dummy
INIT
	// all ports are resistive (aka must connect to a control volume)
	f_suc.is_C = FALSE
	f_econ.is_C = FALSE
	f_dis.is_C = FALSE
CONTINUOUS
	// Suction: gets sucked into the suction chamber based on polynomial flow rate
	P_suc := f_suc.P
	h_suc := f_suc.h
	rho_suc := CRYO_PF_prop_vs_ph(f_suc.fluid,P_suc,h_suc,fprop_density,ier,iex,iey)
	s_suc := CRYO_PF_prop_vs_ph(f_suc.fluid,P_suc,h_suc,fprop_entropy,ier,iex,iey)
	m_suc := 1.07*(mc[1]*P_dis/P_suc + mc[2]) * rho_suc * RPM.signal[1]/RPM_nom * 1e-3 // 1e-3 to convert to kg/s. Curve fit coeff were for g/s
	//m_suc := (RPM.signal[1]/60)*disp*0.9*rho_suc
	
	vDot_suc := m_suc/rho_suc
	f_suc.mf = IF (m_suc>0) m_suc ELSE 0
	f_suc.mb = 0
	mh_suc := m_suc*h_suc
	
	// Intermediate: enthalpy due to polytropic compression of suction refrigerant
	P_int := P_suc * z**y
	h_int := CRYO_PF_prop_vs_ps(f_suc.fluid,P_int,s_suc,fprop_enthalpy,ier,iex,iey)
	rho_int := CRYO_PF_prop_vs_ph(f_suc.fluid,P_int,h_int,fprop_density,ier,iex,iey)
	
	// Economizer Port: reverse flow *can* occur here
	P_econ := f_econ.P
	f_econ.hb = h_int // in case of reverse flow, intermediate enthalpy flows into econ line
	rho_econ := CRYO_PF_prop_vs_ph(f_suc.fluid,P_econ,f_econ.hf,fprop_density,ier,iex,iey)
	
	//m_econ := 1.5e-6 * sqrt(donor_cell(P_econ-P_int,rho_econ,rho_int)*1e5) * regRoot2(P_econ-P_int) // hard-coded flow coefficient
	m_econ1 := Cv_econ * sqrt(rho_econ*1e5) * regRoot2(P_econ-P_int)
	m_econ2 := Cv_econ_backflow * sqrt(rho_int*1e5) * regRoot2(P_econ-P_int)
	m_econ := spliceFunction(m_econ1,m_econ2,P_econ-P_int,0.001)
	//m_econ1 := (-0.802849e-3*(P_suc/P_econ) + 0.597687e-3) * rho_econ
	//m_econ2 := -min(-(-0.802849e-3*(P_suc/P_econ) + 0.597687e-3) * rho_int, m_suc - 1e-7)
	//m_econ := spliceFunction(m_econ1, m_econ2, 0.597687e-3/0.802849e-3 - (P_suc/P_econ), 5e-3)	
	f_econ.mf = IF (m_econ>0) m_econ ELSE 0
	f_econ.mb = IF (m_econ<0) -m_econ ELSE 0
	mh_econ = semiLinear(m_econ,f_econ.hf,h_int)
	
	// Scroll-set mixing of the econ line and compressod vapour from suction chamber
	m_suc*h_int + mh_econ - m_dis*h_mix = 0 // energy balance (calculates h_mix)
	rho_mix := CRYO_PF_prop_vs_ph(f_suc.fluid,P_int,h_mix,fprop_density,ier,iex,iey)
	vDot_mix := m_dis/rho_mix
	
	// Discharge line: witnesses the result of all the aforementioned mixing
	P_dis := f_dis.P
	m_dis := m_suc + m_econ
	m_dis*h_dis = mh_suc + mh_econ + power*etaMot*etaMech // energy balance
	rho_dis := CRYO_PF_prop_vs_ph(f_suc.fluid,P_dis,h_dis,fprop_density,ier,iex,iey)
	f_dis.hf = h_dis
	f_dis.mf = IF (m_dis>0) m_dis ELSE 0
	f_dis.mb = 0
	
	// Power
	power := max(m_suc*(h_int-h_suc), 0.97 * (pc[1]*P_suc*1e5*vDot_suc*((P_int/P_suc)**pc[2]-1) + pc[3]*P_int*1e5*vDot_mix*((P_dis/P_int)**pc[4]-1)))
	motorLoss := power * (1-etaMot*etaMech)
	
	f_suc.fluid = f_dis.fluid
	f_suc.fluid = f_econ.fluid
	f_suc.n_fluid = f_dis.n_fluid
	f_suc.n_fluid = f_econ.n_fluid
	f_suc.hb = h_int
	
	motorLoss = Q_loss.signal[1]
END COMPONENT



COMPONENT EconCompressor (
	BOOLEAN ConstantRefHTC = TRUE "If FALSE uses Cooper Pool Boiling HTC",
	BOOLEAN ConstantAirHTC = TRUE "If FALSE, uses ChurchillChu flat plate HTC",
	BOOLEAN DynamicRefHTC = FALSE "If TRUE adds a low pass filter to HTC calculations"
)
PORTS
	IN fluid f_suc
	IN fluid f_econ
	OUT fluid f_dis
	IN analog_signal(n=1) RPM "Revolution speed"
	IN analog_signal(n=1) T_amb "Ambient air temperature"
DATA
	ENUM Material mat = SS_304 "Compressor shell material"
	REAL D_motor = 0.05 UNITS u_m "Motor diameter"
	REAL Do_suc = 0.1396 UNITS u_m "Suction chamber shell outer diameter"
	REAL Th_suc = 0.005 UNITS u_m "Suction chamber shell thickness"
	REAL L_suc = 0.5 UNITS u_m "Suction chamber height"
	REAL D_i[1] = 0.01 UNITS u_m "Suction chamber inlet port diameter"
	REAL H_i[1] = 0.4 UNITS u_m "Suction chamber inlet port height"
	REAL D_o[1] = 0.01 UNITS u_m "Suction chamber outlet port diameter"
	REAL H_o[1] = 0.4 UNITS u_m "Suction chamber outlet port height"
	REAL Do_dis = 0.1396 UNITS u_m "Discharge chamber shell outer diameter"
	REAL Th_dis = 0.005 UNITS u_m "Discharge chamber shell thickness"
	REAL L_dis = 0.06 UNITS u_m "Discharge chamber height"

	REAL disp = 27.42E-6 UNITS u_m3_s "Displacement"
	REAL RPM_nom = 3500 UNITS u_rpm "Nominal (design) rotational speed"
	REAL etaMech = 0.90 "Mechanical efficiency"
	REAL etaMot = 0.90 "Motor efficiency"
	REAL y = 1.35 "Polytropic index (=cp/cv)"
	REAL z = 1.21 "First stage volume ratio"
	REAL mc[2] = {-0.1396,1.9894} "Mass flow curve fit coefficients"
	REAL pc[4] = {1.9,0.61,1.9,0.61} "Power curve fit coefficients"

	REAL alpha_a_suc = 10 UNITS u_W_m2K "Suction chamber natural-convective HTC"
	REAL alpha_a_dis = 20 UNITS u_W_m2K "Discharge chamber natural-convective HTC"

	REAL alpha_l_dis = 750 UNITS u_W_m2K "Discharge chamber - Liquid heat transfer coefficient"
	REAL alpha_tp_dis = 1500 UNITS u_W_m2K "Discharge chamber - Two phase heat transfer coefficient"
	REAL alpha_g_dis = 500 UNITS u_W_m2K "Discharge chamber - Vapour heat transfer coefficient"
	REAL alpha_suc = 200 UNITS u_W_m2K "Suction chamber heat transfer coefficient"
	REAL m_steady = 0.05 UNITS u_kg_s "Steady state refrigerant mass flow rate at which HTCs were calculated"
	REAL m0 = 0.01 UNITS u_kg_s "Nominal mass flow rate"
	REAL dP0 = 0.02 UNITS u_bar "Nominal pressure drop in discharge chamber"

	ENUM InitialConditionsAccu init_suc = P_h "Initial condition specification for Suction chamber"
	REAL P0_suc = 17 UNITS u_bar "Initial suction pressure"
	REAL h0_suc = 300000 UNITS u_J_kg "Initial suction enthalpy"
	REAL T0_suc = 300 UNITS u_K "Initial suction temperature"
	ENUM InitialConditions init = Ph "Initial condition specification for Econ and Discharge CV"
	REAL P0_dis = 17 UNITS u_bar "Initial discharge pressure"
	REAL h0_dis = 300000 UNITS u_J_kg "Initial discharge enthalpy"
	REAL T0_dis = 300 UNITS u_K "Initial discharge temperature"
	REAL P0_econ = 17 UNITS u_bar "Econ port initial pressure"
	REAL h0_econ = 300000 UNITS u_J_kg "Econ port initial enthalpy"
	REAL T0_econ = 300 UNITS u_K "Econ port initial temperature"
	REAL Tw0_suc = 300 UNITS u_K "Initial suction chamber shell temperature"
	REAL Tw0_dis = 300 UNITS u_K "Initial discharge chamber shell temperature"
TOPOLOGY
	EconCompressorBase EconCompressor (
		disp = disp,
		RPM_nom = RPM_nom,
		etaMech = etaMech,
		etaMot = etaMot,
		y = y,
		z = z,
		mc = mc,
		pc = pc
	)
	SuctionChamberVolume (ConstantRefHTC=ConstantRefHTC,DynamicRefHTC=DynamicRefHTC) SuctionChamber (
		D = Do_suc-2*Th_suc,
		L = L_suc,
		Di = D_motor,
		D_in = D_i,
		z_in = H_i,
		D_out = D_o,
		z_out = H_o,
		alpha_ref = alpha_suc,
		init = init_suc,
		P0 = P0_suc,
		h0 = h0_suc,
		T0 = T0_suc
	)
	AnnularWall_naturalConvection(n=1,ConstantAirHTC=ConstantAirHTC) SuctionChamberWall (
		mat = mat,
		Do = Do_suc,
		Th = Th_suc,
		alpha_a0 = alpha_a_suc,
		Tw0_i = Tw0_suc,
		Tw0_o = Tw0_suc
	)
	CartridgeHeater(n=1) Motor (
		D = D_motor,
		mat = mat,
		L = L_suc,
		Tw0_i = Tw0_suc,
		Tw0_o = Tw0_suc		
	)	
	CV_Ph (ConstantHTC=ConstantRefHTC,DynamicHTC=DynamicRefHTC,voidFractionModel=Homogeneous) DischargeChamber (
		D = Do_dis-2*Th_dis,
		L = L_dis,
		alpha_l = alpha_l_dis,
		alpha_tp = alpha_tp_dis,
		alpha_g = alpha_g_dis,
		m_steady = m_steady,
		init = init,
		P0 = P0_dis,
		h0 = h0_dis,
		T0 = T0_dis
	)
	AnnularWall_naturalConvection(n=1, ConstantAirHTC=ConstantAirHTC, accountForRadiation=TRUE) DischargeChamberWall (
		mat = mat,
		Do = Do_dis,
		Th = Th_dis,
		alpha_a0 = alpha_a_dis,
		Tw0_i = Tw0_dis,
		Tw0_o = Tw0_dis	
	)
	AdiabVol_Ph EconPort (
		D = 0.01,
		L = 0.1, // hard coded
		init = init,
		P0 = P0_econ,
		h0 = h0_econ,
		T0 = T0_econ
	)
	dP dp (
		m0 = m0,
		dP0 = dP0
	)
	CONNECT T_amb TO SuctionChamberWall.T_amb
	CONNECT T_amb TO DischargeChamberWall.T_amb
	CONNECT f_suc TO SuctionChamber.f_in[1]
	CONNECT SuctionChamber.f_out[1] TO EconCompressor.f_suc
	CONNECT f_econ TO EconPort
	CONNECT EconPort TO EconCompressor.f_econ
	CONNECT EconCompressor.f_dis TO DischargeChamber
	CONNECT DischargeChamber TO dp
	CONNECT dp TO f_dis
	CONNECT Motor.tp_out TO SuctionChamber.tp_in
	CONNECT SuctionChamber.tp_out TO SuctionChamberWall.tp_in
	CONNECT EconCompressor.Q_loss TO Motor.Q_in
	CONNECT DischargeChamber.tp_out TO DischargeChamberWall.tp_in
	CONNECT RPM TO EconCompressor.RPM
END COMPONENT






--------------------------------------------------------------------------------
// PUMP MODELS
--------------------------------------------------------------------------------
COMPONENT SimplePump IS_A Compressor
// same as parent component
END COMPONENT



COMPONENT SimplePump_Ctrl IS_A Compressor
PORTS
	IN bool_signal(n=1) start "Start signal"
DATA
	REAL RPM_max = 140 UNITS u_rpm "Maximum RPM"
	REAL minFlow = 0 UNITS u_kg_s "Minimum flow rate through pump"
DECLS
	CLOSE Q_loss
CONTINUOUS
<:m>	m = ZONE(start.signal[1]) max(minFlow,RPM.signal[1]/100 * RPM_max/60 * disp * etaVol * rho) OTHERS minFlow
END COMPONENT



COMPONENT Pump_Lewa IS_A CO2.AbsCompressor (
	INTEGER nHead = 3 "Number of heads"
)
"Lewa diaphragm pump model with variable frequency drive"
// Data for this model can be found in the pump datasheets provided by Lewa
// Default values for LDG3 pump
PORTS
	IN analog_signal(n=1) stroke "Stroke setting in %"
	IN analog_signal(n=1) nuSig "Motor frequency [Hz]"
	IN bool_signal(n=1) on "Start signal"
DATA
	REAL D_plunger = 60 UNITS u_mm "Plunger diameter"
	REAL maxStroke = 60 UNITS u_mm "Maximum stroke length"
	TABLE_1D MotorRPM = {{5,50},{176,1485}} "Frequency vs Motor RPM"
	REAL red = 8.33 "Motor gear reduction ratio 1:red"
	REAL etaVol = 0.99 "Volumetric efficiency"
	REAL etaIse = 0.80 "Isentropic efficiency"
DECLS
	REAL actualRPM UNITS u_rpm
	REAL disp UNITS u_m3 "Displacement volume per stroke"
	REAL strokePerMin "Strokes per minute"
	REAL vDot UNITS u_l_s "Volumetric flow rate"
CONTINUOUS
	actualRPM = linearInterp1D(MotorRPM,nuSig.signal[1]) // number of motor rotations every minute
	strokePerMin = actualRPM/red // number of strokes every minute
	disp = MATH.PI * (D_plunger*1e-3/2)**2 * (max(0,min(stroke.signal[1],100))/100)*(maxStroke*1e-3) // volume displaced per stroke
	m = ZONE (on.signal[1]) strokePerMin/60 * disp * etaVol * rho * nHead OTHERS 0 // mass flow per second
	vDot = m/rho * 1000 // liters per second
<hf> 	f_out.hf = f_in.hf + (hs-f_in.hf) / etaIse
END COMPONENT


COMPONENT Pump_Lewa_pn IS_A CO2.AbsCompressor (
	INTEGER nHead = 3 "Number of heads"
)
"Lewa diaphragm pump model with pump curves"
PORTS
	IN analog_signal(n=1) stroke "Stroke setting in %"
	IN analog_signal(n=1) nuSig "Motor frequency [Hz]"
	IN bool_signal(n=1) on "Start signal"
DATA
	REAL etaIse = 0.9 "Isentropic efficiency"
DECLS
	CONST REAL c0 = -46.604083183870486
	CONST REAL c1 = -0.066853714487060
	CONST REAL c2 = 0.039088368714418
	CONST REAL c3 = 0.365330508675078	
	REAL vDot UNITS u_l_s "Volumetric flow rate"
	REAL dP UNITS u_bar "Pump pressure head"
CONTINUOUS
	dP := f_out.P - f_in.P
	m := ZONE (on.signal[1]) max(0, nHead/3*(c1*dP+c2*rho+c3*stroke.signal[1]+c0) * nuSig.signal[1]*1e-3) OTHERS 0 // [kg/s]
	vDot := m/rho * 1000 // liters per second
<hf> 	f_out.hf := f_in.hf + (hs-f_in.hf) / etaIse
END COMPONENT



COMPONENT Pump_LewaSucn (
	INTEGER nHead = 3 "Number of heads"
)
"Lewa pump with suction chamber"
PORTS
	IN fluid f_in
	OUT fluid f_out
	IN analog_signal(n=1) stroke "Stroke setting in %"
	IN analog_signal(n=1) nuSig "Motor frequency [Hz]"
	IN bool_signal(n=1) on "Start signal"
	IN analog_signal(n=1) T_amb "Ambient temperature in K"
	IN analog_signal(n=1) Q_heat "Heat leak into shell from hydraulics and motor"
DATA
	// Geometry
	REAL Mw = 5 UNITS u_kg "Mass of remote head"
	REAL As = 1 UNITS u_m2 "Outer surface area of remote head"
	REAL D_plunger = 60 UNITS u_mm "Plunger diameter"
	REAL maxStroke = 60 UNITS u_mm "Maximum stroke length"
	TABLE_1D MotorRPM = {{5,50},{176,1485}} "Frequency vs Motor RPM"
	REAL red = 8.33 "Motor gear reduction ratio 1:red"
	REAL etaVol = 0.99 "Volumetric efficiency"
	REAL etaIse = 0.80 "Isentropic efficiency"
	REAL alpha_a = 20 UNITS u_W_m2K "Air-to-shell heat transfer coefficient"
	REAL alpha_l = 500 UNITS u_W_m2K
	REAL alpha_tp = 500 UNITS u_W_m2K
	REAL alpha_g = 500 UNITS u_W_m2K
	ENUM InitialConditions init = Ph "Initial condition specification for Econ and Discharge CV"
	REAL P0 = 20 UNITS u_bar "Initial pressure"
	REAL h0 = 100000 UNITS u_J_kg "Initial enthalpy"
	REAL T0 = 250 UNITS u_K "Initial temperature"
	REAL Tw0 = 293.15 UNITS u_K "Initial wall temperature"
TOPOLOGY
	CV_Ph (ConstantHTC=TRUE,DynamicHTC=FALSE,voidFractionModel=Homogeneous) SuctionChamber (
		D = D_plunger/1000*nHead,
		L = maxStroke/1000/etaVol * nHead,
		alpha_l = alpha_l,
		alpha_tp = alpha_tp,
		alpha_g = alpha_g,
		m_steady = 0.5*nHead,
		init = init,
		P0 = P0,
		h0 = h0,
		T0 = T0
	)
	GenericShell Shell (
		mat = SS_316,
		Mw = Mw,
		As = As,
		alpha_a = alpha_a,
		Tw0 = Tw0
	)
	Pump_Lewa (nHead=nHead) Pump (
		D_plunger = D_plunger,
		maxStroke = maxStroke,
		MotorRPM = MotorRPM,
		red = red,
		etaVol = etaVol,
		etaIse = etaIse
	)
	
	CONNECT f_in TO SuctionChamber.f_in
	CONNECT SuctionChamber.f_out TO Pump.f_in
	CONNECT Pump.f_out TO f_out
	CONNECT stroke TO Pump.stroke
	CONNECT nuSig TO Pump.nuSig
	CONNECT on TO Pump.on
	CONNECT SuctionChamber.tp_out TO Shell.tp_in
	CONNECT T_amb TO Shell.T_amb
	CONNECT Q_heat TO Shell.Q_heat
END COMPONENT