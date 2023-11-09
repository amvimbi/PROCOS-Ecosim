/*-----------------------------------------------------------------------------------------
 LIBRARY: CO2
 FILE: HEATEXCHANGERS
 AUTHOR: Viren Bhanot
 COMPANY: CERN
 DESCRIPTION: Heat Exchanger Models
 CREATION DATE: 21/10/2016
-----------------------------------------------------------------------------------------*/
USE CRYOLIB
USE MATH
USE PORTS_LIB
USE THERMAL


--------------------------------------------------------------------------------
// PIPES WITH CARTRIDGE HEATER INSIDE
--------------------------------------------------------------------------------
COMPONENT DummyLoad (
	INTEGER n = 3 "Number of control volumes",
	BOOLEAN ConstantRefHTC = TRUE "FALSE if using correlations",
	BOOLEAN DynamicHTC = TRUE "If TRUE, applies a low-pass filter on HTC",
	ENUM VoidFractionModel vfModel = Homogeneous "Two-phase flow void fraction model"
)
"Pipe model with a cartridge heater inside"
PORTS
	IN fluid f_in
	OUT fluid f_out
	IN analog_signal(n=1) K "Heater control signal"
DATA
	ENUM PipeMat TubeMaterial = SS_304
	ENUM PipeMat HeaterMaterial = SS_304
	REAL Do = 0.018 UNITS u_m "Pipe outer diameter"
	REAL Th = 0.001 UNITS u_m "Pipe thickness"
	REAL Dh = 0.012 UNITS u_m "Heater diameter"
	REAL L = 1 UNITS u_m "Pipe length"
	REAL z_in = 0 UNITS u_m "Elevation of inlet wrt user-defined base"
	REAL dz = 0 UNITS u_m "Elevation change between inlet and outlet, negative for downward slopes"

	ENUM SinglePhaseHTC choice_1p = UseDittusBoelter //Gnielinski,DittusBoelter,Colburn
	ENUM TwoPhaseHTC choice_2p = UseChen
	REAL tauAlpha = 3 UNITS u_s "Time constant for low-pass HTC filter"
	REAL Q_max = 1000 UNITS u_W "Maximum heater power"
	REAL alpha_l = 1000 UNITS u_W_m2K "Liquid HTC, for Constant HTC"
	REAL alpha_tp = 1000 UNITS u_W_m2K "Two Phase HTC, for Constant HTC"
	REAL alpha_g = 1000 UNITS u_W_m2K "Vapour HTC, for Constant HTC"
	REAL m_steady = 0.05 UNITS u_kg_s "Mass flow rate at which constant HTCs calculated"
	REAL m0 = 0.005 UNITS u_kg_s "Nominal (design) mass flow rate"
	REAL dP0 = 0.01 UNITS u_bar "Nominal (design) pressure drop"

	ENUM InitialConditions init = Ph "Specify initial condition variables"
	REAL P0_i = 20 UNITS u_bar "Initial inlet pressure"
	REAL h0_i = 250000 UNITS u_J_kg "Initial inlet specific enthalpy"	
	REAL T0_i = 255 UNITS u_K "Initial inlet temperature"
	REAL x0_i = 0.5 "Initial inlet vapour quality"
	REAL P0_o = 20 UNITS u_bar "Initial outlet pressure"
	REAL h0_o = 250000 UNITS u_J_kg "Initial outlet specific enthalpy"	
	REAL T0_o = 255 UNITS u_K "Initial outlet temperature"
	REAL x0_o = 0.5 "Initial outlet vapour quality"
	REAL Tw0_i = 245 UNITS u_K "Initial wall temperature - first segment" // Same initial temperatures for heater and wall
	REAL Tw0_o = 245 UNITS u_K "Initial wall temperature - last segment"
TOPOLOGY
	AnnularStream_Ph (
		n = n,
		ConstantRefHTC = ConstantRefHTC,
		DynamicHTC = DynamicHTC,
		vfModel = vfModel
	) 
	Stream (
		Do = Do-2*Th,
		Di = Dh,
		L = L,
		z_in = z_in,
		dz = dz,
		choice_1p = choice_1p,
		choice_2p = choice_2p,
		tauAlpha = tauAlpha,
		alpha_l = alpha_l,
		alpha_tp = alpha_tp,
		alpha_g = alpha_g,
		m_steady = m_steady,
		m0 = m0,
		dP0 = dP0,
		init = init,
		P0_i = P0_i,
		h0_i = h0_i,
		T0_i = T0_i,
		x0_i = x0_i,
		P0_o = P0_o,
		h0_o = h0_o,
		T0_o = T0_o,
		x0_o = x0_o
	)
	CartridgeHeaterCtrl (n=n) Heater (
		mat = HeaterMaterial,
		D = Dh,
		L = L,
		Q_max = Q_max,
		Tw0_i = Tw0_i,
		Tw0_o = Tw0_o
	)
	AnnularWall_onePort (n=n) Wall (
		mat = TubeMaterial,
		Do = Do,
		Th = Th,
		L = L,
		Tw0_i = Tw0_i,
		Tw0_o = Tw0_o
	)
	CONNECT K TO Heater.K
	CONNECT Heater.tp_out TO Stream.tp_in
	CONNECT Stream.tp_out TO Wall.tp_in
	CONNECT f_in TO Stream.f_in
	CONNECT Stream.f_out TO f_out
END COMPONENT



COMPONENT DummyLoadTout IS_A DummyLoad
"Adds ports that output average heater and wall temperatures"
PORTS
	OUT analog_signal(n=1) T_heater
	OUT analog_signal(n=1) T_wall
DATA
	BOOLEAN ToutInCelsius = TRUE "If false, outputs Kelvin"
CONTINUOUS
	T_heater.signal[1] = ZONE(ToutInCelsius) SUM(i IN 1,n; Heater.Tw[i])/n - 273.15 OTHERS SUM(i IN 1,n ; Heater.Tw[i])/n
	T_wall.signal[1] = ZONE(ToutInCelsius) SUM(i IN 1,n; Wall.Tw[i])/n - 273.15 OTHERS SUM(i IN 1,n ; Wall.Tw[i])/n
END COMPONENT


COMPONENT DummyLoad_Tamb (
	INTEGER n = 3 "Number of control volumes",
	BOOLEAN ConstantRefHTC = TRUE "FALSE if using correlations",
	BOOLEAN DynamicHTC = TRUE "If TRUE, applies a low-pass filter on HTC",
	ENUM VoidFractionModel vfModel = Homogeneous "Two-phase flow void fraction model"
)
"Pipe model with a cartridge heater inside"
PORTS
	IN fluid f_in
	OUT fluid f_out
	IN analog_signal(n=1) T_amb "Ambient Temperature (K)"
	IN analog_signal(n=1) K "Heater control signal"
	OUT analog_signal(n=1) T_heater
	OUT analog_signal(n=1) T_wall
DATA
	BOOLEAN ToutInCelsius = TRUE "If false, outputs Kelvin"

	ENUM PipeMat TubeMaterial = SS_304
	ENUM PipeMat HeaterMaterial = SS_304
	REAL Do = 0.018 UNITS u_m "Pipe outer diameter"
	REAL Th = 0.001 UNITS u_m "Pipe thickness"
	REAL Dh = 0.012 UNITS u_m "Heater diameter"
	REAL L = 1 UNITS u_m "Pipe length"
	REAL z_in = 0 UNITS u_m "Elevation of inlet wrt user-defined base"
	REAL dz = 0 UNITS u_m "Elevation change between inlet and outlet, negative for downward slopes"

	ENUM SinglePhaseHTC choice_1p = UseDittusBoelter //Gnielinski,DittusBoelter,Colburn
	ENUM TwoPhaseHTC choice_2p = UseChen
	REAL tauAlpha = 3 UNITS u_s "Time constant for low-pass HTC filter"
	REAL Q_max = 1000 UNITS u_W "Maximum heater power"
	REAL alpha_a0 = 10 UNITS u_W_m2K
	REAL alpha_l = 1000 UNITS u_W_m2K "Liquid HTC, for Constant HTC"
	REAL alpha_tp = 1000 UNITS u_W_m2K "Two Phase HTC, for Constant HTC"
	REAL alpha_g = 1000 UNITS u_W_m2K "Vapour HTC, for Constant HTC"
	REAL m_steady = 0.05 UNITS u_kg_s "Mass flow rate at which constant HTCs calculated"
	REAL m0 = 0.005 UNITS u_kg_s "Nominal (design) mass flow rate"
	REAL dP0 = 0.01 UNITS u_bar "Nominal (design) pressure drop"

	ENUM InitialConditions init = Ph "Specify initial condition variables"
	REAL P0_i = 20 UNITS u_bar "Initial inlet pressure"
	REAL h0_i = 250000 UNITS u_J_kg "Initial inlet specific enthalpy"	
	REAL T0_i = 255 UNITS u_K "Initial inlet temperature"
	REAL x0_i = 0.5 "Initial inlet vapour quality"
	REAL P0_o = 20 UNITS u_bar "Initial outlet pressure"
	REAL h0_o = 250000 UNITS u_J_kg "Initial outlet specific enthalpy"	
	REAL T0_o = 255 UNITS u_K "Initial outlet temperature"
	REAL x0_o = 0.5 "Initial outlet vapour quality"
	REAL Tw0_i = 245 UNITS u_K "Initial wall temperature - first segment" // Same initial temperatures for heater and wall
	REAL Tw0_o = 245 UNITS u_K "Initial wall temperature - last segment"
TOPOLOGY
	AnnularStream_Ph (
		n = n,
		ConstantRefHTC = ConstantRefHTC,
		DynamicHTC = DynamicHTC,
		vfModel = vfModel
	) Stream (
		Do = Do-2*Th,
		Di = Dh,
		L = L,
		z_in = z_in,
		dz = dz,
		choice_1p = choice_1p,
		choice_2p = choice_2p,
		tauAlpha = tauAlpha,
		alpha_l = alpha_l,
		alpha_tp = alpha_tp,
		alpha_g = alpha_g,
		m_steady = m_steady,
		m0 = m0,
		dP0 = dP0,
		init = init,
		P0_i = P0_i,
		h0_i = h0_i,
		T0_i = T0_i,
		x0_i = x0_i,
		P0_o = P0_o,
		h0_o = h0_o,
		T0_o = T0_o,
		x0_o = x0_o
	)
	CartridgeHeaterCtrl (n=n) Heater (
		mat = HeaterMaterial,
		D = Dh,
		L = L,
		Q_max = Q_max,
		Tw0_i = Tw0_i,
		Tw0_o = Tw0_o
	)
	AnnularWall_naturalConvection (
		n = n,
		accountForRadiation = FALSE,
		ConstantAirHTC = TRUE
	) Wall (
		mat = TubeMaterial,
		Do = Do,
		Th = Th,
		L = L,
		alpha_a0 = alpha_a0,
		Tw0_i = Tw0_i,
		Tw0_o = Tw0_o
	)
	CONNECT T_amb TO Wall.T_amb
	CONNECT K TO Heater.K
	CONNECT Heater.tp_out TO Stream.tp_in
	CONNECT Stream.tp_out TO Wall.tp_in
	CONNECT f_in TO Stream.f_in
	CONNECT Stream.f_out TO f_out
CONTINUOUS
	T_heater.signal[1] = ZONE(ToutInCelsius) SUM(i IN 1,n; Heater.Tw[i])/n - 273.15 OTHERS SUM(i IN 1,n ; Heater.Tw[i])/n
	T_wall.signal[1] = ZONE(ToutInCelsius) SUM(i IN 1,n; Wall.Tw[i])/n - 273.15 OTHERS SUM(i IN 1,n ; Wall.Tw[i])/n
END COMPONENT




COMPONENT DummyLoad_ru_Tamb (
	INTEGER n = 3 "Number of control volumes",
	BOOLEAN ConstantRefHTC = TRUE "FALSE if using correlations",
	BOOLEAN DynamicHTC = TRUE "If TRUE, applies a low-pass filter on HTC",
	ENUM VoidFractionModel vfModel = Homogeneous "Two-phase flow void fraction model"
)
"Pipe model with a cartridge heater inside"
PORTS
	IN fluid f_in
	OUT fluid f_out
	IN analog_signal(n=1) T_amb "Ambient Temperature (K)"
	IN analog_signal(n=1) K "Heater control signal"
	OUT analog_signal(n=1) T_heater
	OUT analog_signal(n=1) T_wall
DATA
	BOOLEAN ToutInCelsius = TRUE "If false, outputs Kelvin"

	ENUM PipeMat TubeMaterial = SS_304
	ENUM PipeMat HeaterMaterial = SS_304
	REAL Do = 0.018 UNITS u_m "Pipe outer diameter"
	REAL Th = 0.001 UNITS u_m "Pipe thickness"
	REAL Dh = 0.012 UNITS u_m "Heater diameter"
	REAL L = 1 UNITS u_m "Pipe length"
	REAL z_in = 0 UNITS u_m "Elevation of inlet wrt user-defined base"
	REAL dz = 0 UNITS u_m "Elevation change between inlet and outlet, negative for downward slopes"

	ENUM SinglePhaseHTC choice_1p = UseDittusBoelter //Gnielinski,DittusBoelter,Colburn
	ENUM TwoPhaseHTC choice_2p = UseChen
	REAL tauAlpha = 3 UNITS u_s "Time constant for low-pass HTC filter"
	REAL Q_max = 1000 UNITS u_W "Maximum heater power"
	REAL alpha_a0 = 10 UNITS u_W_m2K
	REAL alpha_l = 1000 UNITS u_W_m2K "Liquid HTC, for Constant HTC"
	REAL alpha_tp = 1000 UNITS u_W_m2K "Two Phase HTC, for Constant HTC"
	REAL alpha_g = 1000 UNITS u_W_m2K "Vapour HTC, for Constant HTC"
	REAL m_steady = 0.05 UNITS u_kg_s "Mass flow rate at which constant HTCs calculated"
	REAL m0 = 0.005 UNITS u_kg_s "Nominal (design) mass flow rate"
	REAL dP0 = 0.01 UNITS u_bar "Nominal (design) pressure drop"

	ENUM InitialConditions init = Ph "Specify initial condition variables"
	REAL P0_i = 20 UNITS u_bar "Initial inlet pressure"
	REAL h0_i = 250000 UNITS u_J_kg "Initial inlet specific enthalpy"	
	REAL T0_i = 255 UNITS u_K "Initial inlet temperature"
	REAL x0_i = 0.5 "Initial inlet vapour quality"
	REAL P0_o = 20 UNITS u_bar "Initial outlet pressure"
	REAL h0_o = 250000 UNITS u_J_kg "Initial outlet specific enthalpy"	
	REAL T0_o = 255 UNITS u_K "Initial outlet temperature"
	REAL x0_o = 0.5 "Initial outlet vapour quality"
	REAL Tw0_i = 245 UNITS u_K "Initial wall temperature - first segment" // Same initial temperatures for heater and wall
	REAL Tw0_o = 245 UNITS u_K "Initial wall temperature - last segment"
TOPOLOGY
	AnnularStream_ru (
		n = n,
		ConstantRefHTC = ConstantRefHTC,
		DynamicHTC = DynamicHTC,
		vfModel = vfModel
	) Stream (
		Do = Do-2*Th,
		Di = Dh,
		L = L,
		z_in = z_in,
		dz = dz,
		choice_1p = choice_1p,
		choice_2p = choice_2p,
		tauAlpha = tauAlpha,
		alpha_l = alpha_l,
		alpha_tp = alpha_tp,
		alpha_g = alpha_g,
		m_steady = m_steady,
		m0 = m0,
		dP0 = dP0,
		init = init,
		P0_i = P0_i,
		h0_i = h0_i,
		T0_i = T0_i,
		x0_i = x0_i,
		P0_o = P0_o,
		h0_o = h0_o,
		T0_o = T0_o,
		x0_o = x0_o
	)
	CartridgeHeaterCtrl (n=n) Heater (
		mat = HeaterMaterial,
		D = Dh,
		L = L,
		Q_max = Q_max,
		Tw0_i = Tw0_i,
		Tw0_o = Tw0_o
	)
	AnnularWall_naturalConvection (
		n = n,
		accountForRadiation = FALSE,
		ConstantAirHTC = TRUE
	) Wall (
		mat = TubeMaterial,
		Do = Do,
		Th = Th,
		L = L,
		alpha_a0 = alpha_a0,
		Tw0_i = Tw0_i,
		Tw0_o = Tw0_o
	)
	CONNECT T_amb TO Wall.T_amb
	CONNECT K TO Heater.K
	CONNECT Heater.tp_out TO Stream.tp_in
	CONNECT Stream.tp_out TO Wall.tp_in
	CONNECT f_in TO Stream.f_in
	CONNECT Stream.f_out TO f_out
CONTINUOUS
	T_heater.signal[1] = ZONE(ToutInCelsius) SUM(i IN 1,n; Heater.Tw[i])/n - 273.15 OTHERS SUM(i IN 1,n ; Heater.Tw[i])/n
	T_wall.signal[1] = ZONE(ToutInCelsius) SUM(i IN 1,n; Wall.Tw[i])/n - 273.15 OTHERS SUM(i IN 1,n ; Wall.Tw[i])/n
END COMPONENT


--------------------------------------------------------------------------------
// Concentric Tube Heat Exchangers
--------------------------------------------------------------------------------
ABSTRACT COMPONENT CTHX (
	INTEGER n = 3 "Number of control volumes",
	BOOLEAN ConstantRefHTC = TRUE "False if using correlations",
	BOOLEAN DynamicHTC = TRUE "If TRUE applies low-pass filter to HTC",
	BOOLEAN Countercurrent = TRUE "False for parallel flow",
	ENUM VoidFractionModel vfModel = Homogeneous "Two-phase flow slip-ration void fraction model"
)
"Base model for concentric tube heat exchangers	"
// Does not contain fluid streams, only defines internal wall
PORTS
	IN fluid fi_in // inner inlet
	OUT fluid fi_out // inner outlet
	IN fluid fo_in // outer inlet
	OUT fluid fo_out // outer outlet
DATA
	ENUM PipeMat mat = SS_304 "Wall material"
	REAL Do_o = 0.008 UNITS u_m "Outer-tube outer diameter"
	REAL Th_o = 0.001 UNITS u_m "Outer-tube wall thickness"
	REAL Do_i = 0.003 UNITS u_m "Inner-tube outer diameter"
	REAL Th_i = 0.001 UNITS u_m "Inner-tube wall thickness"
	REAL L = 10 UNITS u_m "Tube length"
	REAL z_in = 0 UNITS u_m "Elevation of inner-tube inlet wrt user-defined base"
	REAL dz = 0 UNITS u_m "Elevation change between inner-tube inlet and outlet, negative for downward slopes"

	ENUM SinglePhaseHTC choice_1p = UseDittusBoelter
	ENUM TwoPhaseHTC choice_2p = UseChen
	REAL tauAlpha = 3 UNITS u_s "Time constant for low-pass HTC filter"
	REAL alpha_l_o = 1000 UNITS u_W_m2K "Outer stream liquid HTC"
	REAL alpha_tp_o = 1000 UNITS u_W_m2K "Outer stream two-phase HTC"
	REAL alpha_g_o = 1000 UNITS u_W_m2K "Outer stream vapour HTC"
	REAL m_steady_o = 0.015 UNITS u_kg_s "Mass flow rate at which outer stream HTCs calculated"
	REAL alpha_l_i = 1000 UNITS u_W_m2K "Inner stream liquid HTC"
	REAL alpha_tp_i = 1000 UNITS u_W_m2K "Inner stream two-phase HTC"
	REAL alpha_g_i = 1000 UNITS u_W_m2K "Inner stream vapour HTC"
	REAL m_steady_i = 0.02 UNITS u_kg_s "Mass flow rate at which inner stream HTCs calculated"
	REAL m0_o = 0.005 UNITS u_kg_s "Outer tube - Nominal (design) mass flow rate"
	REAL dP0_o = 0.01 UNITS u_bar "Outer tube - Nominal (design) pressure drop"
	REAL m0_i = 0.005 UNITS u_kg_s "Inner tube - Nominal (design) mass flow rate"
	REAL dP0_i = 0.01 UNITS u_bar "Inner tube - Nominal (design) pressure drop"

	ENUM InitialConditions init_i = Ph "Initial condition specification for inner stream"
	REAL Pi0_i = 20 UNITS u_bar "Inner tube - Initial inlet pressure"
	REAL hi0_i = 250000 UNITS u_J_kg "Inner tube - Initial inlet specific enthalpy"
	REAL Ti0_i = 255 UNITS u_K "Inner tube - Initial inlet temperature"
	REAL xi0_i = 0.5 "Inner tube - Initial inlet vapour quality"
	REAL Pi0_o = 20 UNITS u_bar "Inner tube - Initial outlet pressure"
	REAL hi0_o = 250000 UNITS u_J_kg "Inner tube - Initial outlet specific enthalpy"
	REAL Ti0_o = 255 UNITS u_K "Inner tube - Initial outlet temperature"
	REAL xi0_o = 0.5 "Inner tube - Initial outlet vapour quality"
	REAL Twi0_i = 245 UNITS u_K "Inner tube wall - Initial temperature - first segment"
	REAL Twi0_o = 245 UNITS u_K "Inner tube wall - Initial temperature - last segment"

	ENUM InitialConditions init_o = Ph "Initial condition specification for outer stream"
	REAL Po0_i = 20 UNITS u_bar "Outer tube - Initial inlet pressure"
	REAL ho0_i = 250000 UNITS u_J_kg "Outer tube - Initial inlet specific enthalpy"
	REAL To0_i = 255 UNITS u_K "Outer tube - Initial inlet temperature"
	REAL xo0_i = 0.5 "Outer tube - Initial inlet vapour quality"
	REAL Po0_o = 20 UNITS u_bar "Outer tube - Initial outlet pressure"
	REAL ho0_o = 250000 UNITS u_J_kg "Outer tube - Initial outlet specific enthalpy"
	REAL To0_o = 255 UNITS u_K "Outer tube - Initial outlet temperature"
	REAL xo0_o = 0.5 "Outer tube - Initial outlet vapour quality"
	REAL Two0_i = 245 UNITS u_K "Outer tube wall - Initial temperature - first segment"
	REAL Two0_o = 245 UNITS u_K "Outer tube wall - Initial temperature - last segment"
DECLS
	// Cannot account for height changes mid-simulation:
	DISCR REAL z_in_outer UNITS u_m // discr variable calculated only once
	DISCR REAL dz_outer UNITS u_m
TOPOLOGY
	AnnularWall_twoPort (n=n, Countercurrent=Countercurrent) InnerWall (
		mat = mat,
		Do = Do_i,
		Th = Th_i,
		L = L,
		Tw0_i = Twi0_i,
		Tw0_o = Twi0_o
	)
INIT
	IF (Countercurrent) THEN
		z_in_outer = z_in + dz
		dz_outer = -dz
	ELSE
		z_in_outer = z_in
		dz_outer = dz
	END IF
DISCRETE
	ASSERT (abs(dz)<=L) FATAL "CTHX elevation change entered improperly, abs(dz) must always be <= L"
END COMPONENT



ABSTRACT COMPONENT CTHX_Insulated IS_A CTHX
"Assumes outer wall is perfectly insulated against ambient"
TOPOLOGY
	AnnularWall_onePort (n=n) OuterWall (
		mat = mat,
		Do = Do_o,
		Th = Th_o,
		L = L,
		Tw0_i = Two0_i,
		Tw0_o = Two0_o
	)
END COMPONENT



COMPONENT CTHX_Ph IS_A CTHX_Insulated
"Concentric Tube Heat Exchanger with Ph state variables"
TOPOLOGY
	FluidStream_Ph (
		n = n,
		ConstantRefHTC = ConstantRefHTC,
		DynamicHTC = DynamicHTC,
		vfModel = vfModel
	) InnerStream (
		D = Do_i-2*Th_i,
		L = L,
		z_in = z_in,
		dz = dz,
		choice_1p = choice_1p,
		choice_2p = choice_2p,
		tauAlpha = tauAlpha,
		alpha_l = alpha_l_i,
		alpha_tp = alpha_tp_i,
		alpha_g = alpha_g_i,
		m_steady = m_steady_i,
		m0 = m0_i,
		dP0 = dP0_i,
		init = init_i,
		P0_i = Pi0_i,
		h0_i = hi0_i,
		T0_i = Ti0_i,
		x0_i = xi0_i,
		P0_o = Pi0_o,
		h0_o = hi0_o,
		T0_o = Ti0_o,
		x0_o = xi0_o
	)
	AnnularStream_Ph (
		n = n,
		ConstantRefHTC = ConstantRefHTC,
		DynamicHTC = DynamicHTC,
		vfModel = vfModel
	) OuterStream (
		Do = Do_o-2*Th_o,
		Di = Do_i,
		L = L,
		z_in = z_in_outer,
		dz = dz_outer,
		choice_1p = choice_1p,
		choice_2p = choice_2p,
		tauAlpha = tauAlpha,
		alpha_l = alpha_l_o,
		alpha_tp = alpha_tp_o,
		alpha_g = alpha_g_o,
		m_steady = m_steady_o,
		m0 = m0_o,
		dP0 = dP0_o,
		init = init_o,
		P0_i = Po0_i,
		h0_i = ho0_i,
		T0_i = To0_i,
		x0_i = xo0_i,
		P0_o = Po0_o,
		h0_o = ho0_o,
		T0_o = To0_o,
		x0_o = xo0_o
	)
	CONNECT fi_in TO InnerStream
	CONNECT InnerStream TO fi_out 
	CONNECT fo_in TO OuterStream
	CONNECT OuterStream TO fo_out
	CONNECT InnerStream.tp_out TO InnerWall.tp_in
	CONNECT InnerWall.tp_out TO OuterStream.tp_in
	CONNECT OuterStream.tp_out TO OuterWall.tp_in
END COMPONENT



COMPONENT CTHX_ru IS_A CTHX_Insulated
"Concentric Tube Heat Exchanger with Ph state variables"
TOPOLOGY
	FluidStream_ru (
		n = n,
		ConstantRefHTC = ConstantRefHTC,
		DynamicHTC = DynamicHTC,
		vfModel = vfModel
	) InnerStream (
		D = Do_i-2*Th_i,
		L = L,
		z_in = z_in,
		dz = dz,
		choice_1p = choice_1p,
		choice_2p = choice_2p,
		tauAlpha = tauAlpha,
		alpha_l = alpha_l_i,
		alpha_tp = alpha_tp_i,
		alpha_g = alpha_g_i,
		m_steady = m_steady_i,
		m0 = m0_i,
		dP0 = dP0_i,
		init = init_i,
		P0_i = Pi0_i,
		h0_i = hi0_i,
		T0_i = Ti0_i,
		x0_i = xi0_i,
		P0_o = Pi0_o,
		h0_o = hi0_o,
		T0_o = Ti0_o,
		x0_o = xi0_o
	)
	AnnularStream_ru (
		n = n,
		ConstantRefHTC = ConstantRefHTC,
		DynamicHTC = DynamicHTC,
		vfModel = vfModel
	) OuterStream (
		Do = Do_o-2*Th_o,
		Di = Do_i,
		L = L,
		z_in = z_in_outer,
		dz = dz_outer,
		choice_1p = choice_1p,
		choice_2p = choice_2p,
		tauAlpha = tauAlpha,
		alpha_l = alpha_l_o,
		alpha_tp = alpha_tp_o,
		alpha_g = alpha_g_o,
		m_steady = m_steady_o,
		m0 = m0_o,
		dP0 = dP0_o,
		init = init_o,
		P0_i = Po0_i,
		h0_i = ho0_i,
		T0_i = To0_i,
		x0_i = xo0_i,
		P0_o = Po0_o,
		h0_o = ho0_o,
		T0_o = To0_o,
		x0_o = xo0_o
	)
	CONNECT fi_in TO InnerStream
	CONNECT InnerStream TO fi_out 
	CONNECT fo_in TO OuterStream
	CONNECT OuterStream TO fo_out
	CONNECT InnerStream.tp_out TO InnerWall.tp_in
	CONNECT InnerWall.tp_out TO OuterStream.tp_in
	CONNECT OuterStream.tp_out TO OuterWall.tp_in
END COMPONENT



--------------------------------------------------------------------------------
// BRAZED PLATE HEAT EXCHANGERS
--------------------------------------------------------------------------------
COMPONENT BPHXStream (
	INTEGER n = 5 "Number of segments",
	ENUM VoidFractionModel vfModel = Smith "Slip ratio based void fraction model"
)
"Generic stream with simplified geometry used for BPHX models"
PORTS
	IN fluid f_in
	OUT fluid f_out
	OUT thermal(n=n) tp
DATA
	REAL D_port = 0.024 UNITS u_m "Port diameter"
	REAL V = 1.16E-3 UNITS u_m3 "Hold up volume"
	REAL As = 1.56 UNITS u_m2 "Surface area of contact"

	REAL alpha_l = 500 UNITS u_W_m2K "Liquid HTC, for ConstantHTC==TRUE"
	REAL alpha_tp = 1000 UNITS u_W_m2K "2Phase HTC, for ConstantHTC==TRUE"
	REAL alpha_g = 400 UNITS u_W_m2K "Vapour HTC, for ConstantHTC==TRUE"
	REAL m_steady = 0.05 UNITS u_kg_s "Steady state mass flow rate for which HTC calculated"
	REAL m0 = 0.01 UNITS u_kg_s "Nominal mass flow rate"
	REAL dP0 = 0.001 UNITS u_bar "Nominal pressure drop"
	REAL dP_small = 1 UNITS u_bar "Range for smoothing dP correlation, leave default if unsure"
	ENUM InitialConditions init = Ph "Specify initial condition variables"
	REAL P0_i = 20 UNITS u_bar "Initial inlet pressure"
	REAL h0_i = 250000 UNITS u_J_kg "Initial inlet specific enthalpy"	
	REAL T0_i = 255 UNITS u_K "Initial inlet temperature"
	REAL x0_i = 0.5 "Initial inlet vapour quality"
	REAL P0_o = 20 UNITS u_bar "Initial outlet pressure"
	REAL h0_o = 250000 UNITS u_J_kg "Initial outlet specific enthalpy"	
	REAL T0_o = 255 UNITS u_K "Initial outlet temperature"
	REAL x0_o = 0.5 "Initial outlet vapour quality"
DECLS
	DISCR REAL Vol UNITS u_dm3 "Internal volume (Liters)"
	REAL drho_dh[n] "Partial density derivative wrt enthalpy"
	REAL drho_dP[n] "Partial density derivative wrt pressure"
	REAL h[n] UNITS u_J_kg "Specific density-weighted enthalpy"
	REAL M UNITS u_kg "Refrigerant charge"
	REAL m[n+1] UNITS u_kg_s "Mass flow rates between thermal cells"
	REAL mh[n+1] UNITS u_W "Enthalpy flow rate"
	REAL P[n] UNITS u_bar "Pressure"
	REAL rho[n] UNITS u_kg_m3 "Density"
	REAL T[n] UNITS u_K "Temperature"
	REAL Tsat[n] UNITS u_K "Saturated temperature"
	REAL delT[n] UNITS u_K "Superheat/Subcooling (positive is superheat)"
	REAL u[n] UNITS u_J "Specific internal energy"
	REAL x[n] "Vapour quality"
	PRIVATE REAL cp[n] UNITS u_J_kgK "Specific heat capacity"
	PRIVATE REAL gamma[n] "Void fraction"
	PRIVATE REAL k[n] UNITS u_W_mK "Thermal conductivity"
	PRIVATE REAL mu[n] UNITS u_Pas "Dynamic viscosity"
	PRIVATE REAL Pr[n] "Prandtl number"
	PRIVATE REAL sigma[n] UNITS u_N_m // surface tension
	PRIVATE REAL vsound[n]
	ENUM Phase phase[n]
	INTEGER ier,ipx,ipy // error codes
	DISCR REAL Pcrit UNITS u_bar
	EXPL REAL alpha[n] UNITS u_W_m2K
	EXPL REAL Q[n] UNITS u_W
	PRIVATE REAL P_Pa[n] UNITS u_Pa
	PRIVATE REAL m_avg[n] UNITS u_kg_s
	//PRIVATE REAL Pred[n] // reduced pressure
	REAL h_corr[n]
	REAL h_flow[n] // flow enthalpy
	REAL Qr_total UNITS u_W
OBJECTS
	SaturationProperties sat[n]
	RefGeometryRecord geo[n]
INIT
	f_out.is_C = TRUE
	Vol = V * 1000

	FOR (i IN 1,n)
		// set unused variables to zero to avoid confusion
		geo[i].EndSurfaces = 0
		geo[i].Do = 0
		geo[i].Di = 0
		geo[i].Th = 0
		geo[i].L = 0
		geo[i].As_i = 0
		geo[i].As_o = 0
		geo[i].V = V/n
		geo[i].Ap = (MATH.PI * (D_port/2)**2)
		geo[i].As = As/n
	END FOR
	linspace(P0_i*1E5,P0_o*1E5,n,P_Pa,1)
	assignInitDiscretized(n,init,f_in.fluid,P0_i,P0_o,h0_i,h0_o,T0_i,T0_o,x0_i,x0_o,P,h,rho,u,x)
	Pcrit = 73.77 // CRYO_PF_CritProp(f_in.fluid,fprop_pressure,ier)	
CONTINUOUS
	M = SUM(i IN 1,n; geo[i].V*rho[i])
<m_i> 	m[1] = f_in.m
<m_o> 	m[n+1] = f_out.m
<m>	EXPAND (i IN 2,n) m[i] = m0*regRoot2((P[i-1]-P[i]),dP_small)/sqrt((dP0/n))
	EXPAND (i IN 1,n) delT[i] = T[i] - Tsat[i]
	// Heat transfer calculations
<Q>	EXPAND_BLOCK (i IN 1,n)
		m_avg[i] = 0.5*(m[i] + m[i+1])
		alpha[i] = max(10,getSmoothHTC(alpha_g,alpha_l,alpha_tp,m_avg[i],m_steady,x[i]))
		Q[i] = alpha[i]*geo[i].As*(T[i]-tp.Tk[i])
		tp.q[i] = Q[i]
	END EXPAND_BLOCK
	Qr_total = SUM(i IN 1,n; Q[i])
	// Slip ratio
	EXPAND_BLOCK (i IN 1,n)
		h_corr[i] = getCorrectionEnthalpy(vfModel,h[i],rho[i],sat[i])
		h_flow[i] = h[i] + h_corr[i]
	END EXPAND_BLOCK
	// Enthalpy flow rate
<mhi>	mh[1] = m[1]*donor_cell(m[1],f_in.hf,h_flow[1])
<mho>	mh[n+1] = m[n+1]*donor_cell(m[n+1],h_flow[n],f_out.hb)
<mh>	EXPAND(i IN 2,n) mh[i] = m[i]*donor_cell(m[i],h_flow[i-1],h_flow[i])
	// Governing equations
<gov>	EXPAND_BLOCK (i IN 1,n)
		P[i] = P_Pa[i]*1E-5
		geo[i].V*(drho_dP[i]*P_Pa[i]'+drho_dh[i]*h[i]') = m[i] - m[i+1]
		geo[i].V*((h[i]*drho_dP[i]-1)*P_Pa[i]' + (h[i]*drho_dh[i]+rho[i])*h[i]') = mh[i] - mh[i+1] - Q[i]
		getStatePh(f_in.fluid,P[i],h[i],cp[i],drho_dP[i],drho_dh[i],
			gamma[i],k[i],mu[i],phase[i],Pr[i],rho[i],sigma[i],
			T[i],Tsat[i],u[i],vsound[i],x[i],sat[i],ier,ipx,ipy)
	END EXPAND_BLOCK
	// Port calculations
	f_in.fluid = f_out.fluid
	f_in.n_fluid = f_out.n_fluid		
	f_in.Gcrit = CRYO_Gcrit_fun(f_in.fluid,P[1]*1E5,rho[1],T[1],0,f_in.v,vsound[1],ier,ipy)
<hb> 	f_in.hb = h[1]
	f_in.rho = rho[1]
<Pi> 	f_in.P = P[1]	
	f_out.Gcrit = CRYO_Gcrit_fun(f_in.fluid,P[n]*1E5,rho[n],T[n],0,f_out.v,vsound[n],ier,ipy)
<hf> 	f_out.hf = h[n]
	f_out.rho = rho[n]
<Po> 	f_out.P = P[n]
	f_in.A = geo[1].Ap
<Pia> 	f_in.P_aux = P[1] + geo[1].V*rho[1]*vsound[1]*1E-5/geo[1].As
	f_in.I = 0.5*geo[1].V/geo[1].As**2
	f_in.v = f_in.m/(rho[1]*geo[1].Ap)
	f_out.A = geo[n].Ap
<Poa> 	f_out.P_aux = P[n] + geo[n].V*rho[n]*vsound[n]*1E-5/geo[n].As
	f_out.I = 0.5*geo[n].V/geo[n].As**2
	f_out.v = f_out.m/(rho[n]*geo[n].Ap)
END COMPONENT



COMPONENT BPHX_SWEP (
	//INTEGER nPlate = 40 "Number of plates",
	//INTEGER nPass = 1 "Number of passes",
	INTEGER n = 5 "Number of segments to discretize the HX into",
	ENUM VoidFractionModel vfModel = Smith "Slip ratio based void fraction model",
	BOOLEAN Countercurrent = FALSE "If TRUE, outlet of chiller matches inlet of secondary"
)
"Component model of SWEP B16DWx30/1P-SC-U Brazed Plate Heat Exchanger used in CORA"
// Inputs derived from results of SWEP's SSP G7 software
PORTS
	IN fluid ch_in
	OUT fluid ch_out
	IN fluid f_in
	OUT fluid f_out
	OUT analog_signal(n=1) Q_ch
DATA
	// Geometry
	ENUM PipeMat mat = SS_316 "Material of parts in contact with fluid"
	REAL D_port = 0.024 UNITS u_m "Port diameter"
	REAL As = 1.56 UNITS u_m2 "Total heat transfer surface area"
	REAL V_i = 1.16E-3 UNITS u_m3 "Hold-up volume, inner circuit"
	REAL V_o = 1.22E-3 UNITS u_m3 "Hold up volume, outer circuit"
	//REAL HX_W = 119.5E-3 UNITS u_m "Heat exchanger width - Dimension B" // used for outer surface area calculations
	//REAL HX_H = 329E-3 UNITS u_m "Heat exchanger height - Dimension C"
	//REAL HX_D = 92E-3 UNITS u_m "Heat exchanger depth - Dimension F"
	REAL Mw = 14.1 UNITS u_kg "Heat exchanger weight"
	// Chiller side
	REAL alpha_l_ch = 750 UNITS u_W_m2K
	REAL alpha_tp_ch = 2500 UNITS u_W_m2K
	REAL alpha_g_ch = 500 UNITS u_W_m2K
	REAL m_steady_ch = 0.0208 UNITS u_kg_s
	REAL m0_ch = 0.02 UNITS u_kg_s 
	REAL dP0_ch = 0.01 UNITS u_bar
	REAL dP_small_ch = 1 UNITS u_bar "Range for low mass flow smoothing, leave default if unsure"
	// Primary side
	REAL alpha_l_sec = 1000 UNITS u_W_m2K
	REAL alpha_tp_sec = 10000 UNITS u_W_m2K
	REAL alpha_g_sec = 800 UNITS u_W_m2K
	REAL m_steady_sec = 0.015 UNITS u_kg_s
	REAL m0_sec = 0.015 UNITS u_kg_s
	REAL dP0_sec = 0.1 UNITS u_bar
	REAL dP_small_sec = 1 UNITS u_bar "Range for low mass flow smoothing, leave default if unsure"
	// Initial conditions
	ENUM InitialConditions init_ch = Ph "Specify initial condition variables for primary side"
	REAL P0_i_ch = 10.908 UNITS u_bar "Chiller side initial inlet pressure"
	REAL h0_i_ch = 301750 UNITS u_J_kg "Chiller side initial inlet enthalpy"
	REAL T0_i_ch = 293.15 UNITS u_K "Chiller side initial inlet temperature"
	REAL x0_i_ch = 0.5 "Chiller side initial inlet vapour quality"
	REAL P0_o_ch = 10.908 UNITS u_bar "Chiller side initial outlet pressure"
	REAL h0_o_ch = 301750 UNITS u_J_kg "Chiller side initial outlet enthalpy"
	REAL T0_o_ch = 293.15 UNITS u_K "Chiller side initial outlet temperature"
	REAL x0_o_ch = 0.5 "Chiller side initial outlet vapour quality"
	ENUM InitialConditions init_sec = Ph "Specify initial condition variables for secondary side"
	REAL P0_i_sec = 57.291 UNITS u_bar "Secondary side initial inlet pressure"
	REAL h0_i_sec = 331870 UNITS u_J_kg "Secondary side initial inlet enthalpy"
	REAL T0_i_sec = 293.15 UNITS u_K "Secondary side initial inlet temperature"
	REAL x0_i_sec = 0.5 "Secondary side initial inlet vapour quality"
	REAL P0_o_sec = 57.291 UNITS u_bar "Secondary side initial outlet pressure"
	REAL h0_o_sec = 331870 UNITS u_J_kg "Secondary side initial outlet enthalpy"
	REAL T0_o_sec = 293.15 UNITS u_K "Secondary side initial outlet temperature"
	REAL x0_o_sec = 0.5 "Secondary side initial outlet vapour quality"
	REAL Tw0_i = 255 UNITS u_K "Initial wall temperature - first segment"
	REAL Tw0_o = 255 UNITS u_K "Initial wall temperature - last segment"
TOPOLOGY
	BPHXStream (n=n, vfModel=vfModel) Chiller (
		D_port = D_port,
		V = V_i,
		As = As,
		alpha_l = alpha_l_ch,
		alpha_tp = alpha_tp_ch,
		alpha_g = alpha_g_ch,
		m_steady = m_steady_ch,
		m0 = m0_ch,
		dP0 = dP0_ch,
		dP_small = dP_small_ch,
		init = init_ch,
		P0_i = P0_i_ch,
		h0_i = h0_i_ch,
		T0_i = T0_i_ch,
		x0_i = x0_i_ch,
		P0_o = P0_o_ch,
		h0_o = h0_o_ch,
		T0_o = T0_o_ch,
		x0_o = x0_o_ch
	)
	GenericWall (n=n) Wall (
		mat = mat,
		Mw = Mw,
		Tw0_i = Tw0_i,
		Tw0_o = Tw0_o
	)
	BPHXStream (n=n, vfModel=vfModel) Secondary (
		D_port = D_port,
		V = V_o,
		As = As,
		alpha_l = alpha_l_sec,
		alpha_tp = alpha_tp_sec,
		alpha_g = alpha_g_sec,
		m_steady = m_steady_sec,
		m0 = m0_sec,
		dP0 = dP0_sec,
		dP_small = dP_small_sec,
		init = init_sec,
		P0_i = P0_i_sec,
		h0_i = h0_i_sec,
		T0_i = T0_i_sec,
		x0_i = x0_i_sec,
		P0_o = P0_o_sec,
		h0_o = h0_o_sec,
		T0_o = T0_o_sec,
		x0_o = x0_o_sec
	)
	CONNECT f_in TO Secondary.f_in
	CONNECT f_out TO Secondary.f_out
	CONNECT ch_in TO Chiller.f_in
	CONNECT ch_out TO Chiller.f_out
	CONNECT Chiller.tp TO Wall.tp_in
	CONNECT Wall.tp_out TO Secondary.tp
CONTINUOUS
	Q_ch.signal[1] = -Chiller.Qr_total
END COMPONENT