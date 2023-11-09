/*-----------------------------------------------------------------------------------------
 LIBRARY: CO2
 FILE: BOUNDARIES
 AUTHOR: Viren Bhanot
 COMPANY: CERN
 DESCRIPTION: Component boundary conditions
 CREATION DATE: 10/07/2017
-----------------------------------------------------------------------------------------*/
USE CRYOLIB
USE MATH
USE PORTS_LIB



/*
 * INDEX
 * 
 * 1. ABSTRACT COMPONENTS
 * 	1. SourceSink
 * 	2. Source
 * 	3. Sink
 * 2. REFRIGERANT SOURCES
 * 	1. RefSource
 *	2. RefSource_Const
 * 3. REFRIGERANT SINKS
 * 	1. RefSink
 *	2. RefSink_Const
 * 4. AIR BOUNDARIES
 *	1. AirSource
 *	2. AirSource_Const
 *	3. AirSink
 *	4. Fan
 */



ENUM BoundaryVariable = {Enthalpy,Temperature,Quality}



--------------------------------------------------------------------------------
// ABSTRACT COMPONENTS
--------------------------------------------------------------------------------
ABSTRACT COMPONENT SourceSink (
	BOOLEAN CapacitiveBoundary = TRUE "TRUE if Pressure source/sink. False if mass flow",
	ENUM BoundaryVariable bound = Enthalpy "Specify second boundary variable"
)
"Base class for refrigerant source and sink components"
DECLS
	ENUM ChemName fld
	REAL P_val UNITS u_bar "Pressure value"
	REAL h_val UNITS u_J_kg "Enthalpy value"
	REAL T_val UNITS u_K "Temperature value"
	REAL x_val "Vapour quality value"
	REAL rho_val UNITS u_kg_m3 "Density value"
	REAL m_val UNITS u_kg_s "Mass flow rate"
	REAL vsound UNITS u_m_s
	INTEGER ier,jx,jy,jz
CONTINUOUS
	IF (bound==Enthalpy) INSERT
		vsound = CRYO_PF_prop_vs_ph(fld,P_val,h_val,fprop_vsound,ier,jx,jy)
		rho_val = CRYO_PF_prop_vs_ph(fld,P_val,h_val,fprop_density,ier,jx,jy)
		T_val = CRYO_PF_prop_vs_ph(fld,P_val,h_val,fprop_temperature,ier,jx,jy)
		x_val = CRYO_PF_prop_vs_ph(fld,P_val,h_val,fprop_quality,ier,jx,jy)
	ELSEIF (bound==Quality) INSERT
		vsound = CRYO_PF_prop_vs_Px(fld,P_val,x_val,fprop_vsound,ier,jx,jy)
		h_val = CRYO_PF_prop_vs_Px(fld,P_val,x_val,fprop_enthalpy,ier,jx,jy)
		rho_val = CRYO_PF_prop_vs_Px(fld,P_val,x_val,fprop_density,ier,jx,jy)
		T_val = CRYO_PF_prop_vs_Px(fld,P_val,x_val,fprop_temperature,ier,jx,jy)
	ELSE
		vsound = CRYO_PF_prop_vs_pT(fld,P_val,T_val,fprop_vsound,ier,jx,jy)
		h_val = CRYO_PF_prop_vs_pT(fld,P_val,T_val,fprop_enthalpy,ier,jx,jy)
		rho_val = CRYO_PF_prop_vs_pT(fld,P_val,T_val,fprop_density,ier,jx,jy)
		x_val = CRYO_PF_prop_vs_pT(fld,P_val,T_val,fprop_quality,ier,jx,jy)
	END IF
END COMPONENT



ABSTRACT COMPONENT Source IS_A SourceSink
"Base class for refrigerant Sources. Adds port assignments to SourceSink"
PORTS
	OUT fluid f_out
INIT
	f_out.is_C = CapacitiveBoundary
CONTINUOUS
	fld = f_out.fluid

	f_out.hf = h_val
	f_out.P = P_val
	IF (CapacitiveBoundary) INSERT
		f_out.P_aux = P_val
		f_out.rho = rho_val
		f_out.A = 1E5
		f_out.I = 0
		f_out.v = 0
		f_out.Gcrit = CRYO_Gcrit_fun(fld,P_val*1E5,rho_val,T_val,0,0,vsound,ier,jz)
	ELSE
		f_out.mb = ZONE (m_val<0) m_val OTHERS 0
		f_out.mf = ZONE (m_val>=0) m_val OTHERS 0
	END IF
END COMPONENT



ABSTRACT COMPONENT Sink IS_A SourceSink
"Base class for Sinks. Adds port assignments to SourceSink"
PORTS
	IN fluid f_in
INIT
	f_in.is_C = CapacitiveBoundary
CONTINUOUS
	fld = f_in.fluid

	f_in.hb = h_val
	f_in.P = P_val
	IF (CapacitiveBoundary) INSERT
		f_in.P = P_val
		f_in.P_aux = P_val
		f_in.rho = rho_val
		f_in.A = 1E5
		f_in.I = 0
		f_in.v = 0
		f_in.Gcrit = CRYO_Gcrit_fun(fld,P_val*1E5,rho_val,T_val,0,0,vsound,ier,jz)
	ELSE
		f_in.mb = ZONE (m_val<0) -m_val OTHERS 0
		f_in.mf = ZONE (m_val>=0) m_val OTHERS 0
	END IF
END COMPONENT



--------------------------------------------------------------------------------
// REFRIGERANT SOURCES
--------------------------------------------------------------------------------
COMPONENT RefSource IS_A Source
"Refrigerant source with port connections for boundary variables"
// Dynamically changing variables can be assigned using Ramps etc.
PORTS
	IN analog_signal(n=1) var1
	IN analog_signal(n=1) var2
CONTINUOUS
	// First variable
	IF (CapacitiveBoundary) INSERT
		P_val = var1.signal[1]
		m_val = 0
	ELSE
		m_val = var1.signal[1]
	END IF
	// Second variable
	IF (bound==Enthalpy) INSERT
		h_val = var2.signal[1]
	ELSEIF (bound==Temperature) INSERT
		T_val = var2.signal[1]
	ELSE
		x_val = var2.signal[1]
	END IF
END COMPONENT



COMPONENT RefSource_Const IS_A Source
"Refrigerant source for constant boundary conditions"
DATA
	REAL P = 20 UNITS u_bar "Pressure"
	REAL m = 0.05 UNITS u_kg_s "Mass flow rate"
	REAL h = 250000 UNITS u_J_kg "Enthalpy"
	REAL T = 255 UNITS u_K "Temperature"
	REAL x = 0.1 "Vapour quality"
CONTINUOUS
	// First variable
	IF (CapacitiveBoundary) INSERT
		P_val = P
		m_val = 0
	ELSE
		m_val = m
	END IF
	// Second variable
	IF (bound==Enthalpy) INSERT
		h_val = h
	ELSEIF (bound==Temperature) INSERT
		T_val = T
	ELSE
		x_val = x
	END IF
END COMPONENT
	


--------------------------------------------------------------------------------
// REFRIGERANT SINKS
--------------------------------------------------------------------------------
COMPONENT RefSink IS_A Sink
"Refrigerant sink with port connections for boundary variables"
PORTS
	IN analog_signal(n=1) var1
	IN analog_signal(n=1) var2
CONTINUOUS
	// First variable
	IF (CapacitiveBoundary) INSERT
		P_val = var1.signal[1]
		m_val = 0
	ELSE
		m_val = var1.signal[1]
	END IF
	// Second variable
	IF (bound==Enthalpy) INSERT
		h_val = var2.signal[1]
	ELSEIF (bound==Temperature) INSERT
		T_val = var2.signal[1]
	ELSE
		x_val = var2.signal[1]
	END IF
END COMPONENT



COMPONENT RefSink_Const IS_A Sink
"Refrigerant sink for constant boundary conditions"
DATA
	REAL P = 20 UNITS u_bar "Pressure"
	REAL m = 0.05 UNITS u_kg_s "Mass flow rate"
	REAL h = 250000 UNITS u_J_kg "Enthalpy"
	REAL T = 255 UNITS u_K "Temperature"
	REAL x = 0.1 "Vapour quality"
CONTINUOUS
	// First variable
	IF (CapacitiveBoundary) INSERT
		P_val = P
		m_val = 0
	ELSE
		m_val = m
	END IF
	// Second variable
	IF (bound==Enthalpy) INSERT
		h_val = h
	ELSEIF (bound==Temperature) INSERT
		T_val = T
	ELSE
		x_val = x
	END IF
END COMPONENT



--------------------------------------------------------------------------------
// AIR BOUNDARIES
--------------------------------------------------------------------------------
ENUM HumidityVariable = {RelativeHumidity,WetBulbTemperature}



COMPONENT AirSource (
	ENUM HumidityVariable choice = RelativeHumidity
)
PORTS
	IN analog_signal(n=1) m
	IN analog_signal(n=1) T_dry
	IN analog_signal(n=1) bdry
	OUT air a_out
DECLS
	CONST REAL m_small = 0.01
CONTINUOUS
	a_out.m = max(m.signal[1],m_small)
	a_out.T_dry = T_dry.signal[1]
	a_out.P = 101325
	IF (choice==RelativeHumidity) INSERT
		a_out.w = W_DBRHP(T_dry.signal[1],bdry.signal[1],101325) // external C++ function
	ELSE
		a_out.w = W_DBWBP(T_dry.signal[1],bdry.signal[1],101325)
	END IF
END COMPONENT



COMPONENT AirSource_Const (
	ENUM HumidityVariable choice = RelativeHumidity

)
"Air source with constant values"
PORTS
	OUT air a_out
DATA
	REAL m = 2 UNITS u_kg_s
	REAL T_dry = 300 UNITS u_K
	REAL bdry = 50
CONTINUOUS
	a_out.m = m
	a_out.T_dry = T_dry
	a_out.P = 101325
	IF (choice==RelativeHumidity) INSERT
		a_out.w = W_DBRHP(T_dry,bdry,101325)
	ELSE
		a_out.w = W_DBWBP(T_dry,bdry,101325)
	END IF
END COMPONENT



COMPONENT AirSink
"Air sink"
PORTS
	IN air a_out
DECLS
	REAL m UNITS u_kg_s
	REAL T_dry UNITS u_K
	REAL w "Humidity ratio"
	REAL P UNITS u_Pa
CONTINUOUS
	a_out.m = m
	a_out.T_dry = T_dry
	a_out.w = w
	a_out.P = P
END COMPONENT



COMPONENT Fan
"Simplified fan component. Merely imposes a mass flow rate on air stream"
PORTS
	IN air a_in
	OUT air a_out
	IN analog_signal(n=1) mDot
CONTINUOUS
	a_in.m = mDot.signal[1]
	a_out.m = a_in.m
	a_out.P = 101325
	a_in.T_dry = a_out.T_dry
	a_in.w = a_out.w
	a_in.P = a_out.P
END COMPONENT



COMPONENT Fan_speed
"Fan with prescribed mass flow rate"
PORTS
	OUT air a_out
	IN analog_signal(n=1) speed
DATA
	REAL mDot_max = 2 UNITS u_kg_s "Maximum mass flow rate"
	REAL T_dry = 300 UNITS u_K "Dry bulb temperature"
	REAL RH = 50 "Relative humidity"
CONTINUOUS
	a_out.m = mDot_max * max(speed.signal[1]/100,0)
	a_out.T_dry = T_dry
	a_out.w = W_DBRHP(T_dry,RH,101325)
	a_out.P = 101325
END COMPONENT