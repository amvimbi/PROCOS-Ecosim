/*-----------------------------------------------------------------------------------------
 LIBRARY: CO2
 FILE: SENSORS
 AUTHOR: Viren Bhanot
 COMPANY: CERN
 DESCRIPTION: Fluid sensors
 CREATION DATE: 25/09/2017
-----------------------------------------------------------------------------------------*/
USE PORTS_LIB
USE MATH
USE PLC
USE CRYOLIB
USE CONTROL
USE THERMO_TABLE_INTERP


/*
 * INDEX
 * 1. SaturationSensor
 * 2. SuperheatSensor
 * 3. TemperatureSensor_Celsius
 * 4. PressureTemperature
 * 5. AirTemperatureSensor
 * 6. EnthalpySensor
 * 7. SaturationTemperature
 * 8. Subcooling
 */



COMPONENT SaturationSensor IS_A CRYOLIB.FluidSensor (
	BOOLEAN OutputTemperature = TRUE "If FALSE, outputs saturation pressure",
	BOOLEAN OutputCelsius = TRUE "If FALSE, Kelvin"
)
"Output saturation temperature or pressure"
DECLS
	REAL T UNITS u_K
	INTEGER jx,jy// error code
CONTINUOUS
	IF (OutputTemperature) INSERT
		IF (OutputCelsius) INSERT
			v = CRYO_PF_Tsat_vs_p(f_in.fluid, f_in.P, ier, jx) - 273.15
		ELSE
			v = CRYO_PF_Tsat_vs_p(f_in.fluid, f_in.P, ier, jx)
		END IF
		T = 0
	ELSE
		T = CRYO_PF_prop_vs_ph(f_in.fluid,f_in.P,f_in.h,fprop_temperature,ier,jx,jy)
		v = CRYO_PF_psat_vs_T(f_in.fluid, T, ier, jx)
	END IF
END COMPONENT



COMPONENT SuperheatSensor IS_A FluidSensor
"Superheating is positive, subcooling is negative"
DECLS
	REAL T UNITS u_K
	REAL Tsat UNITS u_K
	INTEGER jx, jy
CONTINUOUS
	T = CRYO_PF_prop_vs_ph(f_in.fluid, f_in.P, f_in.hf, fprop_temperature,ier,jx,jy)
	Tsat = CRYO_PF_Tsat_vs_p(f_in.fluid, f_in.P, ier, jx)
	v = T - Tsat // positive for superheating
END COMPONENT



ENUM SuperheatSubcool = {Superheat, Subcool} // ENUMs are considered better practice than Booleans

COMPONENT DeltaTSensor IS_A FluidSensor (
	ENUM SuperheatSubcool type = Superheat "Whether superheat should be positive or subcooling"
)
"Superheating or Subcooling sensor"
DECLS
	REAL T UNITS u_K
	REAL Tsat UNITS u_K
	INTEGER jx, jy
CONTINUOUS
	T = CRYO_PF_prop_vs_ph(f_in.fluid, f_in.P, f_in.hf, fprop_temperature,ier,jx,jy)
	Tsat = CRYO_PF_Tsat_vs_p(f_in.fluid, f_in.P, ier, jx)
	IF (type==Superheat) INSERT
		v = T-Tsat 
	ELSE
		v = Tsat-T
	END IF
END COMPONENT


COMPONENT TemperatureSensor_Celsius IS_A FluidSensor (BOOLEAN Threshold=FALSE)
"Temperature sensor (°C)"
DECLS
	INTEGER jx,jy
CONTINUOUS
	IF (Threshold == TRUE) INSERT
		v = IF (f_in.mb>TH) CRYO_PF_prop_vs_ph(f_in.fluid,  f_in.P, f_in.hb, fprop_temperature,ier2, jx,jy) - 273.15
			ELSE CRYO_PF_prop_vs_ph(f_in.fluid, f_in.P, f_in.hf, fprop_temperature,ier2, jx,jy) - 273.15
	ELSE 
		v = CRYO_PF_prop_vs_ph(f_in.fluid, f_in.P, f_in.hf, fprop_temperature,ier2, jx,jy) - 273.15
	END IF
END COMPONENT



COMPONENT PressureTemperature IS_A CRYOLIB.FluidChannel (
	BOOLEAN OutputCelsius = TRUE "FALSE: Kelvin. Note: Both Tsat and T will be in C",
	BOOLEAN SuperheatPositive = TRUE "If FALSE: Subcooling will be positive"
)
"Pressure/Temperature cross. Also outputs Subcooling and Saturation Temperature"
// Note: Regardless of reverse flow, always takes inlet port values
PORTS
	OUT analog_signal(n=1) TT "Temperature"
	OUT analog_signal(n=1) PT "Pressure"
	OUT analog_signal(n=1) ST "Saturation temperature"
	OUT analog_signal(n=1) SC "Subcooling (subcool positive, superheat negative)"
DECLS
	INTEGER jx,jy
CONTINUOUS
	PT.signal[1] = f_in.P
	IF OutputCelsius INSERT
		TT.signal[1] = CRYO_PF_prop_vs_ph(f_in.fluid,f_in.P,f_in.hf,fprop_temperature,ier,jx,jy) - 273.15
		ST.signal[1] = CRYO_PF_prop_vs_Px(f_in.fluid,f_in.P,0,fprop_temperature,ier,jx,jy) - 273.15 
	ELSE
		TT.signal[1] = CRYO_PF_prop_vs_ph(f_in.fluid,f_in.P,f_in.hf,fprop_temperature,ier,jx,jy)
		ST.signal[1] = CRYO_PF_prop_vs_Px(f_in.fluid,f_in.P,0,fprop_temperature,ier,jx,jy)
	END IF
	IF SuperheatPositive INSERT
		SC.signal[1] = TT.signal[1] - ST.signal[1]
	ELSE
		SC.signal[1] = ST.signal[1] - TT.signal[1]
	END IF
	f_in.Gcrit = f_out.Gcrit
	f_in.mf = f_out.mf
	f_in.mb = f_out.mb
	f_in.rho = f_out.rho	
	f_out.P = f_in.P
	f_out.hf = f_in.hf
	f_out.hb = f_in.hb
	f_in.is_C = f_out.is_C	
	f_in.v = f_out.v
	f_in.I = f_out.I	
	f_in.A = f_out.A	
	f_in.P_aux = f_out.P_aux
END COMPONENT



COMPONENT AirTemperatureSensor
"Returns dry bulb temperature in an air stream"
PORTS
	IN air a_in
	OUT air a_out
	OUT analog_signal(n=1) T_dry
DECLS
	REAL T_degC
CONTINUOUS
	a_in.P = a_out.P
	a_in.m = a_out.m
	a_in.T_dry = a_out.T_dry
	a_in.w = a_out.w
	T_dry.signal[1] = a_in.T_dry
	T_degC = a_in.T_dry - 273.15
END COMPONENT



COMPONENT EnthalpySensor IS_A FluidSensor
CONTINUOUS
	v = ZONE (f_in.m>=0) f_in.hf OTHERS f_in.hb
END COMPONENT



COMPONENT VapourQualitySensor IS_A FluidSensor
DECLS
	REAL h
	INTEGER ier,jx,jy
CONTINUOUS
	h = ZONE (f_in.m>=0) f_in.hf OTHERS f_in.hb
	v = CRYO_PF_prop_vs_ph(f_in.fluid,f_in.P,h,fprop_quality,ier,jx,jy)
END COMPONENT