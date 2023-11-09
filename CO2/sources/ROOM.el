/*-----------------------------------------------------------------------------------------
 LIBRARY: CO2
 FILE: ROOM
 AUTHOR: Viren Bhanot
 COMPANY: CERN
 DESCRIPTION: Confined air space model
 CREATION DATE: 06/11/2017
-----------------------------------------------------------------------------------------*/
USE THERMAL
USE CRYOLIB
USE MATH
USE PLC
USE PORTS_LIB



COMPONENT LumpedRoomModel (
	INTEGER nSource = 2 "Number of sources of heat"
)
"Lumped volume with multiple internal heat sources"
PORTS
	IN air a_in
	OUT air a_out
	IN analog_signal(n=nSource) Qsource
	OUT analog_signal(n=1) Temp
DATA
	REAL L = 0.54 UNITS u_m "Length"
	REAL W = 0.54 UNITS u_m "Width"
	REAL H = 0.54 UNITS u_m "Height"
	REAL T0 = 299 UNITS u_K "Air initial temperature"
DECLS
	DISCR REAL V UNITS u_m3 "Air volume"
	REAL Cp_air UNITS u_J_kgK "Air specific heat"
	REAL hin
	REAL hout
	REAL k_air
	REAL mu_air
	REAL vDot "Volumetric flow rate"
	REAL ACH "Air changes per hour"
	REAL Q_int UNITS u_W "Total internal heat generation"
	REAL rho_air UNITS u_kg_m3 "Air density"
	REAL T UNITS u_K "Air temperature"
	INTEGER ier,iex,jx,jy "Error handling codes"
INIT
	V = L*W*H
	T = T0
CONTINUOUS
	hin = CRYO_PF_prop_vs_pT(CRYOLIB.Air,a_in.P*1E-5,a_in.T_dry,fprop_enthalpy,iex,jx,jy) // treat as dry air
	hout = CRYO_PF_prop_vs_pT(CRYOLIB.Air,a_out.P*1E-5,T,fprop_enthalpy,iex,jx,jy)
	Fluid_prop_vs_pT(THERMAL.Air,a_in.P,T,Cp_air,mu_air,rho_air,k_air,ier)
	Q_int = SUM(i IN 1,nSource; Qsource.signal[i])
	(V*rho_air*Cp_air) * T' = a_in.m*Cp_air*(a_in.T_dry - T) + Q_int
	vDot = a_in.m/rho_air // m3/s
	ACH = vDot*3600/V // 1/hour
	
	
	// Port assignments
	Temp.signal[1] = T
	a_in.m = a_out.m
	a_in.P = a_out.P
	a_in.w = a_out.w
	a_out.T_dry = T
END COMPONENT