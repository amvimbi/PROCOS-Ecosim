/*-----------------------------------------------------------------------------------------
 LIBRARY: CO2
 FILE: TEST_CODE
 CREATION DATE: 12/07/2017
-----------------------------------------------------------------------------------------*/
USE CRYOLIB
USE MATH
USE PORTS_LIB
USE THERMAL
USE THERMO_TABLE_INTERP


COMPONENT testCPP
"Test to confirm that the EcosimProRefprop library is consistent with Matlab refprop calls"
// This test confirms that properties for the three refprop functions are correct
// Verified against Refprop. Partial derivatives: verified against matlab code
DECLS
	REAL P = 20e5
	REAL h
	REAL rho
	REAL ddph, ddhp
	CONST REAL hhf = 180000
	CONST REAL hhg = 400000
	REAL arhof,arhog,ahfs,ahgs,addphf,addhpf,addphg,addhpg
CONTINUOUS
	h = ZONE (TIME<1) 100e3 ZONE(TIME>=1 AND TIME<2) 250e3 OTHERS 500e3
	rho = rp_rhoPh(P,h)
	rp_partialDers(P,h,ddph,ddhp)
	rp_all(P,hhf,hhg,arhof,arhog,ahfs,ahgs,addphf,addhpf,addphg,addhpg)
END COMPONENT



FUNCTION NO_TYPE getPartialDerFD (
	IN ENUM ChemName chem "Working chemical constituent",
	IN REAL p_bar "Pressure",
	IN REAL h "Enthalpy",
	OUT REAL drho_dP "drho/dp at constant h",
	OUT REAL drho_dh "drho/dh at constant P",
	OUT INTEGER ier,
	OUT INTEGER ipx,
	OUT INTEGER ipy
)
"Calculate partial derivatives using the forward difference method"
// Very low accuracy
DECLS
	REAL dP, dh
	REAL P_Pa
	REAL eps = 1e-6
BODY
	P_Pa = p_bar*1e5
	dP = P_Pa*eps
	dh = h*eps
	drho_dP = (CRYO_PF_prop_vs_ph(chem,(P_Pa+dP)*1e-5,h,fprop_density,ier,ipx,ipy) - CRYO_PF_prop_vs_ph(chem,p_bar,h,fprop_density,ier,ipx,ipy)) / dP
	drho_dh = (CRYO_PF_prop_vs_ph(chem,p_bar,h+dh,fprop_density,ier,ipx,ipy) - CRYO_PF_prop_vs_ph(chem,p_bar,h,fprop_density,ier,ipx,ipy)) / dh
END FUNCTION



COMPONENT test_PartialDers
"Compared partial derivatives: Refprop vs Lookup Table vs Forward Difference"
PORTS
	IN analog_signal(n=1) pr
	IN analog_signal(n=1) enth
DECLS
	ENUM ChemName fluid = Carbon_Dioxide
	REAL P
	REAL ddph_eqn, ddph_fd
	REAL ddhp_eqn, ddhp_fd
	REAL ddph_rp, ddhp_rp
	REAL h
	REAL a,b // dummy
	INTEGER ier,xx,yy
	REAL x
CONTINUOUS
	P = pr.signal[1]
	h = enth.signal[1]
	x = CRYO_PF_prop_vs_ph(fluid,P,h,fprop_quality,ier,xx,yy)
	getPartialDers(fluid,P,h,ddph_eqn,ddhp_eqn,a,ier,xx,yy)
	getPartialDerFD(fluid,P,h,ddph_fd,ddhp_fd,ier,xx,yy)
	rp_partialDers(P*1e5,h,ddph_rp,ddhp_rp)
END COMPONENT