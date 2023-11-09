/*-----------------------------------------------------------------------------------------
 LIBRARY: CO2
 FILE: PORTS
 CREATION DATE: 04/08/2017
-----------------------------------------------------------------------------------------*/
USE MATH
USE CRYOLIB



/*
 * moistair.lib:
 * Minimal psychrometrics function library written in C++. Uses ASHRAE's simplified
 * empirical correlations to evaluate relationships between humidity parameters
 */
"C++" FUNCTION REAL Pws_T(IN REAL v1) IN "moistairlib.lib"
"C++" FUNCTION REAL h_DBRHP(IN REAL DB, IN REAL RH, IN REAL P) IN "moistairlib.lib"
"C++" FUNCTION REAL h_DBTdP(IN REAL DB, IN REAL Td, IN REAL P) IN "moistairlib.lib"
"C++" FUNCTION REAL h_DBWBP(IN REAL DB, IN REAL WB, IN REAL P) IN "moistairlib.lib"
"C++" FUNCTION REAL Mu_DBRHP(IN REAL DB, IN REAL RH, IN REAL P) IN "moistairlib.lib"
"C++" FUNCTION REAL Mu_DBTdP(IN REAL DB, IN REAL Td, IN REAL P) IN "moistairlib.lib"
"C++" FUNCTION REAL Mu_DBWBP(IN REAL DB, IN REAL WB, IN REAL P) IN "moistairlib.lib"
"C++" FUNCTION REAL Pw_DBRHP(IN REAL DB, IN REAL RH, IN REAL P) IN "moistairlib.lib"
"C++" FUNCTION REAL Pw_DBWBP(IN REAL DB, IN REAL WB, IN REAL P) IN "moistairlib.lib"
"C++" FUNCTION REAL RH_DBTdP(IN REAL DB, IN REAL Td, IN REAL P) IN "moistairlib.lib"
"C++" FUNCTION REAL RH_DBWBP(IN REAL DB, IN REAL WB, IN REAL P) IN "moistairlib.lib"
"C++" FUNCTION REAL Td_DBRHP(IN REAL DB, IN REAL RH, IN REAL P) IN "moistairlib.lib"
"C++" FUNCTION REAL Td_DBWBP(IN REAL DB, IN REAL WB, IN REAL P) IN "moistairlib.lib"
"C++" FUNCTION REAL rho_DBRHP(IN REAL DB, IN REAL RH, IN REAL P) IN "moistairlib.lib"
"C++" FUNCTION REAL rho_DBTdP(IN REAL DB, IN REAL Td, IN REAL P) IN "moistairlib.lib"
"C++" FUNCTION REAL rho_DBWBP(IN REAL DB, IN REAL WB, IN REAL P) IN "moistairlib.lib"
"C++" FUNCTION REAL v_DBRHP(IN REAL DB, IN REAL RH, IN REAL P) IN "moistairlib.lib"
"C++" FUNCTION REAL v_DBTdP(IN REAL DB, IN REAL Td, IN REAL P) IN "moistairlib.lib"
"C++" FUNCTION REAL v_DBWBP(IN REAL DB, IN REAL WB, IN REAL P) IN "moistairlib.lib"
"C++" FUNCTION REAL W_DBRHP(IN REAL DB, IN REAL RH, IN REAL P) IN "moistairlib.lib"
"C++" FUNCTION REAL W_DBWBP(IN REAL DB, IN REAL WB, IN REAL P) IN "moistairlib.lib"
"C++" FUNCTION REAL WB_DBRHP(IN REAL DB, IN REAL RH, IN REAL P) IN "moistairlib.lib"
"C++" FUNCTION REAL WB_DBTdP(IN REAL DB, IN REAL Td, IN REAL P) IN "moistairlib.lib"
"C++" FUNCTION REAL Ws_TP(IN REAL T, IN REAL P) IN "moistairlib.lib"



PORT air
"Moist air properties port"
// Uses external C++ library for calculating moist air properties.
// Uses ASHRAE handbook equations for approximate values.
	SUM REAL m = 0.0 UNITS u_kg_s RANGE -1000, 1000 "Air mass flow rate"
	EQUAL REAL T_dry = 273.15 UNITS u_K "Dry bulb temperature"
	EQUAL REAL w = 0.5 "Humidity ratio"
	EQUAL REAL P = 101325 UNITS u_Pa "Air pressure"
END PORT