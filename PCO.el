/*-----------------------------------------------------------------------------------------
 LIBRARY: DEMO
 FILE: PCO
 AUTHOR: Viren Bhanot
 COMPANY: CERN
 DESCRIPTION: Controller logic for Demo plant and accumulators
 CREATION DATE: 02/06/2023
-----------------------------------------------------------------------------------------*/
USE CO2
USE PORTS_LIB
USE CRYOLIB
USE CONTROL
USE MATH
USE PLC

--------------------------------------------------------------------------------
// CONTROL COMPONENTS
--------------------------------------------------------------------------------

COMPONENT BoolMux
PORTS
	IN bool_signal(n=1) s_in[2]
	OUT bool_signal(n=2) s_out
CONTINUOUS
SEQUENTIAL
	s_out.signal[1] = s_in[1].signal[1]
	s_out.signal[2] = s_in[2].signal[1]
/*
	s_out.signal[3] = s_in[3].signal[1]
	s_out.signal[4] = s_in[4].signal[1]
	s_out.signal[5] = s_in[5].signal[1]
	s_out.signal[5] = s_in[6].signal[1]
	s_out.signal[5] = s_in[7].signal[1]
*/
END SEQUENTIAL
END COMPONENT



COMPONENT BoolDemux (
	INTEGER n_out = 8 "Number of outputs"
)
PORTS
	IN bool_signal (n=n_out) s_in "Inlet port"
	OUT bool_signal (n=1) s_out[n_out] "Outlet ports"
TOPOLOGY
	PATH s_in TO s_out
CONTINUOUS
SEQUENTIAL
	s_out[1].signal[1] = s_in.signal[1]
	s_out[2].signal[1] = s_in.signal[2]
	s_out[3].signal[1] = s_in.signal[3]
	s_out[4].signal[1] = s_in.signal[4]
	s_out[5].signal[1] = s_in.signal[5]
	s_out[6].signal[1] = s_in.signal[6]
	s_out[7].signal[1] = s_in.signal[7]
	s_out[8].signal[1] = s_in.signal[8]
END SEQUENTIAL
END COMPONENT



COMPONENT Mux15 IS_A CONTROL.Mux
DECLS
	CLOSE n_in = 15
END COMPONENT



COMPONENT Demux15 IS_A CONTROL.Demux
DECLS
	CLOSE n_out = 15
END COMPONENT






--------------------------------------------------------------------------------
// PLANT PCO
--------------------------------------------------------------------------------
ENUM Plant_OptionMode = {S2PACL,S2PACL_StB,T2PACL_FastRecover}



COMPONENT DEMO_PCO_Plant (
	INTEGER nHead = 3 "Number of heads of plant"
)
"Process Control Object for the DEMO plant"
// Chiller is controlled by max(SHC5x40,DT5x40)
// TODO: Fix DT5x40 and SH5x40 setpoints → should start with SP:=PV
PORTS
	IN bool_signal(n=2) Bool_in
	IN analog_signal(n=15) Control_in
	OUT bool_signal(n=8) Bool_out
	OUT analog_signal(n=15) Control_out
	IN bool_signal(n=1) TakeOver_Rq "Used only for StB mode"
	IN bool_signal(n=1) SwapRequest "Used only for T2PACL mode"
DATA
	ENUM Plant_OptionMode OptionMode = S2PACL "Choose mode for startup"
	ENUM ChemName fluid = Carbon_Dioxide
	BOOLEAN DTC_Mode = FALSE "If TRUE, use chiller in Superheat control mode"
	REAL DPC1014_SPreq = 5 UNITS u_bar "Set point request for CV1014"
	REAL DPC1s14_SPreq = 10 UNITS u_bar "Set point request for pump startup"
	REAL DP2022_SPreq = 10 UNITS u_bar "Detector dP request"
	REAL emptyLL = 12 UNITS u_pct "Level at which to stop emptying, should be a bit higher than setPoint"
	REAL emptyHH = 40 UNITS u_pct "Level at which maxiumum emptying speed"
	REAL fillLL = 5 UNITS u_pct "Level low-limit for maximum filling speed"
	REAL posHH = 100 UNITS u_pct "Max opening (%) of surface storage valves"
	REAL posLL = 0 UNITS u_pct "Minimum opening (%) of surface storage valves"
	REAL setPoint = 10 UNITS u_pct "Set point for accumulator level"
	REAL SC9x98_TtL = 2 UNITS u_K "Subcooling threshold RESET for pump inlet"
	REAL SC9x98_TtH = 5 UNITS u_K "Subcooling threshold SET for pump inlet"
	REAL TrTime = 120 UNITS u_s "Minimum transition time between steps"
DECLS
	// Inputs
	BOOLEAN RunOrder
	BOOLEAN AccuReq
	INTEGER AccuStp
	REAL DP1014 UNITS u_bar "Pressure difference between plant outlet and surface storage basement"
	REAL DP1s14 UNITS u_bar "Pressure difference between plant outlet and back pressure regulator inlet"
	REAL DP1016 UNITS u_bar "Pressure difference between plant supply and return"
	REAL DP2022 UNITS u_bar "Pressure drop across the transfer line bypass valve" // controlled by CV1a16
	REAL LT4080 UNITS u_pct "Local accumulator liquid level"
	REAL PT1s18 UNITS u_bar "Surface storage pressure measured on plant skid"
	REAL PT9082 UNITS u_bar "Pressure at outlet of surface storage low pressure valve"
	REAL SC9x98 UNITS u_K "Pump inlet Subcooling"
	REAL CV1s14_PosSt UNITS u_pct "CV1s14 position"
	REAL TT1014 UNITS u_C "Temperature at inlet of CV1014"
	REAL P9_DP1016_SP UNITS u_bar "DP1016 of Plant 9, used in S2PACL_StB stepper"
	REAL P9_CV1a16_PosSt UNITS u_pct "Position of Plant 9's CV1a16"
	REAL ST4080 UNITS u_C "Accumulator saturation temperature"
	REAL DL_TTmax UNITS u_C "Maximum 'detector', aka dummy load, temperature"
	// Outputs (Controlled Objects)
	BOOLEAN DPC1014_TR_S
	BOOLEAN DPC1s14_TR_S
	BOOLEAN DPC1016_TR_S
	BOOLEAN DPC2022_TR_S
	BOOLEAN PC9082_TR_S
	BOOLEAN EV9s82
	BOOLEAN DTC5x40_TR_S
	BOOLEAN SHC5x40_TR_S
	REAL LP1002 UNITS u_Hz "Pump frequency request"
	REAL DPC1014_PosR UNITS u_pct
	REAL DPC1014_SP UNITS u_bar "Set point for CV1014"
	REAL DPC1s14_PosR UNITS u_pct
	REAL DPC1s14_SP UNITS u_bar
	REAL DPC1016_PosR UNITS u_pct
	REAL DPC1016_SP UNITS u_bar
	REAL DPC2022_PosR UNITS u_pct
	REAL DPC2022_SP UNITS u_bar
	REAL PC9082_PosR UNITS u_pct
	REAL PC9082_SP UNITS u_bar
	REAL DTC5x40_PosR UNITS u_pct
	REAL SHC5x40_PosR UNITS u_pct
	// States
	// BOOLEAN LP1002_SC_OK "Pump subcooling OK flag" // TODO Add superheat OK flag as output
	BOOLEAN LP1002_start
	// Transitions
	INTEGER Stp "Stepper position"
	INTEGER Stp_prev "Previous step"
	REAL timer UNITS u_s "Time counter"
	REAL time_prev UNITS u_s "Time at which to start counting"
	BOOLEAN T0
	BOOLEAN S2PACL_Tr[6], S2PACL_StB_Tr[5], T2PACL_Tr[5]
	// Local Variables
	REAL LP1002_p1 = 35 UNITS u_Hz "Nominal startup value for pump frequency"
	REAL CV1s14_Rq1
	REAL CV9s82_Rq1
	CONST REAL dP_ctrl = 0.1 UNITS u_bar "Deadband dP for checking stepper conditions"
	CONST REAL SC9x98_p1 = 6 UNITS u_K
	CONST REAL CV1014_p1 = 0 UNITS u_pct
	CONST REAL DP1s14_p1 = 8 UNITS u_bar
	CONST REAL CV1s14_p3 = 1 UNITS u_pct
	CONST REAL CV1a16_p1 = 0 UNITS u_pct
	CONST REAL CV1a16_p2 = 100 UNITS u_pct
	CONST REAL P1_DPC1016_P1 = 2 UNITS u_bar "dP for keeping StB restart pump dP lower than spare plant" // prevent accidental take-over
	CONST REAL CV9s82_p1 = 100 UNITS u_pct
	CONST REAL PC9082_p1 = 0.5 UNITS u_bar
	CONST REAL TT1014_p1 = -30 UNITS u_C 
	CONST REAL P9_CV1a16_pX = 95 UNITS u_pct "CV1a16 position on P9, for use in StB stepper"
	
	PRIVATE INTEGER ier,ipx,ipy // error codes
INIT
	timer = 0
	Stp = 0
CONTINUOUS
SEQUENTIAL
	--------------------------------------------------
	// DEMUX INPUT SIGNALS
	--------------------------------------------------
	RunOrder = Bool_in.signal[1]
	AccuReq = Bool_in.signal[2]
	AccuStp = Control_in.signal[1]
	CV1s14_PosSt = Control_in.signal[2]
	DP1014 = Control_in.signal[3]
	DP1s14 = Control_in.signal[4]
	DP1016 = Control_in.signal[5]
	DP2022 = Control_in.signal[6]
	LT4080 = Control_in.signal[7]
	PT1s18 = Control_in.signal[8]
	PT9082 = Control_in.signal[9]
	SC9x98 = Control_in.signal[10]
	TT1014 = Control_in.signal[11] // Celsius
	P9_DP1016_SP = Control_in.signal[12]
	P9_CV1a16_PosSt = Control_in.signal[13]
	ST4080 = Control_in.signal[14] // Celsius
	DL_TTmax = Control_in.signal[15] // Celsius
	
	--------------------------------------------------
	// VALCALC: CALCULATE VALUES
	--------------------------------------------------
	CV1s14_Rq1 = spliceFunction(posHH,posLL,LT4080-(emptyLL+(emptyHH-emptyLL)/2),(emptyHH-emptyLL)/2)
	CV9s82_Rq1 = spliceFunction(posLL,posHH,LT4080-(fillLL+(setPoint-fillLL)/2),(setPoint-fillLL)/2)
	LP1002_p1 = 35 // TODO: add logic

	--------------------------------------------------
	// STEPPER
	--------------------------------------------------

	// S2PACL
	----------------------------------------
	// In this mode we first equalise with surface storage
	// Then, chiller liquefies, pump turns on in local
	// circulation, and then finally we deliver liquid.
	IF (OptionMode==S2PACL) THEN
		IF (Stp==0) THEN
			LP1002 = 0
			// Valves
			DPC1014_TR_S = TRUE //TRS Assumes PLC drives???
			DPC1014_PosR = 100 // Normally Open (NO)
			DPC1014_SP = DP1014
			DPC1016_TR_S = TRUE
			DPC1016_PosR = 100 // NO
			DPC1016_SP = DP1016
			DPC2022_TR_S = TRUE
			DPC2022_PosR = 100 // NO
			DPC2022_SP = DP2022
			DPC1s14_TR_S = TRUE
			DPC1s14_PosR = 0 // Normally Closed (NC)
			DPC1s14_SP = DP1s14
			PC9082_TR_S = TRUE
			PC9082_PosR = 0 // NC
			EV9s82 = FALSE // NC
			// Chiller
			DTC5x40_TR_S = TRUE
			DTC5x40_PosR = 0 // NC
			SHC5x40_TR_S = TRUE
			SHC5x40_PosR = 0 // NC
		ELSEIF (Stp==1) THEN
		// Equalise with surface storage
			LP1002 = 0
			// Valves
			DPC1014_TR_S = TRUE
			DPC1014_PosR = 100
			DPC1014_SP = DP1014
			DPC1016_TR_S = TRUE
			DPC1016_PosR = CV1a16_p1
			DPC1016_SP = DP1016
			DPC2022_TR_S = TRUE
			DPC2022_PosR = CV1a16_p1
			DPC2022_SP = DP2022
			DPC1s14_TR_S = TRUE
			DPC1s14_PosR = 0
			DPC1s14_SP = DP1s14
			PC9082_TR_S = FALSE
			PC9082_PosR = 0 
			PC9082_SP = PT1s18
			EV9s82 = FALSE
			// Chiller
			DTC5x40_TR_S = TRUE
			DTC5x40_PosR = 0
			SHC5x40_TR_S = TRUE
			SHC5x40_PosR = 0
		ELSEIF (Stp==2) THEN
		// Liquefy pump inlet
			LP1002 = 0
			// Valves
			DPC1014_TR_S = TRUE
			DPC1014_PosR = 0
			DPC1014_SP = DP1014
			DPC1016_TR_S = TRUE
			DPC1016_PosR = CV1a16_p1
			DPC1016_SP = DP1016
			DPC2022_TR_S = TRUE
			DPC2022_PosR = CV1a16_p1
			DPC2022_SP = DP2022
			DPC1s14_TR_S = TRUE
			DPC1s14_PosR = 100
			DPC1s14_SP = DP1s14
			PC9082_TR_S = TRUE
			PC9082_PosR = CV9s82_p1
			EV9s82 = TRUE
			// Chiller
			DTC5x40_TR_S = TRUE
			DTC5x40_PosR = 0
			SHC5x40_TR_S = FALSE
			SHC5x40_PosR = 0
		ELSEIF (Stp==3) THEN
		// Pump run
			LP1002 = LP1002_p1
			// Valves
			DPC1014_TR_S = TRUE
			DPC1014_PosR = 0
			DPC1014_SP = DP1014
			DPC1016_TR_S = TRUE
			DPC1016_PosR = CV1a16_p1
			DPC1016_SP = DP1016
			DPC2022_TR_S = TRUE
			DPC2022_PosR = CV1a16_p1
			DPC2022_SP = DP2022
			IF (TT1014>TT1014_p1) THEN
				DPC1s14_TR_S = TRUE
				DPC1s14_PosR = 100
			ELSE
				DPC1s14_TR_S = FALSE
				DPC1s14_SP = DPC1s14_SPreq
			END IF
			PC9082_TR_S = TRUE
			PC9082_PosR = CV9s82_p1
			EV9s82 = TRUE
			// Chiller
			DTC5x40_TR_S = TRUE
			DTC5x40_PosR = 0
			SHC5x40_TR_S = FALSE
			SHC5x40_PosR = 0
		ELSEIF (Stp==4) THEN
		// Deliver liquid
			LP1002 = LP1002_p1
			// Valves
			DPC1014_TR_S = FALSE
			DPC1014_PosR = 0
			DPC1014_SP = DPC1014_SPreq
			DPC1016_TR_S = TRUE
			DPC1016_PosR = CV1a16_p1
			DPC1016_SP = DP1016
			DPC2022_TR_S = FALSE
			DPC2022_PosR = 0
			DPC2022_SP = DP2022_SPreq
			DPC1s14_TR_S = FALSE
			DPC1s14_PosR = 100
			DPC1s14_SP = DPC1s14_SPreq
			PC9082_PosR = 100
			PC9082_TR_S = TRUE
			EV9s82 = TRUE
			// Chiller
			DTC5x40_TR_S = TRUE
			DTC5x40_PosR = 0
			SHC5x40_TR_S = FALSE
			SHC5x40_PosR = 0
		ELSEIF (Stp==5) THEN
		// Operation S-2PACL
			LP1002 = LP1002_p1
			// Valves
			DPC1014_TR_S = FALSE
			DPC1014_PosR = 0
			DPC1014_SP = DPC1014_SPreq
			DPC1016_TR_S = TRUE
			DPC1016_PosR = CV1a16_p1
			DPC1016_SP = DP1016
			DPC2022_TR_S = FALSE
			DPC2022_PosR = 0
			DPC2022_SP = DP2022_SPreq
			DPC1s14_TR_S = TRUE
			DPC1s14_PosR = CV1s14_Rq1 // Accu emptying position request
			DPC1s14_SP = DP1s14
			PC9082_TR_S = TRUE
			PC9082_PosR = CV9s82_Rq1 // Accu filling position request
			PC9082_SP = PT9082
			EV9s82 = FALSE
			// Chiller
			DTC5x40_TR_S = TRUE
			DTC5x40_PosR = 0
			SHC5x40_TR_S = FALSE
			SHC5x40_PosR = 0	
		END IF
	
	// S2PACL IN STANBDY MODE
	----------------------------------------
	// In this mode, we are already connected to Accumulator
	// No need to equalise against surface storage.
	// Surface storage _is_ available in this mode, just not used.
	ELSEIF (OptionMode==S2PACL_StB) THEN
		IF (Stp==0) THEN
		// Safety position
			LP1002 = 0
			// Valves
			DPC1014_TR_S = TRUE
			DPC1014_PosR = 100 // Normally Open (NO)
			DPC1014_SP = DP1014
			DPC1016_TR_S = TRUE
			DPC1016_PosR = 100 // NO
			DPC1016_SP = DP1016
			DPC2022_TR_S = TRUE
			DPC2022_PosR = 100 // NO
			DPC2022_SP = DP2022
			DPC1s14_TR_S = TRUE
			DPC1s14_PosR = 0 // Normally Closed (NC)
			DPC1s14_SP = DP1s14
			PC9082_TR_S = TRUE
			PC9082_PosR = 0 // NC
			PC9082_SP = PT9082
			EV9s82 = FALSE // NC
			// Chiller
			DTC5x40_TR_S = TRUE
			DTC5x40_PosR = 0 // NC
			SHC5x40_TR_S = TRUE
			SHC5x40_PosR = 0 // NC
		ELSEIF (Stp==1) THEN
		// Liquefy
			LP1002 = 0
			// Valves
			DPC1014_TR_S = TRUE
			DPC1014_PosR = 100
			DPC1014_SP = DP1014
			DPC1016_TR_S = TRUE
			DPC1016_PosR = CV1a16_p1 // close
			DPC1016_SP = DP1016
			DPC2022_TR_S = TRUE
			DPC2022_PosR = CV1a16_p2
			DPC2022_SP = DP2022
			DPC1s14_TR_S = TRUE
			DPC1s14_PosR = 0
			DPC1s14_SP = DP1s14
			PC9082_TR_S = TRUE
			PC9082_PosR = 0
			PC9082_SP = PT9082
			EV9s82 = FALSE
			// Chiller
			DTC5x40_TR_S = TRUE
			DTC5x40_PosR = 0
			SHC5x40_TR_S = FALSE
			SHC5x40_PosR = 0
		ELSEIF (Stp==2) THEN
		// Pump run
			LP1002 = LP1002_p1
			// Valves
			DPC1014_TR_S = FALSE // not in tracking mode now
			DPC1014_PosR = 100
			DPC1014_SP = DPC1014_SPreq
			DPC1016_TR_S = FALSE // not in tracking mode now
			DPC1016_PosR = CV1a16_p2
			DPC1016_SP = P9_DP1016_SP - P1_DPC1016_P1 // stay below P9 DP to not accidentally take over
			DPC2022_TR_S = TRUE // still in tracking mode
			DPC2022_PosR = CV1a16_p1
			DPC2022_SP = DP2022
			DPC1s14_TR_S = TRUE // no s0 use
			DPC1s14_PosR = 0
			DPC1s14_SP = DP1s14
			PC9082_TR_S = TRUE
			PC9082_PosR = 0
			PC9082_SP = PT9082
			EV9s82 = FALSE
			// Chiller
			DTC5x40_TR_S = TRUE
			DTC5x40_PosR = 0
			SHC5x40_TR_S = FALSE
			SHC5x40_PosR = 0
		ELSEIF (Stp==3) THEN
		// Ready for substitution
			LP1002 = LP1002_p1
			// Valves
			DPC1014_TR_S = FALSE
			DPC1014_PosR = 100
			DPC1014_SP = DPC1014_SPreq
			DPC1016_TR_S = FALSE
			DPC1016_PosR = CV1a16_p1
			DPC1016_SP = P9_DP1016_SP // now match P9 DP set point // TODO actual value
			DPC2022_TR_S = TRUE
			DPC2022_PosR = CV1a16_p1
			DPC2022_SP = DP2022
			DPC1s14_TR_S = TRUE
			DPC1s14_PosR = 0
			DPC1s14_SP = DP1s14
			PC9082_TR_S = TRUE
			PC9082_PosR = 0
			PC9082_SP = PT9082
			EV9s82 = FALSE
			// Chiller
			DTC5x40_TR_S = TRUE
			DTC5x40_PosR = 0
			SHC5x40_TR_S = FALSE
			SHC5x40_PosR = 0
		ELSEIF(Stp==4) THEN
		// Take over
			// Pump
			LP1002 = LP1002_p1
			// Valves
			DPC1014_TR_S = FALSE
			DPC1014_PosR = 0
			DPC1014_SP = DPC1014_SPreq
			DPC1016_TR_S = TRUE
			DPC1016_PosR = CV1a16_p1 // = 0 (now close the controller)
			DPC1016_SP = DP1016
			DPC2022_TR_S = FALSE // start this one instead
			DPC2022_PosR = CV1a16_p1
			DPC2022_SP = DP2022_SPreq
			DPC1s14_TR_S = TRUE
			DPC1s14_PosR = 0
			DPC1s14_SP = DP1s14
			PC9082_TR_S = TRUE
			PC9082_PosR = 0
			PC9082_SP = PT9082
			EV9s82 = FALSE
			// Chiller
			DTC5x40_TR_S = TRUE
			DTC5x40_PosR = 0
			SHC5x40_TR_S = FALSE
			SHC5x40_PosR = 0
		ELSEIF (Stp==5) THEN
		// S-2PACL operation
			// Unable to change ENUMs in EcosimPro, so have to copy this from S2PACL
			LP1002 = LP1002_p1
			// Valves
			DPC1014_TR_S = FALSE
			DPC1014_PosR = 0
			DPC1014_SP = DPC1014_SPreq
			DPC1016_TR_S = TRUE
			DPC1016_PosR = CV1a16_p1
			DPC1016_SP = DP1016
			DPC2022_TR_S = FALSE
			DPC2022_PosR = 0
			DPC2022_SP = DP2022_SPreq
			DPC1s14_TR_S = TRUE
			DPC1s14_PosR = CV1s14_Rq1 // Accu emptying position request
			DPC1s14_SP = DP1s14
			PC9082_TR_S = TRUE
			PC9082_PosR = CV9s82_Rq1 // Accu filling position request
			PC9082_SP = PT9082
			EV9s82 = FALSE
			// Chiller
			DTC5x40_TR_S = TRUE
			DTC5x40_PosR = 0
			SHC5x40_TR_S = FALSE
			SHC5x40_PosR = 0
		END IF
	
	// Traditional 2PACL startup
	----------------------------------------
	// Same old same old. Pressurise to liquefy
	// Then cool down and start up. Assumes no
	// Surface storage. S0 valves are untouched
	ELSEIF (OptionMode==T2PACL_FastRecover) THEN
		IF (Stp==0) THEN
			LP1002 = 0
			// Valves
			DPC1014_TR_S = TRUE
			DPC1014_PosR = 100 // Normally Open (NO)
			DPC1014_SP = DP1014
			DPC1016_TR_S = TRUE
			DPC1016_PosR = 100 // NO
			DPC1016_SP = DP1016
			DPC2022_TR_S = TRUE
			DPC2022_PosR = 100 // NO
			DPC2022_SP = DP2022
			DPC1s14_TR_S = TRUE
			DPC1s14_PosR = 0 // Normally Closed (NC)
			DPC1s14_SP = DP1s14
			PC9082_TR_S = TRUE
			PC9082_PosR = 0 // NC
			PC9082_SP = PT9082
			EV9s82 = FALSE // NC
			// Chiller
			DTC5x40_TR_S = TRUE
			DTC5x40_PosR = 0 // NC
			SHC5x40_TR_S = TRUE
			SHC5x40_PosR = 0 // NC
		ELSEIF (Stp==1) THEN
		// 31: Liquefy
			LP1002 = 0
			// Valves
			DPC1014_TR_S = TRUE
			DPC1014_PosR = 100
			DPC1014_SP = DP1014
			DPC1016_TR_S = TRUE
			DPC1016_PosR = CV1a16_p1 // close
			DPC1016_SP = DP1016
			DPC2022_TR_S = TRUE
			DPC2022_PosR = CV1a16_p2 // open
			DPC2022_SP = DP2022
			DPC1s14_TR_S = TRUE
			DPC1s14_PosR = 0
			DPC1s14_SP = DP1s14
			PC9082_TR_S = TRUE
			PC9082_PosR = 0
			PC9082_SP = PT9082
			EV9s82 = FALSE
			// Chiller
			DTC5x40_TR_S = TRUE
			DTC5x40_PosR = 0
			SHC5x40_TR_S = FALSE
			SHC5x40_PosR = 0
		ELSEIF (Stp==2) THEN
		// 32: Pump run
			LP1002 = LP1002_p1
			// Valves
			DPC1014_TR_S = TRUE
			DPC1014_PosR = 100
			DPC1014_SP = DP1014
			DPC1016_TR_S = TRUE
			DPC1016_PosR = CV1a16_p1 // close
			DPC1016_SP = DP1016
			DPC2022_TR_S = FALSE
			IF (AccuStp<4) THEN
				// Detector pressurisation
				DPC2022_SP = 1 // bar
			ELSE
				DPC2022_SP = DP2022_SPreq
			END IF
			DPC1s14_TR_S = TRUE
			DPC1s14_PosR = 0
			DPC1s14_SP = DP1s14
			PC9082_TR_S = TRUE
			PC9082_PosR = 0
			PC9082_SP = PT9082
			EV9s82 = FALSE
			// Chiller
			DTC5x40_TR_S = TRUE
			DTC5x40_PosR = 0
			SHC5x40_TR_S = FALSE
			SHC5x40_PosR = 0
		ELSEIF (Stp==3) THEN
		// 33: Operation 2PACL/AC4080
		// (In this step we can swap to DT control mode)
			LP1002 = LP1002_p1
			// Valves
			DPC1014_TR_S = TRUE
			DPC1014_PosR = 100
			DPC1014_SP = DP1014
			DPC1016_TR_S = TRUE
			DPC1016_PosR = CV1a16_p1 // close
			DPC1016_SP = DP1016
			DPC2022_TR_S = FALSE
			IF (AccuStp<4) THEN
				// Detector pressurisation
				DPC2022_SP = 1 // bar
			ELSE
				DPC2022_SP = DP2022_SPreq
			END IF
			DPC2022_PosR = DP2022
			DPC1s14_TR_S = TRUE
			DPC1s14_PosR = 100
			DPC1s14_SP = DP1s14
			PC9082_TR_S = TRUE
			PC9082_PosR = 0
			PC9082_SP = PT9082
			EV9s82 = FALSE
			// Chiller
			DTC5x40_TR_S = TRUE
			DTC5x40_PosR = 0
			SHC5x40_TR_S = FALSE
			SHC5x40_PosR = 0
		ELSEIF (Stp==4) THEN
		// 34: Preparation to swap
			LP1002 = LP1002_p1
			// Valves
			DPC1014_TR_S = FALSE
			DPC1014_PosR = 100
			DPC1014_SP = DPC1014_SPreq
			DPC1016_TR_S = TRUE
			DPC1016_PosR = CV1a16_p1 // close
			DPC1016_SP = DP1016
			DPC2022_TR_S = FALSE
			IF (AccuStp<4) THEN
				// Detector pressurisation
				DPC2022_SP = 1 // bar
			ELSE
				DPC2022_SP = DP2022_SPreq
			END IF
			DPC2022_PosR = DP2022
			DPC1s14_TR_S = TRUE
			DPC1s14_PosR = 0
			DPC1s14_SP = DP1s14
			PC9082_TR_S = TRUE
			PC9082_PosR = 0
			PC9082_SP = PT9082
			EV9s82 = FALSE
			// Chiller
			DTC5x40_TR_S = TRUE
			DTC5x40_PosR = 0
			SHC5x40_TR_S = FALSE
			SHC5x40_PosR = 0
		END IF
	END IF
		
	--------------------------------------------------
	// TRANSITION LOGIC
	--------------------------------------------------
	timer = TIME - time_prev // Reset everytime step changes
	T0 = NOT RunOrder
	// S2PACL	
	S2PACL_Tr[1] = OptionMode==S2PACL AND RunOrder AND AccuStp>=1
	S2PACL_Tr[2] = PT9082<=PT1s18+dP_ctrl AND PT9082>=PT1s18-dP_ctrl AND timer>=TrTime
	S2PACL_Tr[3] = SC9x98>=SC9x98_p1 AND timer>=TrTime
	S2PACL_Tr[4] = DP1s14>=DP1s14_p1 AND timer>=TrTime
	S2PACL_Tr[5] = CV1s14_PosSt<=CV1s14_p3 AND AccuStp>=4 AND AccuStp<100 // TODO: fix AccuStp<100
	
	// S2PACL_StB
	S2PACL_StB_Tr[1] = OptionMode==S2PACL_StB AND RunOrder AND AccuStp>=1
	S2PACL_StB_Tr[2] = SC9x98>=SC9x98_p1 AND timer>=TrTime
	S2PACL_StB_Tr[3] = TakeOver_Rq.signal[1] AND timer>=TrTime
	S2PACL_StB_Tr[4] = P9_CV1a16_PosSt>=P9_CV1a16_pX AND timer>=TrTime
	S2PACL_StB_Tr[5] = TakeOver_Rq.signal[1]
	
	// T2PACL
	T2PACL_Tr[1] = OptionMode==T2PACL_FastRecover AND RunOrder AND AccuStp>=1
	T2PACL_Tr[2] = SC9x98>=SC9x98_p1 AND ST4080-DL_TTmax>=-1 AND timer>=TrTime
	T2PACL_Tr[3] = AccuStp>=4 AND timer>=TrTime
	T2PACL_Tr[4] = SwapRequest.signal[1] AND timer>=TrTime
	
	// Step increment
	Stp_prev = Stp
	IF (OptionMode==S2PACL) THEN
		IF (Stp<5 AND S2PACL_Tr[Stp+1]) THEN
			Stp = Stp+1
			timer = 0
			time_prev = TIME
		END IF
	ELSEIF (OptionMode==S2PACL_StB) THEN
		IF (Stp<5 AND S2PACL_StB_Tr[Stp+1]) THEN
			Stp = Stp+1
			timer = 0
			time_prev = TIME
		END IF
	ELSEIF (OptionMode==T2PACL_FastRecover) THEN
		IF (Stp<4 AND T2PACL_Tr[Stp+1]) THEN
			Stp = Stp+1
			timer = 0
			time_prev = TIME
		END IF
	
	END IF
	
	--------------------------------------------------
	// DL: DEPENDENT LOGIC
	--------------------------------------------------
	// Boolean signals
	Bool_out.signal[1] = DPC1014_TR_S
	Bool_out.signal[2] = DPC1s14_TR_S
	Bool_out.signal[3] = DPC1016_TR_S
	Bool_out.signal[4] = DPC2022_TR_S
	Bool_out.signal[5] = PC9082_TR_S
	Bool_out.signal[6] = EV9s82
	Bool_out.signal[7] = DTC5x40_TR_S
	Bool_out.signal[8] = SHC5x40_TR_S
	// Output Signals
	Control_out.signal[1] = Stp
	Control_out.signal[2] = LP1002
	Control_out.signal[3] = DPC1014_PosR
	Control_out.signal[4] = DPC1014_SP
	Control_out.signal[5] = DPC1s14_PosR
	Control_out.signal[6] = DPC1s14_SP
	Control_out.signal[7] = DPC1016_PosR
	Control_out.signal[8] = DPC1016_SP
	Control_out.signal[9] = DPC2022_PosR
	Control_out.signal[10] = DPC2022_SP
	Control_out.signal[11] = PC9082_PosR
	Control_out.signal[12] = PC9082_SP
	Control_out.signal[13] = DTC5x40_PosR
	Control_out.signal[14] = SHC5x40_PosR
END SEQUENTIAL
END COMPONENT
