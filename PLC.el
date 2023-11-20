/*-----------------------------------------------------------------------------------------
 LIBRARY: CO2
 FILE: PLC
 AUTHOR: Viren Bhanot
 COMPANY: CERN
 DESCRIPTION: PLC Control Components
 CREATION DATE: 04/07/2018
-----------------------------------------------------------------------------------------*/
USE PORTS_LIB
USE CRYOLIB
USE CONTROL
USE MATH
USE PLC



COMPONENT ANALOG_varLim
"ANALOG controller with variable upper and lower limits"
PORTS
	IN analog_signal(n=1) PosR "Auto Position Request (%)"
	OUT analog_signal(n=1) feedback "feedback position (%)"
	IN analog_signal(n=1) Out_Min "Minimum Output"
	IN analog_signal(n=1) Out_Max "Maximum Output"
DATA
	REAL FoV = 0.0 "Forced value (%)"
	BOOLEAN FoMoSt = FALSE "Forced Mode Set"
	REAL InSpd=20.0	"Increase Speed (%/s)"
	REAL DeSpd=20.0	"Decrease Speed (%/s)"
	REAL PreConstraint=0.0 "Pre constraint on stem position (%)"
	REAL Out_init = 0.0 "Initial position (%)"	
DECLS
	REAL PMinRan "Minimum Output"
	REAL PMaxRan "Maximum Output"
	REAL s_stem "Stem position(%)"
	REAL Pos "Position (%)"
	REAL spd "Current speed (%/s)"	
TOPOLOGY
	PATH PosR TO feedback
INIT
	IF (FoMoSt == TRUE) THEN
		Pos = FoV
	ELSE
		Pos = Out_init
	END IF
	s_stem = Out_init - PreConstraint
	feedback.signal[1] = Out_init
CONTINUOUS
	PMinRan = Out_Min.signal[1]
	PMaxRan = Out_Max.signal[1]
	Pos = ZONE(FoMoSt == TRUE) FoV OTHERS min(max(PMinRan,PosR.signal[1]),PMaxRan)
	feedback.signal[1] = s_stem + PreConstraint
	SEQUENTIAL
		IF(s_stem<= Pos - PreConstraint - InSpd*CINT AND s_stem<PMaxRan) THEN
			spd = InSpd 
		ELSEIF(s_stem<Pos-PreConstraint AND s_stem> Pos-PreConstraint- InSpd*CINT AND s_stem<PMaxRan) THEN
			spd = (Pos-PreConstraint-s_stem)*10/CINT
		ELSEIF(s_stem>= Pos-PreConstraint + DeSpd*CINT AND s_stem > PMinRan)THEN
			spd = -DeSpd 
		ELSEIF(s_stem>Pos-PreConstraint AND s_stem< Pos-PreConstraint + DeSpd*CINT AND s_stem > PMinRan) THEN
			spd = (Pos-PreConstraint-s_stem)*10/CINT
		ELSE
			spd = 0
		END IF	
	END SEQUENTIAL
	s_stem'= spd
END COMPONENT



COMPONENT ANALOG_varLimSpd
"ANALOG controller with variable upper and lower limits"
PORTS
	IN analog_signal(n=1) PosR "Auto Position Request (%)"
	OUT analog_signal(n=1) feedback "feedback position (%)"
	IN analog_signal(n=1) Out_Min "Minimum Output"
	IN analog_signal(n=1) Out_Max "Maximum Output"
	IN analog_signal(n=1) Spd "Speed Signal"
DATA
	REAL  FoV = 0.0 "Forced value (%)"
	BOOLEAN FoMoSt = FALSE "Forced Mode Set"
	REAL PreConstraint=0.0 "Pre constraint on stem position (%)"
	REAL Out_init = 0.0 "Initial position (%)"	
DECLS
	REAL InSpd "Increase Speed (%/s)"
	REAL DeSpd "Decrease Speed (%/s)"
	REAL PMinRan "Minimum Output"
	REAL PMaxRan "Maximum Output"
	REAL s_stem "Stem position(%)"
	REAL Pos "Position (%)"
	REAL spd "Current speed (%/s)"	
TOPOLOGY
	PATH PosR TO feedback // defines default path from PosR to feedback
INIT
	IF (FoMoSt == TRUE) THEN
		Pos = FoV
	ELSE
		Pos = Out_init 
	END IF
	InSpd = Spd.signal[1]
	DeSpd = Spd.signal[1]
	s_stem = Out_init - PreConstraint
	feedback.signal[1] = Out_init
CONTINUOUS
	InSpd = Spd.signal[1]
	DeSpd = Spd.signal[1]
	PMinRan = Out_Min.signal[1]
	PMaxRan = Out_Max.signal[1]
	Pos = ZONE(FoMoSt == TRUE) FoV OTHERS min(PosR.signal[1],PMaxRan)
	feedback.signal[1] = s_stem + PreConstraint
	SEQUENTIAL
		IF(s_stem<= Pos - PreConstraint - InSpd*CINT AND s_stem<PMaxRan) THEN
			spd = InSpd 
		ELSEIF(s_stem<Pos-PreConstraint AND s_stem> Pos-PreConstraint- InSpd*CINT AND s_stem<PMaxRan) THEN
			spd = (Pos-PreConstraint-s_stem)*10
		ELSEIF(s_stem>= Pos-PreConstraint + DeSpd*CINT AND s_stem > PMinRan)THEN
			spd = -DeSpd 
		ELSEIF(s_stem>Pos-PreConstraint AND s_stem< Pos-PreConstraint + DeSpd*CINT AND s_stem > PMinRan) THEN
			spd = (Pos-PreConstraint-s_stem)*10
		ELSE
			spd = 0
		END IF	
	END SEQUENTIAL
	s_stem'= spd
END COMPONENT



/*
 * NOTES ON PLC CONTROL OF CORA:
 * In CORA, the run order is a boolean. There is no longer any Stepper
 * The run order is either TRUE or FALSE. There are, however, transition phases.
 * They are as follows:
 *	- T0: No Run Order
 *	- T1: Run Order
 *	- T2: Run Order AND AC4082_Tsp <= ST4082
 *	- T3: Run Order AND TT1092 < -30°C AND LP1002 is ON
 */
COMPONENT CORA_PCO_CO2
"Process Control Object for CORA. Includes all the CO2 control logic"
PORTS
	IN bool_signal(n=1) Run_Order "Run order from Master"
	IN analog_signal(n=5) Control_in
	OUT analog_signal(n=5) Control_out
DATA
	BOOLEAN TurnOnAtSP = FALSE "If TRUE, heater power only turns on once Usp obtained"
DECLS
	// Inputs
	REAL AC4082_Usp UNITS u_C "Accumulator user set point"
	REAL TT1092 UNITS u_C "Pump inlet temperature"
	REAL ST1092 UNITS u_C "Pump inlet saturated temperature"
	REAL ST4082 UNITS u_C "Accumulator saturated temperature"
	// Outputs
	REAL AC4082_Tsp UNITS u_C "Accumulator temperature set point"
	REAL EH1a52 "Heater power signal in % (0-100)"
	REAL EH4082_OutMax "Accumulator heater power high limit"
	REAL LP1002 "Pump speed signal"
	// States
	REAL AC4082_Asp UNITS u_C "Accumulator auto set point"
	REAL SC1092
	BOOLEAN LP1002_start "Pump start signal"
	BOOLEAN LP1002_SC_OK "Flag to check if subcooling at pump inlet is OK"
	BOOLEAN LP1002_SC_ST "Flag to check if subcooling at pump inlet is OK"
	// Transitions
	BOOLEAN Tr0 "Stand by/Safety position"
	BOOLEAN Tr1 "Start position"
	BOOLEAN Tr2 "Chiller start"
	BOOLEAN Tr3 "Accumulator cool down"
	REAL Phase "Sequencer state"
	// Local Variables
	CONST REAL LP1002_Tt = 0 UNITS u_K "Subcooling threshold for pump inlet"
	CONST REAL deadband = 0.5 UNITS u_K
INIT
	Phase = 0
	LP1002_SC_OK = FALSE
	LP1002_SC_ST = FALSE
CONTINUOUS
SEQUENTIAL
	// DEMUX INPUT SIGNALS
	Phase = Control_in.signal[1] // read old phase first before updating
	AC4082_Usp = Control_in.signal[2]
	ST4082 = Control_in.signal[3]
	ST1092 = Control_in.signal[4]
	TT1092 = Control_in.signal[5]
	
	// VALCALC: CALCULATE VALUES
	SC1092 = ST1092 - TT1092 // Positive if actual temperature lower than saturation temperature
	// accumulator auto set point
	IF (TT1092+6 >= 27) THEN
		AC4082_Asp = 27
	ELSE
		AC4082_Asp = TT1092 + 6
	END IF
	// accumulator set point
	AC4082_Tsp = max(AC4082_Asp, AC4082_Usp)

	// accumulator heater power high limit
	IF (ST4082<=17) THEN
		EH4082_OutMax = 100
	ELSEIF (ST4082<=19) THEN
		EH4082_OutMax = 75
	ELSEIF (ST4082<=22) THEN
		EH4082_OutMax = 50
	ELSE
		EH4082_OutMax = 25
	END IF
	EH4082_OutMax = 100
	// pump subcooling Set-Reset logic
	IF (LP1002_SC_ST==TRUE) THEN // when sc was ok before (>6)
		IF (SC1092 > LP1002_Tt) THEN // LP1002_Tt is enough to keep OK status 
			LP1002_SC_OK = TRUE
			LP1002_SC_ST = TRUE
		ELSE // sc < LP1002_Tt
			LP1002_SC_ST = FALSE
			LP1002_SC_OK = FALSE
		END IF
	ELSE
		IF (SC1092 > LP1002_Tt + 6) THEN
			LP1002_SC_ST = TRUE
			LP1002_SC_OK = TRUE
		ELSE
			LP1002_SC_ST = FALSE
			LP1002_SC_OK = FALSE
		END IF
	END IF
	// pump start logic
	IF (Tr2 OR Tr3) THEN
		IF (LP1002_SC_OK) THEN
			LP1002_start = TRUE
			LP1002 = 100
		ELSE
			LP1002_start = FALSE
			LP1002 = 0
		END IF
	ELSEIF (Tr0 OR Tr1) THEN
		LP1002_start = FALSE
	ELSE
		LP1002_start = FALSE // default
		LP1002 = 0
	END IF
	// heater output signal
	IF (Run_Order.signal[1] AND LP1002_start) THEN
		IF (TurnOnAtSP) THEN
			IF (ST4082>=AC4082_Usp-deadband AND ST4082<=AC4082_Usp+deadband) THEN
				EH1a52 = 100
			END IF
		ELSEIF (TIME>25000) THEN
			EH1a52 = 150
		ELSE
			EH1a52 = 100
		END IF
	ELSE
		EH1a52 = 0
	END IF
	
	// TL: TRANSITION LOGIC
	Tr3 = LogicalAND(Run_Order.signal[1],LP1002_start)
	Tr2 = LogicalAND(Run_Order.signal[1],AC4082_Tsp <= ST4082)
	Tr1 = LogicalAND(Run_Order.signal[1],LP1002_start)
	Tr0 = LogicalNOT(Run_Order.signal[1])
	Phase = TransitionLogic(Tr0,Tr1,Tr2,Tr3)
	
	// DL: DEPENDENT LOGIC
	//Mux Output Signals
	Control_out.signal[1] = Phase // Sequencer phase
	Control_out.signal[2] = AC4082_Tsp
	Control_out.signal[3] = EH4082_OutMax
	Control_out.signal[4] = EH1a52
	Control_out.signal[5] = LP1002
END SEQUENTIAL
END COMPONENT



COMPONENT CORA_PCO_Chiller
"Process control object for CORA chiller. Contains R404A control logic"
PORTS
	IN bool_signal(n=1) Run_Order
	IN analog_signal(n=1) Control_in
	OUT analog_signal(n=1) Control_out
DECLS
	CONST REAL delta = 0.01
	// Inputs
	REAL Phase
	// Outputs
	REAL GP5002
CONTINUOUS
SEQUENTIAL
	Phase = Control_in.signal[1]
	IF TIME >=10000 THEN
		GP5002 = 16
	ELSEIF Run_Order.signal[1] AND Phase >= 2 THEN
		GP5002 = 100
	ELSE
		GP5002 = 0
	END IF
END SEQUENTIAL
	Control_out.signal[1] = GP5002
END COMPONENT



COMPONENT CORA_PCO_Chiller_pid
"CORA Chiller PCO for PID compressor control"
PORTS
	IN bool_signal(n=1) Run_Order
	IN analog_signal(n=1) Phase
	OUT bool_signal(n=1) TR_S
	OUT analog_signal(n=1) PosR
	OUT analog_signal(n=1) SP
CONTINUOUS
SEQUENTIAL
	PosR.signal[1] = 0
	SP.signal[1] = 1 // 1 bar set point (hard coded)
	IF (Phase.signal[1] >= 2 AND Run_Order.signal[1]) THEN
		TR_S.signal[1] = FALSE
	ELSE
		TR_S.signal[1] = TRUE
	END IF
END SEQUENTIAL
END COMPONENT



COMPONENT LevelController
"Controls level of accumulator using Set-Reset type logic"
/* When emptying, the rate is limited by how much we allow the surface storage
 * valve to steal from the bypass valve. Some flow must always be allowed through
 * the bypass valve, otherwise we lose ability to control detector flow rate
 *
 * When filling the plant, we have the highest flow rate when the accumulator
 * level is the lowest. Then, as we get closer to the cut-off point, we ease
 * up on the opening level. A linear equation is used for this.
 */
PORTS
	IN analog_signal(n=1) Phase
	IN analog_signal(n=1) level
	IN analog_signal(n=1) bypassFlow
	
	OUT analog_signal(n=1) Dis
	OUT analog_signal(n=1) Suc
	OUT analog_signal(n=1) Bypass
DATA
	REAL lowLimit = 10 "Level below which to start filling accu"
	REAL highLimit = 80 "Upper limit at which to start emptying accu"
	REAL cutOff = 25 "Cut off at which to stop filling accumulator"
	REAL lowFlow = 50 "Minimum value for bypass flow" // turn valve off at this point
	REAL highFlow = 150 "Value at which to start limiting ss flow"
	REAL highPos = 100 "High position for valves"
	REAL lowPos = 0 "Low pasition for valves"
DECLS
	BOOLEAN fill
	BOOLEAN empty
	DISCR REAL epsilon = 5
	REAL emptyPos "Position of surface storage high pressure valve"
	REAL fillPos "Position of surface storage low pressure valve"
DISCRETE
	WHEN (level.signal[1]<lowLimit) THEN
		fill = TRUE
	END WHEN
	WHEN (level.signal[1]>cutOff) THEN
		fill = FALSE
	END WHEN
	WHEN (level.signal[1]>highLimit) THEN
		empty = TRUE
	END WHEN
	WHEN (level.signal[1]<lowLimit+epsilon) THEN
		empty = FALSE
	END WHEN
CONTINUOUS
	// make sure we don't steal too much flow from bypass:
	// (see Test_RateLimitingAlgo.m matlab file for plot)
	emptyPos = spliceFunction(highPos,lowPos,bypassFlow.signal[1]-(highFlow+lowFlow)/2,(highFlow-lowFlow)/2)
	fillPos = spliceFunction(lowPos,highPos,level.signal[1]-cutOff/2,cutOff/2)
	SEQUENTIAL
		IF (Phase.signal[1]==0) THEN
			Suc.signal[1] = 0
			Dis.signal[1] = 0
			Bypass.signal[1] = 0
		ELSEIF (Phase.signal[1]>=1 AND Phase.signal[1]<=3) THEN
			Suc.signal[1] = 100
			Dis.signal[1] = 100
			Bypass.signal[1] = 100
		ELSEIF(Phase.signal[1]==4 OR Phase.signal[1]==5) THEN
			Dis.signal[1] = 50
		ELSEIF(Phase.signal[1]>5) THEN
			// Enable local accumulator level control using surface storage
			Bypass.signal[1] = 0
			IF (fill) THEN
				Dis.signal[1] = 0
				Suc.signal[1] = fillPos
			ELSEIF (empty) THEN
				Dis.signal[1] = emptyPos
				Suc.signal[1] = 0
			ELSE
				Suc.signal[1] = 0
				Dis.signal[1] = 0
			END IF
		END IF
	END SEQUENTIAL
END COMPONENT



COMPONENT LevelController_splice
"Level control without boolean on/offs. Instead, continuous control to maintain constant level"
// also removes the phase input
PORTS
	IN analog_signal(n=1) level "Accumulator level"
	IN analog_signal(n=1) bypassFlow "Flow through detector bypass valve"
	OUT analog_signal(n=1) emptyAccu "Signal to low-pressure-side valve"
	OUT analog_signal(n=1) fillAccu "Signal to high-pressure-side valve"
DATA
	REAL fillLL = 5 UNITS u_pct "Level low-limit for maximum filling speed"
	REAL setPoint = 10 UNITS u_pct "Level set point"
	REAL emptyLL = 12 UNITS u_pct "Level at which to stop emptying, should be a bit higher than setPoint"
	REAL emptyHH = 40 UNITS u_pct "Level at which maxiumum emptying speed"
	REAL bypassHH= 0.15 UNITS u_kg_s "Value of bypass flow at which to start limiting surface storage valve"
	REAL bypassLL = 0.05 UNITS u_kg_s "Minimum value for bypass flow" // turn valve off at this point
	REAL posHH = 100 UNITS u_pct "Max opening (%) of valves"
	REAL posLL = 0 UNITS u_pct "Minimum opening (%) of valves"
DECLS
	REAL ePos "Position of plant emptying valve"
	REAL fPos "Position of plant filling valve"
	REAL maxPos "Maximum allowed opening for emptying valve based on bypass flow limits"
	REAL lvl, bp // variables to aid code readability
CONTINUOUS
	// make equations more readable:
	lvl = level.signal[1]
	bp = bypassFlow.signal[1]
	// make sure we don't steal too much flow from bypass:
	maxPos = spliceFunction(posHH,posLL,bp-(bypassHH+bypassLL)/2,(bypassHH-bypassLL)/2)
	// interpolate valve positions:
	ePos = spliceFunction(maxPos,posLL,lvl-(emptyLL+(emptyHH-emptyLL)/2),(emptyHH-emptyLL)/2)
	fPos = spliceFunction(posLL,posHH,lvl-(fillLL+(setPoint-fillLL)/2),(setPoint-fillLL)/2)
	// assign signals:
	ePos = emptyAccu.signal[1]
	fPos = fillAccu.signal[1]
END COMPONENT



--------------------------------------------------------------------------------
// OBSOLETE COMPONENTS
--------------------------------------------------------------------------------
COMPONENT DEMO_PCO_CO2
"Process Control Object for the DEMO plant [OBSOLETE]"
// Accumulator and its related components have a separate PCO
PORTS
	IN bool_signal(n=1) Run_Order
	IN analog_signal(n=20) Control_in
	OUT analog_signal(n=14) Control_out
	OUT bool_signal(n=5) Bool_out
DATA
	REAL TrTime = 120 UNITS u_s "Minimum transition time between steps"
DECLS
	// Inputs
	REAL TT3004 "Actual temperature at detector inlet"
	REAL ST3004 "Saturation temperature at detector inlet"
	REAL CV7c08_feedback "Actual position of valve"
	REAL WT4080 "Local accumulator liquid level"
	REAL PT6030 "Pressure at the base of surface storage"
	REAL PT4c80 "Pressure at outlet of surface storage low pressure valve"
	REAL SC1a98 "Subcooling at pump inlet"
	REAL DP7a10 "Pressure drop across detector volume"
	REAL ST4080 "Accumulator saturation temperature"
	REAL AC4080_Usp "User set point for accumulator"
	REAL PT7008 "Pump outlet pressure"
	REAL FT1a16 "Detector bypass mass flow rate"
	
	// Outputs
	REAL LP1a02 "Pump speed"
	REAL GP5002 "Compressor speed"
	REAL LC4080_PV "Accumulator level controller PID Process Value"
	REAL LC4080_SP "Accumulator level controller PID Set Point"
	REAL LC4080_AtOutL "Accumulator level controller output lower limit"
	REAL LC4080_AtOutH "Accumulator level controller output higher limit"
	BOOLEAN CV7c08_TR_S "Surface storage high pressure valve"
	REAL CV7c08_PosR
	BOOLEAN CV4c80_TR_S
	REAL CV4c80_PosR "Surface storage low pressure valve"
	REAL CV7008_SP "Pump outlet pressure control valve"
	BOOLEAN CV7008_TR_S
	REAL CV7008_PosR
	REAL CV7a10_SP
	BOOLEAN CV7a10_TR_S "Detector flow rate control valve"
	REAL CV7a10_PosR
	BOOLEAN PV7c10
	REAL EH3006 "Heater power"
	
	// States
	BOOLEAN LP1a02_SC_OK "Pump subcooling OK flag"
	BOOLEAN LP1a02_start
	REAL SC3004 "Subcooling at detector inlet"
	REAL DP7c08 "Pressure drop across SS high-pressure valve"
	REAL ST6010
	
	// Transitions
	INTEGER Stp "Stepper position"
	INTEGER Stp_prev "Previous step"
	REAL timer UNITS u_s "Time counter"
	REAL time_prev UNITS u_s "Time at which to start counting"
	BOOLEAN Tr0, Tr1, Tr2, Tr3, Tr4, Tr5, Tr6, Tr7, Tr8 // conditions
	BOOLEAN Tr[8]
	
	// Local Variables
	CONST REAL CV7a10_SPval = 12 UNITS u_bar "Set point value for detector bypass valve"
	CONST REAL SC1a98_TtL = 2 UNITS u_K "Subcooling threshold RESET for pump inlet"
	CONST REAL SC1a98_TtH = 5 UNITS u_K "Subcooling threshold SET for pump inlet"
	CONST REAL dP_ctrl = 0.1 UNITS u_bar "Deadband dP for checking stepper conditions"
	CONST REAL epsilon = 2 "Deadband"
	CONST ENUM ChemName fluid = Carbon_Dioxide
	PRIVATE INTEGER ier,ipx,ipy // error codes
	
	// Level Controller-related
	CONST REAL fillLL = 5 UNITS u_pct "Level low-limit for maximum filling speed"
	CONST REAL setPoint = 10 UNITS u_pct "Set point for accumulator level"
	DISCR REAL emptyLL = 12 UNITS u_pct "Level at which to stop emptying, should be a bit higher than setPoint"
	CONST REAL emptyHH = 40 UNITS u_pct "Level at which maxiumum emptying speed"
	CONST REAL bypassHH = 0.1 UNITS u_kg_s "Value of bypass flow at which to start limiting surface storage valve"
	CONST REAL bypassLL = 0.01 UNITS u_kg_s "Minimum value for bypass flow" // turn valve off at this point
	CONST REAL posHH = 100 UNITS u_pct "Max opening (%) of surface storage valves"
	CONST REAL posLL = 0 UNITS u_pct "Minimum opening (%) of surface storage valves"
	REAL maxPos "limiter on HH of emptying valve" // ensures we don't steal too much from bypass flow
	REAL ePos "Emptying valve position"
	REAL fPos "Filling valve position"
	
INIT
	Stp = 0
	timer = 0
	LP1a02_SC_OK = FALSE
	CV7c08_feedback = 0

DISCRETE
	// Safety position
	WHEN (Tr0) THEN
		Stp = 0
	END WHEN
	// Pump subcooling
	WHEN(SC1a98 > SC1a98_TtH) THEN
		LP1a02_SC_OK = TRUE
	END WHEN
	WHEN(SC1a98 < SC1a98_TtL) THEN
		LP1a02_SC_OK = FALSE
	END WHEN

CONTINUOUS
SEQUENTIAL
	//DEMUX INPUT SIGNALS
	TT3004 = Control_in.signal[1]
	ST3004 = Control_in.signal[2]
	WT4080 = Control_in.signal[3]
	PT6030 = Control_in.signal[4]
	PT4c80 = Control_in.signal[5]
	CV7c08_feedback = Control_in.signal[6]
	SC1a98 = Control_in.signal[7]
	DP7a10 = Control_in.signal[8]
	ST4080 = Control_in.signal[9]
	AC4080_Usp = Control_in.signal[10]
	PT7008 = Control_in.signal[11]	
	FT1a16 = Control_in.signal[12]
	
	// VALCALC: CALCULATE VALUES
	ST6010 = CRYO_PF_prop_vs_Px(fluid,PT6030-5,0,fprop_temperature,ier,ipx,ipy) - 273.15
	DP7c08 = PT7008 - PT6030 // dP between pump discharge and surface storage
	SC3004 = ST3004 - TT3004 // detector inlet subcooling
	
	// Chiller
	IF (Stp>1) THEN
		GP5002 = 100
	ELSE
		GP5002 = 0
	END IF
	// Pump
	IF (Stp>=3 AND LP1a02_SC_OK) THEN
		LP1a02_start = TRUE
		LP1a02 = 100
	ELSE
		LP1a02_start = FALSE
		LP1a02 = 0
	END IF
	// Detector
	IF (Stp==7) THEN
		EH3006 = 15 // establish small vapour quality
	ELSEIF (Stp==8) THEN
		EH3006 = 100 // detector on
	ELSE
		EH3006 = 0 // detector off
	END IF
	
	// CONTROL VALVES
	/*
	 * NOTE on surface storage valves: 
	 * There are two places where the PID needs to _regulate_
	 * instead of being in a forced position. 
	 * 
	 * The first is at Step 1 where we equalize surface storage 
	 * and plant pressures. Only 4c80 is involved in this step.
	 * The _more_ it opens, the _lower_ the dP. Thus, this is
	 * reverse action
	 * 
	 * The second is at Step 3 and 4 where 7c08 regulates the 
	 * dP between pump outlet and low pressure side. In this case
	 * the more the valve is open, the lower the dP. Thus, once
	 * again reverse action.
	 * 
	 * However, we still need one TR_S for each valve because
	 * at some point BOTH valves need to be open.
	 * 
	 * Next thing to note: A _negative_ output of the PID will
	 * operate the Low Pressure valve (4c80), while a positive
	 * output will open the High Pressure valve (7c08).
	 */
	IF (Stp==0) THEN
	// SAFETY POSITION
		// Accumulator Level PID Controller:
		LC4080_PV = DP7c08 // process value for PID controller
		LC4080_SP = 10 // set point for PID controller
		LC4080_AtOutH = 100 // maximum output of PID
		LC4080_AtOutL = 0 // minimum output of PID (both set to garbage high numbers since forced mode active)
		// Surface storage valves
		CV7c08_TR_S = TRUE // aka forced mode active
		CV7c08_PosR = 0 // disconnect surface storage from plant in safety position
		CV4c80_TR_S = TRUE
		CV4c80_PosR = 0 // disconnect surface storage from plant in safety position
		// Pump discharge pressure control valve
		CV7008_SP = 0
		CV7008_TR_S = TRUE
		CV7008_PosR = 100
		// Detector bypass valve
		CV7a10_SP = 0
		CV7a10_TR_S = TRUE
		CV7a10_PosR = 0
		// Transfer line bypass valve
		PV7c10 = FALSE

	ELSEIF(Stp==1) THEN
	// EQUALIZING PRESSURES
		LC4080_PV = DP7c08
		LC4080_SP = 10
		LC4080_AtOutH = 100 // dummy values for controller
		LC4080_AtOutL = 0 // dummy values for controller
		CV7c08_TR_S = TRUE
		CV7c08_PosR = 0
		CV4c80_TR_S = TRUE
		CV4c80_PosR = 100 // small opening to have slow pressure changes
		CV7008_SP = 0
		CV7008_TR_S = TRUE
		CV7008_PosR = 100
		CV7a10_SP = 0
		CV7a10_TR_S = TRUE
		CV7a10_PosR = 0
		PV7c10 = FALSE
	
	ELSEIF(Stp==2) THEN
	// LIQUEFYING PUMP (chiller turns on and cools)
		LC4080_PV = DP7c08
		LC4080_SP = 10
		LC4080_AtOutH = 100
		LC4080_AtOutL = 0
		CV7c08_TR_S = TRUE
		CV7c08_PosR = 100
		CV4c80_TR_S = TRUE
		CV4c80_PosR = 100
		CV7008_SP = 0
		CV7008_TR_S = TRUE
		CV7008_PosR = 0
		CV7a10_SP = 0
		CV7a10_TR_S = TRUE
		CV7a10_PosR = 0
		PV7c10 = TRUE

	ELSEIF(Stp==3) THEN
	// PUMP START UP
		LC4080_PV = DP7c08
		LC4080_SP = 10
		LC4080_AtOutH = 100
		LC4080_AtOutL = 0
		CV7c08_TR_S = FALSE // move to PID control for controlling pump pressure drop
		CV7c08_PosR = 0
		CV4c80_TR_S = TRUE
		CV4c80_PosR = 100 // fully open to make surface storage act as accu
		CV7008_SP = 0
		CV7008_TR_S = TRUE
		CV7008_PosR = 0 // off for the moment
		CV7a10_SP = 0
		CV7a10_TR_S = TRUE
		CV7a10_PosR = 0
		PV7c10 = TRUE
		
	ELSEIF(Stp==4) THEN
	// PRESSURIZE DETECTOR
		// CV7c08: REGULATING towards eventual OFF
		LC4080_AtOutH = 100
		LC4080_AtOutL = 0
		LC4080_PV = DP7c08
		LC4080_SP = 18
		CV7c08_PosR = 0
		/*IF (CV7c08_feedback==0) THEN
			CV7c08_TR_S = TRUE
		ELSE
			CV7c08_TR_S = FALSE // move to PID control
		END IF*/
		CV4c80_TR_S = TRUE // surface storage still acting as accumulator
		CV4c80_PosR = 100 // therefore, keep valve fully open
		CV7008_SP = 5 // Lower than CV7c08, causing flow to go through this path instead
		CV7008_TR_S = FALSE // moved to PID controller
		CV7008_PosR = 100 // dummy value
		CV7a10_SP = 0
		CV7a10_TR_S = TRUE
		CV7a10_PosR = 0
		PV7c10 = TRUE
	
	ELSEIF(Stp==5) THEN
	// DETECTOR LIQUID CIRCULATION
		LC4080_PV = 0
		LC4080_SP = 0
		LC4080_AtOutH = 0
		LC4080_AtOutL = -100
		CV7c08_TR_S = TRUE
		CV7c08_PosR = 0
		CV4c80_TR_S = TRUE
		CV4c80_PosR = 100
		CV7008_SP = 5
		CV7008_TR_S = FALSE
		CV7008_PosR = 100
		CV7a10_SP = CV7a10_SPval
		CV7a10_TR_S = FALSE
		CV7a10_PosR = 0
		PV7c10 = FALSE
	
	ELSEIF(Stp==6) THEN
	// ACCUMULATOR CONTROL
		LC4080_PV = 0
		LC4080_SP = 0
		LC4080_AtOutH = 100
		LC4080_AtOutL = -100
		CV7c08_TR_S = TRUE
		CV4c80_TR_S = TRUE
		CV7008_SP = 5
		CV7008_TR_S = FALSE
		CV7008_PosR = 100
		CV7a10_SP = CV7a10_SPval
		CV7a10_TR_S = FALSE
		CV7a10_PosR = 0
		PV7c10 = FALSE
		
	ELSEIF(Stp==7) THEN
	// DETECTOR COOL-DOWN
		LC4080_PV = 0
		LC4080_SP = 0
		LC4080_AtOutH = 100
		LC4080_AtOutL = -100
		CV7c08_TR_S = TRUE
		CV4c80_TR_S = TRUE
		CV7008_SP = 5
		CV7008_TR_S = FALSE
		CV7008_PosR = 100
		CV7a10_SP = CV7a10_SPval
		CV7a10_TR_S = FALSE
		CV7a10_PosR = 0
		PV7c10 = FALSE
	
	ELSEIF(Stp==8) THEN
	// DETECTOR POWER ALLOW
		// SURFACE STORAGE REGULATION
		LC4080_PV = 0
		LC4080_SP = 0
		LC4080_AtOutH = 100
		LC4080_AtOutL = -100
		CV7c08_TR_S = TRUE
		CV4c80_TR_S = TRUE
		CV7008_SP = 5
		CV7008_TR_S = FALSE
		CV7008_PosR = 100
		CV7a10_SP = CV7a10_SPval
		CV7a10_TR_S = FALSE
		CV7a10_PosR = 0
		PV7c10 = FALSE	
	ELSE
		ASSERT (TRUE) FATAL "Incorrect value for Stp"
	END IF
	
	// ACCUMULATOR LEVEL CONTROLLER
	// make sure we don't steal too much flow from bypass:
	maxPos = spliceFunction(posHH,posLL,FT1a16-(bypassHH+bypassLL)/2,(bypassHH-bypassLL)/2)
	ePos = spliceFunction(maxPos,posLL,WT4080-(emptyLL+(emptyHH-emptyLL)/2),(emptyHH-emptyLL)/2)
	fPos = spliceFunction(posLL,posHH,WT4080-(fillLL+(setPoint-fillLL)/2),(setPoint-fillLL)/2)
	IF (Stp>=6) THEN
		CV4c80_PosR = fPos
		CV7c08_PosR = ePos
	END IF
	
	// TL: TRANSITION LOGIC
	timer = TIME - time_prev // Reset everytime the step changes
	Tr0 = NOT Run_Order.signal[1]
	Tr1 = Run_Order.signal[1]
	Tr2 = PT4c80<=PT6030+dP_ctrl AND PT4c80>=PT6030-dP_ctrl
	Tr3 = LP1a02_SC_OK AND CV7c08_feedback==100 AND timer>=TrTime
	Tr4 = DP7c08>=10 AND DP7c08<=20 AND timer>=TrTime
	Tr5 = SC3004>3 AND timer>=TrTime
	Tr6 = DP7a10>=CV7a10_SPval-dP_ctrl AND DP7a10<=CV7a10_SPval+dP_ctrl AND timer>=TrTime
	Tr7 = ST4080<=ST6010+dP_ctrl AND timer>=TrTime
	Tr8 = ST4080<=AC4080_Usp+0.2 AND timer>=TrTime
	Tr[1] = Tr1
	Tr[2] = Tr2
	Tr[3] = Tr3
	Tr[4] = Tr4
	Tr[5] = Tr5
	Tr[6] = Tr6
	Tr[7] = Tr7
	Tr[8] = Tr8

	// The case of safety position switch activation is in the DISCRETE block.
	Stp_prev = Stp
	IF (Stp<8 AND Tr[Stp+1]) THEN
		Stp = Stp+1
		timer = 0
		time_prev = TIME
	END IF
			
	// DL: DEPENDENT LOGIC
	// Boolean signals
	Bool_out.signal[1] = PV7c10 // Bypass-to-transfer-line path
	Bool_out.signal[2] = CV7c08_TR_S
	Bool_out.signal[3] = CV4c80_TR_S
	Bool_out.signal[4] = CV7008_TR_S
	Bool_out.signal[5] = CV7a10_TR_S
	
	// Output Signals
	Control_out.signal[1] = Stp
	Control_out.signal[2] = GP5002
	Control_out.signal[3] = LP1a02
	Control_out.signal[4] = LC4080_PV
	Control_out.signal[5] = LC4080_SP
	Control_out.signal[6] = LC4080_AtOutH
	Control_out.signal[7] = LC4080_AtOutL
	Control_out.signal[8] = CV7c08_PosR
	Control_out.signal[9] = CV4c80_PosR
	Control_out.signal[10] = CV7008_SP
	Control_out.signal[11] = CV7008_PosR
	Control_out.signal[12] = CV7a10_SP
	Control_out.signal[13] = CV7a10_PosR
	Control_out.signal[14] = EH3006
END SEQUENTIAL
END COMPONENT



COMPONENT DEMO_PCO_Accu
"PCO for the Accumulator of DEMO [OBSOLETE]"
/*
 * NOTE: for controlling the Accumulator heating and cooling
 * there are two modes: damper and set point. In damper mode
 * only the heater works, while for set-point, both heating
 * and cooling have to work. I handle this by setting the 
 * OutMin and OutMax values of the controller, instead of 
 * switching between two different controllers. However, this
 * means only one set of PID params for both cases = not ideal.
 */
PORTS
	IN bool_signal(n=1) Run_Order
	IN analog_signal(n=7) Control_in
	OUT analog_signal(n=4) Control_out
	OUT bool_signal(n=5) Bool_out
DECLS
	// Inputs
	INTEGER Stp "Stepper position"
	REAL ST4080 "Accumulator saturation temperature"
	REAL PT6010 "Surface storage basement pressure"
	REAL ST7076 "Saturation temperature at inlet of BPR"
	REAL TT3004 "Detector inlet temperature"
	REAL AC4080_Usp "User-defined set-point"
	REAL TT1a98 "Pump inlet temperature"
	// Outputs
	BOOLEAN PV7a12 "Detector bypass valve"
	BOOLEAN PV7078 "Back pressure regulator bypass valve"
	REAL CV7078_SP "Back pressure regulator set-point"
	BOOLEAN CV7078_TR_S "BPR Tracking Mode Set"
	REAL CV7078_PosR "BPR Position request"
	BOOLEAN AC4080_SetPtCtrl "Regulator selection" // TRUE: Set-point control, FALSE: Damper mode
	// States
	REAL AC4080_Tsp
	REAL AC4080_Asp
	// Local Variables
	CONST ENUM ChemName fluid = Carbon_Dioxide
	REAL ST6010 "Saturation temperature at 5 bar below surface storage" // used as set point for PV7078
	PRIVATE INTEGER ier,ipx,ipy // error codes
	REAL STC4080_PosR
	BOOLEAN STC4080_TR_S
CONTINUOUS
SEQUENTIAL
	// DEMUX INPUT SIGNAL
	Stp = Control_in.signal[1]
	ST4080 = Control_in.signal[2]
	PT6010 = Control_in.signal[3]
	ST7076 = Control_in.signal[4]
	TT3004 = Control_in.signal[5]
	AC4080_Usp = Control_in.signal[6]
	TT1a98 = Control_in.signal[7]
	ST6010 = CRYO_PF_prop_vs_Px(fluid,PT6010-5,0,fprop_temperature,ier,ipx,ipy)-273.15

	// VALCALC: CALCULATE VALUES
	AC4080_Asp = TT1a98 + 6
	IF (Stp==0) THEN
	// Safety Position
		PV7a12 = TRUE
		PV7078 = TRUE
		CV7078_TR_S = TRUE
		CV7078_PosR = 100
		CV7078_SP = 0
		AC4080_SetPtCtrl = FALSE
		AC4080_Tsp = 0
	ELSEIF (Stp>=1 AND Stp<=3) THEN
	// Damper Mode
		PV7a12 = TRUE
		PV7078 = TRUE
		CV7078_TR_S = TRUE
		CV7078_PosR = 100
		CV7078_SP = 0
		AC4080_SetPtCtrl = FALSE
		AC4080_Tsp = 80 // heater temperature in damper mode
	ELSEIF (Stp==4) THEN
		PV7a12 = TRUE
		PV7078 = FALSE
		CV7078_TR_S = FALSE
		CV7078_PosR = 100
		CV7078_SP = TT3004 + 5
		AC4080_SetPtCtrl = FALSE
		AC4080_Tsp = 80 // heater temperature
	ELSEIF (Stp==5) THEN
		PV7a12 = FALSE
		PV7078 = FALSE
		CV7078_TR_S = FALSE
		CV7078_PosR = 100
		CV7078_SP = TT3004 + 5
		AC4080_SetPtCtrl = FALSE
		AC4080_Tsp = 80 // heater temperature
	ELSE
	// SET-POINT
		PV7a12 = FALSE
		IF (ST4080<=ST6010-5) THEN // AND CV7078 IS CLOSED
			PV7078 = TRUE
		ELSE
			PV7078 = FALSE // open up if set-point in Accu control range
		END IF
		IF (Stp==6) THEN
			CV7078_SP = TT3004 + 5
		ELSE
			CV7078_SP = AC4080_Usp
		END IF
		CV7078_TR_S = FALSE
		CV7078_PosR = 100
		AC4080_SetPtCtrl = TRUE
		AC4080_Tsp = min(max(AC4080_Asp,AC4080_Usp),ST6010)
	END IF
	
	IF (Stp<=5) THEN
		STC4080_TR_S = TRUE
	ELSE
		STC4080_TR_S = FALSE
	END IF
	STC4080_PosR = 0

	// OUTPUT SIGNALS
	Control_out.signal[1] = CV7078_SP
	Control_out.signal[2] = CV7078_PosR
	Control_out.signal[3] = AC4080_Tsp
	Control_out.signal[4] = STC4080_PosR
	
	Bool_out.signal[1] = PV7a12
	Bool_out.signal[2] = PV7078
	Bool_out.signal[3] = CV7078_TR_S
	Bool_out.signal[4] = AC4080_SetPtCtrl
	Bool_out.signal[5] = STC4080_TR_S
END SEQUENTIAL
END COMPONENT









COMPONENT DEMO_PCO_Plant
"Process Control Object for the DEMO plant"
PORTS
	IN bool_signal(n=1) Run_Order
	IN analog_signal(n=10) Control_in
	OUT bool_signal(n=5) Bool_out
	OUT analog_signal(n=10) Control_out
	IN analog_signal(n=1) AccuStp
	OUT analog_signal(n=1) PlantStp
	IN bool_signal(n=1) AccuReq // Only used if acts as spare plant
DATA
	BOOLEAN forceStart = FALSE "If TRUE, skips straight to Step 4"
	BOOLEAN isBackup = FALSE "If TRUE, acts as spare plant and stays connected to surface storage till needed"
	ENUM ChemName fluid = Carbon_Dioxide
	REAL bypassHH = 0.1 UNITS u_kg_s "Value of bypass flow at which to start limiting surface storage valve"
	REAL bypassLL = 0.01 UNITS u_kg_s "Minimum value for bypass flow" // turn valve off at this point
	REAL DP2022_SPreq = 10 UNITS u_bar "Detector dP request"
	REAL emptyLL = 12 UNITS u_pct "Level at which to stop emptying, should be a bit higher than setPoint"
	REAL emptyHH = 40 UNITS u_pct "Level at which maxiumum emptying speed"
	REAL fillLL = 5 UNITS u_pct "Level low-limit for maximum filling speed"
	REAL posHH = 100 UNITS u_pct "Max opening (%) of surface storage valves"
	REAL posLL = 0 UNITS u_pct "Minimum opening (%) of surface storage valves"
	REAL setPoint = 10 UNITS u_pct "Set point for accumulator level"
	REAL SC9a98_TtL = 2 UNITS u_K "Subcooling threshold RESET for pump inlet"
	REAL SC9a98_TtH = 5 UNITS u_K "Subcooling threshold SET for pump inlet"
	REAL TrTime = 120 UNITS u_s "Minimum transition time between steps"
DECLS
	// Inputs
	REAL DP1014 "Pressure difference between plant outlet and surface storage basement"
	REAL DP1s14 "Pressure difference between plant outlet and back pressure regulator inlet"
	REAL DP2022 "Pressure drop across the transfer line bypass valve" // controlled by CV1a16
	REAL FT1a16 "Flow through plant bypass" // used in level control rate limiter logic
	REAL LT4080 "Local accumulator liquid level"
	REAL PT6080 "Pressure at the base of surface storage"
	REAL PT9082 "Pressure at outlet of surface storage low pressure valve"
	REAL SC9a98 "Subcooling at pump inlet"
	REAL CV1s14_feedback
	REAL freq UNITS u_Hz "Pump frequency request"
	
	// Outputs
	REAL LP1a02
	BOOLEAN CV1014_TR_S
	REAL CV1014_PosR
	REAL DPC1014_SP
	BOOLEAN CV1s14_TR_S
	REAL CV1s14_PosR
	REAL DPC1s14_SP
	BOOLEAN CV1a16_TR_S
	REAL CV1a16_PosR
	REAL DPC2022_SP
	REAL GP5002
	BOOLEAN CV9s82_TR_S
	REAL CV9s82_PosR
	REAL DPC9s82_SP
	BOOLEAN EV9s82
	
	// States
	BOOLEAN LP1a02_SC_OK "Pump subcooling OK flag"
	BOOLEAN LP1a02_start
	
	// Transitions
	INTEGER Stp "Stepper position"
	INTEGER Stp_prev "Previous step"
	REAL timer UNITS u_s "Time counter"
	REAL time_prev UNITS u_s "Time at which to start counting"
	BOOLEAN Tr0, Tr1, Tr2, Tr3, Tr4 // Transition conditions
	BOOLEAN Tr[4]
	
	// Local Variables
	CONST REAL dP_ctrl = 0.1 UNITS u_bar "Deadband dP for checking stepper conditions"
	PRIVATE INTEGER ier,ipx,ipy // error codes
	
	// Level Controller-related
	REAL maxPos "limiter on HH of emptying valve" // ensures we don't steal too much from bypass flow
	REAL ePos "Emptying valve position"
	REAL fPos "Filling valve position"
	
INIT
	IF (forceStart) THEN
		Stp = 4
	ELSE
		Stp = 0
	END IF
	timer = 0
	LP1a02_SC_OK = FALSE
	
DISCRETE
	// Safety position
	WHEN (Tr0) THEN
		Stp = 0
	END WHEN
	// Pump subcooling
	WHEN(SC9a98 > SC9a98_TtH) THEN
		LP1a02_SC_OK = TRUE
	END WHEN
	WHEN(SC9a98 < SC9a98_TtL) THEN
		LP1a02_SC_OK = FALSE
	END WHEN

CONTINUOUS
SEQUENTIAL
	// DEMUX INPUT SIGNALS ########################
	DP1014 = Control_in.signal[1]
	DP1s14 = Control_in.signal[2]
	FT1a16 = Control_in.signal[3]
	DP2022 = Control_in.signal[4]
	LT4080 = Control_in.signal[5]
	PT6080 = Control_in.signal[6]
	PT9082 = Control_in.signal[7]
	SC9a98 = Control_in.signal[8]
	CV1s14_feedback = Control_in.signal[9]
	freq = Control_in.signal[10]
	
	// VALCALC: CALCULATE VALUES ########################
	// Chiller
	IF (Stp>1) THEN
		GP5002 = 100
	ELSE
		GP5002 = 0
	END IF
	
	// Pump
	IF (Stp>=3 AND LP1a02_SC_OK) THEN
		LP1a02_start = TRUE
		LP1a02 = freq
	ELSE
		LP1a02_start = FALSE
		LP1a02 = 0
	END IF
	
	// Level controller
	maxPos = spliceFunction(posHH,posLL,FT1a16-(bypassHH+bypassLL)/2,(bypassHH-bypassLL)/2) // (make sure we don't steal too much flow from bypass)
	ePos = spliceFunction(maxPos,posLL,LT4080-(emptyLL+(emptyHH-emptyLL)/2),(emptyHH-emptyLL)/2)
	fPos = spliceFunction(posLL,posHH,LT4080-(fillLL+(setPoint-fillLL)/2),(setPoint-fillLL)/2)
		
	// Control Valves
	IF (Stp==0) THEN
	// Safety position
		CV1014_TR_S = TRUE
		CV1014_PosR = 100 // NO (Normally Open)
		DPC1014_SP = DP1014
		CV1a16_TR_S = TRUE
		CV1a16_PosR = 0 // NO
		DPC2022_SP = DP2022
		
		CV1s14_TR_S = TRUE
		CV1s14_PosR = 0 // NC
		DPC1s14_SP = DP1s14
		CV9s82_TR_S = TRUE
		CV9s82_PosR = 0 // NC
		DPC9s82_SP = PT9082
		EV9s82 = FALSE // NC
		
	ELSEIF(Stp==1) THEN
	// Equalise with surface storage pressure
		CV1014_TR_S = TRUE
		CV1014_PosR = 100
		DPC1014_SP = DP1014
		CV1a16_TR_S = TRUE
		CV1a16_PosR = 0
		DPC2022_SP = DP2022
		
		CV1s14_TR_S = TRUE
		CV1s14_PosR = 0
		DPC1s14_SP = DP1s14
		CV9s82_TR_S = FALSE
		CV9s82_PosR = 0 // dummy
		DPC9s82_SP = PT6080
		EV9s82 = FALSE
		
	ELSEIF(Stp==2) THEN
	// Liquefy pump
		CV1014_TR_S = TRUE
		CV1014_PosR = 0
		DPC1014_SP = DP1014
		CV1a16_TR_S = TRUE
		CV1a16_PosR = 0
		DPC2022_SP = DP2022
		
		CV1s14_TR_S = TRUE
		CV1s14_PosR = 100
		DPC1s14_SP = DP1014
		CV9s82_TR_S = TRUE
		CV9s82_PosR = 100
		DPC9s82_SP = PT9082 // dummy (set-point = process value)
		EV9s82 = TRUE
		
	ELSEIF(Stp==3) THEN
	// Pump startup
		CV1014_TR_S = TRUE
		CV1014_PosR = 0
		DPC1014_SP = DP1014
		CV1a16_TR_S = TRUE
		CV1a16_PosR = 0
		DPC2022_SP = DP2022
		
		CV1s14_TR_S = FALSE
		CV1s14_PosR = 100 // dummy value
		DPC1s14_SP = 10
		CV9s82_TR_S = TRUE
		CV9s82_PosR = 100
		DPC9s82_SP = PT9082 // dummy (set-point = process value)
		EV9s82 = TRUE
		
	ELSEIF(Stp==4) THEN
	// Provide liquid
		CV1014_TR_S = FALSE
		CV1014_PosR = 0 // dummy value
		DPC1014_SP = 5 // wants a lower dP than CV1s14, thus it will steal all the flow slowly
		CV1a16_TR_S = FALSE // start regulating detector pressure drop
		CV1a16_PosR = 0 // dummy value
		IF (AccuStp.signal[1] < 3) THEN
			DPC2022_SP = 1 // lower set-point for detector pressurisation
		ELSE
			DPC2022_SP = DP2022_SPreq // higher set-point for detector flow control
		END IF
		IF (isBackup AND AccuReq.signal[1]==FALSE) THEN
			CV9s82_TR_S = TRUE
			CV9s82_PosR = 100
			CV1a16_TR_S = TRUE // Keep CV1a16 at fixed position for backup
			CV1a16_PosR = 80
		ELSEIF (AccuStp.signal[1]<4) THEN
			CV1s14_TR_S = FALSE
			CV1s14_PosR = CV1s14_feedback
			DPC1s14_SP = 10 // higher than CV1014 set-point
			CV9s82_TR_S = TRUE
			CV9s82_PosR = 100
			DPC9s82_SP = PT9082 // dummy (set-point = process value)
			EV9s82 = TRUE
		ELSE
			CV1s14_TR_S = TRUE // forced mode
			CV1s14_PosR = ePos
			DPC1s14_SP = DP1s14 // make set-point = process-value
			CV9s82_TR_S = TRUE
			CV9s82_PosR = fPos // level control
			EV9s82 = FALSE // close this valve, no longer needed.
		END IF
	END IF
	
	// TL: TRANSITION LOGIC ########################
	timer = TIME - time_prev // Reset everytime the step changes
	
	Tr0 = NOT Run_Order.signal[1]
	Tr1 = Run_Order.signal[1]
	Tr2 = PT9082 <= PT6080 + dP_ctrl AND PT9082 >= PT6080 - dP_ctrl
	Tr3 = LP1a02_SC_OK AND timer >= TrTime
	Tr4 = DP1s14 >= 8 AND timer >= TrTime // set-point hard-coded atm
	
	Tr[1] = Tr1
	Tr[2] = Tr2
	Tr[3] = Tr3
	Tr[4] = Tr4

	// The case of safety position switch activation is in the DISCRETE block.
	Stp_prev = Stp
	IF (forceStart) THEN
		Stp = 4
	ELSEIF (Stp<4 AND Tr[Stp+1]) THEN
		Stp = Stp+1
		timer = 0
		time_prev = TIME
	END IF
			
	// DL: DEPENDENT LOGIC ########################
	// Boolean signals
	Bool_out.signal[1] = CV1014_TR_S
	Bool_out.signal[2] = CV1s14_TR_S
	Bool_out.signal[3] = CV1a16_TR_S
	Bool_out.signal[4] = CV9s82_TR_S
	Bool_out.signal[5] = EV9s82
	
	// Output Signals
	
	PlantStp.signal[1] = Stp
	
	Control_out.signal[1] = GP5002
	Control_out.signal[2] = LP1a02
	Control_out.signal[3] = CV1014_PosR
	Control_out.signal[4] = DPC1014_SP
	Control_out.signal[5] = CV1s14_PosR
	Control_out.signal[6] = DPC1s14_SP
	Control_out.signal[7] = CV1a16_PosR
	Control_out.signal[8] = DPC2022_SP
	Control_out.signal[9] = CV9s82_PosR
	Control_out.signal[10] = DPC9s82_SP
END SEQUENTIAL
END COMPONENT



COMPONENT DEMO_PCO_Accumulator
"Process Control Object for the DEMO Accumulator"
PORTS
	IN bool_signal(n=1) Run_Order
	IN analog_signal(n=10) Control_in
	OUT bool_signal(n=5) Bool_out
	OUT analog_signal(n=5) Control_out
	IN bool_signal(n=1) Plant_Run_Order
	OUT bool_signal(n=1) Spare
	
DATA
	INTEGER nHead = 3 "Number of pump heads"
	BOOLEAN forceStart = FALSE "If TRUE, skips straight to last step"
	REAL TrTime = 120 UNITS u_s "Minimum time between transitions"
	REAL DPC2022_SP = 10 UNITS u_bar "Pressure drop set point for CV1a16"
	REAL ST8074_Det_SP = 290 UNITS u_bar "BPR set point for Stp<5"
DECLS
	// INPUTS
	REAL System_Usp UNITS u_C "User-defined temperature set-point"
	REAL FT1016 UNITS u_kg_s "plant outlet mass flow rate"
	REAL DP2022 UNITS u_bar "Pressure drop across transfer line bypass"
	REAL SC3032 UNITS u_C "Detector refrigerant subcooling"
	REAL ST3032 UNITS u_C "Detector refrigerant saturation temperature"
	REAL TT3032 UNITS u_C "Detector wall temperature"
	REAL ST4080 UNITS u_C "Accumulator saturation temperature"
	REAL PT6080 UNITS u_bar "Surface storage pressure, including static"
	REAL CV8076_feedback UNITS u_pct "Back pressure regulator valve position"
	REAL ST8074 UNITS u_C "BPR inlet saturation temperature"
	// OUTPUTS
	BOOLEAN DamperMode
	BOOLEAN DetectorPowerAllow
	BOOLEAN PV2a22
	REAL STC4080_ActSP UNITS u_C "Active set point for accumulator pressure controller"
	REAL STC4080_SP UNITS u_C "Output to STC4080" // when TR_S is TRUE, this is set to ST4080 to avoid PID jumps during transiton
	REAL STC8074_ActSP UNITS u_C "Active set point, calculated, for back pressure regulator"
	REAL STC8074_SP UNITS u_C "Output to STC8074" // when TR_S is TRUE, this is set to ST8074 to avoid PID jumps during transition
	REAL CV8076_PosR
	BOOLEAN CV8076_TR_S
	BOOLEAN EV8076
	REAL PreVapourisationSignal
	REAL STC4080_AtOutH
	REAL STC4080_AtOutL
	// STATES
	REAL System_Asp "System automatic set-point"
	REAL System_ActSP "System active set-point"
	REAL ST6080 UNITS u_C "Calculated as: ST(PT6020 - 5 bar)"
	// CONSTANTS	
	CONST REAL dP_ctrl = 0.2 UNITS u_bar
	CONST ENUM ChemName fluid = Carbon_Dioxide
	CONST REAL SC_Detector_SP = 3 UNITS u_K "Minimum detector liquefication subcooling"
	DISCR REAL FT1016_Threshold
	// STEPPER
	BOOLEAN Plant_request "Request flow from associated cooling plant"
	INTEGER Stp "Accumulator stepper position"
	INTEGER Stp_prev
	BOOLEAN Tr0, Tr1, Tr2, Tr3, Tr4, Tr5, Tr6	
	BOOLEAN Tr[6]
	REAL timer UNITS u_s "Timer for PCO step transition"
	REAL startTime UNITS u_s "Step transition start time"
	BOOLEAN flowTimerFlag "Flag to start flowTimer"
	REAL flowTimer UNITS u_s "Timer to wait for stable flow"
	REAL flowTimerStartTime UNITS u_s "Time at which flowTimer begins counting"
	BOOLEAN powerTimerFlag "Flag for DetectorPowerAllow"
	REAL powerTimer UNITS u_s "Timer for DetectorPowerAllow"
	REAL powerTimerStartTime UNITS u_s "Time at which powerTimer begins counting"
	INTEGER ier,jx,jy // error codes
INIT
	FT1016_Threshold = 0.2*nHead
	IF (forceStart) THEN
		Stp = 6
	ELSE
		Stp = 0
	END IF

DISCRETE
	// Set-Reset for BPR bypass
	WHEN Stp>2 AND CV8076_feedback>99 THEN
		EV8076 = TRUE
	END WHEN
	WHEN Stp>2 AND CV8076_feedback<97 THEN
		EV8076 = FALSE
	END WHEN

CONTINUOUS
SEQUENTIAL
	// DEMUX INPUT SIGNALS ################
	System_Usp = Control_in.signal[1] // [°C]
	FT1016 = Control_in.signal[2] // flow from pump module
	DP2022 = Control_in.signal[3] // pressure drop over transfer lines
	ST3032 = Control_in.signal[5] // detector refrigerant saturation temperature
	TT3032 = Control_in.signal[6] // detector wall temperature
	ST4080 = Control_in.signal[7]
	PT6080 = Control_in.signal[8]
	ST8074 = Control_in.signal[9]
	CV8076_feedback = Control_in.signal[10]

	// VALCALC: CALCULATE VALUES ################
	ST6080 = CRYO_PF_prop_vs_Px(fluid,PT6080-5,0,fprop_temperature,ier,jx,jy) - 273.15 // already 5 bar lowered
	SC3032 = ST3032 - TT3032
	
	// Temperature Set Point
	IF(Stp<=4) THEN
		System_Asp = TT3032 + 5
	ELSE
		System_Asp = System_Usp
	END IF
	System_ActSP = max(System_Asp,System_Usp)
	IF(System_ActSP<=ST6080) THEN 
		STC8074_ActSP = -60
		STC4080_ActSP = System_ActSP
	ELSE
		STC8074_ActSP = System_ActSP
		STC4080_ActSP = ST6080
	END IF
	
	IF (Stp==0) THEN
	// SAFETY POSITION
		DamperMode = TRUE
		STC4080_SP = ST4080
		PV2a22 = TRUE
		CV8076_TR_S = TRUE
		CV8076_PosR = 100
		STC8074_SP = ST8074
	ELSEIF (Stp==1) THEN
	// PLANT START
		DamperMode = TRUE
		STC4080_SP = ST4080
		PV2a22 = TRUE
		CV8076_TR_S = TRUE
		CV8076_PosR = 100
		STC8074_SP = ST8074
	ELSEIF (Stp==2) THEN
	// DETECTOR PRESSURISATION
		DamperMode = TRUE
		STC4080_SP = ST4080
		PV2a22 = TRUE
		CV8076_TR_S = FALSE
		CV8076_PosR = 100
		STC8074_SP = ST8074_Det_SP
	ELSEIF (Stp==3) THEN
	// DETECTOR CIRCULATION
		DamperMode = TRUE
		STC4080_SP = ST4080
		PV2a22 = FALSE
		CV8076_TR_S = FALSE
		STC8074_SP = ST8074_Det_SP
	ELSEIF (Stp==4) THEN
	// ACCUMULATOR CONTROL
		DamperMode = FALSE
		STC4080_SP = STC4080_ActSP
		PV2a22 = FALSE
		CV8076_TR_S = FALSE
		STC8074_SP = ST8074_Det_SP
	ELSEIF (Stp==5) THEN
	// DETECTOR COOLDOWN
		DamperMode = FALSE
		STC4080_SP = STC4080_ActSP
		PV2a22 = FALSE
		CV8076_TR_S = FALSE
		STC8074_SP = STC8074_ActSP
	ELSEIF (Stp==6) THEN
	// DETECTOR POWER ALLOW
		DamperMode = FALSE
		STC4080_SP = STC4080_ActSP
		PV2a22 = FALSE
		CV8076_TR_S = FALSE
		STC8074_SP = STC8074_ActSP
	END IF
	
	// Timers
	IF (FT1016 > FT1016_Threshold AND flowTimerFlag==FALSE) THEN
		flowTimerFlag = TRUE
		flowTimerStartTime = TIME
	END IF
	IF (Stp==6 AND powerTimerFlag==FALSE) THEN
		powerTimerFlag = TRUE
		powerTimerStartTime = TIME
	END IF
	IF (powerTimerFlag) THEN
		powerTimer = TIME - powerTimerStartTime
	END IF
	IF (Stp==6 AND ((powerTimer>3600 AND TIME<20000) OR (TIME>25000 AND TIME<40000) OR (TIME>45000))) THEN
		DetectorPowerAllow = TRUE
	ELSE
		DetectorPowerAllow = FALSE
	END IF
	
	// Pre-vapourisation
	IF (Stp==6 AND timer>1200 AND DetectorPowerAllow==FALSE) THEN
		PreVapourisationSignal = 100
	ELSE
		PreVapourisationSignal = 0
	END IF
	
	// Detector power
	STC4080_AtOutH = 100
	STC4080_AtOutL = -100
	
	// Spare Plant
	IF (Plant_Run_Order.signal[1] == FALSE) THEN
		Spare.signal[1] = TRUE
	ELSE
		Spare.signal[1] = FALSE
	END IF
	
	
	// TL: TRANSITION LOGIC ################	
	timer = TIME - startTime
	IF (flowTimerFlag) THEN
		flowTimer = TIME - flowTimerStartTime
	END IF
	
	Tr0 = NOT Run_Order.signal[1]
	Tr1 = Run_Order.signal[1]
	Tr2 = FT1016 >= FT1016_Threshold AND flowTimer >= 1800
	Tr3 = SC3032 > SC_Detector_SP
	Tr4 = DP2022 > DPC2022_SP-2 AND timer >= TrTime
	Tr5 = ST4080 <= ST6080 + dP_ctrl AND timer >= TrTime
	Tr6 = ST4080 >= System_Usp - 0.2 AND ST4080 <= System_Usp + 0.2 AND timer >= TrTime
	
	Tr[1] = Tr1
	Tr[2] = Tr2
	Tr[3] = Tr3
	Tr[4] = Tr4
	Tr[5] = Tr5
	Tr[6] = Tr6
	
	Stp_prev = Stp
	IF (forceStart) THEN
		Stp = 6
	ELSEIF (Stp<6 AND Tr[Stp+1]) THEN
		Stp = Stp+1
		timer = 0
		startTime = TIME
	END IF

	// DL: DEPENDENT LOGIC ################
	Bool_out.signal[1] = PV2a22 // Transfer line bypass valve position request
	Bool_out.signal[2] = EV8076 // Back pressure regulator bypass valve position request
	Bool_out.signal[3] = CV8076_TR_S // Back pressure regulator controller tracking mode set
	Bool_out.signal[4] = DetectorPowerAllow // Boolean for detector power OK
	Bool_out.signal[5] = DamperMode // Accumulator to operate in damper mode

	Control_out.signal[1] = Stp // Accumulator stepper step
	Control_out.signal[2] = STC4080_SP // Accumulator temperature set point
	Control_out.signal[3] = CV8076_PosR // Back pressure regulator tracking mode position request
	Control_out.signal[4] = STC8074_SP // Back pressure regulator controller set-point
	Control_out.signal[5] = PreVapourisationSignal // Pre-vapourisation on/off
END SEQUENTIAL
END COMPONENT



COMPONENT LHCb_JDB_PCO
PORTS
	IN bool_signal(n=1) Run_Order
	IN analog_signal(n=10) Control_in
	OUT bool_signal(n=5) Control_out
	OUT analog_signal(n=1) stepper
DATA
	REAL timerTh = 300 UNITS u_s "Threshold for transition"
	REAL deadband = 0.5 UNITS u_bar "Deadband for pressure equalisation"
	REAL SC_setpoint = 0.5 UNITS u_K "Detector subcooling setpoint"
DECLS
	// Input Signals
	REAL PV3010_feedback
	REAL PV3B14_feedback
	REAL PV3C14_feedback
	REAL PV3036_feedback
	REAL PT4060 UNITS u_bar "Accumulator pressure"
	REAL PT3036 UNITS u_bar "Detector pressure"
	REAL SC3C20 UNITS u_K "Subcooling in detector"
	REAL ST4060 UNITS u_C "Accumulator Tsat"
	REAL Usp UNITS u_C "Set-point temperature"
	
	// Output Signals
	BOOLEAN PV3010
	BOOLEAN PV3B14
	BOOLEAN PV3C14
	BOOLEAN PV3036	
	// Transitions
	INTEGER Stp "Stepper position"
	INTEGER Stp_prev "Previous step"
	REAL timer UNITS u_s "Time counter"
	REAL startTime UNITS u_s "Time at which to start counting"
	BOOLEAN Tr0, Tr1, Tr2, Tr3, Tr4, Tr5, Tr6 // Transition conditions
	BOOLEAN Tr[6]
CONTINUOUS
SEQUENTIAL
	// INPUT SIGNALS
	PV3010_feedback = Control_in.signal[1]
	PV3B14_feedback = Control_in.signal[2]
	PV3C14_feedback = Control_in.signal[3]
	PV3036_feedback = Control_in.signal[4]
	PT4060 = Control_in.signal[5]
	PT3036 = Control_in.signal[6]
	SC3C20 = Control_in.signal[7]
	ST4060 = Control_in.signal[8]
	Usp = Control_in.signal[9]
	
	// VALCALC: CALCULATE VALUES
	IF (Stp==0) THEN
	// SAFETY POSITION
		PV3010 = FALSE
		PV3B14 = TRUE
		PV3036 = FALSE
		PV3C14 = FALSE
	ELSEIF (Stp==1) THEN
	// EQUALIZING PRESSURE
		PV3010 = FALSE
		PV3B14 = TRUE
		PV3036 = FALSE
		PV3C14 = FALSE
	ELSEIF (Stp==2) THEN
	// CONNECTING JB
		PV3010 = FALSE
		PV3B14 = TRUE
		PV3036 = TRUE
		PV3C14 = FALSE
	ELSEIF (Stp==3) THEN
	// LIQUEFYING AND CIRCULATING JB
		PV3010 = TRUE
		PV3B14 = TRUE
		PV3036 = TRUE
		PV3C14 = FALSE
	ELSEIF (Stp==4) THEN
	// DETECTOR CIRCULATION
		PV3010 = TRUE
		PV3B14 = TRUE
		PV3036 = TRUE
		PV3C14 = TRUE
	ELSEIF (Stp==5) THEN
	// CLOSE BYPASS
		PV3010 = TRUE
		PV3B14 = FALSE
		PV3036 = TRUE
		PV3C14 = TRUE
	ELSEIF (Stp==6) THEN
	// OPERATION
		PV3010 = TRUE
		PV3B14 = FALSE
		PV3036 = TRUE
		PV3C14 = TRUE
	END IF
	
	// TL: TRANSITION LOGIC ####################
	timer = TIME - startTime
	// Transition Conditions
	Tr0 = NOT Run_Order.signal[1]
	Tr1 = Run_Order.signal[1] AND PV3010_feedback<=1 AND PV3C14_feedback<=1 AND PV3036_feedback<=1
	Tr2 = PT4060>=PT3036-deadband AND PT4060<=PT3036+deadband AND timer>timerTh
	Tr3 = PV3036_feedback>=99 AND timer>timerTh
	Tr4 = PV3010_feedback>=99 AND SC3C20>=SC_setpoint AND timer>timerTh
	Tr5 = PV3C14_feedback>=99 AND timer>timerTh
	Tr6 = ST4060>=Usp-deadband AND ST4060<=Usp+deadband AND timer>timerTh
	
	// Transition array	
	Tr[1] = Tr1
	Tr[2] = Tr2
	Tr[3] = Tr3
	Tr[4] = Tr4
	Tr[5] = Tr5
	Tr[6] = Tr6
	// Step Timer
	IF (Stp<6 AND Tr[Stp+1]) THEN
		Stp = Stp+1
		timer = 0
		startTime = TIME
	END IF
	
	// OUTPUT SIGNALS
	Control_out.signal[1] = PV3010
	Control_out.signal[2] = PV3B14
	Control_out.signal[3] = PV3C14
	Control_out.signal[4] = PV3036
	Control_out.signal[5] = FALSE // dummy signal to not have to create a new Bool Demux
	stepper.signal[1] = Stp
END SEQUENTIAL
END COMPONENT



COMPONENT LHCb_PCO
PORTS
	IN analog_signal(n=6) Control_in
	OUT analog_signal(n=3) Control_out
DECLS
	REAL Stp
	REAL ST3036
	REAL TT3030
	REAL Asp
	REAL Usp
	REAL SC1002
	REAL ActSp
	REAL TT1002
	REAL PumpSpd
	REAL DetectorPower
CONTINUOUS
SEQUENTIAL
	Stp = Control_in.signal[1]
	TT3030 = Control_in.signal[2] // in °C
	Usp = Control_in.signal[3]
	SC1002 = Control_in.signal[4]
	TT1002 = Control_in.signal[5]
	ST3036 = Control_in.signal[6]

	// Automatic set point
	IF (Stp<3) THEN
		Asp = ST3036
	ELSEIF (Stp==3) THEN
		Asp = TT3030 // equalise with JDB
	ELSE
		Asp = TT1002 + 6 // limit pump subcooling loss
	END IF
	
	ActSp = max(Usp,Asp) // NOT ALWAYS.
	
	IF (SC1002>=5) THEN
		PumpSpd = 100
	ELSE
		PumpSpd = 0
	END IF
	
	IF Stp==6 THEN
		DetectorPower = 100
	ELSE
		DetectorPower = 0
	END IF
	
	Control_out.signal[1] = ActSp
	Control_out.signal[2] = PumpSpd
	Control_out.signal[3] = DetectorPower
END SEQUENTIAL
END COMPONENT