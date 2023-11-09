/*-----------------------------------------------------------------------------------------
 LIBRARY: CO2
 FILE: SplitRangeController
 AUTHOR: Viren Bhanot
 COMPANY: CERN
 DESCRIPTION: Splits signal between 
 CREATION DATE: 05/10/2016
-----------------------------------------------------------------------------------------*/
USE CRYOLIB
USE MATH
USE PORTS_LIB



ENUM OutputLimit = {NoLimit,UpperLimit,LowerLimit} "Limits for the Ramp component"



COMPONENT Constant
PORTS
	OUT analog_signal(n=1) s
DATA
	REAL signal =	0.0
CONTINUOUS
	s.signal[1] = signal
END COMPONENT



COMPONENT BoolToggle
PORTS
	IN analog_signal(n=1) in
	OUT bool_signal(n=1) out
CONTINUOUS
SEQUENTIAL
	IF (in.signal[1]==0) THEN
		out.signal[1] = FALSE 
	ELSE
		out.signal[1] = TRUE
	END IF
END SEQUENTIAL
END COMPONENT



--------------------------------------------------------------------------------
// SPLIT RANGE CONTROLLER
--------------------------------------------------------------------------------
COMPONENT SplitRangeController
"Split range controller for a single input PID"
PORTS
	IN analog_signal(n=1) in
	OUT analog_signal(n=1) out1 "Cooling signal (when input negative)"
	OUT analog_signal(n=1) out2 "Heating signal (when input positive)"
INIT
	out1.signal[1] = 0
	out2.signal[1] = 0
CONTINUOUS
<cool> 	out1.signal[1] = ZONE (in.signal[1]<0) abs(in.signal[1]) OTHERS 0
<heat> 	out2.signal[1] = ZONE (in.signal[1]>0) abs(in.signal[1]) OTHERS 0
END COMPONENT



COMPONENT SplitRangeController_phase IS_A SplitRangeController
"For the CORA model"
PORTS
	IN analog_signal(n=1) phase
CONTINUOUS
<:cool> out1.signal[1] = ZONE (in.signal[1]<0 AND phase.signal[1]>2) abs(in.signal[1]) OTHERS 0
END COMPONENT



COMPONENT SplitRangeController_2PID
"Split Range Controller for 2 input PIDs"
PORTS
	IN analog_signal(n=1) in1				
	IN analog_signal(n=1) in2
	IN analog_signal(n=1) PV
	IN analog_signal(n=1) SP			
	OUT analog_signal(n=1) out1
	OUT analog_signal(n=1) out2	
INIT
	out1.signal[1] = 0
	out2.signal[1] = 0
CONTINUOUS
	out1.signal[1] = ZONE (PV.signal[1]>(SP.signal[1])+5) in1.signal[1] OTHERS 0
	out2.signal[1] = ZONE (PV.signal[1]<(SP.signal[1])-5) in2.signal[1] OTHERS 0
END COMPONENT



COMPONENT Ramp (ENUM OutputLimit limit = UpperLimit)
"Linear ramp up to (or down to) a predefined limit"
PORTS
	OUT analog_signal(n=1) output
DATA
	REAL initial = 0
	REAL slope = 1
	REAL limitValue = 100 "Maximum/Minimum value output by Ramp"
CONTINUOUS
	IF (limit == UpperLimit) INSERT
		output.signal[1] = min(slope*TIME+initial, limitValue)
	ELSEIF (limit == LowerLimit) INSERT
		output.signal[1] = max(slope*TIME+initial, limitValue)
	ELSE
		output.signal[1] = slope*TIME+initial
	END IF
END COMPONENT



COMPONENT StartupRamp
PORTS
	OUT analog_signal(n=1) output
DATA
	REAL initial = 0 "Initial value"
	REAL final = 3500 "Final value"
	REAL timeToFinal = 15 UNITS u_s "Time to go from initial to final value"
DECLS
	DISCR REAL slope
INIT
	slope = (final - initial)/timeToFinal
CONTINUOUS
	output.signal[1] = ZONE(final>=initial) min(final,initial+slope*TIME) OTHERS max(final,initial+slope*TIME)
END COMPONENT



COMPONENT SmoothStartup
"Creates a log curve startup ramp"
PORTS
	OUT analog_signal(n=1) output
DATA
	REAL w_nom = 3 "Steady-state output to ramp up to"
	REAL tau = 10 UNITS u_s "Time constant for startup ramp"
DECLS
	REAL w
INIT
	w = 0
CONTINUOUS
	w' = (w_nom - w) / tau
	output.signal[1] = w
END COMPONENT



COMPONENT SmoothStartupIO
"Creates a log curve startup ramp"
PORTS
	IN analog_signal(n=1) input
	OUT analog_signal(n=1) output
DATA
	REAL tau = 10 UNITS u_s "Time constant for startup ramp"
	REAL w0 = 0 "Initial condition for output signal"
DECLS
	REAL w
INIT
	w = w0
CONTINUOUS
	w' = (input.signal[1] - w) / tau
	output.signal[1] = w
END COMPONENT


COMPONENT polynomial (
	INTEGER order = 1 "Polynomial order"
)
"Polynomial function of input signal"
PORTS
	IN analog_signal(n=1) in
	OUT analog_signal(n=1) out
DATA
	REAL a[order+1] = 1 "Polynomial coefficients → a[1]*x^n + a[2]*x^(n-1)..."
DECLS
	INTEGER i
	REAL x, y
DISCRETE
	ASSERT (order>=0) FATAL "Negative polynomial order not possible"
CONTINUOUS
	x = in.signal[1]
	SEQUENTIAL
		y = a[order+1] // constant term
		FOR (i IN 1,order)
			y = y + a[order+1-i] * x**i
		END FOR
	END SEQUENTIAL
	y = out.signal[1]
END COMPONENT



COMPONENT SetReset
PORTS
	IN analog_signal(n=1) input
	OUT analog_signal(n=1) output
DATA
	REAL highVal = 10
	REAL lowVal = 5
	REAL outputHigh = 1
	REAL outputLow = 0
	REAL outputInit = 1
INIT
	output.signal[1] = outputInit
DISCRETE
	WHEN (input.signal[1]>highVal) THEN
		output.signal[1] = outputHigh
	END WHEN
	WHEN (input.signal[1]<lowVal) THEN
		output.signal[1] = outputLow
	END WHEN
END COMPONENT



--------------------------------------------------------------------------------
// PR CONDITIONERS
--------------------------------------------------------------------------------
COMPONENT PRConditioner
"Converts deg C to 4-20 mA signal for PLC"
PORTS
	IN analog_signal(n=1) T_in
	OUT analog_signal(n=1) AI // Analog Input
DATA
	REAL T_low = -40 UNITS u_K "Lower limit"
	REAL T_high = 70 UNITS u_C "Higher limit"
CONTINUOUS
	AI.signal[1] = ((20-4) / (T_high-T_low)) * (T_in.signal[1]-T_low) + 4
END COMPONENT



COMPONENT PRCond_INT 
PORTS
	IN analog_signal(n=1) T_in
DATA
	REAL RangeMin = 0
	REAL RangeMax = 100
	REAL RawMin = 0
	REAL RawMax = 27648
DECLS
	REAL out
CONTINUOUS
	out = ((RawMax-RawMin) / (RangeMax-RangeMin)) * (T_in.signal[1]-RangeMin) + RawMin
END COMPONENT



--------------------------------------------------------------------------------
// MUXES AND DEMUXES
--------------------------------------------------------------------------------
COMPONENT Mux10 IS_A CONTROL.Mux
DECLS
	CLOSE n_in = 10
END COMPONENT



COMPONENT Mux11 IS_A CONTROL.Mux
DECLS
	CLOSE n_in = 11
END COMPONENT



COMPONENT Mux15 IS_A CONTROL.Mux
DECLS
	CLOSE n_in = 15
END COMPONENT



COMPONENT Mux20 IS_A CONTROL.Mux
DECLS
	CLOSE n_in = 20
END COMPONENT



COMPONENT Demux10 IS_A CONTROL.Demux
DECLS
	CLOSE n_out = 10
END COMPONENT



COMPONENT Demux12 IS_A CONTROL.Demux
DECLS
	CLOSE n_out = 12
END COMPONENT



COMPONENT Demux14 IS_A CONTROL.Demux
DECLS
	CLOSE n_out = 14
END COMPONENT



// BOOLEAN DEMUX
COMPONENT boolDemux5
PORTS
	IN bool_signal(n=5) s_in
	OUT bool_signal(n=1) s_out[5]
CONTINUOUS
	SEQUENTIAL
		s_out[1].signal[1] = s_in.signal[1]
		s_out[2].signal[1] = s_in.signal[2]
		s_out[3].signal[1] = s_in.signal[3]
		s_out[4].signal[1] = s_in.signal[4]
		s_out[5].signal[1] = s_in.signal[5]
	END SEQUENTIAL
END COMPONENT




COMPONENT OnOff (BOOLEAN ReverseAction = FALSE)
PORTS
	IN analog_signal(n=1) a
	OUT analog_signal(n=1) b
DATA
	REAL highlimit = 10
	REAL lowlimit = 0
	REAL outOn = 100
	REAL outOff = 0
DECLS
	BOOLEAN low
	BOOLEAN high
INIT
	low = FALSE
	high = FALSE
DISCRETE
	WHEN (a.signal[1]<lowlimit) THEN
		low = TRUE
		high = FALSE
	END WHEN
	WHEN (a.signal[1]>highlimit) THEN
		low = FALSE
		high = TRUE
	END WHEN
CONTINUOUS
SEQUENTIAL
	IF (ReverseAction) THEN
		IF (low) THEN
			b.signal[1] = outOn
		ELSE
			b.signal[1] = outOff
		END IF
	ELSE
		IF (high) THEN
			b.signal[1] = outOn
		ELSE
			b.signal[1] = outOff
		END IF
	END IF
END SEQUENTIAL
END COMPONENT



COMPONENT FlowSplitter
PORTS
	IN analog_signal(n=1) inlet
	OUT analog_signal(n=1) outlet1
	OUT analog_signal(n=1) outlet2
DATA
	REAL out1 = 0.7 "Portion to outlet 1 as fraction of 1"
CONTINUOUS
	outlet1.signal[1] = out1*inlet.signal[1]
	outlet2.signal[1] = (1-out1)*inlet.signal[1]
END COMPONENT
