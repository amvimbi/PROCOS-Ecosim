/*-----------------------------------------------------------------------------------------
 LIBRARY: CO2
 FILE: THERMAL
 AUTHOR: Viren Bhanot
 COMPANY: CERN
 DESCRIPTION: Thermal components such as wall
 CREATION DATE: 10/05/2017
-----------------------------------------------------------------------------------------*/
USE THERMAL
USE CRYOLIB
USE MATH
USE PLC
USE PORTS_LIB



/*
 * INDEX
 * 1. Abstract Components
 *	1. Wall
 *	2. CylinderWall
 *	3. AnnularWall
 * 	4. GenericWall
 *	5. Wall_axial
 *	6. CylinderWall_axial
 * 2. Cylinder Wall Components
 *	1. CartridgeHeater
 *	2. CartridgeHeaterCtrl
 *	3. CartridgeHeaterPWM
 *	4. CartridgeHeater_axial
 * 3. Annular Wall Components
 *	1. AnnularWall_onePort
 *	2. AnnularWall_twoPort
 *	3. AnnularWall_naturalConvection
 *	4. LumpedShell
 *	5. LumpedShell_insulated
 *	6. SuctionChamberWall
 */



/*
 * NOTE ON SIGN CONVENTION:
 * EcosimPro sign convention for ports is strange: Fluxes entering an IN port or leaving 
 * an OUT port are positive. In addition, if you connect two IN ports to each other, 
 * EcosimPro will INVERT the sign of the heat flux, by itself:
 *	ComponentA.tp_in.q[1] = -ComponentB.tp_in.q[1]
 * This can lead to unwanted behaviour and therefore, the following rules
 * should be kept in mind:
 *
 	1. For thermal components, the energy balance should always be written as:
 *	M*cpw*Tw' = (tp_in.q[1] - tp_out.q[i]), aka (InPort - OutPort)
 *	DO NOT rely on the signs of q to take care of themselves. In Modelica,
 *	for example, you would write: M*cpw*der(Tw) = portA.Qflow + portB.Qflow
 *
 * 	2. For the refrigerant streams, the energy balance should be:
 *	(...) = mhIn - mhOut + tp_in.q - tp_out.q
 *
 * 	3. Lastly, in the refrigerant streams, use these for heat transfer:
 *	tp_in.q = HTC*As*(tp_in.Tk - T) // positive when wall is hotter
 *	tp_out.q = HTC*As*(T - tp_out.Tk) // positive when Tref > Twall
 *
 * The two equations above handle the IN port flux being positive, OUT flux 
 * being negative. If the inlet thermal port has Tk > T
 * then heat will flow INTO refrigerant. Therefore, positive from inlet port
 * On the other hand, if outlet port Tk > T, then heat will still flow into 
 * refrigerant, but EcosimPro will be expecting an outflow of heat. Therefore,
 * a negative sign will result for tp_out.q.
 *
 * (Note: behaviour of connecting two IN ports has not been checked)
 */



--------------------------------------------------------------------------------
// ABSTRACT COMPONENTS
--------------------------------------------------------------------------------
ABSTRACT COMPONENT Wall (INTEGER n = 3 "Number of control volumes")
"Base class for thermal walls. Assigns temperature properties and initial temperatures"
// Does not assign the shape of the wall
// Does not define the number of ports
DATA
	ENUM PipeMat mat = SS_304 "Wall material"
	REAL Tw0_i = 255 UNITS u_K "Initial wall temperature - first segment"
	REAL Tw0_o = 255 UNITS u_K "Initial wall temperature - last segment"
DECLS
	REAL Mw UNITS u_kg "Metal mass"
	REAL Tw[n] UNITS u_K
	PRIVATE INTEGER tmp[n] //dummy
OBJECTS
	WallRecord wallGeo[n] //does not specify shape yet
INIT
	linspace(Tw0_i,Tw0_o,n,Tw,1)
DISCRETE
	ASSERT (mat!=THERMAL.None) WARNING "Material type 'None' isn't supported"
CONTINUOUS
	EXPAND_BLOCK (i IN 1,n)
		wallGeo[i].setTempProps(Tw[i],tmp[i])
	END EXPAND_BLOCK
<Mw> 	Mw = SUM(i IN 1,n; wallGeo[i].M)
END COMPONENT



ABSTRACT COMPONENT CylinderWall IS_A CO2.Wall
"Solid cylindrical metal wall"
DATA
	REAL D = 0.1 UNITS u_m "Diameter"
	REAL L = 1 UNITS u_m
INIT
	FOR (i IN 1,n)
		wallGeo[i].setCylindricalProperties(mat,D,L/n)
	END FOR
END COMPONENT



ABSTRACT COMPONENT AnnularWall IS_A CO2.Wall
"Annular wall base component"
DATA
	REAL Do = 0.1 UNITS u_m "Outer diameter"
	REAL Th = 0.01 UNITS u_m "Thickness"
	REAL L = 1 UNITS u_m
INIT
	FOR (i IN 1,n)
		wallGeo[i].setAnnularProperties(mat,Do,Th,L/n)
	END FOR
END COMPONENT



COMPONENT GenericWall (
	INTEGER n = 5 "Number of segments"
)
"Generic wall component without shape specification"
// No shape required, just wall thermal mass
// for BPHX
PORTS
	IN thermal(n=n) tp_in
	OUT thermal(n=n) tp_out
DATA
	ENUM PipeMat mat = SS_316 "Wall material"
	REAL Mw = 14.1 UNITS u_kg "Wall mass"
	REAL Tw0_i = 255 UNITS u_K "Initial wall temperature - first segment"
	REAL Tw0_o = 255 UNITS u_K "Initial wall temperature - last segment"
DECLS
	REAL Cp[n] UNITS u_J_kgK "Specific heat capacity"
	REAL Tw[n] UNITS u_K
INIT
	linspace(Tw0_i,Tw0_o,n,Tw,1)
CONTINUOUS
	EXPAND_BLOCK (i IN 1,n)
		Cp[i] = FunVarSolidProp(mat,SpecificHeat,Tw[i],1,ier)
		Tw[i]' = (tp_in.q[i] - tp_out.q[n-i+1])/(Mw/n*Cp[i])
		tp_in.Tk[i] = Tw[i]
		tp_out.Tk[n-i+1] = Tw[i] //countercurrent flow
	END EXPAND_BLOCK
END COMPONENT



COMPONENT GenericShell
"Generic shell component without shape specification"
// No shape required, just wall thermal mass
// for sucton chamber
PORTS
	IN thermal(n=1) tp_in
	IN analog_signal(n=1) T_amb "Ambient temperature"
	IN analog_signal(n=1) Q_heat "Heat leak from hydraulics and motor"
DATA
	ENUM PipeMat mat = SS_316 "Wall material"
	REAL As = 1 UNITS u_m2 "Outer surface area"
	REAL Mw = 14.1 UNITS u_kg "Wall mass"
	REAL alpha_a = 10 UNITS u_W_m2K "Convective heat transfer coefficient"
	REAL Tw0 = 255 UNITS u_K "Initial wall temperature - first segment"
DECLS
	REAL Cp UNITS u_J_kgK "Specific heat capacity"
	REAL Tw UNITS u_K
	REAL Qa_conv UNITS u_W "Convective heat transfer"
INIT
	Tw = Tw0
CONTINUOUS
	Qa_conv = alpha_a * As * (Tw - T_amb.signal[1]) // positive if wall hotter
	Cp = FunVarSolidProp(mat,SpecificHeat,Tw,1,ier)
	Tw' = (tp_in.q[1] + Q_heat.signal[1] - Qa_conv) / (Mw*Cp)
	tp_in.Tk[1] = Tw
END COMPONENT



ABSTRACT COMPONENT Wall_axial (
	INTEGER nz = 5 "Number of axial nodes",
	INTEGER nr = 3  "Number of radial nodes"
)
"Cylinder wall with axial discretization"
PORTS
	OUT thermal (n=nz) tp_out
DATA
	REAL L = 1 UNITS u_m "Length"
	ENUM PipeMat mat = SS_304 "Material"
	REAL Tw0 = 255 UNITS u_K "Initial temperature [uniform for simplicity]"
DECLS
	REAL Mw
	REAL Tw[nr,nz] UNITS u_K
	PRIVATE INTEGER tmp[nr,nz]
OBJECTS
	WallRecord wallGeo[nr,nz]
INIT
	FOR (i IN 1,nr)
		FOR (j IN 1,nz)
			Tw[i,j] = Tw0
		END FOR
	END FOR
CONTINUOUS
	EXPAND (i IN 1,nr) EXPAND (j IN 1,nz) wallGeo[i,j].setTempProps(Tw[i,j],tmp[i,j])
<Mw> 	Mw = SUM(i IN 1,nr; SUM (j IN 1,nz; wallGeo[i,j].M))
END COMPONENT



ABSTRACT COMPONENT CylinderWall_axial IS_A Wall_axial
DATA
	REAL D = 0.01 UNITS u_m "Diameter"
DECLS
	DISCR REAL dr UNITS u_m "Radial nodal distance"
	DISCR REAL dz UNITS u_m "Axial nodal distance"
INIT
	dr = (D/2)/nr
	dz = L/nz
	FOR (i IN 1,nr)
		FOR (j IN 1,nz)
			IF (i==1) THEN
				wallGeo[i,j].setCylindricalProperties(mat,2*dr,dz)
			ELSE
				wallGeo[i,j].setAnnularProperties(mat,2*i*dr,dr,dz)
			END IF
		END FOR
	END FOR
CONTINUOUS
END COMPONENT



--------------------------------------------------------------------------------
// CYLINDRICAL WALL COMPONENTS
--------------------------------------------------------------------------------
COMPONENT CartridgeHeater IS_A CylinderWall
"Note: Heat load analog_signal port has n nodes, therefore, the connected \
inputs are already discretized into n volumes."
PORTS
	IN analog_signal(n=1) Q_in
	OUT thermal(n=n) tp_out
DECLS
	REAL T_avg UNITS u_K "Average temperature of heater"
CONTINUOUS
	EXPAND_BLOCK(i IN 1,n)
		Tw[i]' = (Q_in.signal[1]/n - tp_out.q[i])/wallGeo[i].C
		tp_out.Tk[i] = Tw[i]
	END EXPAND_BLOCK
	T_avg = SUM(i IN 1,n ; Tw[i]) / n
END COMPONENT



COMPONENT CartridgeHeaterCtrl IS_A CylinderWall
PORTS
	IN analog_signal(n=1) K "Control signal in % (0-100)"
	OUT thermal(n=n) tp_out
DATA
	REAL Q_max = 1000 UNITS u_W "Maximum heater power"
DECLS
	REAL Q_sig UNITS u_W "Input heater power"
	REAL T_avg UNITS u_K "Average temperature of heater"
CONTINUOUS
	EXPAND_BLOCK(i IN 1,n)
		Tw[i]' = ((K.signal[1]/100)*(Q_max/n) - tp_out.q[i]) / wallGeo[i].C
		tp_out.Tk[i] = Tw[i]
	END EXPAND_BLOCK
	Q_sig = K.signal[1]/100*Q_max
	T_avg = SUM(i IN 1,n ; Tw[i])/n
END COMPONENT


COMPONENT CartridgeHeaterCtrlTout IS_A CartridgeHeaterCtrl (
	BOOLEAN ToutInCelsius = TRUE
)
PORTS
	OUT analog_signal(n=1) Tout
CONTINUOUS
	Tout.signal[1] = ZONE (ToutInCelsius) T_avg - 273.15 OTHERS T_avg
END COMPONENT



COMPONENT CartridgeHeaterPWM IS_A CylinderWall
PORTS
	IN bool_signal(n=1) DO "Control signal (Boolean)"
	OUT thermal(n=n) tp_out
DATA
	REAL Q_max = 1000 UNITS u_W "Maximum heater power"
DECLS
	REAL Q_act UNITS u_W "Actual Heating Power"
	REAL T_avg UNITS u_K "Average heater temperature"
CONTINUOUS
	Q_act = ZONE (DO.signal[1]) Q_max OTHERS 0
	EXPAND_BLOCK(i IN 1,n)
		Tw[i]' = (Q_act/n - tp_out.q[i])/wallGeo[i].C
		tp_out.Tk[i] = Tw[i]
	END EXPAND_BLOCK
	T_avg = SUM(i IN 1,n ; Tw[i]) / n
END COMPONENT



COMPONENT CartridgeHeater_axial IS_A CylinderWall_axial
PORTS
	IN analog_signal(n=1) K "Control signal in % [0-100]"
DATA
	REAL Q_max = 1000 UNITS u_W "Maximum heater power"
DECLS
	REAL Q[nr,nz] UNITS u_W "Heat transfer rate"
	REAL C[nr,nz] UNITS u_J_K "Heat capacity"
CONTINUOUS
	EXPAND_BLOCK (i IN 1,nr)
		EXPAND_BLOCK (j IN 1,nz)
			
		END EXPAND_BLOCK
	END EXPAND_BLOCK
END COMPONENT



--------------------------------------------------------------------------------
// ANNULAR WALL COMPONENTS
--------------------------------------------------------------------------------
COMPONENT AnnularWall_onePort IS_A AnnularWall
"1D wall. Thermal port on inner wall, insulated on the outside"
PORTS
	IN thermal(n=n) tp_in
INIT
	FOR (i IN 1,n) 
		wallGeo[i].EndSurfaces = 0
	END FOR
CONTINUOUS
	EXPAND_BLOCK(i IN 1,n)
		Tw[i]' = tp_in.q[i]/wallGeo[i].C
		tp_in.Tk[i] = Tw[i]
	END EXPAND_BLOCK
END COMPONENT



COMPONENT AnnularWall_twoPort IS_A AnnularWall (
	BOOLEAN Countercurrent = TRUE "FALSE if parallel flow"
)
"1D wall. Thermal ports on both inner and outer wall."
PORTS
	IN thermal(n=n) tp_in
	OUT thermal(n=n) tp_out
CONTINUOUS
	IF (Countercurrent) INSERT
		EXPAND_BLOCK (i IN 1,n)
			Tw[i]' = (tp_in.q[i] - tp_out.q[n-i+1])/wallGeo[i].C
			tp_in.Tk[i] = Tw[i]
			tp_out.Tk[n-i+1] = Tw[i]
		END EXPAND_BLOCK
	ELSE
		EXPAND_BLOCK(i IN 1,n)
			Tw[i]' = (tp_in.q[i] - tp_out.q[i])/wallGeo[i].C
			tp_in.Tk[i] = Tw[i]
			tp_out.Tk[i] = Tw[i]
		END EXPAND_BLOCK
	END IF
END COMPONENT



COMPONENT AnnularWall_naturalConvection IS_A AnnularWall (
	BOOLEAN accountForRadiation = TRUE "Recommended true if high surface temperatures",
	BOOLEAN ConstantAirHTC = TRUE "If FALSE, uses Churchill-Chu 1975 natural convective correlation"
)
"Annular wall with natural convective heat transfer on outer surface"
/*
 * Significance of terms:
 *   - Prandtl compares thickness of momentum vs thermal boundary layer.
 *   - Grashof number gives ratio of buoyancy to viscous forces on a fluid.
 *     (higher Gr means the fluid has an easier time rising)
 *   - Rayleigh number is product of Prandtl and Grashof. Above a value, 
 *     convection dominates. Otherwise, mostly conduction to the air.
 */
PORTS
	IN analog_signal(n=1) T_amb "Ambient air temperature"
	IN thermal(n=n) tp_in
	OUT analog_signal(n=1) Q_air "Total air-side heat transfer rate"
DATA
	REAL alpha_a0 = 10 UNITS u_W_m2K "Constant convective heat transfer coefficient"
	REAL epsilon = 0.85 "Wall emissivity"
DECLS
	CONST ENUM ThFluids medum = Air // THERMAL library has its own, simpler air calculations (FuncsPureFluids.el)
	CONST REAL P_air = 101325 UNITS u_Pa "Air pressure (treated constant)"
	CONST REAL sigma = 5.67E-8 "Stefan-Boltzmann constant"
	REAL Gr_air[n] "Grashof number"
	REAL Ra_air[n] "Rayleigh number"
	REAL alpha_a[n] UNITS u_W_m2K "Air-side heat transfer coefficient"
	REAL Qa_conv[n] UNITS u_W  "Air side convective heat transfer"
	REAL Qa_rad[n] UNITS u_W "Air side radiative heat transfer"
	REAL Qa[n] UNITS u_W
	REAL Qa_total UNITS u_W
	INTEGER ier,iex,iey,dummy
OBJECTS
	AirProperties air
INIT
	FOR (i IN 1,n) 
		wallGeo[i].EndSurfaces = 0 //by default, assumed open at both ends
	END FOR
CONTINUOUS
	air.setProps(P_air,T_amb.signal[1],dummy)
	EXPAND_BLOCK(i IN 1,n)
		Gr_air[i] = 9.81 * air.beta * wallGeo[i].Do**3 * abs(Tw[i]-T_amb.signal[1]) / (air.mu/air.rho)**2 //Grashof for pipes
		Ra_air[i] = Gr_air[i] * air.Pr
		IF (ConstantAirHTC) INSERT
			alpha_a[i] = alpha_a0
		ELSE
			alpha_a[i] = ChurchillChu(Ra_air[i],air.Pr,air.k,wallGeo[i].Do)
		END IF
		Qa_conv[i] = alpha_a[i] * wallGeo[i].As_o * (Tw[i] - T_amb.signal[1])
		IF (accountForRadiation) INSERT
			Qa_rad[i] = sigma * epsilon * wallGeo[i].As_o * (Tw[i]**4 - T_amb.signal[1]**4)
		ELSE
			Qa_rad[i] = 0
		END IF
		Qa[i] = Qa_conv[i] + Qa_rad[i]
		Tw[i]' = (tp_in.q[i] - Qa[i]) / wallGeo[i].C
		tp_in.Tk[i] = Tw[i]
	END EXPAND_BLOCK
	Qa_total = SUM (i IN 1,n; Qa[i])
	Q_air.signal[1] = Qa_total
END COMPONENT




COMPONENT AnnularWall_naturalConvectionTout IS_A AnnularWall_naturalConvection (
	BOOLEAN ToutInCelsius = TRUE
)
"Outputs temperature at the outlet port"
PORTS
	OUT analog_signal(n=1) Tout
CONTINUOUS
	Tout.signal[1] = ZONE (ToutInCelsius) Tw[n] - 273.15 OTHERS Tw[n]
END COMPONENT


COMPONENT LumpedShell IS_A AnnularWall_naturalConvection
"Lumped shell component. Closed at both ends"
DATA
	REAL Tw0 = 255 UNITS u_K "Initial wall temperature"
DECLS
	CLOSE n = 1 //shells always considered lumped
	CLOSE Tw0_i = 0 //we use a single Tw0 instead of these two
	CLOSE Tw0_o = 0 //easier for the user to see one Tw instead of two for a lumped component
INIT
	wallGeo[1].EndSurfaces = 2
	Tw[1] = Tw0
END COMPONENT



COMPONENT LumpedShell_insulated IS_A AnnularWall_onePort
"Lumped shell, closed at both ends and perfectly insulated"
DATA
	REAL Tw0 = 255 UNITS u_K "Initial wall temperature"
DECLS
	CLOSE n = 1
	CLOSE Tw0_i = 0
	CLOSE Tw0_o = 0
INIT
	wallGeo[1].EndSurfaces = 2
	Tw[1] = Tw0
END COMPONENT



COMPONENT SuctionChamberWall IS_A LumpedShell
"Suction chamber wall. Open at one end"
INIT
	wallGeo[1].EndSurfaces = 1
END COMPONENT