/*-----------------------------------------------------------------------------------------
 LIBRARY: CO2
 FILE: CLASSES
 CREATION DATE: 28/07/2017

 DESCRIPTION:
 Classes in this library are employed in the same manner as Modelica Records.
 Specifically, they are used to organize the handling of geometries and 
 material properties. This is good practice since passing these enormous
 number of parameters is simplifed to just passing an instance of the class
 in a function.
 
 The downside of using classes is that they aren't visible in the simulation
 monitor, and thus these variables cannot be plotted or exported. To circumvent
 this, dummy variables must be created in the Component models that are assigned
 the values we want to plot. For things like areas and volumes.
-----------------------------------------------------------------------------------------*/
USE MATH
USE CRYOLIB
USE THERMAL



/*
 * INDEX
 * 1. BASIC
 *	1. SaturationProperties
 *	2. ThermodynamicState
 *	3. AirProperties
 * 2. REFRIGERANT GEOMETRY RECORDS
 *	1. RefrigerantRecord
 *	2. CylinderRefRecord
 *	3. AnnulusRefRecord
 * 3. WALL PROPERTIES RECORDS
 *	1. WallRecord
 *	2. CylinderWallRecord
 *	3. AnnularWallRecord
 * 4. PLATE FIN HEAT EXCHANGER RECORDS
 *	1. AirFlowPassageRecord
 *	2. TubeFinGeometry
 *	3. TubeFinMaterialProperties
 */



--------------------------------------------------------------------------------
// BASIC CLASSES
--------------------------------------------------------------------------------
CLASS SaturationProperties
"Record of thermodynamic saturation properties"
DECLS
	REAL cp_l, cp_g UNITS u_J_kgK "Specific heat"
	REAL h_l, h_g UNITS u_J_kg "Enthalpy"
	REAL k_l, k_g UNITS u_W_mK "Thermal conductivity"
	REAL mu_l, mu_g UNITS u_Pas "Dynamic viscosity"
	REAL rho_l, rho_g UNITS u_kg_m3 "Density"
END CLASS



CLASS ThermodynamicState
"Record of thermodynamic state properties"
// This class is used to calculate in one go all the state properties of the
// fluid in a control volume.
DECLS
	ENUM ChemName fluid
	REAL P UNITS u_bar
	REAL h UNITS u_J_kg
	REAL rho UNITS u_kg_m3
	REAL u UNITS u_J_kg
	REAL s UNITS u_J_K
	REAL T UNITS u_K
	REAL x
	ENUM Phase phase
METHODS
	METHOD NO_TYPE setState_ph (
		IN REAL P_val,
		IN REAL h_val,
		OUT INTEGER tmp=DUMMY_INTEGER
	)
	"Set state parameters for P,h state variables"
	DECLS
		INTEGER tmp1,tmp2,tmp3
	BODY
		P = P_val
		h = h_val
		getThermodynamicState_ph(fluid,P,h,rho,s,T,u,x,phase,tmp1,tmp2,tmp3)
	END METHOD
	
	METHOD NO_TYPE setState_ru(
		IN REAL rho_val,
		IN REAL u_val,
		OUT INTEGER tmp=DUMMY_INTEGER
	)
	"Set state parameters for rho,u state variables"
	DECLS
		INTEGER tmp1,tmp2,tmp3
	BODY
		rho = rho_val
		u = u_val
		getThermodynamicState_ru(fluid,rho,u,P,h,s,T,x,phase,tmp1,tmp2,tmp3)
	END METHOD
END CLASS



CLASS AirProperties
DECLS
	REAL beta "Isothermal compressibilty" 
	REAL Cp UNITS u_J_kgK
	REAL k UNITS u_W_mK
	REAL mu UNITS u_Pas
	REAL rho UNITS u_kg_m3
	REAL Pr "Prandtl"
METHODS
	METHOD NO_TYPE setProps (
		IN REAL P,
		IN REAL T,
		OUT INTEGER tmp=DUMMY_INTEGER
	)
	BODY
		Fluid_prop_vs_pT(THERMAL.Air,P,T,Cp,mu,rho,k,ier)
		beta = Fluid_beta_vs_T(THERMAL.Air,T,ier)
		Pr = Cp * mu / k
	END METHOD
END CLASS



--------------------------------------------------------------------------------
// REFRIGERANT GEOMETRY RECORD
--------------------------------------------------------------------------------
CLASS RefGeometryRecord
"Generic class to hold geometric properties of a control volume"
// Note that the outer diameter refers to the outer edge of the REFRIGERANT volume,
// and not of the metal shell that exists outside. Thus, the INNER diameter of a 
// cylindrical tube would actual be the outer diameter here.
DECLS
	REAL Do UNITS u_m2 "Outer diameter"
	REAL Di UNITS u_m2 "Inner diameter, only non zero for annular geometries"
	REAL Th UNITS u_m2 "Thickness"
	REAL L UNITS u_m2 "Height/Length"
	REAL Ap UNITS u_m2 "Port area"
	REAL Ac UNITS u_m2 "Cross sectional area"
	REAL As_i UNITS u_m2 "Inner surface area"
	REAL As_o UNITS u_m2 "Outer surface area"
	REAL As UNITS u_m2 "Total surface area"
	REAL V UNITS u_m3 "Volume"
	INTEGER EndSurfaces
METHODS
	METHOD NO_TYPE setCylindricalGeometry(IN REAL Do_val, IN REAL L_val)
	BODY
		Do = Do_val
		Di = 0
		Th = 0
		L = L_val
		Ap = MATH.PI*(Do/2)**2
		Ac = Ap // currently, port area assumed to be cross sectional area -> only valid for cylindrical pipes
		As_i = 0
		As_o = MATH.PI*Do*L
		As = As_o + EndSurfaces*Ap
		V = Ap*L
	END METHOD
	
	METHOD NO_TYPE setAnnularGeometry(IN REAL Do_val, IN REAL Th_val, IN REAL L_val)
	BODY
		Do = Do_val
		Th = Th_val
		L = L_val
		Di = Do - 2*Th
		Ap = MATH.PI*((Do/2)**2 - (Di/2)**2)
		Ac = Ap
		As_i = MATH.PI*Di*L
		As_o = MATH.PI*Do*L
		As = As_i + As_o
		V = Ap*L
	END METHOD
	
	METHOD NO_TYPE setAccumulatorGeometry(IN REAL Do_val, IN REAL L_val)
	"Set geometry of an accumulator: Cylindrical tube with ellipsoidal caps of half the tube diameter"
	DECLS
		REAL Vcap // volume of the two caps combined (bottom and top)s
		REAL Acap // area of the two caps combined (bottom and top)
	BODY
		EndSurfaces = 2
		Do = Do_val // refrigerant volume outer diameter
		Di = 0 // Di is only non-zero for annular geometries
		Th = 0 // no concept of thickness for cylindrical geometry
		L = L_val
		Ac = MATH.PI*(Do/2)**2 // cross section area doesn't have much meaning for curved ends.
		Ap = Ac // TODO: should be different to C/S area
		As_i = 0
		Acap = MATH.PI *(((Do*Do)**1.6 + (Do*Do/2)**1.6 + ((Do/2)*Do)**1.6)/3)**(1/1.6)
		As_o = MATH.PI*Do*(L-Do/2) + Acap
		As = As_o
		Vcap = 4/3 * MATH.PI * (Do*Do*Do/2)/8 // ellipsoid volume is 4/3.pi.(a.b.c) where abc are the three radii
		V = Ac*(L-Do/2) + Vcap
		// we have the full ellipsoid, aka half a cap at each end.
	END METHOD
END CLASS



--------------------------------------------------------------------------------
// WALL RECORDS
--------------------------------------------------------------------------------
CLASS WallRecord
"Geometry and material properties for thermal walls"
/*
 * NOTE: The annular wall doesn't mean that it is open on both ends. It can 
 * still be closed at both ends. For instance, an accumulator shell is 
 * annular, in that it is hollow, but is still closed at both ends. By contrast
 * a cylindrical wall will always be closed at both ends, since an open cylinder 
 * is by definition annular (has some metal thickness).
 */
DECLS
	ENUM Material mat
	REAL cpw UNITS u_J_kgK "Specific heat capacity"
	REAL C UNITS u_J_K "Heat capacity"
	REAL k UNITS u_W_mK "Thermal conductivity" // currently unused: all models are 1D
	REAL M UNITS u_kg
	REAL rho UNITS u_kg_m3
	REAL Do UNITS u_m "Outer diameter"
	REAL Di UNITS u_m "Inner diameter"
	REAL Th UNITS u_m "Thickness"
	REAL L UNITS u_m "Height/Length"
	REAL As_i UNITS u_m2 "Inner curved surface area"
	REAL As_o UNITS u_m2 "Outer curved surface area"
	REAL As UNITS u_m2 "Surface area"
	REAL Ap UNITS u_m2 "Cross-sectional (port) area"
	REAL V UNITS u_m3 "Volume"
	INTEGER EndSurfaces "Number of closed end surfaces"
METHODS
	METHOD NO_TYPE setMaterialProps(IN ENUM Material matval)
	"Set material properties: density and mass"
	// You will never yourself need to call this function. This will be 
	// called from the setProps function of the child classes.
	BODY
		mat = matval
		rho = FunConstSolidProp(mat,Density,1)
		M = rho*V
	END METHOD

	METHOD NO_TYPE setTempProps (IN REAL Tw, OUT INTEGER tmp=DUMMY_INTEGER)
	"Set temperature dependent properties: conductivity and specific heat"
	/* 
	 * A function in EcosimPro must always have a RETURN statement in the
	 * CONTINUOUS block. This function doesn't need to return anything
	 * To solve this, EcosimPro has DUMMY_INTEGER as an OUT argument
	 * and I assign it to a tmp var in the Component model.
	 */
	BODY
		k = FunVarSolidProp(mat,Conductivity,Tw,1,ier) // this "ier" is a global variable in CRYOLIB
		cpw = FunVarSolidProp(mat,SpecificHeat,Tw,1,ier)
		C = cpw*M
	END METHOD
	
	METHOD NO_TYPE setCylindricalProperties (IN ENUM Material matval, IN REAL D_val, IN REAL L_val)
	BODY
		EndSurfaces = 2
		Do = D_val
		Di = 0
		Th = Do/2
		L = L_val
		Ap = PI * (Do/2)**2
		As_i = 0
		As_o = PI*Do*L + EndSurfaces*Ap // cylinder will always be closed
		As = As_o
		V = MATH.PI * (Do/2)**2 * L
		setMaterialProps(matval)
	END METHOD
	
	METHOD NO_TYPE setAnnularProperties (IN ENUM Material matval, IN REAL Do_val, IN REAL Th_val, IN REAL L_val)
	BODY
		Do = Do_val
		Th = Th_val
		Di = Do - 2*Th
		L = L_val
		Ap = PI * ((Do/2)**2 - (Di/2)**2)
		As_i = PI*Di*L + EndSurfaces*PI*(Di/2)**2
		As_o = PI*Do*L + EndSurfaces*PI*(Do/2)**2
		As = As_i + As_o
		V = Ap * L
		setMaterialProps(matval)
	END METHOD
END CLASS



--------------------------------------------------------------------------------
// ROUND-TUBE PLATE FIN HEAT EXCHANGER RECORDS
--------------------------------------------------------------------------------
ENUM PlateFinType = {BareTube,WithCollar,WithoutCollar}



CLASS AirFlowPassageRecord
"Record of propreties of air flow passage through a round-tube, plate-fin HX"
DECLS
	INTEGER nTube
	INTEGER nBank
	BOOLEAN staggered
	ENUM PlateFinType ftype
	REAL Do UNITS u_m "Tube outer diameter"
	REAL Dc UNITS u_m "Diameter of fin collar"
	REAL L UNITS u_m "Tube length"
	REAL A_pri UNITS u_m2 "Primary flow area"
	REAL A_sec UNITS u_m2 "Secondary flow area"
	REAL Ao UNITS u_m2 "Effective air flow area"
	REAL A_min UNITS u_m2
	REAL eta_fin "Fin efficiency"
	REAL Pt UNITS u_m "Transverse fin spacing"
	REAL Pl UNITS u_m "Longitudinal fin spacing"
	REAL Th_f UNITS u_m "Fin thickness"
	REAL FPI "Fins per inch"
	REAL FPM "Fins per meter"
	REAL Fp UNITS u_m "Fin spacing"
	REAL mu
	REAL HX_l UNITS u_m "Heat exchanger length"
	REAL HX_d UNITS u_m "Heat exchanger depth"
	REAL HX_h UNITS u_m "Heat exchanger height"
METHODS
	METHOD NO_TYPE setProps (
		IN INTEGER nTube_val,
		IN INTEGER nBank_val,
		IN BOOLEAN staggered_val,
		IN REAL Do_val,
		IN REAL L_val,
		IN REAL Pt_val,
		IN REAL Pl_val,
		IN ENUM PlateFinType ftype_val,
		IN REAL FPI_val,
		IN REAL Th_f_val,
		IN REAL K_f,
		IN REAL alpha_a
	)
	"Sets constant geometry parameters for two dimensional air flow"
	DECLS
		REAL r,Re,phi,m,x,y,z // temporary variables
	BODY
		nTube = nTube_val
		nBank = nBank_val
		staggered = staggered_val
		Do = Do_val
		L = L_val
		Pt = Pt_val
		Pl = Pl_val
		ftype = ftype_val
		FPI = FPI_val
		HX_h = nTube*Pt
		HX_l = nBank*Pl
		HX_d = L
		IF (ftype==BareTube) THEN
			FPM = 0
			Th_f = 0
			Fp = 0
			A_pri = PI*Do*L*nTube*nBank
			A_sec = 0
			IF (staggered) THEN
				x = 0.5*(Pt-Do)
				y = ((Pt/2)**2 + Pl**2)**0.5 - Do
				z = 2*min(x,y)
				A_min = ((HX_h/Pt-1)*z + (Pt-Do))*HX_d
			ELSE
				A_min = (Pt-Do)*L*nTube
			END IF
			eta_fin = 0
		ELSE
			Th_f = Th_f_val
			FPM = FPI/0.0254
			Fp = 1/(FPM + Th_f)
			IF (ftype==WithCollar) THEN
				Dc = Do + 2*Th_f
			ELSE // Without collar
				Dc = Do
			END IF
			A_pri = PI*Dc*L*(1-Th_f*FPM)*nTube*nBank
			A_sec = 2*(HX_l*HX_h - PI*Dc**2*nTube*nBank/4)*FPM*L + 2*HX_h*Th_f*FPM*L
			IF (staggered) THEN
				x = 0.5*(Pt-Dc)*(1-Th_f*FPM)
				y = ((Pt/2)**2 + Pl**2)**0.5 - Dc - (Pt-Dc)*Th_f*FPM
				z = 2*min(x,y)
				A_min = ((HX_h/Pt-1)*z + (Pt-Dc)*(1-Th_f*FPM))*L
			ELSE
				A_min = ((Pt-Dc)*L*(1 - Th_f*FPM))*HX_h/Pt
			END IF
			r = Do/2 + Th_f
			Re = r * 1.28 * Pt/(Do + 2*Th_f) * (Pl/Pt - 0.2)**0.5
			phi = (Re/r - 1) * (1 + 0.35*log(Re/r))
			m = sqrt(2*alpha_a/(K_f*Th_f))
			eta_fin = tanh(m*r*phi) * cos(0.1*m*r*phi) / (m*r*phi)
		END IF
		Ao = A_pri + eta_fin*A_sec
	END METHOD
END CLASS



CLASS TubeFinGeometry
"Record of tube and fin geometry (not material properties)"
// basically needs to return volumes for material properties to use
DECLS
	INTEGER nTube "Number of tubes per bank"
	INTEGER nBank "Number of banks of tubes"
	BOOLEAN staggered "False means inline tubes"
	ENUM PlateFinType ftype
	REAL Do UNITS u_m "Tube outer diameter"
	REAL Di UNITS u_m "Tube inner diameter"
	REAL Dc UNITS u_m "Collar diameter"
	REAL Th_t UNITS u_m "Tube wall thickness"
	REAL L UNITS u_m "Tube length"
	REAL Th_f UNITS u_m "Fin thickness"
	REAL FPI "Fins per inch"
	REAL FPM "Fins per meter"
	REAL Fp UNITS u_m "Fin spacing"
	REAL Pt UNITS u_m "Tube vertical spacing"
	REAL Pl UNITS u_m "Tube horizontal spacing"
	REAL A_pri UNITS u_m2 "Primary area"
	REAL A_sec UNITS u_m2 "Secondary area"
	REAL V_t UNITS  u_m3 "Tube material volume"
	REAL V_f UNITS u_m3 "Fin material volume"
	REAL HX_h UNITS u_m
	REAL HX_d UNITS u_m
	REAL HX_l UNITS u_m
METHODS
	METHOD NO_TYPE setProps (
		IN INTEGER nTube_val,
		IN INTEGER nBank_val,
		IN BOOLEAN staggered_val,
		IN REAL Do_val,
		IN REAL Th_t_val,
		IN REAL L_val,
		IN REAL Pt_val,
		IN REAL Pl_val,
		IN ENUM PlateFinType ftype_val,
		IN REAL FPI_val,
		IN REAL Th_f_val
	)
	"Sets geometry and material properties for the entire heat exchanger"
	BODY
		nTube = nTube_val
		nBank = nBank_val
		staggered = staggered_val
		Do = Do_val
		Th_t = Th_t_val
		L = L_val
		Pt = Pt_val
		Pl = Pl_val
		ftype = ftype_val
		FPI = FPI_val
		Th_f = Th_f_val
		Di = Do - 2*Th_t
		HX_h = nTube*Pt
		HX_l = nBank*Pl
		HX_d = L
		IF (ftype==BareTube) THEN
			FPM = 0
			Th_f = 0
			Fp = 0
			A_pri = PI*Do*L*nTube*nBank
			A_sec = 0
		ELSE
			FPM = FPI/0.0254
			Fp = 1/(FPM + Th_f)
			IF (ftype==WithCollar) THEN
				Dc = Do + 2*Th_f
			ELSE // Without collar
				Dc = Do
			END IF
			A_pri = PI*Dc*L*(1-Th_f*FPM)*nTube*nBank
			A_sec = 2*(HX_l*HX_h - PI*Dc**2*nTube*nBank/4)*FPM*L + 2*HX_h*Th_f*FPM*L
		END IF
		V_t = PI*(Do**2 - Di**2)/4 * L * nTube * nBank 
		V_f = A_sec/2 * Th_f
	END METHOD
END CLASS



CLASS TubeFinMaterialProperties
"Record of tube and fin material properties"
// Each control volume will have its own instance
DECLS
	ENUM Material mat_t = Aluminum
	REAL rho_t UNITS u_kg_m3
	REAL cp_t UNITS u_J_kgK
	REAL k_t UNITS u_W_mK
	REAL M_t UNITS u_kg
	ENUM Material mat_f = Aluminum
	REAL rho_f UNITS u_kg_m3
	REAL cp_f UNITS u_J_kgK
	REAL k_f UNITS u_W_mK
	REAL M_f UNITS u_kg
	REAL M UNITS u_kg
	REAL C UNITS u_J_K
	REAL V UNITS u_m3
METHODS
	METHOD NO_TYPE setConstProps (
		IN ENUM Material mat_t_val,
		IN ENUM Material mat_f_val,
		IN REAL V_t_val,
		IN REAL V_f_val
	)
	BODY
		rho_t = FunConstSolidProp(mat_t,Density,1)
		rho_f = FunConstSolidProp(mat_f,Density,1)
		M_t = V_t_val * rho_t
		M_f = V_f_val * rho_f
		M = M_t + M_f
	END METHOD
	
	METHOD NO_TYPE setTempProps (IN REAL Tw, OUT INTEGER tmp = DUMMY_INTEGER)
	BODY
		k_t = FunVarSolidProp(mat_t,Conductivity,Tw,1,ier)
		cp_t = FunVarSolidProp(mat_t,SpecificHeat,Tw,1,ier)
		k_f = FunVarSolidProp(mat_f,Conductivity,Tw,1,ier)
		cp_f = FunVarSolidProp(mat_f,SpecificHeat,Tw,1,ier)
		C = cp_f*M_f + cp_t*M_t
	END METHOD
END CLASS