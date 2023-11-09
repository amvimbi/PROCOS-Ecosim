/*-----------------------------------------------------------------------------------------
 LIBRARY: CO2
 FILE: ADDED_FUNCTIONALITIES
 CREATION DATE: 07/06/2023
-----------------------------------------------------------------------------------------*/
USE MATH 
USE CRYOLIB
USE THERMAL
USE THERMO_TABLE_INTERP



FUNCTION NO_TYPE characterize_partially_submerged_circular_surface (
	IN REAL zl 	UNITS u_m	"liquid height i.e. liquid level starting from the port bottom",
	IN REAL D 	UNITS u_m	"Diameter of circular cross section",
	OUT REAL Al 	UNITS u_m2	"Area submerged by liquid",
	OUT REAL Av 	UNITS u_m2	"Area in contact with vapor",
	OUT REAL gamma 	"Void surface fraction")
    DECLS
	REAL R		UNITS u_m	"Radius"
	REAL theta 			"Angle in radian between (i) the radius semiline (angle = Pi in trogonometric circle) and (ii) the semiline that crosses the intersection point of the circle and the line defining the liquid level " 
	REAL zv		UNITS u_m	"vapor height"
	REAL Atot 	UNITS u_m2	"Total area of circular port"
    BODY
	R = D/2
	Atot = MATH.PI * (D/2)**2 
	IF (zl <= R) THEN
		theta = asin((R-zl)/R)
		Al = (0.5*MATH.PI-theta)*R**2 - R*cos(theta)*(R-zl)
		Av = Atot - Al
		gamma = Av/Atot
	ELSE 
		zv = D-zl
		theta = asin((R-zv)/R)
		Av = (0.5*MATH.PI-theta)*R**2 - R*cos(theta)*(R-zv)
		Al = Atot - Av
		gamma = Av/Atot
	END IF
END FUNCTION 



FUNCTION REAL CRYO_PF_CritProp_CORR (
        IN ENUM ChemName chem      "Working chemical constituent",
        IN ENUM  FluidProps fprop  "Wanted chem property",
        OUT INTEGER ier)
    DECLS
        REAL vprop            	"Property value (?)"
        REAL prop_crit[20]    	"Fluid properties at critical point"
        INTEGER ig         	"Fluid number"
    BODY 
        IF(ReadPro[chem] == 0) THEN
	        ig = setofPos(ChemName, chem)
		ReadPro[chem] = ig
	ELSE
		ig = ReadPro[chem]
	END IF
	IF(itab[ig] == 0) THEN
		IF(path_prop =="") THEN
	           	path_prop = expandFilePath(Pathprops)
	         	FOR (i IN ChemName)   -- loop to transfer defined CEA labels to FORTRAN
	            		read_label(CEA_GAS_Name[i]) 
	         	END FOR
	           	-- reading of CEA coefficients of all declared ChemName
	           	read_cea(setofSize(ChemName), path_prop, MW_ch, ier)
	        END IF
		WRITE("\n   ******    Reading properties file %s.dat     ******\n",ChemFileName[chem])
           	table_read(" ",ig, path_prop, ChemFileName[chem], ig, ier, fml, UsrFluiddH[chem])
           	ng_tot   = ng_tot + 1
           	itab[ig] = ng_tot
        END IF
        ASSERT (ier < 100) FATAL "CRYO_PF_CritProp : problems reading thermo properties file"
	-------------------------------
	--- Missing code
	IF (Prop_crit[ig,1]<1E-4) THEN
		crit_cond(ig, prop_crit, ier)
      	  	FOR (i IN 1,20)   -- loop to transfer properties to global array
			IF i == 1 THEN
	        		Prop_crit[ig,i] = prop_crit[i]/1e5 -- Pressure conversion from Pa to Bar
			ELSE
				Prop_crit[ig,i] = prop_crit[i]
			END IF
	   	END FOR
	ELSE 
	END IF
	-------------------------------
        IF (fprop == fprop_pressure) THEN
        	vprop = Prop_crit[ig,1]  
        ELSEIF (fprop == fprop_temperature) THEN
        	vprop = Prop_crit[ig,2]
        ELSEIF (fprop == fprop_density) THEN 
           	vprop = Prop_crit[ig,3]
        ELSE 
           	ASSERT (FALSE) FATAL "Property not supported by CRYO_PF_CritProp"
        END IF
        RETURN vprop
END FUNCTION



COMPONENT FP_CritProp_compute 

	DATA
		ENUM ChemName chemical =  Carbon_Dioxide
		
	DECLS
		REAL rho		UNITS u_kg_m3		"Density"
		REAL Txx 		UNITS u_K			"Temperature"
		REAL Pxx		UNITS u_bar			"Pressure"
		REAL h		UNITS u_J_kg		"Enthalpy"
		REAL u		UNITS u_J_kg		"Specific internal energy"
		REAL x		UNITS no_units		"Quality"

		INTEGER ier1, ier2, ier3, ier4, ier5, ier6, ier7, ier8, ier9, ier10, ier11, ier12, ier13, ier14, ier15, ier16, ier17, ier18
	
	CONTINUOUS 
		//rho = CRYO_PF_CritProp_CORR(chemical, fprop_density, ier1)
		//T = CRYO_PF_CritProp_CORR(chemical, fprop_temperature, ier2)
		Pxx = CRYO_PF_CritProp_CORR(chemical, fprop_pressure, ier3)
		
		//rho = CRYO_PF_prop_vs_ph(chemical, 40, 400000, fprop_density, ier1, ier2, ier3)
		Txx = CRYO_PF_prop_vs_ph(chemical, 20, 250000, fprop_temperature, ier4, ier5, ier6)//-273.15
		//P = CRYO_PF_prop_vs_ph(chemical, 40, 400000, fprop_pressure, ier7, ier8, ier9)
		//h = CRYO_PF_prop_vs_ph(chemical, 40, 400000, fprop_enthalpy, ier10, ier11, ier12)
		//u = CRYO_PF_prop_vs_ph(chemical, 40, 400000, fprop_energy, ier13, ier14, ier15)
		//x = CRYO_PF_prop_vs_ph(chemical, 40, 400000, fprop_quality, ier16, ier17, ier18)
END COMPONENT