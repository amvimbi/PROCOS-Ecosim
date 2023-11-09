/*-----------------------------------------------------------------------------------------
 LIBRARY: CO2
 FILE: HEATPUMPS
 AUTHOR: Viren Bhanot
 COMPANY: CERN
 CREATION DATE: 11/08/2017
 DESCRIPTION: Components related to vapour compression systems
-----------------------------------------------------------------------------------------*/
USE CRYOLIB
USE MATH
USE THERMAL
USE PORTS_LIB



COMPONENT PlateFinHeatExchanger (
	ENUM InitialConditions init = Ph "Specify inital conditions",
	INTEGER nTube = 5 "Number of tubes per bank",
	INTEGER nBank = 4 "Number of banks of tubes",
	BOOLEAN staggered = TRUE "False if inline tubes",
	BOOLEAN UseAirHTCCorrelation = FALSE "If TRUE uses the correlation of Wang et al. (1999)",
	ENUM VoidFractionModel vfModel = Smith "Void fraction model for better charge prediction"
)
"Round tube plate fin heat exchanger"
PORTS
	IN fluid f_in
	OUT fluid f_out
	IN air a_in
	OUT air a_out
DATA
	ENUM Material mat_t = Copper "Tube material"
	REAL Do = 0.0095 UNITS u_m "Tube outer diameter"
	REAL Th_t = 0.0008 UNITS u_m "Tube thickness"
	REAL L = 2.16*22/5 UNITS u_m "Tube length"
	REAL Pt = 0.025 UNITS u_m "Tube vertical spacing"
	REAL Pl = 0.025 UNITS u_m "Tube horizontal spacing"

	ENUM PlateFinType ftype = WithCollar "BareTubes means no fins"
	ENUM Material mat_f = Aluminum "Fin material"
	REAL FPI = 21 "Fins uper inch"
	REAL Th_f = 0.0001 UNITS u_m "Fin thickness"

	REAL alpha_l = 750 UNITS u_W_m2K "Liquid heat transfer coefficient"
	REAL alpha_tp = 1000 UNITS u_W_m2K "Two-phase heat transfer coefficient"
	REAL alpha_g = 500 UNITS u_W_m2K "Vapour heat transfer coefficient"
	REAL m0 = 0.01 UNITS u_kg_s "Nominal (design) mass flow rate"
	REAL dP0 = 0.02 UNITS u_bar "Nominal (design) pressure drop"
	REAL m_steady = 0.05 UNITS u_kg_s "Steady state mass flow rate - used for HTC correction"

	REAL alpha_a0 = 50 UNITS u_W_m2K "Constant air heat transfer coefficient"
	REAL m0_air = 2 UNITS u_kg_s "Nominal (design) air mass flow rate"

	REAL P0_i = 17 UNITS u_bar "Initial pressure - first segment"
	REAL h0_i = 300000 UNITS u_J_kg "Initial enthalpy - first segment"
	REAL T0_i = 300 UNITS u_K "Initial temperature - first segment"
	REAL x0_i = 0.8 "Initial vapour quality - first segment"
	REAL P0_o = 17 UNITS u_bar "Initial pressure - last segment"
	REAL h0_o = 300000 UNITS u_J_kg "Initial enthalpy - last segment"
	REAL T0_o = 300 UNITS u_K "Initial temperature - last segment"
	REAL x0_o = 0.2 "Initial vapour quality - last segment"
	REAL Tw0_i = 300 UNITS u_K "Initial wall temperature - first segment"
	REAL Tw0_o = 300 UNITS u_K "Initial wall temperature - last segment"
DECLS
	DISCR REAL Vol UNITS u_dm3 "Overall internal volume in Liters"
	REAL P[nBank,nTube] UNITS u_bar
	REAL P_Pa[nBank,nTube] UNITS u_Pa
	REAL h[nBank,nTube] UNITS u_J_kg
	REAL h_corr[nBank,nTube] UNITS u_J_kg
	REAL h_flow[nBank,nTube] UNITS u_J_kg
	REAL T[nBank,nTube] UNITS u_K
	REAL Tsat[nBank,nTube] UNITS u_K
	REAL delT[nBank,nTube] UNITS u_K
	REAL u[nBank,nTube] UNITS u_J_kg
	REAL rho[nBank,nTube] UNITS u_kg_m3
	REAL x[nBank,nTube]
	ENUM Phase phase[nBank,nTube]
	REAL drho_dP[nBank,nTube]
	REAL drho_dh[nBank,nTube]
	REAL m[nBank,nTube+1] UNITS u_kg_s
	REAL mh[nBank,nTube+1] UNITS u_W
	REAL M UNITS u_kg
	REAL Tw[nBank,nTube] UNITS u_K
	REAL cp[nBank,nTube] UNITS u_J_kgK
	REAL k[nBank,nTube] UNITS u_W_mK
	REAL gamma[nBank,nTube]
	REAL mu[nBank,nTube] UNITS u_Pas
	REAL Pr[nBank,nTube]
	REAL sigma[nBank,nTube] UNITS u_N_m
	REAL vsound[nBank,nTube] UNITS u_m_s
	REAL alpha[nBank,nTube] UNITS u_W_m2K
	REAL alpha_a UNITS u_W_m2K
	REAL Qr[nBank,nTube] UNITS u_W
	REAL Qa[nBank,nTube] UNITS u_W
	REAL Qr_total UNITS u_W
	REAL Qa_total UNITS u_W
	REAL Ta[nBank+1,nTube] UNITS u_K
	REAL m_avg[nBank,nTube] UNITS u_kg_s

	CONST ENUM ThFluids medium = Air
	REAL beta_air "Coefficient of thermal expansion"
	REAL cp_air UNITS u_J_kgK "Air specific heat"
	REAL k_air UNITS u_W_mK "Thermal conductivity"
	REAL mu_air UNITS u_Pas "Air viscosity"
	REAL rho_air UNITS u_kg_m3 "Air density"
	REAL Pr_air "Prandtl number"
	REAL Gc_air
	REAL Re_air
	REAL u_max UNITS u_m_s "Max air velocity"
	INTEGER iex,ipx,ipy
	INTEGER tmp[nBank,nTube] // temporary variables
	
OBJECTS
	RefGeometryRecord geo[nBank,nTube]
	SaturationProperties sat[nBank,nTube]
	TubeFinGeometry wall // single instance for entire heat exchanger geometry
	TubeFinMaterialProperties matProp[nBank,nTube]
	AirFlowPassageRecord air
TOPOLOGY
	PATH f_in TO f_out
INIT
	wall.setProps(nTube,nBank,staggered,Do,Th_t,L,Pt,Pl,ftype,FPI,Th_f)
	air.setProps(nTube,nBank,staggered,Do,L,Pt,Pl,ftype,FPI,Th_f,205,alpha_a0) // aluminium thermal conductivity hard-coded
	assignInit2D(nBank,nTube,init,f_in.fluid,P0_i,P0_o,h0_i,h0_o,T0_i,T0_o,x0_i,x0_o,P,h)
	FOR (i IN 1,nBank)
		FOR (j IN 1,nTube)
			geo[i,j].EndSurfaces = 0 // all tubes open at both ends
			geo[i,j].setCylindricalGeometry(Do-2*Th_t,L)
			Tw[i,j] = Tw0_i + ((i-1)*nTube + j - 1) / (nBank*nTube) * (Tw0_o - Tw0_i)
			matProp[i,j].setConstProps(mat_t,mat_f,wall.V_t/(nBank*nTube),wall.V_f/(nBank*nTube))
			P_Pa[i,j] = P[i,j]*1E5
		END FOR
	END FOR
	f_out.is_C = TRUE
	Vol = SUM(i IN 1,nBank; SUM(j IN 1,nTube; geo[i,j].V)) * 1000
CONTINUOUS
	Qr_total = SUM(i IN 1,nBank; SUM(j IN 1,nTube; Qr[i,j]))
	Qa_total = SUM(i IN 1,nBank; SUM(j IN 1,nTube; Qa[i,j]))
	
	Fluid_prop_vs_pT(medium,101325,a_in.T_dry,cp_air,mu_air,rho_air,k_air,ier)
	beta_air = Fluid_beta_vs_T(medium,a_in.T_dry,ier)
	u_max = (a_in.m)/(rho_air*air.A_min)
	Gc_air = rho_air*u_max
	Re_air = Gc_air * air.Dc / mu_air
	Pr_air = cp_air * mu_air / k_air
	alpha_a = ZONE (UseAirHTCCorrelation) max(5,ChangWang(Re_air,Pr_air,Gc_air,cp_air,air))
		OTHERS (a_in.m/m0_air)**0.8 * alpha_a0
	
	EXPAND_BLOCK (i IN 1,nBank)
		EXPAND_BLOCK (j IN 1,nTube)
			h_corr[i,j] = getCorrectionEnthalpy(vfModel,h[i,j],rho[i,j],sat[i,j])
			h_flow[i,j] = h[i,j] + h_corr[i,j]
		END EXPAND_BLOCK
	END EXPAND_BLOCK
	
	m[1,1] = f_in.m
	m[nBank,nTube+1] = f_out.m
	mh[1,1] = semiLinear(m[1,1],f_in.hf,h_flow[1,1])
	mh[nBank,nTube+1] = semiLinear(f_out.m,h_flow[nBank,nTube],f_out.hb)
	EXPAND_BLOCK (i IN 1,nBank)
		EXPAND_BLOCK (j IN 2,nTube)
			m[i,j] = m0*regRoot2((P[i,j-1]-P[i,j])/(dP0/(nTube*nBank)))
			mh[i,j] = semiLinear(m[i,j],h_flow[i,j-1],h_flow[i,j])
		END EXPAND_BLOCK
	END EXPAND_BLOCK
	
	EXPAND_BLOCK (i IN 1,nBank-1)
		m[i,nTube+1] = m0*regRoot2((P[i,nTube]-P[i+1,1])/(dP0/(nTube*nBank)))
		m[i+1,1] = m[i,nTube+1]
		mh[i,nTube+1] = semiLinear(m[i,nTube+1],h_flow[i,nTube],h_flow[i+1,1])
		mh[i+1,1] = mh[i,nTube+1]
	END EXPAND_BLOCK
	
	EXPAND (i IN 1,nTube) Ta[1,i] = a_in.T_dry // first bank has uniform temperatures
	
	EXPAND_BLOCK (i IN 1,nBank)
		EXPAND_BLOCK (j IN 1,nTube)
			m_avg[i,j] = 0.5 * (abs(m[i,j]) + abs(m[i,j+1])) // used for HTC correction
			P[i,j] = P_Pa[i,j]*1E-5
			getStatePh(f_in.fluid,P[i,j],h[i,j],cp[i,j],drho_dP[i,j],drho_dh[i,j],
				gamma[i,j],k[i,j],mu[i,j],phase[i,j],Pr[i,j],rho[i,j],sigma[i,j],
				T[i,j],Tsat[i,j],u[i,j],vsound[i,j],x[i,j],sat[i,j],ier,ipx,ipy)
			matProp[i,j].setTempProps(Tw[i,j],tmp[i,j])
			alpha[i,j] = max(10,getSmoothHTC(alpha_g,alpha_l,alpha_tp,m_avg[i,j],m_steady,x[i,j]))
			Qr[i,j] = alpha[i,j] * geo[i,j].As * (T[i,j] - Tw[i,j])
			
			// alternative could be to use div_safe for the division by a_in.m, in which case, it would just start remembering the previous
			// value in case of zero mass flow rate. But what if a_in.m is zero at the _first_ time step itself?
			Ta[i+1,j] = ZONE(a_in.m != 0) Ta[i,j] + (Tw[i,j] - Ta[i,j]) * (1 - exp(-alpha_a * (air.Ao/(nTube*nBank)) / (a_in.m/nTube * cp_air))) OTHERS Ta[i,j]
			Qa[i,j] = a_in.m/nTube * cp_air * (Ta[i+1,j] - Ta[i,j])
			
			geo[i,j].V*(drho_dP[i,j]*P_Pa[i,j]'+drho_dh[i,j]*h[i,j]') = m[i,j] - m[i,j+1]
			geo[i,j].V*((h[i,j]*drho_dP[i,j]-1)*P_Pa[i,j]' + (h[i,j]*drho_dh[i,j]+rho[i,j])*h[i,j]') = mh[i,j]-mh[i,j+1] - Qr[i,j]
			
			matProp[i,j].C * Tw[i,j]' = Qr[i,j] - Qa[i,j]
			
			delT[i,j] = T[i,j] - Tsat[i,j]
			
		END EXPAND_BLOCK
	END EXPAND_BLOCK
	
	M = SUM(i IN 1,nBank; SUM(j IN 1,nTube; geo[i,j].V*rho[i,j]))
	
	f_in.fluid = f_out.fluid
	f_in.n_fluid = f_out.n_fluid		
	f_in.Gcrit = CRYO_Gcrit_fun(f_in.fluid,P[1,1]*1E5,rho[1,1],T[1,1],0,f_in.v,vsound[1,1],ier,ipy)
	<hb> f_in.hb = h[1,1]
	f_in.rho = rho[1,1]
	f_in.P = P[1,1]
	f_in.A = geo[1,1].Ap
	f_in.P_aux = P[1,1] + geo[1,1].V*rho[1,1]*vsound[1,1]*1E-5/geo[1,1].As
	f_in.I = 0.5*geo[1,1].V/geo[1,1].As**2
	f_in.v = f_in.m/(rho[1,1]*geo[1,1].Ap)

	f_out.Gcrit = CRYO_Gcrit_fun(f_in.fluid,P[nBank,nTube]*1E5,rho[nBank,nTube],T[nBank,nTube],0,f_out.v,vsound[nBank,nTube],ier,ipy)
	<hf> f_out.hf = h[nBank,nTube]
	f_out.rho = rho[nBank,nTube]
	f_out.P = P[nBank,nTube]
	f_out.A = geo[nBank,nTube].Ap
	f_out.P_aux = P[nBank,nTube] + geo[nBank,nTube].V*rho[nBank,nTube]*vsound[nBank,nTube]*1E-5/geo[nBank,nTube].As
	f_out.I = 0.5*geo[nBank,nTube].V/geo[nBank,nTube].As**2
	f_out.v = f_out.m/(rho[nBank,nTube]*geo[nBank,nTube].Ap)
	
	a_out.m = a_in.m
	a_out.T_dry = SUM(i IN 1,nTube; Ta[nBank+1,i])/nTube // average value
	a_out.w = 0
	a_out.P = a_in.P
END COMPONENT



COMPONENT ACoil (
	ENUM InitialConditions init = Ph "Specify initial conditions",
	INTEGER nCoil = 2 "Number of parallel coils",
	INTEGER nTube = 2 "Number of tubes per bank",
	INTEGER nBank = 4 "Number of banks of tubes",
	BOOLEAN staggered = TRUE "False if inline tubes",
	ENUM VoidFractionModel choice = Smith "Slip-ratio based model to implement",
	BOOLEAN UseAirHTC = FALSE "If TRUE, uses ChangWang correlation"
)
"Round tube plate fin heat exchanger with multiple coils in parallel"
PORTS
	IN fluid f_in
	OUT fluid f_out
	IN air a_in
	OUT air a_out
DATA
	ENUM Material mat_t = Copper "Tube material"
	REAL Do = 0.0095 UNITS u_m "Tube outer diameter"
	REAL Th_t = 0.0008 UNITS u_m "Tube thickness"
	REAL L = 0.503*22/2 UNITS u_m "Tube length"
	REAL Pt = 0.025 UNITS u_m "Tube vertical spacing"
	REAL Pl = 0.025 UNITS u_m "Tube horizontal spacing"

	ENUM PlateFinType ftype = WithCollar "BareTubes means no fins"
	ENUM Material mat_f = Aluminum "Fin material"
	REAL FPI = 21 "Fins per inch"
	REAL Th_f = 0.0001 UNITS u_m "Fin thickness"

	REAL alpha_l = 750 UNITS u_W_m2K "Liquid heat transfer coefficient"
	REAL alpha_tp = 1000 UNITS u_W_m2K "Two-phase heat transfer coefficient"
	REAL alpha_g = 500 UNITS u_W_m2K "Vapour heat transfer coefficient"
	REAL m0 = 0.01 UNITS u_kg_s "Nominal (design) mass flow rate"
	REAL dP0 = 0.02 UNITS u_bar "Nominal (design) pressure drop"
	REAL m_steady = 0.05 UNITS u_kg_s "Steady state mass flow rate - used for HTC correction"

	REAL alpha_a0 = 50 UNITS u_W_m2K "Constant air heat transfer coefficient"
	REAL m0_air = 1.2 UNITS u_kg_s "Nominal (design) air mass flow rate"

	REAL P0_i = 17 UNITS u_bar "Initial pressure - first segment"
	REAL h0_i = 300000 UNITS u_J_kg "Initial enthalpy - first segment"
	REAL T0_i = 300 UNITS u_K "Initial temperature - first segment"
	REAL x0_i = 0.8 "Initial vapour quality - first segment"
	REAL P0_o = 17 UNITS u_bar "Initial pressure - last segment"
	REAL h0_o = 300000 UNITS u_J_kg "Initial enthalpy - last segment"
	REAL T0_o = 300 UNITS u_K "Initial temperature - last segment"
	REAL x0_o = 0.2 "Initial vapour quality - last segment"
	REAL Tw0_i = 300 UNITS u_K "Initial wall temperature - first segment"
	REAL Tw0_o = 300 UNITS u_K "Initial wall temperature - last segment"
DECLS
	DISCR REAL Vol UNITS u_dm3
	REAL P[nBank,nTube] UNITS u_bar
	REAL P_Pa[nBank,nTube] UNITS u_Pa
	REAL h[nBank,nTube] UNITS u_J_kg
	REAL h_corr[nBank,nTube] UNITS u_J_kg
	REAL h_flow[nBank,nTube] UNITS u_J_kg
	REAL T[nBank,nTube] UNITS u_K
	REAL Tsat[nBank,nTube] UNITS u_K
	REAL delT[nBank,nTube] UNITS u_K
	REAL u[nBank,nTube] UNITS u_J_kg
	REAL rho[nBank,nTube] UNITS u_kg_m3
	REAL x[nBank,nTube]
	ENUM Phase phase[nBank,nTube]
	REAL drho_dP[nBank,nTube]
	REAL drho_dh[nBank,nTube]
	REAL m[nBank,nTube+1] UNITS u_kg_s
	REAL mh[nBank,nTube+1] UNITS u_W
	REAL M UNITS u_kg
	REAL Tw[nBank,nTube] UNITS u_K
	REAL cp[nBank,nTube] UNITS u_J_kgK
	REAL k[nBank,nTube] UNITS u_W_mK
	REAL gamma[nBank,nTube]
	REAL mu[nBank,nTube] UNITS u_Pas
	REAL Pr[nBank,nTube]
	REAL sigma[nBank,nTube] UNITS u_N_m
	REAL vsound[nBank,nTube] UNITS u_m_s
	REAL alpha[nBank,nTube] UNITS u_W_m2K
	REAL alpha_a UNITS u_W_m2K
	REAL Qr[nBank,nTube] UNITS u_W
	REAL Qa[nBank,nTube] UNITS u_W
	REAL Qr_total UNITS u_W
	REAL Qa_total UNITS u_W
	REAL Ta[nBank+1,nTube] UNITS u_K
	REAL m_avg[nBank,nTube] UNITS u_kg_s

	CONST ENUM ThFluids medium = Air
	REAL beta_air "Coefficient of thermal expansion"
	REAL cp_air UNITS u_J_kgK "Air specific heat"
	REAL k_air UNITS u_W_mK "Thermal conductivity"
	REAL mu_air UNITS u_Pas "Air viscosity"
	REAL rho_air UNITS u_kg_m3 "Air density"
	REAL Pr_air "Prandtl number"
	REAL Gc_air
	REAL Re_air
	REAL u_max UNITS u_m_s "Max air velocity"
	INTEGER iex,ipx,ipy
	INTEGER tmp[nBank,nTube] // temporary variables
OBJECTS
	RefGeometryRecord geo[nBank,nTube]
	SaturationProperties sat[nBank,nTube]
	TubeFinGeometry wall
	TubeFinMaterialProperties matProp[nBank,nTube]
	AirFlowPassageRecord air
TOPOLOGY
	PATH f_in TO f_out
INIT
	wall.setProps(nTube,nBank,staggered,Do,Th_t,L,Pt,Pl,ftype,FPI,Th_f)
	air.setProps(nTube,nBank,staggered,Do,L,Pt,Pl,ftype,FPI,Th_f,205,alpha_a0) // aluminium thermal conductivity hard coded
	assignInit2D(nBank,nTube,init,f_in.fluid,P0_i,P0_o,h0_i,h0_o,T0_i,T0_o,x0_i,x0_o,P,h)
	FOR (i IN 1,nBank)
		FOR (j IN 1,nTube)
			geo[i,j].EndSurfaces = 0
			geo[i,j].setCylindricalGeometry(Do-2*Th_t,L)
			Tw[i,j] = Tw0_i + ((i-1)*nTube + j - 1) / (nBank*nTube) * (Tw0_o - Tw0_i)
			matProp[i,j].setConstProps(mat_t,mat_f,wall.V_t/(nBank*nTube),wall.V_f/(nBank*nTube))
			P_Pa[i,j] = P[i,j]*1E5
		END FOR
	END FOR
	f_out.is_C = FALSE
	Vol = nCoil * SUM(i IN 1,nBank; SUM(j IN 1,nTube; geo[i,j].V)) * 1000 //L
CONTINUOUS
	Qr_total = nCoil * SUM(i IN 1,nBank; SUM(j IN 1,nTube; Qr[i,j]))
	Qa_total = nCoil * SUM(i IN 1,nBank; SUM(j IN 1,nTube; Qa[i,j]))
	
	Fluid_prop_vs_pT(medium,101325,a_in.T_dry,cp_air,mu_air,rho_air,k_air,ier)
	beta_air = Fluid_beta_vs_T(medium,a_in.T_dry,ier)
	u_max = (a_in.m/nCoil)/(rho_air*air.A_min)
	Gc_air = rho_air*u_max
	Re_air = Gc_air * air.Dc / mu_air
	Pr_air = cp_air * mu_air / k_air
	IF (UseAirHTC) INSERT
		alpha_a = max(5,ChangWang(Re_air,Pr_air,Gc_air,cp_air,air))
	ELSE
		alpha_a = (a_in.m/m0_air)**2 * alpha_a0
	END IF
		
	EXPAND_BLOCK (i IN 1,nBank)
		EXPAND_BLOCK (j IN 1,nTube)
			h_corr[i,j] = getCorrectionEnthalpy(choice,h[i,j],rho[i,j],sat[i,j])
			h_flow[i,j] = h[i,j] + h_corr[i,j]
		END EXPAND_BLOCK
	END EXPAND_BLOCK
	
	m[1,1] = f_in.m/nCoil
	m[nBank,nTube+1] = f_out.m/nCoil
	mh[1,1] = semiLinear(m[1,1],f_in.hf,h_flow[1,1])
	mh[nBank,nTube+1] = semiLinear(m[nBank,nTube+1],h_flow[nBank,nTube],f_out.hb)
	EXPAND_BLOCK (i IN 1,nBank)
		EXPAND_BLOCK (j IN 2,nTube)
			m[i,j] = m0/nCoil * regRoot2((P[i,j-1]-P[i,j])/(dP0/nBank/nTube))
			mh[i,j] = semiLinear(m[i,j],h_flow[i,j-1],h_flow[i,j])
		END EXPAND_BLOCK
	END EXPAND_BLOCK
	EXPAND_BLOCK (i IN 1,nBank-1)
		m[i,nTube+1] = m0/nCoil * regRoot2((P[i,nTube]-P[i+1,1])/(dP0/nBank/nTube))
		m[i+1,1] = m[i,nTube+1]
		mh[i,nTube+1] = semiLinear(m[i,nTube+1],h_flow[i,nTube],h_flow[i+1,1])
		mh[i+1,1] = mh[i,nTube+1]
	END EXPAND_BLOCK
	
	EXPAND (i IN 1,nTube) Ta[1,i] = a_in.T_dry // first bank has uniform temperatures
	
	EXPAND_BLOCK (i IN 1,nBank)
		EXPAND_BLOCK (j IN 1,nTube)
			m_avg[i,j] = 0.5 * (abs(m[i,j]) + abs(m[i,j+1])) // used for HTC correction
			P[i,j] = P_Pa[i,j]*1E-5
			getStatePh(f_in.fluid,P[i,j],h[i,j],cp[i,j],drho_dP[i,j],drho_dh[i,j],
				gamma[i,j],k[i,j],mu[i,j],phase[i,j],Pr[i,j],rho[i,j],sigma[i,j],
				T[i,j],Tsat[i,j],u[i,j],vsound[i,j],x[i,j],sat[i,j],ier,ipx,ipy)
			matProp[i,j].setTempProps(Tw[i,j],tmp[i,j])
			alpha[i,j] = max(10,getSmoothHTC(alpha_g,alpha_l,alpha_tp,m_avg[i,j],m_steady/nCoil,x[i,j]))
			Qr[i,j] = alpha[i,j] * geo[i,j].As * (T[i,j] - Tw[i,j])
			 
			Ta[i+1,j] = Ta[i,j] + (Tw[i,j] - Ta[i,j]) * (1 - exp(-alpha_a * (air.Ao/(nTube*nBank)) / (a_in.m/(nCoil*nTube) * cp_air)))
			Qa[i,j] = a_in.m/(nCoil*nTube) * cp_air * (Ta[i+1,j] - Ta[i,j])

			geo[i,j].V*(drho_dP[i,j]*P_Pa[i,j]'+drho_dh[i,j]*h[i,j]') = m[i,j] - m[i,j+1]
			geo[i,j].V*((h[i,j]*drho_dP[i,j]-1)*P_Pa[i,j]' + (h[i,j]*drho_dh[i,j]+rho[i,j])*h[i,j]') = mh[i,j]-mh[i,j+1] - Qr[i,j]
			
			matProp[i,j].C * Tw[i,j]' = Qr[i,j] - Qa[i,j]
			delT[i,j] = T[i,j] - Tsat[i,j]
		END EXPAND_BLOCK
	END EXPAND_BLOCK
	
	M = nCoil * SUM(i IN 1,nBank; SUM(j IN 1,nTube; geo[i,j].V*rho[i,j]))
	
	f_in.fluid = f_out.fluid
	f_in.n_fluid = f_out.n_fluid		
	f_in.Gcrit = CRYO_Gcrit_fun(f_in.fluid,P[1,1]*1E5,rho[1,1],T[1,1],0,f_in.v,vsound[1,1],ier,ipy)
<hb> 	f_in.hb = h_flow[1,1]
	f_in.rho = rho[1,1]
	f_in.P = P[1,1]
	f_in.A = geo[1,1].Ap
	f_in.P_aux = P[1,1] + geo[1,1].V*rho[1,1]*vsound[1,1]*1E-5/geo[1,1].As
	f_in.I = 0.5*geo[1,1].V/geo[1,1].As**2
	f_in.v = f_in.m/(rho[1,1]*geo[1,1].Ap)

	f_out.Gcrit = CRYO_Gcrit_fun(f_in.fluid,P[nBank,nTube]*1E5,rho[nBank,nTube],T[nBank,nTube],0,f_out.v,vsound[nBank,nTube],ier,ipy)
<hf> 	f_out.hf = h_flow[nBank,nTube]
	f_out.rho = rho[nBank,nTube]
	f_out.P = P[nBank,nTube]
	f_out.A = geo[nBank,nTube].Ap
	f_out.P_aux = P[nBank,nTube] + geo[nBank,nTube].V*rho[nBank,nTube]*vsound[nBank,nTube]*1E-5/geo[nBank,nTube].As
	f_out.I = 0.5*geo[nBank,nTube].V/geo[nBank,nTube].As**2
	f_out.v = f_out.m/(rho[nBank,nTube]*geo[nBank,nTube].Ap)
	
	a_out.m = a_in.m
	a_out.T_dry = SUM(i IN 1,nTube; Ta[nBank+1,i])/nTube // average value
	a_out.w = 0
	a_out.P = a_in.P
END COMPONENT