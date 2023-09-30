var documenterSearchIndex = {"docs":
[{"location":"Photosynthesis/#Module-Photosynthesis","page":"Photosynthesis","title":"Module Photosynthesis","text":"","category":"section"},{"location":"Photosynthesis/","page":"Photosynthesis","title":"Photosynthesis","text":"CurrentModule = Ecophys.Photosynthesis","category":"page"},{"location":"Photosynthesis/#CO2-assimilation-and-stomatal-conductance","page":"Photosynthesis","title":"CO2 assimilation and stomatal conductance","text":"","category":"section"},{"location":"Photosynthesis/","page":"Photosynthesis","title":"Photosynthesis","text":"C3\nC3Q\nC4\nC4Q\nphotosynthesis","category":"page"},{"location":"Photosynthesis/#Ecophys.Photosynthesis.C3","page":"Photosynthesis","title":"Ecophys.Photosynthesis.C3","text":"C3(Sco25 = 2800.0, E_Sco = -24.46e3, Kmc25 = 270.0, E_Kmc = 80.99e3, \n    Kmo25 = 165.0e3, E_Kmo = 23.72e3, Vcmax25 = 120.0, E_Vcmax = 65.33e3, \n    simpleJ = false, k2ll = 0.35, theta = 0.7, Phi2 = 0.82, sigma2 = 0.5, \n    beta = 0.85, fcyc = 0.1, \n    fpseudo = 0.05, Jmax25 = 230.0, E_Jmax = 30.0e3, D_Jmax = 200.0e3, \n    S_Jmax = 650.0, TPU25 = 12.0, E_TPU = 53.1e3, D_TPU = 20.18e3,\n    S_TPU = 650.0, Rd25 = 1.2, E_Rd = 46.39e3, gm25 = 0.4, E_gm = 49.6e3, \n    D_gm = 437.4e3, S_gm = 1400.0, gso = 0.01, a1 = 0.85, b1 = 0.14e-3)\n\nData structure to store all the parameters for the C3 photosynthesis model.\n\nArguments\n\nSco25: Sc/o parameter at 25 C\nE_Sco: Apparent activation energy of Sc/o (J/mol)\nKmc25: Km for CO2 at 25 C (μmol/mol)\nE_Kmc: Activation energy of Kmc (J/mol)\nKmo25: Km for O2 at 25 C (umol/mol)\nE_Kmo: Activation energy of Kmo (J/mol)\nVcmax25: Maximum rate of carboxylation at 25 C (μmol/m2/s)\nE_Vcmax: Activation energy of Vcmax (J/mol)\ntheta: Curvature parameter\nsimpleJ: Use k2ll rather than calculating from other parameters\nk2ll: Low-light use efficiency for electron transport\nPhi2: Low-light PSII quantum yield\nsigma2: Partitioning of excitation between PSII and PSI\nbeta: Leaf absorptance of PAR\nfcyc:  Fraction of electrons at PSI that follow cyclic transport around PSI\nfpseudo: Fraction of electrons at PSI that are used by alternative electron sinks\nJmax25: Maximum rate of electron transport (μmol/m2/s)\nE_Jmax: Activation energy Jmax (J/mol)\nD_Jmax: Deactivation energy of Jmax (J/mol)\nS_Jmax: Entropy term for Jmax (K)\nTPU25: Maximum rate of triose phosphate utilisation (μmol/m2/s)\nE_TPU: Activation energy TPU (J/mol)\nD_TPU: Deactivation energy of TPU (J/mol)\nS_TPU: Entropy term for TPU (K)\nRd25: Respiration rate at 25 C (μmol/m2/s)\nE_Rd: Activation energy of Rd (J/mol)\ngm25: Maximum rate of CO2 assimilation at 25 C (μmol/m2/s)\nE_gm: Activation energy of gm (J/mol)\nD_gm: Deactivation energy of gm (J/mol)\nS_gm: Entropy term for gm (K)\ngso: Minimum stomatal conductance to fluxes of CO2 in darkness (mol/m2/s/Pa)\na1: Empirical parameter in gs formula\nb1: Empirical parameter in gs formula\n\n\n\n\n\n","category":"type"},{"location":"Photosynthesis/#Ecophys.Photosynthesis.C3Q","page":"Photosynthesis","title":"Ecophys.Photosynthesis.C3Q","text":"C3Q(Sco25 = 2800.0, E_Sco = -24.46e3J/mol, Kmc25 = 270.0μmol/mol, E_Kmc = 80.99e3J/mol,\n     Kmo25 = 165.0e3μmol/mol, E_Kmo = 23.72e3J/mol, Vcmax25 = 120.0μmol/m^2/s, E_Vcmax = 65.33e3J/mol,\n     simpleJ = false, k2ll = 0.35, theta = 0.7, Phi2 = 0.82, sigma2 = 0.5, beta = 0.85, fcyc = 0.1, fpseudo = 0.05, \n     Jmax25 = 230.0μmol/m^2/s, E_Jmax = 30.0e3J/mol, D_Jmax = 200.0e3J/mol, S_Jmax = 650.0J/mol/K, \n     TPU25 = 12.0μmol/m^2/s, E_TPU = 53.1e3J/mol, D_TPU = 201.8e3J/mol, S_TPU = 650.0K, \n     Rd25 = 1.2μmol/m^2/s, E_Rd = 46.39e3J/mol, gm25 = 0.4mol/m^2/s, E_gm = 49.6e3J/mol, \n     D_gm = 437.4e3J/mol, S_gm = 1400.0K, gso = 0.01mol/m^2/s, a1 = 0.85, b1 = 0.14e-3/Pa)\n\nData structure to store all the parameters for the C3 photosynthesis model using Quantity objects from Unitful.jl.\n\nArguments\n\nSco25: Sc/o parameter at 25 C\nE_Sco: Apparent activation energy of Sc/o (J/mol)\nKmc25: Km for CO2 at 25 C (μmol/mol)\nE_Kmc: Activation energy of Kmc (J/mol)\nKmo25: Km for O2 at 25 C (umol/mol)\nE_Kmo: Activation energy of Kmo (J/mol)\nVcmax25: Maximum rate of carboxylation at 25 C (μmol/m2/s)\nE_Vcmax: Activation energy of Vcmax (J/mol)\ntheta: Curvature parameter\nsimpleJ: Use k2ll rather than calculating from other parameters\nk2ll: Low-light use efficiency for electron transport\nPhi2: Low-light PSII quantum yield\nsigma2: Partitioning of excitation between PSII and PSI\nbeta: Leaf absorptance of PAR\nfcyc:  Fraction of electrons at PSI that follow cyclic transport around PSI\nfpseudo: Fraction of electrons at PSI that are used by alternative electron sinks\nJmax25: Maximum rate of electron transport (μmol/m2/s)\nE_Jmax: Activation energy Jmax (J/mol)\nD_Jmax: Deactivation energy of Jmax (J/mol)\nS_Jmax: Entropy term for Jmax (J/K/mol)\nTPU25: Maximum rate of triose phosphate utilisation (μmol/m2/s)\nE_TPU: Activation energy TPU (J/mol)\nD_TPU: Deactivation energy of TPU (J/mol)\nS_TPU: Entropy term for TPU (J/K/mol)\nRd25: Respiration rate at 25 C (μmol/m2/s)\nE_Rd: Activation energy of Rd (J/mol)\ngm25: Mesophyll conductance at 25 C (mol/m2/s)\nE_gm: Activation energy of gm (J/mol)\nD_gm: Deactivation energy of gm (J/mol)\nS_gm: Entropy term for gm (J/K/mol)\ngso: Minimum stomatal conductance to fluxes of CO2 in darkness (mol/m2/s)\na1: Empirical parameter in gs formula\nb1: Empirical parameter in gs formula (1/kPa)\n\n\n\n\n\n","category":"type"},{"location":"Photosynthesis/#Ecophys.Photosynthesis.C4","page":"Photosynthesis","title":"Ecophys.Photosynthesis.C4","text":"C4(Sco25 = 2590.0, E_Sco = -24.46e3, Kmc25 = 650.0, E_Kmc = 79.43e3, Kmo25 = 450e3, \n   E_Kmo = 36380.0, Vcmax25 = 120.0, E_Vcmax = 65.33, theta = 0.7, Phi2 = 0.83, sigma2 = 0.5, \n   beta = 0.85, fQ = 1.0, fpseudo = 0.1, h = 4.0, Jmax25 = 230.0, E_Jmax = 48e3, D_Jmax = 200e3, \n   S_Jmax = 630.0, x = 0.4, alpha = 0.1, kp25 = 0.7, E_kp = 46.39e3, gbs = 0.003, Rd25 = 1.2, \n   E_Rd = 46.39e3, gso = 0.01, a1 = 0.9, b1 = 0.15e-3)\n\nData structure to store all the parameters for the C3 photosynthesis model.\n\nArguments\n\nSco25: Sc/o parameter at 25 C\nE_Sco: Apparent activation energy of Sc/o (J/mol)\nKmc25: Km for CO2 at 25 C (μmol/mol)\nE_Kmc3: Activation energy of Kmc (J/mol)\nKmo25: Km for O2 at 25 C (umol/mol)\nE_Kmo: Activation energy of Kmo (J/mol)\nVcmax25: Maximum rate of carboxylation at 25 C (μmol/m2/s)\nE_Vcmax: Activation energy of Vcmax (J/mol)\ntheta: Curvature parameter\nPhi2: Low-light PSII quantum yield\nsigma2: Partitioning of excitation between PSII and PSI\nbeta: Leaf absorptance of PAR\nfQ: Fraction of electrons at reduced plastoquinone that follow the Q-cycle\nfpseudo: Fraction of electrons at PSI that follow cyclic transport around PSI\nh: Number of protons required to produce one ATP\nJmax25: Maximum rate of electron transport (μmol/m2/s)\nE_Jmax: Activation energy Jmax (J/mol)\nD_Jmax: Deactivation energy of Jmax (J/mol)\nS_Jmax: Entropy coefficient of Jmax (J/mol/K)\nx: Fraction of electron transport partitioned to mesophyll cells\nalpha: Fraction of O2 evolution occuring in the bundle sheath\nkp25: Initial carboxylation efficiency of the PEP carboxylase (mol/m2/s)\nE_kp: Activation energy of kp (J/mol)\ngbs: Bundle sheath conductance (mol/m^2/s)\nRd25:: Respiration rate at 25 C (μmol/m2/s)\nE_Rd: Activation energy of Rd (J/mol)\ngso: Minimum stomatal conductance to fluxes of CO2 in darkness (mol/m2/s)\na1: Empirical parameter in gs formula\nb1: Empirical parameter in gs formula (1/kPa)\n\n\n\n\n\n","category":"type"},{"location":"Photosynthesis/#Ecophys.Photosynthesis.C4Q","page":"Photosynthesis","title":"Ecophys.Photosynthesis.C4Q","text":"C4(Sco25 = 2590.0, E_Sco = -24.46e3J/mol, Kmc25 = 650.0μmol/mol, E_Kmc = 79.43e3J/mol,\n   Kmo25 = 450e3μmol/mol, E_Kmo = 36380.0J/mol, Vcmax25 = 120.0μmol/m^2/s, E_Vcmax = 65.33J/mol,\n   theta = 0.7, Phi2 = 0.83, sigma2 = 0.5, beta = 0.85, fQ = 1.0, fpseudo = 0.1, h = 4.0, \n   Jmax25 = 230.0μmol/m^2/s, E_Jmax = 48e3J/mol, D_Jmax = 200e3J/mol, S_Jmax = 630.0J/mol/K, \n   x = 0.4, alpha = 0.1, kp25 = 0.7mol/m^2/s, E_kp = 46.39e3J/mol, gbs = 0.003mol/m^2/s, \n   Rd25 = 1.2μmol/m^2/s, E_Rd = 46.39e3J/mol, gso = 0.01mol/m^2/s, a1 = 0.9, b1 = 0.15e-3/Pa)\n\nData structure to store all the parameters for the C4 photosynthesis model using Quantity objects from Unitful.jl.\n\nArguments\n\nSco25: Sc/o parameter at 25 C\nE_Sco: Apparent activation energy of Sc/o (J/mol)\nKmc25: Km for CO2 at 25 C (μmol/mol)\nE_Kmc3: Activation energy of Kmc (J/mol)\nKmo25: Km for O2 at 25 C (umol/mol)\nE_Kmo: Activation energy of Kmo (J/mol)\nVcmax25: Maximum rate of carboxylation at 25 C (μmol/m2/s)\nE_Vcmax: Activation energy of Vcmax (J/mol)\ntheta: Curvature parameter\nPhi2: Low-light PSII quantum yield\nsigma2: Partitioning of excitation between PSII and PSI\nbeta: Leaf absorptance of PAR\nfQ: Fraction of electrons at reduced plastoquinone that follow the Q-cycle\nfpseudo: Fraction of electrons at PSI that follow cyclic transport around PSI\nh: Number of protons required to produce one ATP\nJmax25: Maximum rate of electron transport (μmol/m2/s)\nE_Jmax: Activation energy Jmax (J/mol)\nD_Jmax: Deactivation energy of Jmax (J/mol)\nS_Jmax: Entropy coefficient of Jmax (J/mol/K)\nx: Fraction of electron transport partitioned to mesophyll cells\nalpha: Fraction of O2 evolution occuring in the bundle sheath\nkp25: Initial carboxylation efficiency of the PEP carboxylase (mol/m2/s)\nE_kp: Activation energy of kp (J/mol)\ngbs: Bundle sheath conductance (mol/m^2/s)\nRd25:: Respiration rate at 25 C (μmol/m2/s)\nE_Rd: Activation energy of Rd (J/mol)\ngso: Minimum stomatal conductance to fluxes of CO2 in darkness (mol/m2/s)\na1: Empirical parameter in gs formula\nb1: Empirical parameter in gs formula (1/kPa)\n\n\n\n\n\n","category":"type"},{"location":"Photosynthesis/#Ecophys.Photosynthesis.photosynthesis","page":"Photosynthesis","title":"Ecophys.Photosynthesis.photosynthesis","text":"photosynthesis(par::C3, PAR = 1000.0, RH = 0.75, Tleaf = 298.0, Ca = 400.0, O2 = 210e3, gb = 0.5, net = true)\nphotosynthesis(par::C4, PAR = 1000.0, RH = 0.75, Tleaf = 298.0, Ca = 400.0, O2 = 210e3, gb = 0.5, net = true)\nphotosynthesis(par::C3Q, PAR = 1000.0μmol/m^2/s, RH = 0.75, Tleaf = 298.0K, Ca = 400.0μmol/mol, O2 = 210e3μmol/mol, gb = 0.5mol/m^2/s, net = true)\nphotosynthesis(par::C4Q, PAR = 1000.0μmol/m^2/s, RH = 0.75, Tleaf = 298.0K, Ca = 400.0μmol/mol, O2 = 210e3μmol/mol, gb = 0.5mol/m^2/s, net = true)\n\nCalculate net or gross CO2 assimilation (umol/m2/s) and stomatal condutance to fluxes of CO2 (mol/m2/s) as a function of  photosynthetically active radiation (PAR, umol/m2/s), relative humidity (RH), leaf temperature (Tleaf, K), air CO2 partial pressure (Ca, μmol/mol), oxygen (O2, μmol/mol) and boundary layer  conductance to CO2 (gb, mol/m2/s). Environmental inputs must be scalar. The argument net indicates whether the net or gross CO2 assimilation should be returned.\n\n\n\n\n\n","category":"function"},{"location":"Photosynthesis/#Boundary-layer-conductance","page":"Photosynthesis","title":"Boundary layer conductance","text":"","category":"section"},{"location":"Photosynthesis/","page":"Photosynthesis","title":"Photosynthesis","text":"gb\nsimplegb\nsimplegbQ\ngbAngle\ngbAngleQ","category":"page"},{"location":"Photosynthesis/#Ecophys.Photosynthesis.gb","page":"Photosynthesis","title":"Ecophys.Photosynthesis.gb","text":"gb(p::gbType, ws, Tleaf, Tair, P)\n\nCompute boundary layer conductance for heat, water vapor and CO2.\n\nArguments\n\np: Model of boundary layer conductance\nws: Wind speed (m/s)\nTleaf: Leaf temperature (K)\nTair: Air temperature (K)\nP: Air pressure (Pa)\n\nReturns\n\ngbh: Boundary layer conductance for heat (W/m²/K)\ngbw: Boundary layer conductance for water vapor (mol/m²/s)\ngbc: Boundary layer conductance for CO2 (mol/m²/s)\n\n\n\n\n\n","category":"function"},{"location":"Photosynthesis/#Ecophys.Photosynthesis.simplegb","page":"Photosynthesis","title":"Ecophys.Photosynthesis.simplegb","text":"simplegb(; d = 0.01)\n\nSimple model of boundary layer conductance.\n\nArguments\n\nd: Characteristic leaf length (m)\n\n\n\n\n\n","category":"type"},{"location":"Photosynthesis/#Ecophys.Photosynthesis.simplegbQ","page":"Photosynthesis","title":"Ecophys.Photosynthesis.simplegbQ","text":"simplegbQ(; d = 0.01m)\n\nSimple model of boundary layer conductance using Quantity from Unitful.jl.\n\nArguments\n\nd: Characteristic leaf length (m)\n\n\n\n\n\n","category":"type"},{"location":"Photosynthesis/#Ecophys.Photosynthesis.gbAngle","page":"Photosynthesis","title":"Ecophys.Photosynthesis.gbAngle","text":"gbAngle(; d = 0.01, ang = 0.0, ar = 1.0, fangm = 1.381, fangk = 0.034, \nα = 2.738, b0_0 = 0.455, d_b0 = 2.625, b0_n = 0.373, b0_KAR = 28.125,  db_0 = 0.085, \nd_db = 0.437, db_n = 5.175, db_KAR = 0.884,  β_0 = 3.362, d_β = 17.664, β_n = 4.727, \nβ_KAR = 0.677)\n\nModel of boundary layer conductance that accounts for inclination angle and leaf aspect ratio (see documentation for details).\n\nArguments\n\nd: Characteristic leaf length (m)\nang: Leaf inclination angle (°)\nar: Leaf aspect ratio (length/width)\nfangm: Maximum enhancement factor due to inclination angle\nfangk: Exponent in response to inclination angle\nα: Effect on back boundary layer conductance due to leaf inclination angle and aspect ratio\nb0_0: Parameter in the effect of aspect ratio (see documentation)\nd_b0: Parameter in the effect of aspect ratio (see documentation)\nb0_n: Parameter in the effect of aspect ratio (see documentation)\nb0_KAR: Parameter in the effect of aspect ratio (see documentation)\ndb_0: Parameter in the effect of aspect ratio (see documentation)\nd_db: Parameter in the effect of aspect ratio (see documentation)\ndb_n: Parameter in the effect of aspect ratio (see documentation)\ndb_KAR: Parameter in the effect of aspect ratio (see documentation)\nβ_0: Parameter in the effect of aspect ratio (see documentation)\nd_β: Parameter in the effect of aspect ratio (see documentation)\nβ_n: Parameter in the effect of aspect ratio (see documentation)\nβ_KAR: Parameter in the effect of aspect ratio (see documentation)\n\n\n\n\n\n","category":"type"},{"location":"Photosynthesis/#Ecophys.Photosynthesis.gbAngleQ","page":"Photosynthesis","title":"Ecophys.Photosynthesis.gbAngleQ","text":"gbAngleQ(; d = 0.01m, ang = 0.0, ar = 1.0, fangm = 1.381, fangk = 0.034, \nα = 2.738, b0_0 = 0.455, d_b0 = 2.625, b0_n = 0.373, b0_KAR = 28.125,  db_0 = 0.085, \nd_db = 0.437, db_n = 5.175, db_KAR = 0.884,  β_0 = 3.362, d_β = 17.664, β_n = 4.727, \nβ_KAR = 0.677)\n\nModel of boundary layer conductance that accounts for inclination angle and leaf aspect ratio (see documentation for details) using Quantity for Unitful.jl.\n\nArguments\n\nd: Characteristic leaf length (m)\nang: Leaf inclination angle (°)\nar: Leaf aspect ratio (length/width)\nfangm: Maximum enhancement factor due to inclination angle\nfangk: Exponent in response to inclination angle\nα: Effect on back boundary layer conductance due to leaf inclination angle and aspect ratio\nb0_0: Parameter in the effect of aspect ratio (see documentation)\nd_b0: Parameter in the effect of aspect ratio (see documentation)\nb0_n: Parameter in the effect of aspect ratio (see documentation)\nb0_KAR: Parameter in the effect of aspect ratio (see documentation)\ndb_0: Parameter in the effect of aspect ratio (see documentation)\nd_db: Parameter in the effect of aspect ratio (see documentation)\ndb_n: Parameter in the effect of aspect ratio (see documentation)\ndb_KAR: Parameter in the effect of aspect ratio (see documentation)\nβ_0: Parameter in the effect of aspect ratio (see documentation)\nd_β: Parameter in the effect of aspect ratio (see documentation)\nβ_n: Parameter in the effect of aspect ratio (see documentation)\nβ_KAR: Parameter in the effect of aspect ratio (see documentation)\n\n\n\n\n\n","category":"type"},{"location":"Photosynthesis/#Energy-balance","page":"Photosynthesis","title":"Energy balance","text":"","category":"section"},{"location":"Photosynthesis/","page":"Photosynthesis","title":"Photosynthesis","text":"SimpleOptical\nenergybalance\nsolve_energy_balance\ntranspiration","category":"page"},{"location":"Photosynthesis/#Ecophys.Photosynthesis.SimpleOptical","page":"Photosynthesis","title":"Ecophys.Photosynthesis.SimpleOptical","text":"SimpleOptical(; αPAR = 0.85, αNIR = 0.20, ϵ = 0.95)\n\nSimple optical properties of a leaf.\n\nArguments\n\nαPAR: Absorption coefficient of PAR\nαNIR: Absorption coefficient of NIR\nϵ: Emissivity for thermal radiation\n\n\n\n\n\n","category":"type"},{"location":"Photosynthesis/#Ecophys.Photosynthesis.energybalance","page":"Photosynthesis","title":"Ecophys.Photosynthesis.energybalance","text":"energybalance(pgb, pAgs, pEb, PAR, NIR, ws, RH, Tair, Ca, P, O2)\n\nCalculate the energy balance of a leaf.\n\nArguments\n\npgb: Boundary layer conductance model\npAgs: Photosynthesis and stomatal conductance model\npEb: Optical properties of the leaf\nPAR: Photosynthetically active radiation (umol/m2/s)\nNIR: Near-infrared radiation (W/m2)\nws: Wind speed (m/s)\nRH: Relative humidity\nTair: Air temperature (K)\nCa: Atmospheric CO2 concentration (μmol/mol)\nP: Air pressure (kPa)\nO2: Atmospheric O2 concentration (μmol/mol)\n\nDetails\n\nInputs maybe be either Real or Quantity types (i.e., with physical units).  If Quantity types are used, the output will be a Quantity type.\n\n\n\n\n\n","category":"function"},{"location":"Photosynthesis/#Ecophys.Photosynthesis.solve_energy_balance","page":"Photosynthesis","title":"Ecophys.Photosynthesis.solve_energy_balance","text":"solve_energy_balance(Ags::Union{C3Q, C4Q}; gb = simplegbQ(), \n                     opt = SimpleOptical(), PAR = 1000.0μmol/m^2/s, \n                     NIR = 250.0W/m^2, ws = 1.0m/s, RH = 0.75, \n                     Tair = 298.0K, Ca = 400.0μmol/mol, P = 101.0kPa, \n                     O2 = 210.0mmol/mol, order = Order2(), xatol = 0.01, \n                     maxfnevals = 100, net = true)\nsolve_energy_balance(Ags::Union{C3, C4}; gb = simplegb(), \n                     opt = SimpleOptical(), PAR = 1000.0, NIR = 250.0, \n                     ws = 1.0, RH = 0.75, Tair = 298.0, Ca = 400.0, \n                     P = 101.0e3, O2 = 210.0e3, order = Order2(), xatol = 0.01, \n                     maxfnevals = 100, net = true)\n\nSolve the leaf energy balance coupled to photosynthesis and transpiration.\n\nArguments\n\nAgs: Photosynthesis and stomatal conductance model\ngb: Boundary layer conductance model\nopt: Optical properties of the leaf\nPAR: Photosynthetically active radiation (umol/m2/s)\nNIR: Near-infrared radiation (W/m2)\nws: Wind speed (m/s)\nRH: Relative humidity\nTair: Air temperature (K)\nCa: Atmospheric CO2 concentration (μmol/mol)\nP: Air pressure (Pa)\nO2: Atmospheric O2 concentration (μmol/mol)\norder: Order of the root solving algorithm that finds leaf temperature         (see Roots.jl package for more information).\nxatol: Absolute tolerance of the root solving algorithm (see Roots.jl package for more information),\nmaxfnevals: Maximum number of function evaluations of the root solving algorithm (see Roots.jl package for more information).\nnet: Whether to return net or gross CO2 assimilation.\n\nDetails\n\nInputs maybe be either Real or Quantity types from Unitful.jl (i.e., with  physical units). If Quantity types are used, the output will be a Quantity  type.\n\nReturns\n\nA named tuple with net CO2 assimilation (An, μmol/m^2/s), gross CO2  assimilation (Ag, μmol/m^2/s), transpiration (Tr, mol/m^2/s) and leaf  temperature (Tleaf, K).\n\n\n\n\n\n","category":"function"},{"location":"Photosynthesis/#Ecophys.Photosynthesis.transpiration","page":"Photosynthesis","title":"Ecophys.Photosynthesis.transpiration","text":"transpiration(;gsw = 0.1, gbw = 1.0, Tleaf = 300.0, Tair = 298.0, P = 101e3,\n               RH = 0.75)\n\nCompute transpiration rate (mol/m^2/s) from conductance to water vapor and  environmental variables.\n\nArguments\n\ngsw: Stomatal conductance to water vapor (mol/m^2/s)\ngbw: Boundary layer conductance to water vapor (mol/m^2/s)\nTleaf: Leaf temperature (K)\nTair: Air temperature (K)\nP: Air pressure (Pa)\nRH: Relative humidity\n\n\n\n\n\n","category":"function"}]
}