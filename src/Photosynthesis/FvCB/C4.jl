
abstract type C4 <: FvCB end

# Data structure to store all the C4 parameters without units
mutable struct C4F{T <: Real} <: C4
    # Rubisco
    Sco25::T # Sc/o parameter at 25 C
    E_Sco::T # Apparent activation energy of Sc/o (J/mol)
    Kmc25::T # Km for CO2 at 25 C (μmol/mol)
    E_Kmc::T # Activation energy of Kmc (J/mol)
    Kmo25::T # Km for O2 at 25 C (umol/mol)
    E_Kmo::T # Activation energy of Kmo (J/mol)
    Vcmax25::T # Maximum rate of carboxylation at 25 C (μmol/m2/s)
    E_Vcmax::T # Activation energy of Vcmax (J/mol)
    # Electron transport
    theta::T # Curvature parameter
    Phi2::T # Low-light PSII quantum yield
    sigma2::T # Partitioning of excitation between PSII and PSI
    beta::T # Leaf absorptance of PAR
    fQ::T # Fraction of electrons at reduced plastoquinone that follow the Q-cycle
    fpseudo::T # Fraction of electrons at PSI that follow cyclic transport around PSI
    h::T # Number of protons required to produce one ATP
    Jmax25::T # Maximum rate of electron transport (μmol/m2/s)
    E_Jmax::T # Activation energy Jmax (J/mol)
    D_Jmax::T # Deactivation energy of Jmax (J/mol)
    Topt_Jmax::T # Entropy coefficient of Jmax (J/mol/K)
    x::T # Fraction of electron transport partitioned to mesophyll cells
    alpha::T # Fraction of O2 evolution occuring in the bundle sheath
    # PEP carboxylation
    kp25::T # Initial carboxylation efficiency of the PEP carboxylase (mol/m2/s)
    E_kp::T # Activation energy of kp (J/mol)
    # Bundle sheath conductance
    gbs::T # Bundle sheath conductance (mol/m^2/s)
    # Respiration
    Rd25::T # Respiration rate at 25 C (μmol/m2/s)
    E_Rd::T # Activation energy of Rd (J/mol)
    # Stomatal conductance
    gso::T # Minimum stomatal conductance to fluxes of CO2 in darkness (mol/m2/s)
    a1::T # Empirical parameter in gs formula
    b1::T # Empirical parameter in gs formula (1/kPa)
end

function C4F(; Sco25 = 2590.0, E_Sco = -24.46e3, Kmc25 = 650.0, E_Kmc = 79.43e3, Kmo25 = 450e3, 
    E_Kmo = 36380.0, Vcmax25 = 120.0, E_Vcmax = 65.33, theta = 0.7, Phi2 = 0.83, sigma2 = 0.5, 
    beta = 0.85, fQ = 1.0, fpseudo = 0.1, h = 4.0, Jmax25 = 230.0, E_Jmax = 48e3, D_Jmax = 200e3, 
    Topt_Jmax = 300.5, x = 0.4, alpha = 0.1, kp25 = 0.7, E_kp = 46.39e3, gbs = 0.003, Rd25 = 1.2, 
    E_Rd = 46.39e3, gso = 0.01, a1 = 0.9, b1 = 0.15)

    C4F(Sco25, E_Sco, Kmc25, E_Kmc, Kmo25, E_Kmo, Vcmax25, E_Vcmax, theta, Phi2, sigma2, beta, 
        fQ, fpseudo, h, Jmax25, E_Jmax, D_Jmax, Topt_Jmax, x, alpha, kp25, E_kp, gbs, Rd25, E_Rd, gso, a1, b1)
end


# Data structure to store all the parameters with units
mutable struct C4Q{T <: Real} <: C4
    # Rubisco
    Sco25::T # Sc/o parameter at 25 C
    E_Sco::Quantity{T, dimension(J/mol)} # Apparent activation energy of Sc/o (J/mol)
    Kmc25::Quantity{T, dimension(μmol/mol)} # Km for CO2 at 25 C (μmol/mol)
    E_Kmc::Quantity{T, dimension(J/mol)} # Activation energy of Kmc (J/mol)
    Kmo25::Quantity{T, dimension(μmol/mol)} # Km for O2 at 25 C (umol/mol)
    E_Kmo::Quantity{T, dimension(J/mol)} # Activation energy of Kmo (J/mol)
    Vcmax25::Quantity{T, dimension(μmol/m^2/s)} # Maximum rate of carboxylation at 25 C (μmol/m2/s)
    E_Vcmax::Quantity{T, dimension(J/mol)} # Activation energy of Vcmax (J/mol)
    # Electron transport
    theta::T # Curvature parameter
    Phi2::T # Low-light PSII quantum yield
    sigma2::T # Partitioning of excitation between PSII and PSI
    beta::T # Leaf absorptance of PAR
    fQ::T # Fraction of electrons at reduced plastoquinone that follow the Q-cycle
    fpseudo::T # Fraction of electrons at PSI that follow cyclic transport around PSI
    h::T # Number of protons required to produce one ATP
    Jmax25::Quantity{T, dimension(μmol/m^2/s)} # Maximum rate of electron transport (μmol/m2/s)
    E_Jmax::Quantity{T, dimension(J/mol)} # Activation energy Jmax (J/mol)
    D_Jmax::Quantity{T, dimension(J/mol)} # Deactivation energy of Jmax (J/mol)
    Topt_Jmax::Quantity{T, dimension(K)} # Optimum temperature of Jmax (J/mol/K)
    x::T # Fraction of electron transport partitioned to mesophyll cells
    alpha::T # Fraction of O2 evolution occuring in the bundle sheath
    # PEP carboxylation
    kp25::Quantity{T, dimension(mol/m^2/s)} # Initial carboxylation efficiency of the PEP carboxylase (mol/m2/s)
    E_kp::Quantity{T, dimension(J/mol)} # Activation energy of kp (J/mol)
    # Bundle sheath conductance
    gbs::Quantity{T, dimension(mol/m^2/s)} # Bundle sheath conductance (mol/m^2/s)
    # Respiration
    Rd25::Quantity{T, dimension(μmol/m^2/s)} # Respiration rate at 25 C (μmol/m2/s)
    E_Rd::Quantity{T, dimension(J/mol)} # Activation energy of Rd (J/mol)
    # Stomatal conductance
    gso::Quantity{T, dimension(mol/m^2/s)} # Minimum stomatal conductance to fluxes of CO2 in darkness (mol/m2/s)
    a1::T # Empirical parameter in gs formula
    b1::Quantity{T, dimension(1/kPa)} # Empirical parameter in gs formula (1/kPa)
end

function C4Q(; Sco25 = 2590.0, E_Sco = -24.46e3J/mol, Kmc25 = 650.0μmol/mol, E_Kmc = 79.43e3J/mol,
            Kmo25 = 450e3μmol/mol, E_Kmo = 36380.0J/mol, Vcmax25 = 120.0μmol/m^2/s, E_Vcmax = 65.33J/mol,
            theta = 0.7, Phi2 = 0.83, sigma2 = 0.5, beta = 0.85, fQ = 1.0, fpseudo = 0.1, h = 4.0, 
            Jmax25 = 230.0μmol/m^2/s, E_Jmax = 48e3J/mol, D_Jmax = 200e3J/mol, Topt_Jmax = 300.5K, 
            x = 0.4, alpha = 0.1, kp25 = 0.7mol/m^2/s, E_kp = 46.39e3J/mol, gbs = 0.003mol/m^2/s, 
            Rd25 = 1.2μmol/m^2/s, E_Rd = 46.39e3J/mol, gso = 0.01mol/m^2/s, a1 = 0.9, b1 = 0.15/kPa)

    C4Q(Sco25, E_Sco, Kmc25, E_Kmc, Kmo25, E_Kmo, Vcmax25, E_Vcmax,
    theta, Phi2, sigma2, beta, fQ, fpseudo, h, Jmax25, E_Jmax, D_Jmax, Topt_Jmax, 
    x, alpha, kp25, E_kp, gbs, Rd25, E_Rd, gso, a1, b1)
end


"""
A_gs(par::C4, PAR, RH, Tleaf, Ca, gb)

Calculate net CO2 assimilation (umol/m2/s), transpiration (mmol/m2/s) and stomatal
condutance to fluxes of CO2 (mol/m2/s) as a function of photosynthetically active
radiation (PAR, umol/m2/s), relative humidity (RH, %), leaf temperature (Tleaf,
K), air CO2 partial pressure (Ca, μmol/mol), oxygen (O2, μmol/mol) and boundary layer 
conductance to CO2 (gb, mol/m2.s). Environmental inputs must be scalar. 
"""
function A_gs(p::C4, PAR, RH, Tleaf, Ca, O2, gb)

    # Check whether we are using plain numbers of physical quantities
    mode = typeof(Ca)

    # Calculate VPD at leaf temperature (kPa)
    vpd = VPD(Tleaf, RH)
    fvpd = max(1.0/(1.0/(p.a1 - p.b1*vpd) - 1.0),0.0)

    # Calculate values of parameters at Tleaf
    Rd    = arrhenius(p.Rd25, p.E_Rd, Tleaf) # μmol/m2/s
    Sco   = arrhenius(p.Sco25, p.E_Sco, Tleaf)
    Vcmax = arrhenius(p.Vcmax25, p.E_Vcmax, Tleaf) # μmol/m2/s
    Kmc   = arrhenius(p.Kmc25, p.E_Kmc, Tleaf) # μmol/mol
    Kmo   = arrhenius(p.Kmo25, p.E_Kmo, Tleaf) # mmol/mol
    Jmax  = peaked(p.Jmax25, p.E_Jmax, p.D_Jmax, p.Topt_Jmax, Tleaf) # μmol/m2/s
    kp    = arrhenius(p.kp25, p.E_kp, Tleaf) # μmol/m2/s

    # Respiration in the mesophyll cells
    Rm = 0.5*Rd # μmol/m2/s

    # CO2 compensation point at leaf surface
    Y = 0.5/Sco
    Cs_star = (p.gbs*Y*O2 - Rd*(1 + (Y*p.alpha)/0.047) + Rm)/(p.gbs + kp)

    # C3 limited by Rubisco - C4 limited by PEP carboxylase
    x1_c1 = Vcmax
    x2_c1 = Kmc/Kmo
    x3_c1 = Kmc
    a_c1 = 1 + kp/p.gbs
    b_c1 = zeroflux(mode)
    Ac1 = solveAC4(fvpd, p.gso, gb, Ca, Cs_star, Y, Rd, p.alpha, x3_c1, x2_c1, x1_c1, p.gbs, Rm, O2, a_c1, b_c1)

    # C3 limited by Rubisco. C4 limited by electron transport
    fcyc = 1 - (4*(1 - p.x)*(1 + p.fQ) + 3*p.h*p.fpseudo)/(3*p.h - 4*(1 - p.x))
    k2ll = p.Phi2*p.sigma2*p.beta*(1 - p.fpseudo/(1 - fcyc))
    J = (k2ll*PAR + Jmax - sqrt((k2ll*PAR + Jmax)^2 - 4*p.theta*k2ll*Jmax*PAR))/(2*p.theta)
    z = (2 + p.fQ - fcyc)/(p.h*(1 - fcyc))
    VpJ2 = p.x*J*z/2
    x1_c2 = Vcmax
    x2_c2 = Kmc/Kmo
    x3_c2 = Kmc
    b_c2  = VpJ2
    Ac2 = solveAC4(fvpd, p.gso, gb, Ca, Cs_star, Y, Rd, p.alpha, x3_c2, x2_c2, x1_c2, p.gbs, Rm, O2, 1.0, b_c2)

    # C3 limited by electron transport. C4 limited by PEPC
    x1_j1 = (1 - p.x)*J*z/3
    x2_j1 = 7*Y/3
    x3_j1 = zeroconc(mode)
    a_j1  = 1 + kp/p.gbs
    b_j1  = zeroflux(mode)
    Aj1   = solveAC4(fvpd, p.gso, gb, Ca, Cs_star, Y, Rd, p.alpha, x3_j1, x2_j1, x1_j1, p.gbs, Rm, O2, a_j1, b_j1)

    # C3 limited by electron transport. C4 limited by electron transport
    x1_j2 = (1 - p.x)*J*z/3
    x2_j2 = 7*Y/3
    x3_j2 = zeroconc(mode)
    b_j2  = VpJ2
    Aj2   = solveAC4(fvpd, p.gso, gb, Ca, Cs_star, Y, Rd, p.alpha, x3_j2, x2_j2, x1_j2, p.gbs, Rm, O2, 1.0, b_j2)

    # Take limiting factor
    Ac = min(Ac1, Ac2)
    Aj = min(Aj1, Aj2)
    A  = min(Ac, Aj)

    # Stomatal conductance
    gs = solvegs(p.gso,  A,  Ca, Cs_star,  Rd,  fvpd, gb) # mol/m2/s

    return (A = A, gs = gs) 
end


# Analytical solution of cubic equation due to coupling of A, gb, gs and gm
function solveAC4(fvpd, gso, gb, Ca, Cs_star, Y, Rd, alpha, x3, x2, x1, gbs, Rm, O2, a, b)
    d = gso*Ca - gso*Cs_star + fvpd*Rd
    m = fvpd - gso/gb
    g = (b - Rm - Y*O2*gbs)*x1*m - ((alpha*Y)/0.047 + 1)*x1*d + a*gbs*x1*(Ca*m - d/gb - Ca + Cs_star)
    f = (b - Rm - Y*O2*gbs)*x1*d + a*gbs*x1*Ca*d
    h = -(((alpha*Y)/0.047 + 1)*x1*m + (a*gbs*x1*(m - 1))/gb)
    i = (b - Rm + gbs*x3 + x2*gbs*O2)*d + a*gbs*Ca*d
    j = (b - Rm + gbs*x3 + x2*gbs*O2)*m + ((alpha*x2)/0.047 - 1)*d + a*gbs*(Ca*m - d/gb - Ca + Cs_star)
    l = ((alpha*x2)/0.047 - 1)*m - (a*gbs*(m - 1))/gb
    p = (j - h + l*Rd)/l
    q = (i + j*Rd - g)/l
    r = -(f - i*Rd)/l
    Q = (p^2 - 3*q)/9
    U = (2*p^3 - 9*p*q + 27*r)/54
    Psi = acos(U/sqrt(Q^3))
    A = -2*sqrt(Q)*cos((Psi + 4*pi)/3) - p/3
end
  