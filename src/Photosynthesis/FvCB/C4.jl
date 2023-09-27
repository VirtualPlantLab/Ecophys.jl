### This file contains public API ###
# C4
# C4Q

abstract type C4Type <: FvCB end

# Data structure to store all the C4 parameters without units

"""
    C4(Sco25 = 2590.0, E_Sco = -24.46e3, Kmc25 = 650.0, E_Kmc = 79.43e3, Kmo25 = 450e3, 
       E_Kmo = 36380.0, Vcmax25 = 120.0, E_Vcmax = 65.33, theta = 0.7, Phi2 = 0.83, sigma2 = 0.5, 
       beta = 0.85, fQ = 1.0, fpseudo = 0.1, h = 4.0, Jmax25 = 230.0, E_Jmax = 48e3, D_Jmax = 200e3, 
       S_Jmax = 630.0, x = 0.4, alpha = 0.1, kp25 = 0.7, E_kp = 46.39e3, gbs = 0.003, Rd25 = 1.2, 
       E_Rd = 46.39e3, gso = 0.01, a1 = 0.9, b1 = 0.15e-3)

Data structure to store all the parameters for the C3 photosynthesis model.

# Arguments
- `Sco25`: Sc/o parameter at 25 C
- `E_Sco`: Apparent activation energy of Sc/o (J/mol)
- `Kmc25`: Km for CO2 at 25 C (μmol/mol)
- `E_Kmc3`: Activation energy of Kmc (J/mol)
- `Kmo25`: Km for O2 at 25 C (umol/mol)
- `E_Kmo`: Activation energy of Kmo (J/mol)
- `Vcmax25`: Maximum rate of carboxylation at 25 C (μmol/m2/s)
- `E_Vcmax`: Activation energy of Vcmax (J/mol)
- `theta`: Curvature parameter
- `Phi2`: Low-light PSII quantum yield
- `sigma2`: Partitioning of excitation between PSII and PSI
- `beta`: Leaf absorptance of PAR
- `fQ`: Fraction of electrons at reduced plastoquinone that follow the Q-cycle
- `fpseudo`: Fraction of electrons at PSI that follow cyclic transport around PSI
- `h`: Number of protons required to produce one ATP
- `Jmax25`: Maximum rate of electron transport (μmol/m2/s)
- `E_Jmax`: Activation energy Jmax (J/mol)
- `D_Jmax`: Deactivation energy of Jmax (J/mol)
- `S_Jmax`: Entropy coefficient of Jmax (J/mol/K)
- `x`: Fraction of electron transport partitioned to mesophyll cells
- `alpha`: Fraction of O2 evolution occuring in the bundle sheath
- `kp25`: Initial carboxylation efficiency of the PEP carboxylase (mol/m2/s)
- `E_kp`: Activation energy of kp (J/mol)
- `gbs`: Bundle sheath conductance (mol/m^2/s)
- `Rd25:`: Respiration rate at 25 C (μmol/m2/s)
- `E_Rd`: Activation energy of Rd (J/mol)
- `gso`: Minimum stomatal conductance to fluxes of CO2 in darkness (mol/m2/s)
- `a1`: Empirical parameter in gs formula
- `b1`: Empirical parameter in gs formula (1/kPa)
"""
Base.@kwdef mutable struct C4{T <: Real} <: C4Type
    # Rubisco
    Sco25::T = 2590.0 # Sc/o parameter at 25 C
    E_Sco::T = -24.46e3 # Apparent activation energy of Sc/o (J/mol)
    Kmc25::T = 650.0 # Km for CO2 at 25 C (μmol/mol)
    E_Kmc::T = 79.43e3 # Activation energy of Kmc (J/mol)
    Kmo25::T = 450e3 # Km for O2 at 25 C (umol/mol)
    E_Kmo::T = 36380.0 # Activation energy of Kmo (J/mol)
    Vcmax25::T = 120.0 # Maximum rate of carboxylation at 25 C (μmol/m2/s)
    E_Vcmax::T = 65.33 # Activation energy of Vcmax (J/mol)
    # Electron transport
    theta::T = 0.7 # Curvature parameter
    Phi2::T = 0.83 # Low-light PSII quantum yield
    sigma2::T = 0.5 # Partitioning of excitation between PSII and PSI
    beta::T = 0.85 # Leaf absorptance of PAR
    fQ::T = 1.0 # Fraction of electrons at reduced plastoquinone that follow the Q-cycle
    fpseudo::T = 0.1 # Fraction of electrons at PSI that follow cyclic transport around PSI
    h::T = 4.0 # Number of protons required to produce one ATP
    Jmax25::T = 230.0 # Maximum rate of electron transport (μmol/m2/s)
    E_Jmax::T = 48e3 # Activation energy Jmax (J/mol)
    D_Jmax::T = 200e3 # Deactivation energy of Jmax (J/mol)
    S_Jmax::T = 630.0 # Entropy coefficient of Jmax (J/mol/K)
    x::T = 0.4 # Fraction of electron transport partitioned to mesophyll cells
    alpha::T = 0.1 # Fraction of O2 evolution occuring in the bundle sheath
    # PEP carboxylation
    kp25::T = 0.7 # Initial carboxylation efficiency of the PEP carboxylase (mol/m2/s)
    E_kp::T = 46.39e3 # Activation energy of kp (J/mol)
    # Bundle sheath conductance
    gbs::T = 0.003 # Bundle sheath conductance (mol/m^2/s)
    # Respiration
    Rd25::T = 1.2 # Respiration rate at 25 C (μmol/m2/s)
    E_Rd::T = 46.39e3 # Activation energy of Rd (J/mol)
    # Stomatal conductance
    gso::T = 0.01 # Minimum stomatal conductance to fluxes of CO2 in darkness (mol/m2/s)
    a1::T = 0.9 # Empirical parameter in gs formula
    b1::T = 0.15e-3 # Empirical parameter in gs formula (1/kPa)
end

"""
    C4(Sco25 = 2590.0, E_Sco = -24.46e3J/mol, Kmc25 = 650.0μmol/mol, E_Kmc = 79.43e3J/mol,
       Kmo25 = 450e3μmol/mol, E_Kmo = 36380.0J/mol, Vcmax25 = 120.0μmol/m^2/s, E_Vcmax = 65.33J/mol,
       theta = 0.7, Phi2 = 0.83, sigma2 = 0.5, beta = 0.85, fQ = 1.0, fpseudo = 0.1, h = 4.0, 
       Jmax25 = 230.0μmol/m^2/s, E_Jmax = 48e3J/mol, D_Jmax = 200e3J/mol, S_Jmax = 630.0J/mol/K, 
       x = 0.4, alpha = 0.1, kp25 = 0.7mol/m^2/s, E_kp = 46.39e3J/mol, gbs = 0.003mol/m^2/s, 
       Rd25 = 1.2μmol/m^2/s, E_Rd = 46.39e3J/mol, gso = 0.01mol/m^2/s, a1 = 0.9, b1 = 0.15e-3/Pa)

Data structure to store all the parameters for the C4 photosynthesis model using
`Quantity` objects from Unitful.jl.

# Arguments
- `Sco25`: Sc/o parameter at 25 C
- `E_Sco`: Apparent activation energy of Sc/o (J/mol)
- `Kmc25`: Km for CO2 at 25 C (μmol/mol)
- `E_Kmc3`: Activation energy of Kmc (J/mol)
- `Kmo25`: Km for O2 at 25 C (umol/mol)
- `E_Kmo`: Activation energy of Kmo (J/mol)
- `Vcmax25`: Maximum rate of carboxylation at 25 C (μmol/m2/s)
- `E_Vcmax`: Activation energy of Vcmax (J/mol)
- `theta`: Curvature parameter
- `Phi2`: Low-light PSII quantum yield
- `sigma2`: Partitioning of excitation between PSII and PSI
- `beta`: Leaf absorptance of PAR
- `fQ`: Fraction of electrons at reduced plastoquinone that follow the Q-cycle
- `fpseudo`: Fraction of electrons at PSI that follow cyclic transport around PSI
- `h`: Number of protons required to produce one ATP
- `Jmax25`: Maximum rate of electron transport (μmol/m2/s)
- `E_Jmax`: Activation energy Jmax (J/mol)
- `D_Jmax`: Deactivation energy of Jmax (J/mol)
- `S_Jmax`: Entropy coefficient of Jmax (J/mol/K)
- `x`: Fraction of electron transport partitioned to mesophyll cells
- `alpha`: Fraction of O2 evolution occuring in the bundle sheath
- `kp25`: Initial carboxylation efficiency of the PEP carboxylase (mol/m2/s)
- `E_kp`: Activation energy of kp (J/mol)
- `gbs`: Bundle sheath conductance (mol/m^2/s)
- `Rd25:`: Respiration rate at 25 C (μmol/m2/s)
- `E_Rd`: Activation energy of Rd (J/mol)
- `gso`: Minimum stomatal conductance to fluxes of CO2 in darkness (mol/m2/s)
- `a1`: Empirical parameter in gs formula
- `b1`: Empirical parameter in gs formula (1/kPa)
"""
Base.@kwdef mutable struct C4Q{T <: Real} <: C4Type
    # Rubisco
    Sco25::T = 2590.0 # Sc/o parameter at 25 C
    E_Sco::Quantity{T, dimension(J / mol)} = -24.46e3J / mol # Apparent activation energy of Sc/o (J/mol)
    Kmc25::Quantity{T, dimension(μmol / mol)} = 650.0μmol / mol # Km for CO2 at 25 C (μmol/mol)
    E_Kmc::Quantity{T, dimension(J / mol)} = 79.43e3J / mol # Activation energy of Kmc (J/mol)
    Kmo25::Quantity{T, dimension(μmol / mol)} = 450e3μmol / mol # Km for O2 at 25 C (umol/mol)
    E_Kmo::Quantity{T, dimension(J / mol)} = 36380.0J / mol # Activation energy of Kmo (J/mol)
    Vcmax25::Quantity{T, dimension(μmol / m^2 / s)} = 120.0μmol / m^2 / s # Maximum rate of carboxylation at 25 C (μmol/m2/s)
    E_Vcmax::Quantity{T, dimension(J / mol)} = 65.33J / mol # Activation energy of Vcmax (J/mol)
    # Electron transport
    theta::T = 0.7 # Curvature parameter
    Phi2::T = 0.83 # Low-light PSII quantum yield
    sigma2::T = 0.5 # Partitioning of excitation between PSII and PSI
    beta::T = 0.85 # Leaf absorptance of PAR
    fQ::T = 1.0 # Fraction of electrons at reduced plastoquinone that follow the Q-cycle
    fpseudo::T = 0.1 # Fraction of electrons at PSI that are used by alternative electron sinks
    h::T = 4.0 # Number of protons required to produce one ATP
    Jmax25::Quantity{T, dimension(μmol / m^2 / s)} = 230.0μmol / m^2 / s # Maximum rate of electron transport (μmol/m2/s)
    E_Jmax::Quantity{T, dimension(J / mol)} = 48e3J / mol # Activation energy Jmax (J/mol)
    D_Jmax::Quantity{T, dimension(J / mol)} = 200e3J / mol # Deactivation energy of Jmax (J/mol)
    S_Jmax::Quantity{T, dimension(J / K / mol)} = 630.0J / K / mol # Optimum temperature of Jmax (J/mol/K)
    x::T = 0.4 # Fraction of electron transport partitioned to mesophyll cells
    alpha::T = 0.1 # Fraction of O2 evolution occuring in the bundle sheath
    # PEP carboxylation
    kp25::Quantity{T, dimension(mol / m^2 / s)} = 0.7mol / m^2 / s # Initial carboxylation efficiency of the PEP carboxylase (mol/m2/s)
    E_kp::Quantity{T, dimension(J / mol)} = 46.39e3J / mol # Activation energy of kp (J/mol)
    # Bundle sheath conductance
    gbs::Quantity{T, dimension(mol / m^2 / s)} = 0.003mol / m^2 / s # Bundle sheath conductance (mol/m^2/s)
    # Respiration
    Rd25::Quantity{T, dimension(μmol / m^2 / s)} = 1.2μmol / m^2 / s # Respiration rate at 25 C (μmol/m2/s)
    E_Rd::Quantity{T, dimension(J / mol)} = 46.39e3J / mol # Activation energy of Rd (J/mol)
    # Stomatal conductance
    gso::Quantity{T, dimension(mol / m^2 / s)} = 0.01mol / m^2 / s # Minimum stomatal conductance to fluxes of CO2 in darkness (mol/m2/s)
    a1::T = 0.9 # Empirical parameter in gs formula
    b1::Quantity{T, dimension(1 / kPa)} = 0.15e-3 / Pa # Empirical parameter in gs formula (1/kPa)
end

function photosynthesis(p::C4;
    PAR = 1000.0,
    RH = 0.75,
    Tleaf = 298.0,
    Ca = 400.0,
    O2 = 210e3,
    gb = 0.5,
    net = true)
    photosynthesis(p, PAR, RH, Tleaf, Ca, O2, gb, net)
end
function photosynthesis(p::C4Q; PAR = 1000.0μmol / m^2 / s, RH = 0.75, Tleaf = 298.0K,
    Ca = 400.0μmol / mol, O2 = 210e3μmol / mol, gb = 0.5mol / m^2 / s, net = true)
    photosynthesis(p, PAR, RH, Tleaf, Ca, O2, gb, net)
end

function photosynthesis(p::C4Type, PAR, RH, Tleaf, Ca, O2, gb, net)

    # Check whether we are using plain numbers of physical quantities
    mode = typeof(Ca)

    # Calculate VPD at leaf temperature (kPa)
    vpd = VPD(Tleaf, RH)
    fvpd = max(1.0 / (1.0 / (p.a1 - p.b1 * vpd) - 1.0), 0.0)

    # Calculate values of parameters at Tleaf
    Rd = arrhenius(p.Rd25, p.E_Rd, Tleaf) # μmol/m2/s
    Sco = arrhenius(p.Sco25, p.E_Sco, Tleaf)
    Vcmax = arrhenius(p.Vcmax25, p.E_Vcmax, Tleaf) # μmol/m2/s
    Kmc = arrhenius(p.Kmc25, p.E_Kmc, Tleaf) # μmol/mol
    Kmo = arrhenius(p.Kmo25, p.E_Kmo, Tleaf) # mmol/mol
    Jmax = peaked(p.Jmax25, p.E_Jmax, p.D_Jmax, p.S_Jmax, Tleaf) # μmol/m2/s
    kp = arrhenius(p.kp25, p.E_kp, Tleaf) # μmol/m2/s

    # Respiration in the mesophyll cells
    Rm = 0.5 * Rd # μmol/m2/s

    # CO2 compensation point at leaf surface
    Y = 0.5 / Sco
    Cs_star = (p.gbs * Y * O2 - Rd * (1 + (Y * p.alpha) / 0.047) + Rm) / (p.gbs + kp)

    # C3 limited by Rubisco - C4 limited by PEP carboxylase
    x1_c1 = Vcmax
    x2_c1 = Kmc / Kmo
    x3_c1 = Kmc
    a_c1 = 1 + kp / p.gbs
    b_c1 = zeroflux(mode)
    Ac1 = solveAC4(fvpd,
        p.gso,
        gb,
        Ca,
        Cs_star,
        Y,
        Rd,
        p.alpha,
        x3_c1,
        x2_c1,
        x1_c1,
        p.gbs,
        Rm,
        O2,
        a_c1,
        b_c1)

    # C3 limited by Rubisco. C4 limited by electron transport
    fcyc = 1 -
           (4 * (1 - p.x) * (1 + p.fQ) + 3 * p.h * p.fpseudo) / (3 * p.h - 4 * (1 - p.x))
    k2ll = p.Phi2 * p.sigma2 * p.beta * (1 - p.fpseudo / (1 - fcyc))
    J = (k2ll * PAR + Jmax -
         sqrt((k2ll * PAR + Jmax)^2 - 4 * p.theta * k2ll * Jmax * PAR)) / (2 * p.theta)
    z = (2 + p.fQ - fcyc) / (p.h * (1 - fcyc))
    VpJ2 = p.x * J * z / 2
    x1_c2 = Vcmax
    x2_c2 = Kmc / Kmo
    x3_c2 = Kmc
    b_c2 = VpJ2
    Ac2 = solveAC4(fvpd,
        p.gso,
        gb,
        Ca,
        Cs_star,
        Y,
        Rd,
        p.alpha,
        x3_c2,
        x2_c2,
        x1_c2,
        p.gbs,
        Rm,
        O2,
        1.0,
        b_c2)

    # C3 limited by electron transport. C4 limited by PEPC
    x1_j1 = (1 - p.x) * J * z / 3
    x2_j1 = 7 * Y / 3
    x3_j1 = zeroconc(mode)
    a_j1 = 1 + kp / p.gbs
    b_j1 = zeroflux(mode)
    Aj1 = solveAC4(fvpd,
        p.gso,
        gb,
        Ca,
        Cs_star,
        Y,
        Rd,
        p.alpha,
        x3_j1,
        x2_j1,
        x1_j1,
        p.gbs,
        Rm,
        O2,
        a_j1,
        b_j1)

    # C3 limited by electron transport. C4 limited by electron transport
    x1_j2 = (1 - p.x) * J * z / 3
    x2_j2 = 7 * Y / 3
    x3_j2 = zeroconc(mode)
    b_j2 = VpJ2
    Aj2 = solveAC4(fvpd,
        p.gso,
        gb,
        Ca,
        Cs_star,
        Y,
        Rd,
        p.alpha,
        x3_j2,
        x2_j2,
        x1_j2,
        p.gbs,
        Rm,
        O2,
        1.0,
        b_j2)

    # Take limiting factor
    Ac = min(Ac1, Ac2)
    Aj = min(Aj1, Aj2)
    An = min(Ac, Aj)

    # Stomatal conductance
    gsc = solvegs(p.gso, An, Ca, Cs_star, Rd, fvpd, gb) # mol/m2/s

    # Choose the right output
    A = net ? An : An + Rd
    return (A = A, gs = gsc)
end

# Analytical solution of cubic equation due to coupling of A, gb, gs and gm
function solveAC4(fvpd, gso, gb, Ca, Cs_star, Y, Rd, alpha, x3, x2, x1, gbs, Rm, O2, a, b)
    d = gso * Ca - gso * Cs_star + fvpd * Rd
    m = fvpd - gso / gb
    g = (b - Rm - Y * O2 * gbs) * x1 * m - ((alpha * Y) / 0.047 + 1) * x1 * d +
        a * gbs * x1 * (Ca * m - d / gb - Ca + Cs_star)
    f = (b - Rm - Y * O2 * gbs) * x1 * d + a * gbs * x1 * Ca * d
    h = -(((alpha * Y) / 0.047 + 1) * x1 * m + (a * gbs * x1 * (m - 1)) / gb)
    i = (b - Rm + gbs * x3 + x2 * gbs * O2) * d + a * gbs * Ca * d
    j = (b - Rm + gbs * x3 + x2 * gbs * O2) * m + ((alpha * x2) / 0.047 - 1) * d +
        a * gbs * (Ca * m - d / gb - Ca + Cs_star)
    l = ((alpha * x2) / 0.047 - 1) * m - (a * gbs * (m - 1)) / gb
    p = (j - h + l * Rd) / l
    q = (i + j * Rd - g) / l
    r = -(f - i * Rd) / l
    Q = (p^2 - 3 * q) / 9
    U = (2 * p^3 - 9 * p * q + 27 * r) / 54
    Psi = acos(U / sqrt(Q^3))
    A = -2 * sqrt(Q) * cos((Psi + 4 * pi) / 3) - p / 3
end
