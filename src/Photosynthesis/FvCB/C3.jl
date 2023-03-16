### This file contains public API ###
# C3
# C3Q
# photosynthesis

abstract type FvCB <: Ags end
abstract type C3Type <: FvCB end

"""
    C3(Sco25 = 2800.0, E_Sco = -24.46e3, Kmc25 = 270.0, E_Kmc = 80.99e3, 
        Kmo25 = 165.0e3, E_Kmo = 23.72e3, Vcmax25 = 120.0, E_Vcmax = 65.33e3, 
        theta = 0.7, Phi2 = 0.82, sigma2 = 0.5, beta = 0.85, fcyc = 0.1, 
        fpseudo = 0.05, Jmax25 = 230.0, E_Jmax = 30.0e3, D_Jmax = 200.0e3, 
        Topt_Jmax = 300.5, TPU25 = 12.0, E_TPU = 53.1e3, D_TPU = 20.18e3)

Data structure to store all the parameters for the C3 photosynthesis model.

# Arguments
- `Sco25`: Sc/o parameter at 25 C
- `E_Sco`: Apparent activation energy of Sc/o (J/mol)
- `Kmc25`: Km for CO2 at 25 C (μmol/mol)
- `E_Kmc`: Activation energy of Kmc (J/mol)
- `Kmo25`: Km for O2 at 25 C (umol/mol)
- `E_Kmo`: Activation energy of Kmo (J/mol)
- `Vcmax25`: Maximum rate of carboxylation at 25 C (μmol/m2/s)
- `E_Vcmax`: Activation energy of Vcmax (J/mol)
- `theta`: Curvature parameter
- `Phi2`: Low-light PSII quantum yield
- `sigma2`: Partitioning of excitation between PSII and PSI
- `beta`: Leaf absorptance of PAR
- `fcyc`:  Fraction of electrons at PSI that follow cyclic transport around PSI
- `fpseudo`: Fraction of electrons at PSI that are used by alternative electron sinks
- `Jmax25`: Maximum rate of electron transport (μmol/m2/s)
- `E_Jmax`: Activation energy Jmax (J/mol)
- `D_Jmax`: Deactivation energy of Jmax (J/mol)
- `Topt_Jmax`: Optimal temperature for Jmax (K)
- `TPU25`: Maximum rate of triose phosphate utilisation (μmol/m2/s)
- `E_TPU`: Activation energy TPU (J/mol)
- `D_TPU`: Deactivation energy of TPU (J/mol)
"""
Base.@kwdef mutable struct C3{T <: Real} <: C3Type
    # Rubisco CO2/O2 specificity
    Sco25::T = 2800.0 # Sc/o parameter at 25 C
    E_Sco::T = -24.46e3 # Apparent activation energy of Sc/o (J/mol)
    # Michaelis-Menten constants carboxylation (J/mol/K)
    Kmc25::T = 270.0 # Km for CO2 at 25 C (μmol/mol)
    E_Kmc::T = 80.99e3 # Activation energy of Kmc (J/mol)
    # Michaelis-Menten constants oxygenation
    Kmo25::T = 165.0e3 # Km for O2 at 25 C (umol/mol)
    E_Kmo::T = 23.72e3 # Activation energy of Kmo (J/mol)
    # Rubisco activity
    Vcmax25::T = 120.0 # Maximum rate of carboxylation at 25 C (μmol/m2/s)
    E_Vcmax::T = 65.33e3 # Activation energy of Vcmax (J/mol)
    # Electron transport
    theta::T  = 0.7 # Curvature parameter
    Phi2::T   = 0.82 # Low-light PSII quantum yield
    sigma2::T = 0.5# Partitioning of excitation between PSII and PSI
    beta::T = 0.85 # Leaf absorptance of PAR
    fcyc::T = 0.1 #  Fraction of electrons at PSI that follow cyclic transport around PSI
    fpseudo::T = 0.05 # Fraction of electrons at PSI that are used by alternative electron sinks
    Jmax25::T = 230.0 # Maximum rate of electron transport (μmol/m2/s)
    E_Jmax::T = 30.0e3 # Activation energy Jmax (J/mol)
    D_Jmax::T = 200.0e3 # Deactivation energy of Jmax (J/mol)
    Topt_Jmax::T = 300.5 # Optimal temperature for Jmax (K)
    # Triose phosphate utilisation
    TPU25::T = 12.0 # Maximum rate of triose phosphate utilisation (μmol/m2/s)
    E_TPU::T = 53.1e3 # Activation energy TPU (J/mol)
    D_TPU::T = 20.18e3 # Deactivation energy of TPU (J/mol)
    Topt_TPU::T = 306.5 # Optimal temperature for TOU (K)
    # Respiration
    Rd25::T = 1.2 # Respiration rate at 25 C (μmol/m2/s)
    E_Rd::T = 46.39e3 # Activation energy of Rd (J/mol)
    # Mesophyll conductance
    gm25::T = 0.4 # Mesophyll conductance at 25 C (mol/m2/s)
    E_gm::T = 49.6e3 # Activation energy of gm (J/mol)
    D_gm::T = 437.4e3 # Deactivation energy of gm (J/mol)
    Topt_gm::T = 308.6 # Optimal temperature for gm (K)
    # Stomatal conductance
    gso::T = 0.01 # Minimum stomatal conductance to fluxes of CO2 in darkness (mol/m2/s)
    a1::T = 0.85 # Empirical parameter in gs formula
    b1::T = 0.14e-3 # Empirical parameter in gs formula (1/kPa)
end

"""
    C3Q(Sco25 = 2800.0, E_Sco = -24.46e3J/mol, Kmc25 = 270.0μmol/mol, E_Kmc = 80.99e3J/mol,
         Kmo25 = 165.0e3μmol/mol, E_Kmo = 23.72e3J/mol, Vcmax25 = 120.0μmol/m^2/s, E_Vcmax = 65.33e3J/mol,
         theta = 0.7, Phi2 = 0.82, sigma2 = 0.5, beta = 0.85, fcyc = 0.1, fpseudo = 0.05, 
         Jmax25 = 230.0μmol/m^2/s, E_Jmax = 30.0e3J/mol, D_Jmax = 200.0e3J/mol, Topt_Jmax = 300.5K, 
         TPU25 = 12.0μmol/m^2/s, E_TPU = 53.1e3J/mol, D_TPU = 201.8e3J/mol, Topt_TPU = 306.5K, 
         Rd25 = 1.2μmol/m^2/s, E_Rd = 46.39e3J/mol, gm25 = 0.4mol/m^2/s, E_gm = 49.6e3J/mol, 
         D_gm = 437.4e3J/mol, Topt_gm = 308.6K, gso = 0.01mol/m^2/s, a1 = 0.85, b1 = 0.14e-3/Pa)

Data structure to store all the parameters for the C3 photosynthesis model using
`Quantity` objects from Unitful.jl.

# Arguments
- `Sco25`: Sc/o parameter at 25 C
- `E_Sco`: Apparent activation energy of Sc/o (J/mol)
- `Kmc25`: Km for CO2 at 25 C (μmol/mol)
- `E_Kmc`: Activation energy of Kmc (J/mol)
- `Kmo25`: Km for O2 at 25 C (umol/mol)
- `E_Kmo`: Activation energy of Kmo (J/mol)
- `Vcmax25`: Maximum rate of carboxylation at 25 C (μmol/m2/s)
- `E_Vcmax`: Activation energy of Vcmax (J/mol)
- `theta`: Curvature parameter
- `Phi2`: Low-light PSII quantum yield
- `sigma2`: Partitioning of excitation between PSII and PSI
- `beta`: Leaf absorptance of PAR
- `fcyc`:  Fraction of electrons at PSI that follow cyclic transport around PSI
- `fpseudo`: Fraction of electrons at PSI that are used by alternative electron sinks
- `Jmax25`: Maximum rate of electron transport (μmol/m2/s)
- `E_Jmax`: Activation energy Jmax (J/mol)
- `D_Jmax`: Deactivation energy of Jmax (J/mol)
- `Topt_Jmax`: Optimal temperature for Jmax (K)
- `TPU25`: Maximum rate of triose phosphate utilisation (μmol/m2/s)
- `E_TPU`: Activation energy TPU (J/mol)
- `D_TPU`: Deactivation energy of TPU (J/mol)
"""
Base.@kwdef mutable struct C3Q{T <: Real} <: C3Type
    # Rubisco CO2/O2 specificity
    Sco25::T = 2800.0 # Sc/o parameter at 25 C
    E_Sco::Quantity{T, dimension(J/mol)} = -24.46e3J/mol # Apparent activation energy of Sc/o (J/mol)
    # Michaelis-Menten constants carboxylation (J/mol/K)
    Kmc25::Quantity{T, dimension(μmol/mol)} = 270.0μmol/mol # Km for CO2 at 25 C (μmol/mol)
    E_Kmc::Quantity{T, dimension(J/mol)} = 80.99e3J/mol # Activation energy of Kmc (J/mol)
    # Michaelis-Menten constants oxygenation
    Kmo25::Quantity{T, dimension(μmol/mol)} = 165.0e3μmol/mol # Km for O2 at 25 C (umol/mol)
    E_Kmo::Quantity{T, dimension(J/mol)} = 23.72e3J/mol # Activation energy of Kmo (J/mol)
    # Rubisco activity
    Vcmax25::Quantity{T, dimension(μmol/m^2/s)} = 120.0μmol/m^2/s # Maximum rate of carboxylation at 25 C (μmol/m2/s)
    E_Vcmax::Quantity{T, dimension(J/mol)} = 65.33e3J/mol # Activation energy of Vcmax (J/mol)
    # Electron transport
    theta::T = 0.7 # Curvature parameter
    Phi2::T = 0.82 # Low-light PSII quantum yield
    sigma2::T = 0.5 # Partitioning of excitation between PSII and PSI
    beta::T = 0.85 # Leaf absorptance of PAR
    fcyc::T = 0.1 #  Fraction of electrons at PSI that follow cyclic transport around PSI
    fpseudo::T = 0.05 # Fraction of electrons at PSI that are used by alternative electron sinks
    Jmax25::Quantity{T, dimension(μmol/m^2/s)} = 230.0μmol/m^2/s # Maximum rate of electron transport (μmol/m2/s)
    E_Jmax::Quantity{T, dimension(J/mol)} = 30.0e3J/mol # Activation energy Jmax (J/mol)
    D_Jmax::Quantity{T, dimension(J/mol)} = 200.0e3J/mol # Deactivation energy of Jmax (J/mol)
    Topt_Jmax::Quantity{T, dimension(K)} = 300.5K # Optimal temperature for Jmax (K)
    # Triose phosphate utilisation
    TPU25::Quantity{T, dimension(μmol/m^2/s)} = 12.0μmol/m^2/s # Maximum rate of triose phosphate utilisation (μmol/m2/s)
    E_TPU::Quantity{T, dimension(J/mol)} = 53.1e3J/mol # Activation energy TPU (J/mol)
    D_TPU::Quantity{T, dimension(J/mol)} = 201.8e3J/mol # Deactivation energy of TPU (J/mol)
    Topt_TPU::Quantity{T, dimension(K)} = 306.5K #  Optimal temperature for TPU (K)
    # Respiration
    Rd25::Quantity{T, dimension(μmol/m^2/s)} = 1.2μmol/m^2/s # Respiration rate at 25 C (μmol/m2/s)
    E_Rd::Quantity{T, dimension(J/mol)} = 46.39e3J/mol # Activation energy of Rd (J/mol)
    # Mesophyll conductance
    gm25::Quantity{T, dimension(mol/m^2/s)} = 0.4mol/m^2/s # Mesophyll conductance at 25 C (mol/m2/s)
    E_gm::Quantity{T, dimension(J/mol)} = 49.6e3J/mol # Activation energy of gm (J/mol)
    D_gm::Quantity{T, dimension(J/mol)} = 437.4e3J/mol # Deactivation energy of gm (J/mol)
    Topt_gm::Quantity{T, dimension(K)} = 308.6K # Optimal temperature for gm (K)
    # Stomatal conductance
    gso::Quantity{T, dimension(mol/m^2/s)} = 0.01mol/m^2/s # Minimum stomatal conductance to fluxes of CO2 in darkness (mol/m2/s)
    a1::T = 0.85 # Empirical parameter in gs formula
    b1::Quantity{T, dimension(1/kPa)} = 0.14e-3/Pa # Empirical parameter in gs formula (1/kPa)
end

"""
    photosynthesis(par::C3, PAR = 1000.0, RH = 0.75, Tleaf = 298.0, Ca = 400.0, O2 = 210e3, gb = 0.5)
    photosynthesis(par::C4, PAR = 1000.0, RH = 0.75, Tleaf = 298.0, Ca = 400.0, O2 = 210e3, gb = 0.5)
    photosynthesis(par::C3Q, PAR = 1000.0μmol/m^2/s, RH = 0.75, Tleaf = 298.0K, Ca = 400.0μmol/mol, O2 = 210e3μmol/mol, gb = 0.5mol/m^2/s)
    photosynthesis(par::C4Q, PAR = 1000.0μmol/m^2/s, RH = 0.75, Tleaf = 298.0K, Ca = 400.0μmol/mol, O2 = 210e3μmol/mol, gb = 0.5mol/m^2/s)

Calculate net CO2 assimilation (umol/m2/s), transpiration (mmol/m2/s) and stomatal
condutance to fluxes of CO2 (mol/m2/s) as a function of photosynthetically active
radiation (PAR, umol/m2/s), relative humidity (RH), leaf temperature (Tleaf,
K), air CO2 partial pressure (Ca, μmol/mol), oxygen (O2, μmol/mol) and boundary layer 
conductance to CO2 (gb, mol/m2/s). Environmental inputs must be scalar. 
"""
function photosynthesis(p::C3; PAR = 1000.0, RH = 0.75, Tleaf = 298.0, Ca = 400.0, O2 = 210e3, gb = 0.5)
    photosynthesis(p, PAR, RH, Tleaf, Ca, O2, gb)
end
function photosynthesis(p::C3Q; PAR = 1000.0μmol/m^2/s, RH = 0.75, Tleaf = 298.0K, 
              Ca = 400.0μmol/mol, O2 = 210e3μmol/mol, gb = 0.5mol/m^2/s)
    photosynthesis(p, PAR, RH, Tleaf, Ca, O2, gb)
end

function photosynthesis(p::C3Type, PAR, RH, Tleaf, Ca, O2, gb)
    # Calculate VPD at leaf temperature (kPa)
    vpd = VPD(Tleaf, RH)
    fvpd = max(1.0/(1.0/(p.a1 - p.b1*vpd) - 1.0),0.0)

    # Calculate values of parameters at Tleaf
    Rd = arrhenius(p.Rd25, p.E_Rd, Tleaf) # μmol/m2/s
    Sco = arrhenius(p.Sco25, p.E_Sco, Tleaf)
    Vcmax = arrhenius(p.Vcmax25, p.E_Vcmax, Tleaf) # μmol/m2/s
    Kmc = arrhenius(p.Kmc25, p.E_Kmc, Tleaf) # μmol/mol
    Kmo = arrhenius(p.Kmo25, p.E_Kmo, Tleaf) # mmol/mol
    Jmax = peaked(p.Jmax25, p.E_Jmax, p.D_Jmax, p.Topt_Jmax, Tleaf) # μmol/m2/s
    TPU = peaked(p.TPU25, p.E_TPU, p.D_TPU, p.Topt_TPU, Tleaf) # μmol/m2/s
    gm = peaked(p.gm25, p.E_gm, p.D_gm, p.Topt_gm, Tleaf) # mol/m2/s

    # Parameters for CO2 diffusion
    gamma_star = 0.5*O2/Sco # μmol/mol
    Ci_star = gamma_star - Rd/gm # μmol/mol

    # Photosynthesis limited by Rubisco
    x1_c = Vcmax # μmol/m2/s
    x2_c = Kmc*(1 + O2/Kmo) # μmol/mol
    Ac = solveAC3(gm, gb, p.gso, fvpd, x2_c, x1_c, gamma_star, Rd, Ca) # μmol/m2/s

    # Limitation by electron transport
    alpha = p.Phi2*p.sigma2*p.beta*(1 - p.fpseudo/(1 - p.fcyc))
    J = (alpha*PAR + Jmax - sqrt((alpha*PAR + Jmax)^2.0 - 4.0*p.theta*alpha*Jmax*PAR))/(2*p.theta) # μmol/m2/s
    x1_j = J/4.0 # μmol/m2/s
    x2_j = 2*gamma_star # μmol/mol
    Aj = solveAC3(gm, gb, p.gso, fvpd, x2_j, x1_j, gamma_star, Rd, Ca) # μmol/m2/s

    # Limitation by TPU
    x1_p = 3.0*TPU # μmol/m2/s
    x2_p = -gamma_star # μmol/mol
    Ap = 3.0*TPU - Rd # μmol/m2/s

    # Minimum of three potential rates
    A = min(Ac,min(Aj,Ap)) # μmol/m2/s

    # Stomatal conductance
    gsc = solvegs(p.gso,  A,  Ca, Ci_star,  Rd,  fvpd, gb) # mol/m2/s

    return (A = A, gs = gsc) 

end


# Analytical solution of cubic equation due to coupling of A, gb, gs and gm
function solveAC3(gm, gb, gso, fvpd, x2, x1, gamma_star, Rd, Ca)
    a = gso*(x2 + gamma_star) + (gso/gm + fvpd)*(x1 - Rd)
    b = Ca*(x1 - Rd) - gamma_star*x1 - Rd*x2
    c = Ca + x2 + (x1 - Rd)*(1/gm + 1/gb)
    d = x2 + gamma_star + (x1 - Rd)/gm
    m = 1/gm + (gso/gm + fvpd)*(1/gm + 1/gb)
    p = -(d + (x1 - Rd)/gm + a*(1/gm + 1/gb) + (gso/gm + fvpd)*c)/m
    q = (d*(x1 - Rd) + a*c + (gso/gm + fvpd)*b)/m
    r = (-a*b)/m
    U = (2*p^3 - 9*p*q + 27*r)/54
    Q = (p^2 - 3*q)/9
    psi = acos(min(1.0, U/sqrt(Q^3)))
    A = -2*sqrt(Q)*cos(psi/3) - p/3
  end
  
# Calculate gs once A is known
function solvegs(gso, A, Ca, Ci_star, Rd, fvpd, gb)
    a = Ca - A/gb - Ci_star
    b = -A - Ca*gso + gso*Ci_star - (A + Rd)*fvpd
    c = A*gso
    A = (-b - sqrt(b*b - 4*a*c))/(2*a)
  end
  