###In development


"""
    compute_ASRQ(f::Vector{Float64})

Compute assimilate requirement, CO2 production factor and carbon content according to Penning de Vries et al. 1982.
These values can be used to compute growth respiration. Different tissues (species/organ) have different composition.
Growth costs are computed given the fraction of dry matter allocated to different carbon pools according to specific composition of the tissue:
    f = [carbohydrates, proteins, lipids, lignin, organic acids, minerals]

    The assimilate requirement (ASRQ, g) is the amount of assimilates required to produce 1 g of dry matter.
    The CO2 production factor (CO2PFF, g) is the amount of CO2 produced per g of dry matter.
    The carbon fraction (CF, %) is the fraction of dry matter that is carbon.

Citations:
Penning de Vries, F. W. T., & van Laar, H. H. (1982). Simulation of growth processes and the model BACROS. In F. W. T. Penning de Vries, & H. H. van Laar (Eds.), Simulation of plant growth and crop production (pp. 114-135). (Simulation monographs). Pudoc. https://edepot.wur.nl/172216
"""
function compute_ASRQ(f::Vector{Float64} = [0.53, 0.25, 0.05, 0.05, 0.06, 0.06])
    CF = @SVector [0.4504, 0.5321, 0.7733, 0.6899, 0.3746, 0.0]
    ASRQ = @SVector [1.242, 2.7, 3.106, 2.174, 0.929, 0.05]
    CO2PFF = @SVector [0.170, 2.009, 1.720, 0.659, -0.011, 0.073]
    CF_val = sum(f.*CF)
    ASRQ_val = sum(f.*ASRQ)
    CO2PFF_val = sum(f.*CO2PFF)

    return (ASRQ = ASRQ_val, CO2PFF = CO2PFF_val, CF = CF_val)
end


