# import the pysb module and all its methods and functions
from pysb import *
from pysb.macros import *
from pysb.bng import *

def catalyze_monosubstrate(enz, sub, product, klist):
    """Alias for pysb.macros.catalyze_two_step_reversible with default binding site 'b1'.
    """
    site_name = 'b1'
    return catalyze_two_step_reversible(enz, site_name, sub, site_name, product, site_name, klist)

def catalyze_bisubstrate(enz, sub1, sub2, product, klist):
    """Alias for pysb.macros.catalyze_ordered_bisubstrate_reversible with default binding site 'b'.
    """
    enz_site1_name = 'b1'
    enz_site2_name = 'b2'
    subs_site_name = 'b1'

    return catalyze_ordered_bisubstrate_reversible(enz, enz_site1_name, enz_site2_name, sub1, subs_site_name, sub2, subs_site_name, product, subs_site_name, klist)
    
def define_parameters(parameter_table):
    
    for rn in parameter_table.keys():
        Parameter('k%sf' % rn, parameter_table[rn][0])
        Parameter('k%sb' % rn, parameter_table[rn][1])

def CrtEBI_Kinetic_Model(TIR_crtE, TIR_crtB, TIR_crtI):

    # instantiate a model
    Model()

    #Reference Pathway Variant
    ref_TIR_crtE = 72268.0
    ref_TIR_crtB = 20496.0
    ref_TIR_crtI = 203462.0
    ref_crtEBI_productivity = 196.0 
    Parameter('relative_TIR_crtE', TIR_crtE / ref_TIR_crtE)
    Parameter('relative_TIR_crtB', TIR_crtB / ref_TIR_crtB)
    Parameter('relative_TIR_crtI', TIR_crtI / ref_TIR_crtI)
    Parameter('one_constant', 1.0)
    Parameter('zero_constant', 0.0)
    
    # declare monomers
    #enzymes first.  b1 = 'binding site'.   
    Monomer('idi',  sites=['b1'])
    Monomer('ispA', sites=['b1','b2'])
    Monomer('crtE', sites=['b1','b2'])
    Monomer('crtB', sites=['b1','b2'])
    Monomer('crtI', sites=['b1'])

    #free substrates second.  b1 = 'binding site'
    Monomer('IPP', sites=['b1'])
    Monomer('DMAPP', sites=['b1'])
    Monomer('GPP', sites=['b1'])
    Monomer('FPP', sites=['b1'])
    Monomer('GGPP', sites=['b1'])
    Monomer('prephytoene2P', sites=['b1'])
    Monomer('phytoene', sites=['b1'])
    Monomer('phytofluene', sites=['b1'])
    Monomer('zeta_carotene', sites=['b1'])
    Monomer('neurosporene', sites=['b1'])

    Initial(crtE(b1=None,b2=None), relative_TIR_crtE)
    Initial(crtB(b1=None, b2=None), relative_TIR_crtB)
    Initial(crtI(b1=None), relative_TIR_crtI)
    Initial(idi(b1=None), one_constant)
    Initial(ispA(b1=None, b2=None), one_constant)

    for species in (IPP,DMAPP,GPP,FPP,GGPP,prephytoene2P,phytoene,phytofluene,zeta_carotene, neurosporene):
        Initial(species(b1=None), zero_constant)
        
    #intermediate complexes will be auto-generated!
    #kinetic parameter definitions
    f_IPP = 6.0 / 7.0
    f_DMAPP = 1.0 / 7.0
    kinetic_parameters = {  1 : [8.67, 0.025],
                            2 : [0.84, 0.072],
                            3 :	[183, 609],
                            4 :	[690, 79],
                            5 : [474, 59],
                            6 : [549, 55],
                            7 : [776, 89],
                            8 : [2.68, 133],
                            9 : [283, 0.19],
                           10 : [0.43, 0.42],
                           11 : [2.66, 84.39],
                           12 : [10560, 30.24],
                           13 : [4, 6.09],
                           14 : [87, 2970],
                           15 : [2.21, 84],
                           16 : [2690, 34],
                           17 : [3.21, 23],
                           18 : [20, 0.76],
                           19 : [1.90, 194],
                           20 : [1138, 3.81],
                           21 : [11.39, 372],
                           22 : [411, 10],
                           23 : [11.39, 290],
                           24 : [106, 2.73]
    }

    #Add Parameters objects for kinetic parameters
    define_parameters(kinetic_parameters)

    synthesize(IPP(b1=None), f_IPP)
    synthesize(DMAPP(b1=None), f_DMAPP)
    catalyze_monosubstrate(idi, IPP, DMAPP, [k1f,k1b,k2f,k2b])
    catalyze_bisubstrate(ispA, DMAPP, IPP, GPP, [k3f,k3b,k4f,k4b,k5f,k5b])
    catalyze_bisubstrate(ispA, GPP, IPP, FPP, [k6f,k6b,k7f,k7b,k8f,k8b])
    catalyze_bisubstrate(crtE, GPP, IPP, FPP, [k9f,k9b,k10f,k10b,k11f,k11b])
    catalyze_bisubstrate(crtE, FPP, IPP, GGPP, [k12f,k12b,k13f,k13b,k14f,k14b])
    catalyze_bisubstrate(crtB, GGPP, GGPP, prephytoene2P, [k15f,k15b,k15f,k15b,k16f,k16b])
    catalyze_monosubstrate(crtB, prephytoene2P, phytoene, [k17f,k17b,k18f,k18b])
    catalyze_monosubstrate(crtI, phytoene, phytofluene, [k19f,k19b,k20f,k20b])
    catalyze_monosubstrate(crtI, phytofluene, zeta_carotene, [k21f,k21b,k22f,k22b])
    catalyze_monosubstrate(crtI, zeta_carotene, neurosporene, [k23f,k23b,k24f,k24b])

    Observable('neurosporene_out', neurosporene(b1=None))
    
    t = range(0, 7*3600, 60)
    from pysb.integrate import Solver
    solver = Solver(model, t)
    print model.species
    return solver
    
    
    
    
if __name__ == "__main__":
    ref_TIR_crtE = 72268.0
    ref_TIR_crtB = 20496.0
    ref_TIR_crtI = 203462.0
    ref_crtEBI_productivity = 196.0 
    solver = CrtEBI_Kinetic_Model(ref_TIR_crtE, ref_TIR_crtB, ref_TIR_crtI)
    solver.run()
    print solver.yobs

    TL = len(solver.yobs)
    sim_neurosporene = solver.yobs[TL-1][0]
    sim_ref_ratio = sim_neurosporene / ref_crtEBI_productivity
    
    solver = CrtEBI_Kinetic_Model(1000.0, 40000.0, 200000.0)
    solver.run()
    print solver.yobs
    sim_neurosporene = solver.yobs[TL-1][0]
    print sim_neurosporene / sim_ref_ratio, " ug/gDCW/hour"
    

