import GA_mod.run_GA as ga                                                              
import GA_mod.crossover as co                                                           
from GA_mod.measure_fitness import FitnessFunction                                      
                                                                                 
FILEPATH = "../extra_files/"
fcidump = FILEPATH + "FCIDUMP_20site_Heisenberg"

ref_csf = [1,2,1,2,1,2,1,2,1,2,
           1,2,1,2,1,2,1,2,1,2]
ref_ordering = tuple(range(1,21))

ga.perform_GA(
        FitnessFunction.DIAG_ELEM_SMS_MAPPING,
        num_chroms=100,
        restricted_ordering_len=20,
        elite_size=10,
        mutation_rates=[0.01],
        generations=10000,
        co_function=co.order_co,
        fcidump=fcidump,
        norb=20,
        sms_ref_csf=ref_csf,
        tMinimize=True,
        sms_ref_ordering=ref_ordering
        )
