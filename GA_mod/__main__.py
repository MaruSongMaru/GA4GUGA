import time
import matplotlib.pyplot as plt
import run_GA as ga
import numpy as np
import process_df
import pandas as pd
import crossover as co

def main():
    df = pd.read_pickle('../Dataset/H8_PN-cluster.pkl')

    fitness_dict_tmp = process_df.gen_fitness_ht(df, 'Permutation', 'DominantCSFenergy1')
    fitness_dict = {key: 1/abs(value) for key, value in fitness_dict_tmp.items()}

    reference_dict = process_df.gen_fitness_ht(df, 'Permutation', 'L4Norm')

    co_fn = co.order_co

    start_time = time.time()
    fit_arr, ref_arr = ga.GA_ensembles(num_ensembles=100, pop_size=10,
                                       ordering_len=8, elite_size=2,
                                       mutation_rate=0.05, generations=100,
                                       co_function=co_fn,
                                       fitness_dict=fitness_dict,
                                       reference_dict=reference_dict)
    end_time = time.time()

    execution_time = end_time - start_time
    print("Execution time:", execution_time, "seconds")

    df_L4sort = df.sort_values(by='L4Norm')
    plt.figure()
    plt.plot(df_L4sort['L4Norm'].values)
    plt.savefig('L4_trend.png')

    plt.figure()
    plt.plot(df_L4sort['DominantCSFenergy1'].values)
    plt.savefig('DiagE_trend.png')

    plot_step = 10

    fit_arr = 1/fit_arr
    fig, axs = plt.subplots(1, 10, figsize=(15, 3))
    xstart = np.min(fit_arr)
    xend = np.max(fit_arr)
    data_range = xend - xstart
    num_bins = 20
    bin_width = data_range / num_bins
    bin_edges = np.arange(xstart, xend + bin_width, bin_width)
    for i, ax in enumerate(axs):
        ax.hist(fit_arr[:, i*plot_step], bins=bin_edges)
        ax.set_xlim(xstart-0.1*data_range, xend+0.1*data_range)
        ax.set_ylim(0, 100)
        if i == 0:
            ax.set_ylabel('Frequency')
        else:
            ax.set_yticklabels([])
        ax.set_title(f'Gen.={i*plot_step}')
    plt.savefig('GA_fitness.png')

    fig, axs = plt.subplots(1, 10, figsize=(15, 3))
    xstart = np.min(ref_arr)
    xend = np.max(ref_arr)
    data_range = xend - xstart
    num_bins = 20
    bin_width = data_range / num_bins
    bin_edges = np.arange(xstart, xend + bin_width, bin_width)
    for i, ax in enumerate(axs):
        ax.hist(ref_arr[:, i*plot_step], bins=bin_edges)
        ax.set_xlim(xstart-0.1*data_range, xend+0.1*data_range)
        ax.set_ylim(0, 100)
        if i == 0:
            ax.set_ylabel('Frequency')
        else:
            ax.set_yticklabels([])
        ax.set_title(f'Gen.={i*plot_step}')
    plt.savefig('GA_reference.png')

if __name__ == "__main__":
    main()
