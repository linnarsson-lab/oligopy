import pandas as pd
import numpy as np
from math import factorial
import random
import matplotlib.pyplot as plt
from dask import delayed, compute
from scipy.stats import rankdata
from datetime import datetime
import glob
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import pickle as pkl
 
def optimize_gene_barcode_random(genes, codebook, df_exp, cycles=None, trials=1e6, lean=False, plot=True, save_folder=None):
    """
    Randomly optimizes the assignment of genes to barcodes to minimize maximum expression overlap across cycles.

    This function evaluates many random permutations of gene-to-barcode assignments to find the configuration 
    that minimizes the maximum summed expression across all cell types and barcode cycles. It is useful for 
    optimizing gene panels in spatial or multiplexed transcriptomics where co-labeling should be minimized.

    Args:
        genes (list): List of gene names. Include dummy entries for empty barcodes (e.g., 'Empty_barcode_X').
        codebook (np.ndarray): Boolean matrix of shape (n_barcodes, n_cycles) defining active cycles per barcode.
        df_exp (pd.DataFrame or list of pd.DataFrame): Expression matrix (genes × cell types); supports a list of datasets.
        cycles (int, optional): Number of barcode cycles; if None, inferred from codebook. Default is None.
        trials (int): Number of random permutations to evaluate. Default is 1e6.
        lean (bool): If True, discard intermediate results to reduce memory usage. Default is False.
        plot (bool): If True, plot optimization progress and results. Default is True.
            If lean is True, the plots will only show the distribution of the best solutions of each worker.
            This will thus misrepresent the actual random distribution
        save_folder (str, optional): Path to save plots and results; if None, results are not saved. Default is None.

    Returns:
        best_gene_order (np.ndarray): Gene ordering corresponding to the best (lowest) max expression score.
        results (list): List of result dictionaries for each permutation, each containing:
            - 'gene_order': the shuffled gene list
            - 'max_result': maximum value across cycles per dataset
            - 'full_results' (optional): per-cycle expression values (if lean=False)
        best_row_idx (int): Index of the best result in the results list.
    """

    n_genes = len(genes)
    n_barcodes = codebook.shape[0]
    if n_barcodes > n_genes:
        blanks = [f'Empty_barcode_{i}' for i in range(n_barcodes-n_genes)]
        genes = np.hstack([genes, blanks])

    if cycles is None:
        cycles = codebook.shape[1]
    df_exp = [df_exp] if not isinstance(df_exp, list) else df_exp

    #Get input in optimal formats
    df_gene_to_idx = [{gene: i for i, gene in enumerate(df.index)} for df in df_exp]
    df_exp_numpy = [df.to_numpy() for df in df_exp]
    codebook = codebook.astype(bool)

    #Delay all input
    genes_delayed = delayed(genes)
    df_exp_delayed = [delayed(df) for df in df_exp_numpy]
    codebook_delayed = delayed(codebook)
    df_gene_to_idx_delayed = delayed(df_gene_to_idx)


    def worker(genes, df_exp_numpy, df_gene_to_idx, cycles, codebook, lean=False):

        shuffled_genes = np.random.permutation(genes)
        r = {'gene_order': shuffled_genes,
            'full_results' : np.zeros((len(df_exp_numpy), cycles))}

        for cycle in range(cycles):
            positive_genes = shuffled_genes[codebook[:,cycle]]
            for i, df in enumerate(df_exp_numpy):
                positive_genes_indexes = [df_gene_to_idx[i][g] for g in positive_genes if g in df_gene_to_idx[i]]
                max_sum_exp = df_exp_numpy[i][positive_genes_indexes].sum(axis=0).max()
                r['full_results'][i, cycle] = max_sum_exp
        r['max_result'] = r['full_results'].max(axis=1)

        if lean: #Save memory
            del r['full_results']
        return r

    def batch_worker(genes, df_exp_numpy, df_gene_to_idx, cycles, codebook, batch_size=100, lean=False):   

        #Unpack delayed input
        genes = genes.compute() if hasattr(genes, 'compute') else genes
        df_exp_numpy = [df.compute() if hasattr(df, 'compute') else df for df in df_exp_numpy]
        df_gene_to_idx = df_gene_to_idx.compute() if hasattr(df_gene_to_idx, 'compute') else df_gene_to_idx
        codebook = codebook.compute() if hasattr(codebook, 'compute') else codebook

        if lean: 
            to_beat = np.array([np.inf for _ in df_exp_numpy])
            result = [None]
            for i in range(batch_size):
                r = worker(genes, df_exp_numpy, df_gene_to_idx, cycles, codebook, lean=False)
                if np.all(r['max_result'] < to_beat):
                    result[0] = r
            return result
        else:
            return [worker(genes, df_exp_numpy, df_gene_to_idx, cycles, codebook, lean=False) for _ in range(batch_size)]

    #Compute 
    ncores = len(client.cluster.workers)
    if ncores < 1:
        ncorese = 1
    batch_size = int(np.ceil(trials / ncores))
    print(f'Computing {ncores} batches of {batch_size} trials, total: {ncores * batch_size} trials')

    delayed_results = [delayed(batch_worker)(genes_delayed, df_exp_delayed, df_gene_to_idx_delayed, cycles, codebook_delayed, batch_size=batch_size, lean=lean) for _ in range(ncores)]
    batch_results = compute(*delayed_results)
    
    #Compile results
    results = [item for batch in batch_results for item in batch]
    compiled_results = np.vstack([r['max_result'] for r in results])

    #Rank each column (lower values = lower rank) and find gene shuffle with lowest rank
    ranks = np.vstack([rankdata(compiled_results[:, col], method='ordinal') for col in range(compiled_results.shape[1])]).T
    rank_sums = ranks.sum(axis=1)
    best_row_idx = np.argmin(rank_sums)
    best_row = compiled_results[best_row_idx]
    best_gene_order = results[best_row_idx]['gene_order']
    best_order_results = results[best_row_idx]['full_results']
    
    mean = np.mean(compiled_results, axis=0)
    std = np.std(compiled_results, axis=0)
    delta = (mean - best_row) / std
    delta2 = (mean - compiled_results.min(axis=0)) / std
    print(f'The solution is {delta} standard deviations below the mean of a random gene order (respectively for each dataset).')
    print('Depending on your data expect your solution to be at least ~2std below the mean. If not, try more itterations.')
    print(f'The best solution per dataset that was found {delta2} standard deviations below the mean.')
    print('Likely reaching the optimal solution unattainable for multiple datasets, but if your results are far away from the individual optimum, you could try more itterations.')
    
    #Save
    if save_folder is not None:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        result_to_save = {'best_gene_order': best_gene_order,
                          'best_results': results[best_row_idx]}
        pkl.dump(result_to_save, open(f'{save_folder}/{timestamp}_optimal_gene_order_{trials}_trials.pkl', 'wb'))    
        save_figures = True
    else:
        save_figures = False
    
    
    if plot:
        plt.figure()
        plt.violinplot(compiled_results)
        plt.scatter(np.arange(compiled_results.shape[1])+1, best_row, c='r', label='chosen gene order')
        plt.ylabel('Maximum summed expression of gene permutation')
        plt.legend()
        plt.title('Max expression values')
        if save_figures:
            plt.savefig(f'{save_folder}/{timestamp}_Maximum_plot_{trials}trials.pdf', bbox_inches='tight')
        
        fig, axes = plt.subplots(figsize=(15,5), ncols=compiled_results.shape[1]+1)
        for j in range(compiled_results.shape[1]):
            ax = axes[j]
            y = [compiled_results[0,j]]
            x = [0]
            for i, r in enumerate(compiled_results[:,j]):
                if r < y[-1]:
                    y.append(r)
                    x.append(i)
            y.append(y[-1])
            x.append(len(compiled_results[:,j]))
            ax.plot(x, y, label='Minimum')
            ax.scatter(x[-1], best_row[j], c='g', label='Cross dataset optimum')
            ax.hlines(np.mean(compiled_results[:,j]), 0, x[-1], colors='r', label='Mean value')
            ax.set_xlabel('Iteration with improvement')
            ax.set_ylabel('Lowest max value')
            ax.set_title(f'Dataset {j} optimization')
        ax.legend()

        ax = axes[-1]
        y = [rank_sums[0]]
        x = [0]
        for i, r in enumerate(rank_sums):
            if r < y[-1]:
                y.append(r)
                x.append(i)
        y.append(y[-1])
        x.append(len(rank_sums))
        ax.plot(x, y, c='g')
        ax.set_title('Rank decrease')
        ax.set_xlabel('Iteration with improvement')
        plt.tight_layout()
        if save_figures:
            plt.savefig(f'{save_folder}/{timestamp}_Iteration_plot_{trials}trials.pdf', bbox_inches='tight')
    
    #Can only plot the distribution over the cycles if we have all the data. 
    if lean==False:
        mean_results = []
        for i in range(len(df_exp)):
            mean_results.append(np.vstack([r['full_results'][i] for r in results]))

        fig, axes = plt.subplots(nrows=3, sharex=True, figsize=(12,8))
        for ax, mean_r, best_r in zip(axes, mean_results, best_order_results):
            ax.violinplot(mean_r)
            x = np.linspace(1, len(best_r), len(best_r))
            ax.scatter(x, best_r, c='r')
            #ax.scatter(x, random_r, c='k')
            title_values = [np.max(best_r), np.mean(best_r), np.std(best_r), np.mean(mean_r.ravel()), np.std(mean_r.ravel())]
            title_values = [round(i, 2) for i in title_values]
            ax.set_title(f'Best: max{title_values[0]} mean{title_values[1]}±{title_values[2]}, Mean: {title_values[3]}±{title_values[4]}')
            ax.set_ylabel('Summed expression per cycle')
        ax.set_xlabel('Cycle')
        if save_figures:
            plt.savefig(f'{save_folder}/{timestamp}_Max_per_cycle_plot_{trials}trials.pdf', bbox_inches='tight')
            
    return best_gene_order, results, best_row_idx
                                                     

def evaluate_order(genes, cycles, df_exp, codebook, df_gene_to_idx):
    """
    Computes the summed expression for a given gene order and codebook assignment.

    For each cycle, it identifies the genes that are labeled in that cycle based on the codebook,
    retrieves their expression values from each dataset, and calculates the maximum summed expression
    across cell types.

    Args:
        genes (list): Ordered list of genes (must match order used in codebook).
        cycles (int): Number of barcode cycles.
        df_exp (list of np.ndarray): List of gene expression matrices (genes × cell types).
        codebook (np.ndarray): Boolean matrix indicating gene labeling per cycle.
        df_gene_to_idx (list of dict): Maps gene names to row indices in df_exp per dataset.

    Returns:
        np.ndarray: Array of shape (n_datasets, n_cycles) with maximum summed expression per cycle.
    """
    r =  np.zeros((len(df_exp), cycles))

    for cycle in range(cycles):
        positive_genes = genes[codebook[:,cycle]]
        for i, df in enumerate(df_exp):
            positive_genes_indexes = [df_gene_to_idx[i][g] for g in positive_genes if g in df_gene_to_idx[i]]
            max_sum_exp = df_exp[i][positive_genes_indexes].sum(axis=0).max()
            r[i, cycle] = max_sum_exp
    return r

def swap_value(arr, value):
    """
    Randomly swaps a given value with another entry in an array.

    Useful for mutating a gene order during optimization.

    Args:
        arr (np.ndarray): Array of gene names.
        value (str): Gene name to be swapped with another.

    Returns:
        np.ndarray: New array with the value swapped at a random position.
    """

    # Find index of the value
    matches = np.where(arr == value)[0]
    if len(matches) == 0:
        raise ValueError(f"Value '{value}' not found in array.")
    i = matches[0]

    # Pick a different index to swap with
    n = len(arr)
    while True:
        j = random.randint(0, n - 1)
        if j != i:
            break

    # Perform the swap
    arr = arr.copy()  # avoid modifying the original array
    arr[i], arr[j] = arr[j], arr[i]
    return arr

def optimize_gene_barcode_evolution(genes, codebook, df_exp, cycles=None,
                                    depth=1, maxiter=10, max_plateau=2,
                                    plot=True, save_folder=None):
    """
    Evolutionary optimization of gene-to-barcode assignment using cycle-specific mutation.

    Iteratively improves gene-to-barcode assignments by mutating gene orders to reduce
    the maximum summed expression across barcode cycles and datasets. At each iteration,
    the worst-performing cycle is identified and targeted for mutation.

    The last iteration stored in the dictionary will contain the best result.

    Args:
        genes (list): List of gene names. Include 'Empty_barcode_X' entries as needed.
        codebook (np.ndarray): Boolean array (n_barcodes × n_cycles) defining barcode design.
        df_exp (pd.DataFrame or list of pd.DataFrame): Expression matrix/matrices (genes × cell types).
        cycles (int, optional): Number of barcode cycles. Defaults to codebook shape if None.
        depth (int): Number of mutations per gene per iteration. Default is 1.
        maxiter (int): Maximum number of optimization iterations. Default is 10.
        max_plateau (int): Stop if no improvement is found in this many consecutive iterations. Default is 2.
        plot (bool): If True, plots optimization progress. Default is True.
        save_folder (str, optional): Folder to save results and figures. Required for persistence.
            In the folder a plot will be constantly updated to monitor progress.

    Returns:
        results (dict): Dictionary keyed by iteration index. Each value contains:
            - 'gene_order': Optimized gene list
            - 'max_values': Max per-dataset expression for that iteration
            - 'Cycle_to_optimize': Cycle targeted for mutation
            - 'genes_per_cycle': Number of labeled genes per cycle
    """

    n_genes = len(genes)
    n_barcodes = codebook.shape[0]
    if n_barcodes > n_genes:
        blanks = [f'Empty_barcode_{i}' for i in range(n_barcodes-n_genes)]
        genes = np.hstack([genes, blanks])

    if cycles is None:
        cycles = codebook.shape[1]
    df_exp = [df_exp] if not isinstance(df_exp, list) else df_exp

    #Get input in optimal formats
    df_gene_to_idx = [{gene: i for i, gene in enumerate(df.index)} for df in df_exp]
    df_exp_numpy = [df.to_numpy() for df in df_exp]
    codebook = codebook.astype(bool)
    
    #Track genes per cycle
    hamming_weight = codebook.sum(axis=1)[0]
    n_genes_per_cycle = []
    total_valid_genes = len([i for i in genes if not i.startswith("Empty_barcode_")])
    mean_n_genes = total_valid_genes * hamming_weight / cycles
    n_better_solutions = []

    #Delay all input
    df_exp_delayed = [delayed(df) for df in df_exp_numpy]
    codebook_delayed = delayed(codebook)
    df_gene_to_idx_delayed = delayed(df_gene_to_idx)
    
    #df_exp_delayed = [delayed(df).persist() for df in df_exp_numpy]
    #codebook_d = delayed(codebook).persist()
    #df_gene_to_idx_d = delayed(df_gene_to_idx).persist()

    plateau_counter = 0    
    results = {}
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    r = evaluate_order(genes, cycles,  df_exp_numpy, codebook, df_gene_to_idx)
    r_original_max = r.max(axis=1)

    for iteration in range(maxiter):

        #Evaluate gene order, find worst cycle and the genes that are positive
        r = evaluate_order(genes, cycles,  df_exp_numpy, codebook, df_gene_to_idx)
        current_max = r.max(axis=1) #Save the current best result
        rank_sums = rankdata(r, axis=1).sum(axis=0)
        worst_cycle = np.argmax(rank_sums)
        positive_genes = genes[codebook[:,worst_cycle]]
        print(f'Iteration: {iteration}, optimizing_cycle: {worst_cycle}', end='\r')

        #Mutate the positive genes
        to_test = []
        for g in positive_genes:
            for i in range(depth):
                to_test.append(swap_value(genes, g))

        #Evaluate the new gene orders   
        delayed_results = [delayed(evaluate_order)(gene_order, cycles, df_exp_delayed, codebook_delayed, df_gene_to_idx_delayed) for gene_order in to_test]
        #delayed_results = [evaluate_on_permutation(gene_order) for gene_order in to_test]
        #delayed_results = [delayed(evaluate_order)(gene_order, cycles, df_exp_numpy_future, codebook_future, df_gene_to_idx_future) for gene_order in to_test]        
        
        batch_results = compute(*delayed_results)

        #Select the best gene order
        compiled_results = []
        compiled_indexes = []
        for i, rr in enumerate(batch_results):
            new_max = rr.max(axis=1)
            #Save if new max is below the current best solution
            if np.all(new_max <= current_max):
                compiled_results.append(new_max)
                compiled_indexes.append(i)
        #Pick the best solution of the solutions that are better than the previously best result
        if len(compiled_indexes) != 0:
            n_better_solutions.append(len(compiled_indexes))
            rank_sums = rankdata(np.stack(compiled_results), axis=0)
            new_best = np.argmin(rank_sums.sum(axis=1))
            new_best_index = compiled_indexes[new_best]
            genes = to_test[new_best_index]

            #Calcualte simultaneous genes / cycle
            gpc = []
            for cycle in range(cycles):
                positive_genes = genes[codebook[:,cycle]]
                valid_genes = [i for i in positive_genes if not i.startswith("Empty_barcode_")]
                gpc.append(len(valid_genes))
            n_genes_per_cycle.append(np.array(gpc))
            
            #Save results
            results[iteration] = {'gene_order': genes,
                                 'max_values': compiled_results[new_best],
                                 'Cycle_to_optimize': worst_cycle,
                                 'genes_per_cycle': np.array(gpc)}            
            pkl.dump(results, open(f'{save_folder}/{timestamp}_Evolution_optimization.pkl', 'wb'))
            
            #Plot
            if plot: 
                #Plot improvement relative to absolute values
                to_plot = np.stack([results[k]['max_values'] for k in sorted(list(results.keys()))])
                fig, axes = plt.subplots(ncols=2, nrows=2, figsize=(10,10))
                ax0, ax1, ax2, ax3 = axes.ravel()
                x = np.arange(to_plot.shape[0]+1)
                for ci in range(to_plot.shape[1]):
                    y = to_plot[:,ci]
                    y = [1] + list(y/r_original_max[ci]) #Relative to original gene order best maximum
                    ax0.plot(x, y, label=f'Dataset: {ci}, Start: {r_original_max[ci]}')
                    ax0.text(x[-1], y[-1], round(to_plot[:,ci][-1], 3), verticalalignment='center' )
                ax0.set_ylim(0,1)
                ax0.spines[['right', 'top']].set_visible(False)
                ax0.legend(loc='lower left')
                ax0.set_xlabel('Iteration')
                ax0.set_ylabel('Relative difference to original max values')
                ax0.set_title('Relative loss')
                
                #Plot improvements but now zoom in
                for ci in range(to_plot.shape[1]):
                    y = to_plot[:,ci]
                    y = y/r_original_max[ci]
                    ax1.plot(x[1:], y, label=f'Dataset: {ci}')
                ax1.spines[['right', 'top']].set_visible(False)
                ax1.legend(loc='lower left')
                ax1.set_xlabel('Iteration')
                ax1.set_ylabel('Relative difference to original max values')
                ax1.set_title('Zoom')
                
                #Plot how number of genes per cycle changes over the itterations
                for i, data in enumerate(n_genes_per_cycle):
                    color = plt.cm.gnuplot(i/len(n_genes_per_cycle))
                    ax2.plot(data, c=color)
                # Add colorbar
                cmap = plt.cm.gnuplot
                norm = mcolors.Normalize(vmin=0, vmax=1)
                sm = cm.ScalarMappable(cmap=cmap, norm=norm)
                sm.set_array([])  # required for colorbar in some versions
                cbar = plt.colorbar(sm, ax=ax2, orientation='vertical', fraction=0.05, pad=0.05)
                cbar.set_ticks([0, 1])
                cbar.set_ticklabels(['Iteration 0', f'Iteration {len(n_genes_per_cycle)}'])
                cbar.ax.tick_params(labelsize=8)
                x0, x1 = ax2.get_xlim()
                ax2.hlines(mean_n_genes, x0, x1, colors='g')
                ax2.spines[['right', 'top']].set_visible(False)
                ax2.set_xlabel('Cycle')
                ax2.set_ylabel('genes/cycle')
                ax2.set_title('Genes simultaneously labeled per cycle')
                    
                ax3.plot(n_better_solutions)
                #x0, x1 = ax3.get_xlim()
                #ax3.hlines(len(to_test), x0, x1, colors='r', label='Tries')
                #ax3.legend(loc='lower left')
                ax3.spines[['right', 'top']].set_visible(False)
                ax3.set_xlabel('Iteration')
                ax3.set_ylabel('Number of better solutions')
                ax3.set_title(f'Better solutions out of {len(to_test)} tries')
                
                
                plt.tight_layout()
                fig.savefig(f'{save_folder}/{timestamp}_Evolution_optimization.png')
                plt.close(fig)
        else:
            #No improvement
            plateau_counter += 1 
            if plateau_counter > max_plateau:
                break
                
    print(f'Results can also be found in: {save_folder}/{timestamp}_Evolution_optimization.pkl')
    return results
