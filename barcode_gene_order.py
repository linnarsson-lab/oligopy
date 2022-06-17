import pandas as pd
import numpy as np
from math import factorial
import random
import matplotlib.pyplot as plt


def optimize_gene_barcode(genes, codebook, df_exp, cycles = 15, trials=50, plot=True):
    """
    Function to optimize how genes are divided over the possible barcodes.
    Given a binary code book, the function will randomly shuffle the genes and
    with the df_exp expression matrix it will sum the expression level of the
    genes that are co-labeled in the same cycle. For the number of trials it
    will return the permutation that has the lowest max expression in a cell
    type over all cycles. 
    
    Args:
    genes (list): List of genes.
    codebook (array): Binary array containing the barcodes. Number of rows
        equals the number of genes, colums is the number of cycles.
    df_exp (dataframe): Dataframe containing the genes as rows and cell
        types as columns. Values could be mean expression or max expression.
    cycles (int): Number of cycles to optimize for. Usually same as barcode
        length.
    trials (int): Number of permutations of the genes to try.
    plot (bool): Plots summary statistics of all trials and best trial.
    
    Returns:
    best_gene_order (list): Order of genes that gives the lowest max
        expression over all cell types and cycles for teh given barcodes.
    results_dict(dict): Dictionary containing all results.
        Structure: results_dict[permutaion_number]:
        ['gene_order', 'sum_counts', 'cycle_max', 'max']
        'gene_order': Contains the tested gene order,
        'sum_count': Contains a dataframe with the summed counts of the 
            co-labeled genes for all rounds.
        'cycle_max': Max of all cycles.
        'max': Max of full table.
    best_permutation (int) Number of the best permutation. 
        Access data by: results_dict[best_permutation]
    
    """
    n_genes = len(genes)
    try:
        print('With your gene list of {} genes, there are {:.2E} possible premutations, you are trying: {} of them.'.format(n_genes, factorial(n_genes), trials))
    except OverflowError:
        print('With your gene list of {} genes, there are >10E256 possible premutations, you are trying: {} of them.'.format(n_genes, trials))
    
    results_dict = {}

    for i in range(trials):
        results_dict[i] = {}

        #make dataframe to store the expression level for a round
        results = pd.DataFrame(np.zeros((cycles, df_exp.shape[1])), columns = df_exp.columns)

        #Shuffle the genes randomly and make it into an array
        shuffle_gene = np.array(random.sample(genes, len(genes)))
        results_dict[i]['gene_order'] = shuffle_gene

        #cycle over the cycles that have different genes in them
        for cycle in range(cycles):

            #genes that are simultaneously labeled in this round
            positive_genes = np.array(shuffle_gene)[codebook[:,cycle] == 1]

            #Get the sum of the expression
            sum_exp = df_exp.loc[positive_genes].sum()

            #Add it to the results df
            results.loc[cycle] = sum_exp

        #Add results to results_dict
        results_dict[i]['sum_counts'] = results

        #Find the maxima of all cycles and take the maximum of that
        maxima = results.max(axis=1).max()
        #Add the cycle maxima (this is the maximum summed expression in a single cell type of all genes labeled in this round)
        results_dict[i]['cycle_max'] = results.max(axis=1)

        results_dict[i]['max'] = results.max(axis=1).max()
    
    #Select permutation with lowest max expression
    maxima_data = {results_dict[i]['max'] : i for i in results_dict.keys()}
    best_max = min(maxima_data.keys())
    best_permutation = maxima_data[best_max]
    best_gene_order = results_dict[best_permutation]['gene_order']
    print('Found a gene order where the max expression over all cell types and cycles is: {} '.format(best_max))
    
    if plot == True:
        fig = plt.figure(constrained_layout=True, figsize=(15,5))
        gs = fig.add_gridspec(1, 5)
        ax1 = fig.add_subplot(gs[:, 0])
        ax2 = fig.add_subplot(gs[:, 1:])


        ax1.boxplot(maxima_data.keys())
        ax1.scatter(1, best_max)
        ax1.set_title('Iteration results')
        ax1.set_ylabel('Max expression in dataset')
        ax1.set_xlabel('Blue dot is chosen permutation')

        #find lowest max expression
        ax2.bar(np.arange(1,16,1), results_dict[best_permutation]['cycle_max'], color='grey')
        ax2.boxplot(results_dict[best_permutation]['sum_counts'])
        ax2.set_title('Best permutation {}: Max count over cycles'.format(best_permutation))
        ax2.set_ylabel('Max expression over cell types')
        ax2.set_xlabel('Barcoding cycle')
    
    return best_gene_order, results_dict, best_permutation