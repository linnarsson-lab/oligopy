def select_genes(df, min_exp, max_exp, bad_genes):
    """
    Select genes from cytograph cluster aggregate.
    
    Args:
    df (dataframe): Dataframe containing the clusterNames and their MarkerGenes
    min_exp (int): Minimal required expression.
    max_exp (int): Maximal allowed expression.
    bad_genes (list): List of genes that are not allowed. Use if the probe
        sequence generation program has difficulty designing probes for 
        the gene of question. 
        
    Returns:
    genes (set): Unique list of genes.
    
    """
    genes = {}
    #make max expresion per gene dic
    max_exp_dict = dict(zip(ds.ra['Gene'], ds[:,:].max(axis=1)))

    for i in df.index:
        CN = df.loc[i, 'ClusterName']
        ok=False
        list_expression = np.array([])
        list_genes = []

        for j in range(5):
            #Get gene name
            g = df.loc[i, 'MarkerGenes'].split(' ')[j]
            #Get expression in target cluster
            clust_exp = ds[np.where(ds.ra['Gene'] == g)[0][0], np.where(ds.ca['ClusterName'] == CN)[0][0]]
            if g not in bad_genes:
                list_expression = np.append(list_expression, clust_exp)
                list_genes.append(g)
            #Get max expression of all clusters
            max_exp_all = max_exp_dict[g]
            #Compare minimal expression with the cluster expression
            #Compare max expression with max expression for all clusters
            if min_exp < clust_exp and max_exp_all < max_exp and g not in bad_genes:
                ok = True
                break

        if ok == False:
            #Selection failed, pick gene that best matches the criteria.
            #Get gene that is closest to the middle of the min and max expression.
            diff = np.absolute(list_expression - ((max_exp - min_exp) /2))
            index_best = np.where(diff == diff.min())
            #Select best gene
            g = list_genes[index_best[0][0]]
            print('No marker for {}, best: {} with expression: {}, index {} --> {}'.format(CN, g, ds[np.where(ds.ra['Gene'] == g)[0][0], np.where(ds.ca['ClusterName'] == CN)[0][0]], index_best[0][0], np.round(list_expression,3)))

        #Add to gene set
        genes[CN] = g

    print('\nNumber of unique selected genes: {}  {} options left'.format(len(set(genes.values())), 168-len(set(genes.values()))))    
    return genes
