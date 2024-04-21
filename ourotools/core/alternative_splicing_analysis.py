# # %% alternative_splicing_analysis %%
# # import internal modules

"""
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
This part should be uncommented in core.py
"""

from . import biobookshelf as bk

"""
||||||||||||||||||||||||||||||||
"""

# from ourotools import bk
# bk.Wide( 100 )

"""
This part should be uncommented in jupyter notebook
<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
"""

from typing import Union, List, Literal, Dict
import os
import pandas as pd
import numpy as np
import multiprocessing as mp
import math
import logging
from copy import copy, deepcopy
import pickle
import time
import glob
import gzip  # to handle gzip file
import scanpy as sc
import anndata as ad
import scipy

''' useful variables '''
dict_kwargs_filter_ft_for_retrieving_isoform_only = {
    'flag_exclude_exon' : True,
    'flag_exclude_sj' : True,
    'flag_exclude_genomic_bin' : True,
    'flag_exclude_repeat_element' : True,
    'flag_exclude_reg_element' : True,
    'flag_exclude_not_uniquely_assigned_to_tx' : True,
    'flag_exclude_intron_retention' : True,
    'flag_exclude_realigned' : True,
    'flag_exclude_flag_gene_strand_specific' : True,
}
dict_kwargs_filter_ft_for_retrieving_isoform_exon_sj_only = {
    'flag_exclude_exon' : False,
    'flag_exclude_sj' : False,
    'flag_exclude_genomic_bin' : True,
    'flag_exclude_repeat_element' : True,
    'flag_exclude_reg_element' : True,
    'flag_exclude_not_uniquely_assigned_to_tx' : True,
    'flag_exclude_intron_retention' : True,
    'flag_exclude_realigned' : True,
    'flag_exclude_flag_gene_strand_specific' : True,
}

''' ourotools pipeline-specific functions '''
def exclude_or_include_internal_polyA_primed( adata, flag_include_internal_polyA_primed : bool = False, inplace : bool = True ) :
    ''' 
    exclude internal polyA primed reads during counting or include all reads for counting
    
    # 2023-12-08 21:17:36 
    '''
    if not inplace :
        adata = adata.copy( )
    
    suffix_excluding_internal_polyA_primed = '|excluding_internal_polyA_primed'
    if flag_include_internal_polyA_primed :
        adata_processed = adata[ :, bk.Search_list_of_strings_with_multiple_query( adata.var.index.values, '-' + suffix_excluding_internal_polyA_primed ) ]
    else :
        adata_processed = adata[ :, bk.Search_list_of_strings_with_multiple_query( adata.var.index.values, suffix_excluding_internal_polyA_primed ) ]
        arr_ft = list( e.rsplit( suffix_excluding_internal_polyA_primed, 1 )[ 0 ] for e in adata_processed.var.index.values ) # discard 'suffix_excluding_internal_polyA_primed' suffix from the feature names
        adata_processed.var_names = arr_ft # update var names
        adata_processed.var.index = arr_ft
    return adata_processed
def use_all_reads_or_full_length_reads_only( adata, flag_include_all_reads : bool = False, inplace : bool = True ) :
    ''' 
    include only full-length raeds or include all reads for counting
    
    ===Usage Example===
    >>> adata = use_all_reads_or_full_length_reads_only( adata, flag_include_all_reads = False ) # use full-length reads only
    
    # 2024-01-08 15:06:36 
    '''
    if not inplace :
        adata = adata.copy( )
    
    suffix_including_full_length_reads_only = '|full_length_with_valid_3p_and_5p_ends'
    if flag_include_all_reads :
        adata_processed = adata[ :, bk.Search_list_of_strings_with_multiple_query( adata.var.index.values, '-' + suffix_including_full_length_reads_only ) ]
    else :
        adata_processed = adata[ :, bk.Search_list_of_strings_with_multiple_query( adata.var.index.values, suffix_including_full_length_reads_only ) ]
        arr_ft = list( e.rsplit( suffix_including_full_length_reads_only, 1 )[ 0 ] for e in adata_processed.var.index.values ) # discard 'suffix_including_full_length_reads_only' suffix from the feature names
        adata_processed.var_names = arr_ft # update var names
        adata_processed.var.index = arr_ft
    return adata_processed
def summarize_count( adata ) :
    ''' 
    retrieve quality metrics 
    
    ===Usage Example===
    >>> summarize_count( adata )

    # 2023-12-10 17:14:21 
    '''
    sc.pp.filter_genes( adata, min_cells = 0 )
    sc.pp.filter_genes( adata, min_counts = 0 )
    sc.pp.filter_cells( adata, min_genes = 0 )
    sc.pp.filter_cells( adata, min_counts = 0 )
def classify_and_sort_feature( 
    adata, 
    inplace : bool = True,
    path_folder_ref : Union[ None, str ] = None,
    flag_convert_gene_id_to_gene_name : bool = True,
    flag_assign_features_with_position_tag_to_overlapping_gene : bool = True,
) :
    """
    classify and sort features 
    an important behavior of this function is placing a gene feature before all the sub-gene features belonging to the gene, which is the required sorting order for many of the downstream analysis steps.

    inplace : bool = True,
    path_folder_ref : Union[ None, str ] = None, # path of the Ouro-Tools reference object used for quantification
    flag_convert_gene_id_to_gene_name : bool = True, # if True, convert 'gene_id' values to 'gene_name' values using the annotations stored in 'path_folder_ref'
    flag_assign_features_with_position_tag_to_overlapping_gene : bool = True, # if True, search for features with '|pos=' tag and assign these features to overlapping genes (if it was not already assigned), and its values in 'name_gene' column will be updated. If the features are overlapped with more than one genes, the feature will not be assigned.
    
    ===Usage Example===
    >>> adata = classify_and_sort_feature( adata, path_folder_ref = '/home/project/Single_Cell_Full_Length_Atlas/ourotools.index/Mus_musculus.GRCm38.102.v0.2.4/', flag_convert_gene_id_to_gene_name = True, flag_assign_features_with_position_tag_to_overlapping_gene = True )

    # 2024-03-22 14:54:03 
    """
    if not inplace :
        adata = adata.copy( ) # create a copy of the data

    """ classify features """
    # retrieve list of features
    arr_ft = adata.var.index.values 
    # strand specific
    adata.var[ 'flag_short_read' ] = bk.Search_list_of_strings_Return_mask( arr_ft, '|short_read' )
    # strand specific
    adata.var[ 'flag_strand_specific' ] = bk.Search_list_of_strings_Return_mask( arr_ft, '|strand' )
    # genomic bin mask
    adata.var[ 'flag_genomic_bin' ] = bk.Search_list_of_strings_Return_mask( arr_ft, 'genomic_region|' )
    # repeat element mask
    mask = ~ adata.var[ 'flag_genomic_bin' ].values
    mask &= ( ~ adata.var[ 'flag_short_read' ].values ) # exclude short_read features for classification
    adata.var[ 'flag_repeat_element' ] = False # initialize the column
    adata.var.loc[ mask, 'flag_repeat_element' ] = bk.Search_list_of_strings_Return_mask( arr_ft[ mask ], 'repeatmasker_ucsc|' )
    # reg. element mask
    mask &= ( ~ adata.var[ 'flag_repeat_element' ].values )
    adata.var[ 'flag_reg_element' ] = False # initialize the column
    adata.var.loc[ mask, 'flag_reg_element' ] = bk.Search_list_of_strings_Return_mask( arr_ft[ mask ], 'regulatory_element|' )
    # gene 
    mask &= ( ~ adata.var[ 'flag_reg_element' ].values )
    adata.var[ 'flag_gene' ] = False # initialize the column
    adata.var.loc[ mask, 'flag_gene' ] = bk.Search_list_of_strings_Return_mask( arr_ft[ mask ], '|', is_negative_query = True )
    # sub-gene
    mask &= ( ~ adata.var[ 'flag_gene' ].values )
    adata.var[ 'flag_sub_gene' ] = mask 

    """ pre-process sub-gene feature names """
    mask_subgene = mask
    # convert gene_id in sub-gene feature name to gene_name
    if flag_convert_gene_id_to_gene_name and isinstance( path_folder_ref, str ) and os.path.exists( f'{path_folder_ref}gtf.gene.genome.gtf' ) :
        df_gtf_gene = bk.GTF_Read( f'{path_folder_ref}gtf.gene.genome.gtf' )
        # retrieve the mapping
        dict_gene_id_to_gene_name = df_gtf_gene.set_index( 'gene_id' ).gene_name.to_dict( )
        # retrieve new ft names, where id_gene is replaced by name_gene
        arr_ft_new = list( ( dict_gene_id_to_gene_name[ name_ft.split( '|', 1 )[ 0 ] ] + '|' + name_ft.split( '|', 1 )[ 1 ] if name_ft.split( '|', 1 )[ 0 ] in dict_gene_id_to_gene_name else name_ft ) if flag_sub_gene else name_ft for name_ft, flag_sub_gene in zip( adata.var_names, mask_subgene ) ) 
        # set new feature names
        adata.var_names = arr_ft_new
        adata.var.index = arr_ft_new
        del arr_ft_new

    # retrieve list of features again
    arr_ft = adata.var.index.values 

    ''' set name_gene for sub-gene and gene features '''
    # retrieve name of the gene for sub-gene features
    adata.var.loc[ mask_subgene, 'name_gene' ] = list( e.split( '|', 1 )[ 0 ] for e in arr_ft[ mask_subgene ] )
    # retrieve name of the gene for gene features
    adata.var.loc[ adata.var[ 'flag_gene' ].values, 'name_gene' ] = arr_ft[ adata.var[ 'flag_gene' ].values ]

    """ further classify sub-gene features """
    # intron retention
    adata.var[ 'flag_intron_retention' ] = False # initialize the column
    adata.var.loc[ mask, 'flag_intron_retention' ] = bk.Search_list_of_strings_Return_mask( arr_ft[ mask ], '|intron_retention@splice_donor_and_acceptor_id' )
    # realigned - shared by SJ and Exon
    mask &= ( ~ adata.var[ 'flag_intron_retention' ].values )
    adata.var[ 'flag_realigned' ] = False # initialize the column
    adata.var.loc[ mask, 'flag_realigned' ] = bk.Search_list_of_strings_Return_mask( arr_ft[ mask ], '|realigned' )
    # SJ
    adata.var[ 'flag_sj' ] = False # initialize the column
    adata.var.loc[ mask, 'flag_sj' ] = bk.Search_list_of_strings_Return_mask( arr_ft[ mask ], '|sj_id' )
    # exon
    mask &= ( ~ adata.var[ 'flag_sj' ].values )
    adata.var[ 'flag_exon' ] = False # initialize the column
    adata.var.loc[ mask, 'flag_exon' ] = bk.Search_list_of_strings_Return_mask( arr_ft[ mask ], '|exon_id' )
    # transcript (gene isoform)
    mask &= ( ~ adata.var[ 'flag_exon' ].values )
    adata.var[ 'flag_tx' ] = False # initialize the column
    adata.var.loc[ mask, 'flag_tx' ] = bk.Search_list_of_strings_Return_mask( arr_ft[ mask ], '|tx_id' )
    # not uniquely assigned to tx (ambiguous reads not assigned to transcripts)
    mask &= ( ~ adata.var[ 'flag_tx' ].values )
    adata.var[ 'flag_not_uniquely_assigned_to_tx' ] = False # initialize the column
    adata.var.loc[ mask, 'flag_not_uniquely_assigned_to_tx' ] = bk.Search_list_of_strings_Return_mask( arr_ft[ mask ], '|not_uniquely_assigned_to_tx' )
    # gene (strand specific)
    mask &= ( ~ adata.var[ 'flag_not_uniquely_assigned_to_tx' ].values )
    adata.var[ 'flag_gene_strand_specific' ] = mask 

    """ assign features with the |pos= tag to genes """
    if flag_assign_features_with_position_tag_to_overlapping_gene :
        # retrieve mask for the entries to be searched
        mask_no_gene_name = pd.isnull( adata.var.name_gene ).values # retrieve a mask for features with no gene names
        mask_no_gene_name_with_pos_tag = np.zeros( len( arr_ft ), dtype = bool ) # initialize a mask with no gene name but with position tag # initialize with False
        mask_no_gene_name_with_pos_tag[ mask_no_gene_name ] = bk.Search_list_of_strings_Return_mask( arr_ft[ mask_no_gene_name ], '|pos=' )

        # retrieve interval tree of gene annotations
        dict_it_gene = bk.GTF_Interval_Tree( f'{path_folder_ref}gtf.gene.genome.gtf', feature = [ 'gene' ], value = 'gene_name' )

        ''' perform gene annotation search and assign 'name_gene' to features with |pos= tag '''
        int_num_entries_to_search = mask_no_gene_name_with_pos_tag.sum( ) # retrieve the number of feature entries to search
        dict_pos_tag_content_to_name_gene_assigned = dict( ) # store search result to avoid redundant searching
        arr_name_gene_assigned = np.full( int_num_entries_to_search, np.nan, dtype = object ) # initialize the array
        for idx, name_ft in enumerate( arr_ft[ mask_no_gene_name_with_pos_tag ] ) :
            # retrieve the content of the position tag
            pos_tag_content = name_ft.split( '|pos=', 1 )[ 1 ].split( '|', 1 )[ 0 ]
            # if the search result is already available for the position tag, skip the search process and retrieve the result.
            if pos_tag_content in dict_pos_tag_content_to_name_gene_assigned : 
                arr_name_gene_assigned[ idx ] = dict_pos_tag_content_to_name_gene_assigned[ pos_tag_content ]
                continue

            # parse the content of the position tag 
            name_chr, str_st_and_en = pos_tag_content.split( ':' )
            st, en = str_st_and_en.split( '-' )
            st, en = int( st ) - 1, int( en ) # convert position to integers, and conver to 0-based coordinates
            if name_chr not in dict_it_gene : # skip if 'name_chr' is invalid
                dict_pos_tag_content_to_name_gene_assigned[ pos_tag_content ] = np.nan # store the result
                continue
            l_itv = list( dict_it_gene[ name_chr ][ st : en ] ) # retrieve list of overlapped intervals
            int_n_gene_overlapped = len( l_itv ) # retrieve the number of genes overlapped with the region described in the |pos= tag of the entry
            if int_n_gene_overlapped == 0 : # if no overlaps were found, do not assign to any gene
                dict_pos_tag_content_to_name_gene_assigned[ pos_tag_content ] = np.nan # store the result
                continue
            if int_n_gene_overlapped == 1 : # if the feature was overlapped with a single gene, assign to the gene
                name_gene_assigned = l_itv[ 0 ][ 2 ] # retrieve name_gene of the gene
                dict_pos_tag_content_to_name_gene_assigned[ pos_tag_content ] = name_gene_assigned # store the result
                arr_name_gene_assigned[ idx ] = name_gene_assigned # save the result
            else : # if the feature was overlapped with more than a single gene, do not assign to any gene
                dict_pos_tag_content_to_name_gene_assigned[ pos_tag_content ] = np.nan # store the result
                continue    
        del dict_pos_tag_content_to_name_gene_assigned # delete the dictionary

        adata.var.loc[ mask_no_gene_name_with_pos_tag, 'name_gene' ] = arr_name_gene_assigned # update the 'name_gene' column

    ''' sort data '''
    # sort features based on classifications and name_gene (for gene and sub-gene features)
    adata = adata[ :, adata.var.sort_values( [ 'name_gene', 'flag_gene', 'flag_intron_retention', 'flag_tx', 'flag_exon', 'flag_sj', 'flag_not_uniquely_assigned_to_tx', 'flag_gene_strand_specific', 'flag_realigned', 'flag_genomic_bin', 'flag_repeat_element', 'flag_reg_element', 'flag_strand_specific' ], ascending = False ).index.values ].copy( ) # sort by 'name_gene', followed by 'flag_gene' with ascending = False will ensure gene feature is placed before other features that have been assigned to genes (sub-gene features, meaning the features with high resolution than 'gene')
    return adata # return the result
def normalize( adata, flag_normalize_using_gene_count_only : bool = True, target_sum : int = 5000 ) :
    ''' 
    normalize data based on the total count of gene features 
    
    flag_normalize_using_gene_count_only : bool = True # if True, normalize using total counts in gene features only. if False, total counts of reads uniquely aligned to the genome will be used (calculated from the genomic bin features)
    target_sum : int = 5000,
    
    ===Usage Example===
    >>> normalize( adata, flag_normalize_using_gene_count_only = False, target_sum = 5000 ) # normalize data using n_counts from the entire genome

    # 2023-12-10 17:14:02 
    '''
    # calculate gene-only counts
    adata_gene = adata[ :, adata.var[ 'flag_gene' ].values ] 
    sc.pp.filter_cells( adata_gene, min_counts = 0 )
    adata.obs[ 'n_counts_gene_only' ] = adata_gene.obs[ 'n_counts' ] # retrieve total count values of the genes
    bk.MPL_1D_Sort_Plot( adata.obs.n_counts_gene_only, y_scale ='log' ) # plot distribution of UMI count values
    
    # calculate genome-aligned reads count
    adata_genomic_region = adata[ :, ( adata.var[ 'flag_genomic_bin' ].values & ( ~ adata.var[ 'flag_strand_specific' ].values ) ) ]
    sc.pp.filter_cells( adata_genomic_region, min_counts = 0 )
    adata.obs[ 'n_counts_entire_genome' ] = adata_genomic_region.obs[ 'n_counts' ] # retrieve total count values of the genes
    bk.MPL_1D_Sort_Plot( adata.obs.n_counts_entire_genome, y_scale ='log' ) # plot distribution of UMI count values
    
    # normalize data
    adata.X = scipy.sparse.csr_matrix( ( ( adata.X.T / adata.obs[ 'n_counts_gene_only' if flag_normalize_using_gene_count_only else 'n_counts_entire_genome' ].values ) * target_sum ).T ) # convert to CSR matrix
    return adata_gene, adata_genomic_region
def filter_ft( 
    df_var,
    flag_exclude_gene : bool = False,
    flag_exclude_exon : bool = False,
    flag_exclude_sj : bool = False,
    flag_exclude_tx : bool = False,
    flag_exclude_strand_specific : bool = False,
    flag_exclude_not_uniquely_assigned_to_tx : bool = False,
    flag_exclude_intron_retention : bool = False,
    flag_exclude_realigned : bool = False,
    flag_exclude_flag_gene_strand_specific : bool = False,
    flag_exclude_genomic_bin : bool = False,
    flag_exclude_repeat_element : bool = False,
    flag_exclude_reg_element : bool = False,
    flag_exclude_feature_not_assigned_to_gene : bool = False,
) :
    """ 
    filter features of a given AnnData object or the AnnData.var DataFrame object 
    (Note: for AnnData object, a view will be returned)
    
    df_var # AnnData.var DataFrame or AnnData object itself
    flag_exclude_gene : bool = False,
    flag_exclude_exon : bool = False,
    flag_exclude_sj : bool = False,
    flag_exclude_tx : bool = False,
    flag_exclude_strand_specific : bool = False,
    flag_exclude_not_uniquely_assigned_to_tx : bool = False,
    flag_exclude_intron_retention : bool = False,
    flag_exclude_realigned : bool = False,
    flag_exclude_flag_gene_strand_specific : bool = False,
    flag_exclude_genomic_bin : bool = False,
    flag_exclude_repeat_element : bool = False,
    flag_exclude_reg_element : bool = False,
    flag_exclude_feature_not_assigned_to_gene : bool = False,
    
    ===Usage Example===
    >>> adata = filter_ft( adata, flag_exclude_strand_specific = True, flag_exclude_genomic_bin = True ) # filter features
    
    # 2023-12-14 22:46:39 
    """
    ''' retrieve AnnData.var DataFrame '''
    flag_anndata_input = isinstance( df_var, sc.AnnData )
    if flag_anndata_input :
        adata = df_var
        df_var = adata.var
        
    if not isinstance( df_var, pd.DataFrame ) :
        RuntimeError( f'given df_var is not a pandas.DataFrame object' )
        return
    
    ''' Perform filtering - using columns containing boolean data '''
    mask = np.ones( len( df_var ), dtype = bool ) # initialize the mask 
    for name_col, flag_val in zip( [ 'flag_gene', 'flag_exon', 'flag_sj', 'flag_tx', 'flag_strand_specific', 'flag_not_uniquely_assigned_to_tx', 'flag_intron_retention', 'flag_realigned', 'flag_gene_strand_specific', 'flag_genomic_bin', 'flag_repeat_element', 'flag_reg_element' ], [ flag_exclude_gene, flag_exclude_exon, flag_exclude_sj, flag_exclude_tx, flag_exclude_strand_specific, flag_exclude_not_uniquely_assigned_to_tx, flag_exclude_intron_retention, flag_exclude_realigned, flag_exclude_flag_gene_strand_specific, flag_exclude_genomic_bin, flag_exclude_repeat_element, flag_exclude_reg_element ] ) : 
        if not flag_val : # include the features
            continue
        mask &= ~ df_var[ name_col ].values # exclude the features
    df_var = df_var[ mask ] # perform filtering operation
    
    ''' Perform filtering - using columns containing non-boolean data '''
    if flag_exclude_feature_not_assigned_to_gene :
        mask = ~ pd.isnull( df_var[ 'name_gene' ].values )
        df_var = df_var[ mask ] # perform filtering operation

    ''' return output '''
    # return filtered features
    if flag_anndata_input :
        return adata[ :, df_var.index.values ] # does not return a copy
    else :
        return df_var 
def search_gene( 
    adata, 
    name_gene : str, 
    ** kwargs
) :
    """ 
    search features and return a subset of 'adata.var' of the features matching the given gene name and the filter options
    
    name_gene : str, 
    ** kwargs # will be given to 'filter_ft'

    ===Usage Example===
    >>> df_search_res = search_gene( adata, 'Rps24', flag_exclude_intron_retention = False )
    
    # 2023-12-10 17:44:07 
    """
    df_var = adata.var # retrieve 'df_var'
    df_var = bk.PD_Select( df_var, name_gene = name_gene ) # retrieve features of a specific gene_name
    return filter_ft( df_var, ** kwargs )
def find_markers_at_sub_gene_level( 
    adata, 
    name_col_clus : str,
    l_name_clus: Union[list, np.ndarray, tuple, set, None] = None,
    method_pval: Union[ None, Literal[ "wilcoxon", "t-test", "mann-whitney-u" ] ] = "wilcoxon",
    flag_calculate_auroc : bool = True,
    float_min_gene_count_for_each_pseudobulk : float = 30,
    random_seed : Union[ int, None ] = 42,
    int_num_cpus : int = 32,
    name_col_var_batch_weight : str = 'n_cells',
    int_max_num_genes_in_a_batch : int = 100,
    float_min_total_weight_for_a_batch : float = 1,
    float_max_total_weight_for_a_batch : float = 1000,
    dict_kwargs_filter_ft : dict = {
        'flag_exclude_exon' : True,
        'flag_exclude_sj' : True,
        'flag_exclude_genomic_bin' : True,
        'flag_exclude_repeat_element' : True,
        'flag_exclude_reg_element' : True,
        'flag_exclude_not_uniquely_assigned_to_tx' : True,
        'flag_exclude_intron_retention' : True,
        'flag_exclude_realigned' : True,
        'flag_exclude_flag_gene_strand_specific' : True,
    },
) :
    """
    find marker features at sub-gene level (isoforms, exons, splice junctions, etc.) for each cluster label by calculating a AUROC metric, log2FC, and Wilcoxon (or alternatively, t-test or Mann-Whitney-U rank test).


    name_col_clus : str = name_col_clus
    l_name_clus: Union[list, np.ndarray, tuple, set, None] =  None
    method_pval: Union[ None, Literal[ "wilcoxon", "t-test", "mann-whitney-u" ] ] = "wilcoxon"
    flag_calculate_auroc : bool = True
    float_min_gene_count_for_each_pseudobulk : float = 10
    random_seed : Union[ str, None ] = 42
    int_num_cpus : int = 16
    name_col_var_batch_weight : str = 'n_cells' # weight should be zero or positive values, and should not be negative values.
    int_max_num_genes_in_a_batch : int = 100,
    float_min_total_weight_for_a_batch : float = 1,
    float_max_total_weight_for_a_batch : float = 1000,
    dict_kwargs_filter_ft : dict = {
        'flag_exclude_exon' : True,
        'flag_exclude_sj' : True,
        'flag_exclude_genomic_bin' : True,
        'flag_exclude_repeat_element' : True,
        'flag_exclude_reg_element' : True,
        'flag_exclude_not_uniquely_assigned_to_tx' : True,
        'flag_exclude_intron_retention' : True,
        'flag_exclude_realigned' : True,
        'flag_exclude_flag_gene_strand_specific' : True,
    } # kwargs that will be given to 'filter_ft'.
    
    ===Usage Example===
    >>> df_res = find_markers_at_sub_gene_level( adata, name_col_clus )

    # 2024-03-05 16:27:59 
    """
    from tqdm import tqdm as progress_bar  # for progress bar
    from sklearn.metrics import roc_auc_score
    import scipy.stats

    ''' retrieve function for testing p-value '''
    test = None
    flag_calculate_pval = method_pval is not None
    if flag_calculate_pval :
        if method_pval not in {"wilcoxon", "t-test", "mann-whitney-u"}:
            raise RuntimeError( f"'method_pval' {method_pval} is invalid, exiting" )
        if method_pval == "t-test":
            test = scipy.stats.ttest_ind
        elif method_pval == "wilcoxon":
            test = scipy.stats.ranksums
        elif method_pval == "mann-whitney-u":
            test = scipy.stats.mannwhitneyu

    ''' filter features '''
    dict_kwargs_filter_ft[ 'flag_exclude_gene' ] = False # gene feature should be included
    dict_kwargs_filter_ft[ 'flag_exclude_feature_not_assigned_to_gene' ] = True # features that were not assigned to gene should be excluded
    adata = filter_ft( adata, ** dict_kwargs_filter_ft ) # filter features

    ''' randomly shuffle the barcodes before making pseudobulk samples '''
    arr_idx_bc = np.arange( len( adata.obs ) )
    np.random.seed( random_seed ) # initialize the random function
    np.random.shuffle( arr_idx_bc )
    adata = adata[ arr_idx_bc ].copy( ) # randomly shuffle the barcodes (to avoid barcodes from sample samples being groupped into the same pseudobulk samples)
    
    ''' drop barcodes with np.nan values as name_clus '''
    adata = adata[ ~ pd.isnull( adata.obs[ name_col_clus ] ) ].copy( ) # randomly shuffle the barcodes (to avoid barcodes from sample samples being groupped into the same pseudobulk samples)
    
    # retrieve cluster labels
    col = adata.obs[ name_col_clus ]    
    arr_clus = col.values
    l_name_clus_all = list( sorted( e for e in col.values.unique( ) if isinstance( e, str ) ) ) # retrieve all unique clus labels # exclude float values
    # retrieve name_clus -> index mapping
    dict_name_clus_to_idx_clus = dict( ( e, i ) for i, e in enumerate( l_name_clus_all ) )

    ''' subset adata to include only specific cluster labels '''
    if l_name_clus is not None :
        set_name_clus = set( l_name_clus_all ).intersection( l_name_clus )
        adata = adata[ list( e in set_name_clus for e in arr_clus ) ].copy( ) # subset the adata
        # retrieve cluster labels again
        col = adata.obs[ name_col_clus ]
        arr_clus = col.values

    ''' perform the analysis '''
    arr_name_gene = adata.var.name_gene.values # retrieve list of name_gene 
    dict_name_clus_to_int_num_cells = bk.COUNTER( arr_clus ) # count the number of cells for each cluster label
    def _analyze_sub_gene_features_for_each_gene( p_i, p_o ) :
        '''
        # 2023-12-10 19:31:57 
        '''
        while True :
            ins = p_i.recv( )
            if ins is None : # exit the loop, since the work list is empty
                break
            l_name_gene = ins # parse input
            dict_res = dict( ) # initialize the output 
            int_num_gene_analyzed = len( l_name_gene ) # retrieve the number of genes analyzed (will be used for progress bar status update)
            for name_gene in l_name_gene : # for each gene
                # subset the AnnData
                adata_subset = adata[ :, arr_name_gene == name_gene ]

                arr_X = adata_subset.X.toarray( ) # retrieve expression count
                arr_gene = arr_X[ :, 0 ] # retrieve gene counts

                '''
                calculate the proportion of cells expressing the gene for each cluster label
                '''
                # count the number of cells expressing the gene for each cluster label
                dict_name_clus_to_prop_expression = dict( )
                for name_clus in arr_clus[ arr_gene > 0 ] : # retrieve cluster names of cells expressing the gene
                    if name_clus not in dict_name_clus_to_prop_expression :
                        dict_name_clus_to_prop_expression[ name_clus ] = 1
                    else :
                        dict_name_clus_to_prop_expression[ name_clus ] += 1
                # calculate the proportions
                for name_clus in dict_name_clus_to_int_num_cells :
                    if name_clus not in dict_name_clus_to_prop_expression :
                        dict_name_clus_to_prop_expression[ name_clus ] = 0
                    else :
                        dict_name_clus_to_prop_expression[ name_clus ] = dict_name_clus_to_prop_expression[ name_clus ] / dict_name_clus_to_int_num_cells[ name_clus ] # calculate the proportion

                '''
                assign cells to pseudobulk samples
                '''
                dict_name_clus_to_pseudobulks = dict( ) # initialize a dictionary that will contain pseudobulk assignment results
                dict_name_clus_to_pseudobulk_under_construction = dict( )
                def _add_pseudobulk( name_clus, l_idx ) : 
                    if name_clus not in dict_name_clus_to_pseudobulks :
                        dict_name_clus_to_pseudobulks[ name_clus ] = [ ]
                    dict_name_clus_to_pseudobulks[ name_clus ].append( l_idx )
                def _update_pseudobulk_under_construction( name_clus, idx, count ) :
                    if name_clus not in dict_name_clus_to_pseudobulk_under_construction :
                        dict_name_clus_to_pseudobulk_under_construction[ name_clus ] = {
                            'total_count' : 0,
                            'l_idx' : [ ],
                        }
                    dict_data = dict_name_clus_to_pseudobulk_under_construction[ name_clus ]
                    dict_data[ 'total_count' ] += count
                    dict_data[ 'l_idx' ].append( idx )
                    if dict_data[ 'total_count' ] >= float_min_gene_count_for_each_pseudobulk : # if the total count meets the pseudobulk count threshold, add the current batch as a pseudobulk
                        _add_pseudobulk( name_clus, dict_data[ 'l_idx' ] )
                        # reset the container
                        dict_data[ 'total_count' ] = 0
                        dict_data[ 'l_idx' ] = [ ]
                def _flush_pseudobulk_under_construction( ) :
                    for name_clus in dict_name_clus_to_pseudobulk_under_construction :
                        dict_data = dict_name_clus_to_pseudobulk_under_construction[ name_clus ]
                        if len( dict_data[ 'l_idx' ] ) > 0 : # if the container is not empty
                            if name_clus in dict_name_clus_to_pseudobulks : # if a pseudobulk for the name_clus exists
                                dict_name_clus_to_pseudobulks[ name_clus ][ -1 ].extend( dict_data[ 'l_idx' ] ) # add the indices to the one of the pseudobulk
                                # empty the container
                                dict_data[ 'l_idx' ] = [ ] 
                                dict_data[ 'total_count' ] = 0

                for i in range( len( arr_clus ) ) :
                    name_clus, count = arr_clus[ i ], arr_gene[ i ] # retrieve name_clus and count
                    if count >= float_min_gene_count_for_each_pseudobulk : # if the single cell alone exceed the pseudobulk threshold, convert it to a pseudobulk sample
                        _add_pseudobulk( name_clus, [ i ] )
                    else :
                        _update_pseudobulk_under_construction( name_clus, i, count )
                _flush_pseudobulk_under_construction( )

                '''
                construct pseudobulk count matrix
                '''
                # retrieve list of name_clus in pseudobulks
                dict_name_clus_to_int_num_pseudobulks = dict( ( name_clus, len( dict_name_clus_to_pseudobulks[ name_clus ] ) ) for name_clus in dict_name_clus_to_pseudobulks ) # retrieve name_clus -> num pseudobulk samples for samples with pseudobulks
                l_name_clus_in_pseudobulks = list( dict_name_clus_to_int_num_pseudobulks ) # retrieve list of name_clus in pseudobulks
                int_num_name_clus_in_pseudobulks = len( l_name_clus_in_pseudobulks )

                ''' %% handle when no pseudobulk samples were generated %%  '''
                if int_num_name_clus_in_pseudobulks <= 1 : # if one or less cluster labels have pseudobulk samples, skip the analysis of the gene (since at least two cluster labels are required for the analysis)
                    continue

                # compose 'arr_cluster_label'
                l_cluster_label = list( )
                for idx_clus_pseudobulks, name_clus in enumerate( l_name_clus_in_pseudobulks ) :
                    l_cluster_label.extend( [ idx_clus_pseudobulks ] * dict_name_clus_to_int_num_pseudobulks[ name_clus ] )
                arr_cluster_label_pseudobulks = np.array( l_cluster_label, dtype = int )

                # compose pseudobulk count matrix
                int_num_pseudobulks = len( arr_cluster_label_pseudobulks ) # count the number of pseudobulk samples
                arr_X_pseudobulks = np.zeros( ( int_num_pseudobulks, arr_X.shape[ 1 ] ) ) # initialize the count matrix for pseudobulks
                idx_pseudobulk = 0
                for name_clus in l_name_clus_in_pseudobulks:
                    for l_idx in dict_name_clus_to_pseudobulks[ name_clus ] :
                        arr_X_pseudobulks[ idx_pseudobulk, : ] = arr_X[ l_idx ].sum( axis = 0 )
                        idx_pseudobulk += 1

                '''
                calculate relative ratio to gene
                '''
                # calculate proportion to gene
                arr_X_pseudobulks_gene_count = arr_X_pseudobulks[ :, 0 ] # retrieve gene count of pseudobulk samples
                arr_X_pseudobulks_porportion_to_gene = ( arr_X_pseudobulks[ :, 1 : ].T / arr_X_pseudobulks_gene_count ).T # calculate proportion to gene
                l_name_ft_subgene = adata_subset.var.index.values[ 1 : ] # retrieve subgene feature names (excluding the gene, which should be the first feature name)

                ''' survey the average gene counts in pseudobulk samples for each cluster label (to identify possible issues arising from sequencing coverage differences between pseudobulk samples) '''
                arr_avg_gene_count_in_pseudobulks = np.zeros( int_num_name_clus_in_pseudobulks, dtype = float )
                for idx_clus_pseudobulks, name_clus in enumerate( l_name_clus_in_pseudobulks ) :
                    # retrieve the average gene expression value of the pseudobulk samples of the cluster 
                    arr_avg_gene_count_in_pseudobulks[ idx_clus_pseudobulks ] = arr_X_pseudobulks_gene_count[arr_cluster_label_pseudobulks == idx_clus_pseudobulks].mean( )

                '''
                perform diagnostical/statistical tests (Log2FC / AUROC / p-value)
                '''
                # initialize the output values
                arr_output = np.full( ( 4, len( l_name_ft_subgene ), len( l_name_clus_in_pseudobulks ) ), np.nan, dtype = float ) # [ 'log2fc', 'auroc', 'pval' ]
                # for each sub-gene feature
                for idx_ft, arr_prop_to_gene in enumerate( arr_X_pseudobulks_porportion_to_gene.T.copy( ) ) :
                    name_ft = l_name_ft_subgene[ idx_ft ] # retrieve name of the feature
                    # for each cluster
                    for idx_clus_pseudobulks, name_clus in enumerate( l_name_clus_in_pseudobulks ) :
                        # retrieve relative expression values of cluster and the rest of the barcodes
                        mask = arr_cluster_label_pseudobulks == idx_clus_pseudobulks
                        arr_prop_to_gene_clus = arr_prop_to_gene[mask]
                        arr_prop_to_gene_rest = arr_prop_to_gene[~mask]

                        # calculate log2fc values
                        mean_arr_prop_to_gene_rest = arr_prop_to_gene_rest.mean()
                        mean_arr_prop_to_gene_clus = arr_prop_to_gene_clus.mean()
                        arr_output[ 0, idx_ft, idx_clus_pseudobulks ] = mean_arr_prop_to_gene_clus # collect the average proportion of the cluster
                        if mean_arr_prop_to_gene_rest != 0:
                            try:
                                arr_output[ 1, idx_ft, idx_clus_pseudobulks ] = math.log2( mean_arr_prop_to_gene_clus / mean_arr_prop_to_gene_rest )
                            except ValueError:  # catch math.log2 domain error
                                pass

                        # calculate auroc
                        if flag_calculate_auroc:
                            arr_output[ 2, idx_ft, idx_clus_pseudobulks ] = roc_auc_score( mask, arr_prop_to_gene )

                        # calculate p-value after t-test or similar tests
                        if flag_calculate_pval :
                            arr_output[ 3, idx_ft, idx_clus_pseudobulks ] = test( arr_prop_to_gene_clus, arr_prop_to_gene_rest ).pvalue

                dict_res[ name_gene ] = ( l_name_ft_subgene, l_name_clus_in_pseudobulks, arr_output, dict_name_clus_to_prop_expression, dict_name_clus_to_int_num_pseudobulks, arr_avg_gene_count_in_pseudobulks ) # collect the output 
            res = ( dict_res, int_num_gene_analyzed ) # compose result
            p_o.send( res ) # send result
        p_o.send( 'completed' ) # notify the all works have been completed.

    def _create_batch( ) :
        # retrieve name_gene - weight series
        s_weight = adata.var[ [ 'name_gene', name_col_var_batch_weight ] ].rename( columns = { name_col_var_batch_weight : 'weight' } ).groupby( 'name_gene' ).sum( ).weight # .sort_values( ascending = False )
        s_weight = s_weight[ s_weight >= float_min_total_weight_for_a_batch ] # filter entries with low weights

        # create batchs
        l_batch = [ ] # initialize list of batches
        batch_under_construction = [ ]
        total_weight = 0
        for name_gene, weight in zip( s_weight.index.values, s_weight.values ) :
            if weight >= float_max_total_weight_for_a_batch :
                l_batch.append( [ name_gene ] )
            else :
                # add the current gene to the batch
                batch_under_construction.append( name_gene )
                total_weight += weight
                if total_weight >= float_max_total_weight_for_a_batch or len( batch_under_construction ) >= int_max_num_genes_in_a_batch : # make a decision to flush a batch
                    l_batch.append( batch_under_construction )
                    # reset the container
                    batch_under_construction = [ ]
                    total_weight = 0
        # flush remaining batch
        if len( batch_under_construction ) > 0 : 
            l_batch.append( batch_under_construction )
        return l_batch

    l_batch = _create_batch( ) # compose batches
    int_num_genes_analyzed = int( np.sum( list( len( e ) for e in l_batch ) ) ) # retrieve the number of genes analyzed.
    pbar = progress_bar( desc = f"analyzing sub-gene features of {int_num_genes_analyzed} genes", total = int_num_genes_analyzed, )  # initialize the progress bar
    dict_res_all = dict( ) # initialize the output container
    def _collect_res( res ) :
        dict_res, int_num_gene_analyzed = res # parse result
        pbar.update( int_num_gene_analyzed )  # update the progress bar (the number of genes analyzed)
        dict_res_all.update( dict_res ) # update the container

    bk.Multiprocessing_Batch_Generator_and_Workers(
        gen_batch = l_batch.__iter__( ),
        process_batch = _analyze_sub_gene_features_for_each_gene,
        post_process_batch = _collect_res,
        int_num_threads = int_num_cpus + 1, # one for batch generation
    )
    pbar.close()  # close the progress bar

    '''
    Compose the output result as a DataFrame
    '''
    # survey the number of rows in the result
    # also, compose the name_gene column of the output matrix
    l_n_rows_res = [ ]
    l_name_gene = list( dict_res_all ) # retrieve list of name_gene in the output
    for name_gene in l_name_gene :
        n_ft, n_clus = dict_res_all[ name_gene ][ 2 ].shape[ 1 : ]
        n_rows_res = n_ft * n_clus
        l_n_rows_res.append( n_rows_res )
    int_num_rows_res = int( np.sum( l_n_rows_res ) ) # retrieve the number of rows in the result
    # compose the 'codes' of pandas categorical data column for the 'name_gene' column
    arr_gene_res = np.zeros( int_num_rows_res, dtype = 'i4' )
    st = 0
    for i, n in enumerate( l_n_rows_res ) :
        arr_gene_res[ st : st + n ] = i
        st += n

    # compose the remaining columns of the result dataframe
    l_name_ft = [ ]
    arr_ft_res = np.zeros( int_num_rows_res, dtype = 'i4' )
    arr_clus_res = np.zeros( int_num_rows_res, dtype = 'i4' )
    arr_prop_expr_res = np.zeros( int_num_rows_res, dtype = float )
    arr_avg_gene_count_res = np.zeros( int_num_rows_res, dtype = float )
    arr_avg_prop_to_gene_res = np.zeros( int_num_rows_res, dtype = float )
    arr_num_pseudobulks_res = np.zeros( int_num_rows_res, dtype = float )
    arr_log2fc_res = np.full( int_num_rows_res, np.nan, dtype = float )
    arr_auroc_res = np.full( int_num_rows_res, np.nan, dtype = float )
    arr_pval_res = np.full( int_num_rows_res, np.nan, dtype = float )
    idx_ft, st_ft, st_clus = 0, 0, 0
    for name_gene in l_name_gene :
        # parse result of a gene
        l_ft, l_clus, arr_res, dict_name_clus_to_prop_expression, dict_name_clus_to_int_num_pseudobulks, arr_avg_gene_count_in_pseudobulks = dict_res_all[ name_gene ]
        n_ft, n_clus = len( l_ft ), len( l_clus )
        # update the feature name column and calculation results
        for i, ft in enumerate( l_ft ) :
            sl = slice( st_ft, st_ft + n_clus ) # retrieve slice for the current feature
            arr_ft_res[ sl ] = idx_ft # update the column
            idx_ft += 1 # update index of the feature name
            st_ft += n_clus # update the position of the feature name column

            # update the calculation results
            arr_avg_prop_to_gene_res[ sl ] = arr_res[ 0 ][ i ]
            arr_log2fc_res[ sl ] = arr_res[ 1 ][ i ]
            arr_auroc_res[ sl ] = arr_res[ 2 ][ i ]
            arr_pval_res[ sl ] = arr_res[ 3 ][ i ]

        l_name_ft.extend( l_ft ) # collect categories for feature names column
        # update the cluster name column and 'proportion of cells expressing the gene' and 'average gene counts of pseudobulks for each cluster' columns
        n_rows_res = n_ft * n_clus
        for i, name_clus in enumerate( l_clus ) :
            sl = slice( st_clus + i, st_clus + n_rows_res, n_clus ) # retrieve the slice for the current cluster name
            arr_clus_res[ sl ] = dict_name_clus_to_idx_clus[ name_clus ] # update the cluster name column
            arr_prop_expr_res[ sl ] = dict_name_clus_to_prop_expression[ name_clus ] # update the 'proportion of cells expressing the gene' column
            arr_avg_gene_count_res[ sl ] = arr_avg_gene_count_in_pseudobulks[ i ] # update the 'average gene counts of pseudobulks for each cluster' column
            arr_num_pseudobulks_res[ sl ] = dict_name_clus_to_int_num_pseudobulks[ name_clus ] # update the number of pseudobulks column
        st_clus += n_rows_res # update the position of the cluster name column

    # compose the result dataframe
    df_res = pd.DataFrame( {
        'name_gene' : pd.Categorical.from_codes(codes=arr_gene_res, dtype=pd.CategoricalDtype(l_name_gene, ordered=True)),
        'name_feature' : pd.Categorical.from_codes(codes=arr_ft_res, dtype=pd.CategoricalDtype(l_name_ft, ordered=True)),
        'name_clus' : pd.Categorical.from_codes(codes=arr_clus_res, dtype=pd.CategoricalDtype(l_name_clus_all, ordered=True)),
        'prop_cell_expressing_gene' : arr_prop_expr_res,
        'num_pseudobulks' : arr_num_pseudobulks_res,
        'gene_pseudobulk_avg' : arr_avg_gene_count_res,
        'prop_to_gene_pseudobulk_avg' : arr_avg_prop_to_gene_res,
        'log2fc' : arr_log2fc_res,
        'auroc' : arr_auroc_res,
        'pval' : arr_pval_res,
    } )
    return df_res # return result
def combine_long_and_short_read_data( 
    adata_sr, adata_lr, path_folder_ref : str
) :
    """
    Combine long and short-read data.
    Assumes short-read data were normalized, log-transformed data, while long-read data is normalized data.
    Cell type annotations (adata_sr.obs)  will be transferred from short-read data object to long-read data object
    'obsm' will be copied to combined data object in the following order : short_read, long_read (long-read obsm overwriting short-read obsm)
    
    adata_sr, # short-read AnnData object
    adata_lr, # long-raed AnnData object
    path_folder_ref : str, # reference folder
    
    # 2024-03-22 14:59:51 
    """
    import anndata as ad
    '''
    combine short-read and long-read data
    '''
    # retrieve copy of raw short read counts (normalized) of target cells
    adata_sr_copy = adata_sr.copy( )
    # change varnames so that var names become uniuque
    adata_sr_copy.var.index = list( f"{e}|short_read" for e in adata_sr_copy.var.index.values )

    # retrieve copy of raw long read counts (normalized) of target cells
    adata_lr_copy = adata_lr.copy( )
    sc.pp.log1p(adata_lr_copy) # perform log-transformation

    # retrieve list of shared barcodes
    l_id_bc_shared = np.sort( list( set( adata_sr_copy.obs.index.values ).intersection( adata_lr_copy.obs.index.values ) ) )

    adata_sr_copy, adata_lr_copy = adata_sr_copy[ l_id_bc_shared, : ].copy( ), adata_lr_copy[ l_id_bc_shared, : ].copy( ) # align obs
    adata_combined = ad.concat( [ adata_sr_copy, adata_lr_copy ], axis = 1 ) # combine short-read and long-read data
    
    adata_combined.obs = adata_sr_copy.obs # transfer obs
    
    ''' copy obsm '''
    for a in [ adata_sr_copy, adata_lr_copy ] :
        for k in a.obsm :
            adata_combined.obsm[ k ] = a.obsm[ k ]

    adata_combined = classify_and_sort_feature( adata_combined, path_folder_ref = path_folder_ref, flag_convert_gene_id_to_gene_name = True, flag_assign_features_with_position_tag_to_overlapping_gene = True )
    return adata_combined

''' utility functions that are not specific to the ourotools pipeline '''
def identify_batch_specific_features( 
    adata,
    name_col_batch : str,
    name_col_cluster : str = None,
    l_name_cluster = None,
    min_prop_expr = 0.15,
    min_score_ratio = 3,
) :
    """
    Identify features that are consistently differently expressed in each sample for multiple clusters
    Note)
    The primary aim of this function is for identifying batch specific features in an entire dataset or in a set of clusters to improve UMAP embedding/clustering without/with the help of other batch-correction algorithms.
    Empirically, batch effect can be significantly reduced (with some loss of information) simply by excluding a set of features contributing to the batch effects, which is often sufficient for clustering analysis.
    # 2024-02-26 22:34:55 
    """
    def _filter_marker_feature( df_unique_marker_feature ) :
        return bk.PD_Threshold( df_unique_marker_feature, prop_expr_1sta = min_prop_expr, score_ratioa = min_score_ratio )

    if l_name_cluster is None :
        set_name_feature_to_exclude = set( _filter_marker_feature( search_uniquely_expressed_marker_features( adata, name_col_batch ) ).feature_name.values )
    else :
        l_l_name_feature = list( _filter_marker_feature( search_uniquely_expressed_marker_features( adata[ adata.obs[ name_col_cluster ] == name_cluster ].copy( ), name_col_batch ) ).feature_name.values for name_cluster in l_name_cluster ) # retrieve batch specific features for each cluster
        # identify batch specific features shared between all clusters
        set_name_feature_to_exclude = set( l_l_name_feature[ 0 ] ) # initialize 'set_name_feature_to_exclude'
        for l_name_feature in l_l_name_feature[ 1 : ] :
            set_name_feature_to_exclude = set_name_feature_to_exclude.intersection( l_name_feature )
    return set_name_feature_to_exclude
def summarize_expression_for_each_clus( 
    adata,
    name_col_cluster : str,
) :
    """
    summarize expression of each cluster for the given adata
    # 2024-02-26 21:01:16 
    """
    l_df = [ ]
    def _parse_array( arr ) :
        if len( arr.shape ) == 1 :
            return arr
        else :
            return arr[ 0 ]
    for name_cluster in set( adata.obs[ name_col_cluster ].values ) :
        adata_subset = adata[ adata.obs[ name_col_cluster ] == name_cluster ] # retrieve 
        arr_num_cell_with_expr = _parse_array( np.array( ( adata_subset.X > 0 ).sum( axis = 0 ) ) )
        arr_avg_expr = _parse_array( np.array( adata_subset.X.sum( axis = 0 ) ) ) / arr_num_cell_with_expr # calculate average expression in cells expressing the feature
        arr_avg_expr[ np.isnan( arr_avg_expr ) ] = 0 # when number of cells expressing is zero, set expression values as 0
        arr_prop_expr = arr_num_cell_with_expr / len( adata_subset.obs )
        _df = pd.DataFrame( { 'avg_expr' : arr_avg_expr, 'prop_expr' : arr_prop_expr } )
        _df[ name_col_cluster ] = name_cluster
        _df[ 'feature_name' ] = adata_subset.var.index.values
        l_df.append( _df )
    df_feature_expr = pd.concat( l_df )
    df_feature_expr[ 'score' ] = df_feature_expr.avg_expr * df_feature_expr.prop_expr # calculate the score
    return df_feature_expr
def search_marker_features_unique_to_a_single_cluster( 
    adata,  
    name_col_cluster : str,
    float_max_score : float = 1000,
) :
    """
    find unique markers for each clusters.
    Note)
    score is calculated as the product of the proportion-expressed and the average (log-normalized) expression values.

    name_col_cluster : str, # name of the column in the 'adata.obs' containing cluster labels 
    float_max_score : 1000, # 'infinite' score ratios will be replaced by this value

    # 2024-03-22 19:06:46 
    """
    '''
    raise an error when there is only single cluster label available
    '''
    if len( set( adata.obs[ name_col_cluster ].values ) ) == 1 : # check the number of labels, and if only a single label available, raise an error
        RuntimeError( f"only single label available for '{name_col_cluster}'" )
    
    '''
    survey proportion expressed and avg expression
    '''
    df_feature_expr = summarize_expression_for_each_clus( adata, name_col_cluster )

    '''
    identify marker features uniquely expressed in a single cluster
    '''
    l_l = [ ]
    for feature_name, _df in df_feature_expr.groupby( 'feature_name' ) :
        # retrieve values
        arr_avg_expr, arr_prop_expr, arr_str_int_cell_type_subclustered_temp, _, arr_score = _df.values.T
        # sort by score
        arr_score_argsort = arr_score.argsort( )
        arr_avg_expr, arr_prop_expr, arr_str_int_cell_type_subclustered_temp, arr_score = arr_avg_expr[ arr_score_argsort ], arr_prop_expr[ arr_score_argsort ], arr_str_int_cell_type_subclustered_temp[ arr_score_argsort ], arr_score[ arr_score_argsort ]

        avg_expr_highest, prop_expr_highest = arr_avg_expr.max( ), arr_prop_expr.max( )
        avg_expr_2nd, avg_expr_1st = arr_avg_expr[ -2 : ]
        prop_expr_2nd, prop_expr_1st = arr_prop_expr[ -2 : ]
        str_int_cell_type_subclustered_temp_2nd, str_int_cell_type_subclustered_temp_1st = arr_str_int_cell_type_subclustered_temp[ -2 : ]
        score_2nd, score_1st = arr_score[ -2 : ]
        l_l.append( [ feature_name, str_int_cell_type_subclustered_temp_1st, avg_expr_1st, prop_expr_1st, score_1st, avg_expr_2nd, prop_expr_2nd, score_2nd, avg_expr_highest, prop_expr_highest ] ) # add a record
    df_unique_marker_feature = pd.DataFrame( l_l, columns = [ 'feature_name', 'str_int_cell_type_subclustered_temp_1st', 'avg_expr_1st', 'prop_expr_1st', 'score_1st', 'avg_expr_2nd', 'prop_expr_2nd', 'score_2nd', 'avg_expr_highest', 'prop_expr_highest' ] )
    arr_score_ratio = df_unique_marker_feature.score_1st.values / df_unique_marker_feature.score_2nd.values
    arr_score_ratio[ np.isinf( arr_score_ratio ) ] = float_max_score
    df_unique_marker_feature[ 'score_ratio' ] = arr_score_ratio
    return df_unique_marker_feature