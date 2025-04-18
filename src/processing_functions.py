# Import Libraries
import os
import glob
import scanpy as sc
import logging
import datetime
import numpy as np
import scvelo as scv
import pandas as pd
import pyviper
import copy
import gc
from cellrank.kernels import CytoTRACEKernel

my_random_seed = 666
n_cells_to_subsample = 3e3

def create_h5ad_for_samples(dataset_path):
    """Process all samples and create h5ad files in a central directory if not already created"""

    # Create central directory for h5ad files
    output_dir = os.path.join(dataset_path, "samples_h5ad_files")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created h5ad directory: {output_dir}")
    
    # Find all sample directories
    sample_dirs = glob.glob(os.path.join(dataset_path, "MJ*"))
    
    for sample_dir in sample_dirs:
        sample_name = os.path.basename(sample_dir)
        h5ad_file = os.path.join(output_dir, f"{sample_name}.h5ad")
        
        # Skip if h5ad file already exists
        if os.path.exists(h5ad_file):
            print(f"h5ad file already exists for {sample_name}")
            continue
        
        # Path to 10x data
        filtered_feature_path = os.path.join(sample_dir, "filtered_feature_bc_matrix")
        
        try:
            # Read 10x data and save as h5ad
            adata = sc.read_10x_mtx(filtered_feature_path)
            adata.obs['sample'] = sample_name
            adata.write(h5ad_file)
            print(f"Created {h5ad_file}")
        except Exception as e:
            try:
                adata = sc.read_10x_mtx(sample_dir)
                adata.obs['sample'] = sample_name
                adata.write(h5ad_file)
                print(f"Created {h5ad_file}")
            except Exception as sample_dir_error:
                print(f"Error processing {sample_name}: {str(e)} AND {str(sample_dir_error)}")

def set_up_logger(logging_path, experiment_name):

    # Create timestamp for logger files
    timestamp = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")

    # Create the logging directory if it doesn't exist
    if not os.path.exists(logging_path):
        os.makedirs(logging_path)
        print(f"Created logging directory: {logging_path}")

    file = os.path.join(logging_path , f"{experiment_name}-log-{timestamp}.txt")
    with open(file, 'w') as f:
        sc.logging.print_versions()
        
    ## Setting up the logger
    # create logger
    logger = logging.getLogger('nepc-analysis')
    logger.setLevel(logging.DEBUG)
    # create formatter
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S')
    # create console handler and set level to debug
    # ch = logging.StreamHandler()
    fh = logging.FileHandler( os.path.join(logging_path, experiment_name + "-log.txt") , mode='w')
    # ch.setLevel(logging.DEBUG)
    fh.setLevel(logging.DEBUG)
    # ch.setFormatter(formatter)
    fh.setFormatter(formatter)
    # logger.addHandler(ch)
    logger.addHandler(fh)

    logger.info("-------------- [START OF ANALYSIS] --------------")
    return logger

def get_sample_h5ad_dict(samples_dataset_dir):
    """
    get all h5ad files from the samples_h5ad_files directory
    and return a dictionary mapping sample names to teh h5ad paths
    """
    # Path to the h5ad files directory
    h5ad_files_path = os.path.join(samples_dataset_dir, "samples_h5ad_files")
    
    # Create a dictionary to store sample name -> AnnData object
    sample_file_dict = {}
    
    # Find all h5ad files
    h5ad_files = glob.glob(os.path.join(h5ad_files_path, "MJ*.h5ad"))
    
    for h5ad_file in h5ad_files:
        # Extract sample name from filename
        sample_name = os.path.basename(h5ad_file).replace(".h5ad", "")
        sample_file_dict[sample_name] = h5ad_file

    print(f"sample_file_dict contains {len(sample_file_dict)} samples")
    return sample_file_dict

def load_concat_adata(sample_file_dict, samples_to_pick, logger, dataset_path, subset_name):
    """
    Load and concatenate AnnData objects from samples, then save if h5ad file does not exist
    
    Args:
        sample_file_dict: Dictionary mapping sample names to h5ad file paths
        samples_to_pick: List of sample names to include
        logger: Logger object for recording progress
        output_path: Path to save the concatenated h5ad file
        
    Returns:
        Concatenated AnnData object
    """

    output_dir = os.path.join(dataset_path, "adata_h5ad_outputs")
    os.makedirs(output_dir, exist_ok=True)

    output_path = os.path.join(output_dir, subset_name)

    if os.path.exists(output_path):
        print(f"Output file {output_path} already exists. Loading existing file...")
        return sc.read_h5ad(output_path)

    adatas = {}

    for a_sample_name, h5ad_path in sample_file_dict.items():

        if a_sample_name not in samples_to_pick:
            continue

        print(">>> Loading sample " , a_sample_name )
        sample_adata = sc.read_h5ad(h5ad_path)
        sample_adata.var_names_make_unique()
        sample_adata.obs["sample_id"] = a_sample_name
        sample_adata.obs.set_index(a_sample_name + "_" + sample_adata.obs.index.astype(str) , inplace=True )

        sample_adata.var['mt'] = sample_adata.var_names.str.startswith('mt-') | sample_adata.var_names.str.startswith('Mt-') | sample_adata.var_names.str.startswith('MT-') # annotate the group of mitochondrial genes as 'mt'
        sc.pp.calculate_qc_metrics(
            sample_adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
        ) 

        print(f"Sample {a_sample_name}'s adata shape before filtering: {sample_adata.shape}")       
        sc.pp.filter_cells(sample_adata, min_genes=1000)
        sc.pp.filter_cells(sample_adata, max_genes=10000)
        sc.pp.filter_cells(sample_adata, min_counts=1000)
        sc.pp.filter_cells(sample_adata, max_counts=50000)
        sc.pp.filter_genes(sample_adata, min_cells=10)
        
        sample_adata = sample_adata[sample_adata.obs['pct_counts_mt'] < 20 , :]
        print(f"Sample {a_sample_name}'s adata shape after filtering: {sample_adata.shape}")
        
        logger.info("Running Scrublet on " + a_sample_name)
        sc.external.pp.scrublet(sample_adata)
        sample_adata.obs['doublet_info'] = sample_adata.obs["predicted_doublet"].astype(str)
        logger.info("\nSample " + a_sample_name + " \n")     
        logger.info(sample_adata.obs['predicted_doublet'].value_counts())    

        if n_cells_to_subsample > 0:
            sc.pp.sample(sample_adata, n=min(int(n_cells_to_subsample), len(sample_adata)) , rng=my_random_seed)
        
        adatas[a_sample_name] = sample_adata
    
    # Get a list of sample names and AnnData objects
    sample_names = list(adatas.keys())
    adata_list = list(adatas.values())

    # Concatenate using the first AnnData object's concatenate method
    adata = adata_list[0].concatenate(
        adata_list[1:],
        batch_categories=sample_names,
        batch_key="sample_id",
        join='outer',    
        fill_value=0
    )

    # Add technology column
    sc_samples = ['MJ002','MJ004','MJ005','MJ007','MJ008','MJ014','MJ015']
    adata.obs["technology"] = [ "scRNASeq" if x in sc_samples else "snRNASeq" for x in adata.obs["sample_id"] ]

    # Normallize and log1p
    adata.raw = adata
    sc.pp.normalize_total(adata, inplace=True,target_sum=1e4)
    sc.pp.log1p(adata)

    # hvg annotation
    sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=3000, inplace=True , subset = False)
    print(f"There are {np.sum(adata.var['highly_variable'])} highly variable genes. ")

    # use scVelo's `moments` function for imputation - note that hack we're using here:
    # we're copying our `.X` matrix into the layers because that's where `scv.tl.moments`
    # expects to find counts for imputation
    adata.layers["spliced"] = adata.raw.X
    adata.layers["unspliced"] = adata.raw.X
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

    # Compute the cytotrace score
    ctk = CytoTRACEKernel(adata).compute_cytotrace()

    logger.info(">>> >> Data Scaling ...")
    sc.pp.scale(adata)

    logger.info(f"Final concatenated dataset shape: {adata.shape}")
    
    # Save the concatenated object
    logger.info(f"Saving concatenated AnnData to {output_path}")
    # Ensure directory exists
    output_dir = os.path.dirname(output_path)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Make sure all sata types are str for saving to h5ad
    for col in adata.var.select_dtypes(include="object").columns:
        adata.var[col] = adata.var[col].astype(str)
    for col in adata.raw.var.select_dtypes(include="object").columns:
        adata.raw.var[col] = adata.raw.var[col].astype(str)
        
    adata.write(output_path)
    print(f"Successfully saved concatenated data to {output_path}")

    return adata

def get_protein_activity(adata, network, save_dir, logger, num_cores=1):
    logger.info("Running VIPER with only TFs and coTFs")

    vp_data = pyviper.viper(gex_data=adata, # gene expression signature
                                interactome=network, # gene regulatory network
                                enrichment = "area",
                                output_as_anndata=True,
                                njobs=num_cores,
                                verbose=False)
    