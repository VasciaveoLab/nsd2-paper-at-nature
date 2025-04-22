# Import Libraries
import os
import glob
import scanpy as sc
import logging
import datetime
import numpy as np
import scvelo as scv
import pandas as pd
from scipy.stats import norm
import pyviper
import copy
import shutil
import gzip
from cellrank.kernels import CytoTRACEKernel

my_random_seed = 666
n_cells_to_subsample = 3e3

def filter_peaks(sample_adata, sample_path):

    # Try to locate the features file. It is usually named "features.tsv" or "genes.tsv"
    features_file = os.path.join(sample_path, "features.tsv.gz")
    if os.path.exists(features_file):
        df_features = pd.read_csv(features_file, sep="\t", header=None, compression="gzip")
    else:
        features_file = os.path.join(sample_path, "genes.tsv")
        df_features = pd.read_csv(features_file, sep="\t", header=None)
    
    # Filter out non-Gene Expression features (like "Peaks")
    if df_features.shape[1] >= 3 and "Gene Expression" in df_features[2].unique():
        
        df_filtered = df_features[df_features[2] == "Gene Expression"].copy()
        before = sample_adata.n_vars

        # Only filter if "Peaks" are actually present in the data
        if "Peaks" in df_features[2].unique():

            # Create set of valid genes from the features file with various normalization options
            valid_genes = set()
            valid_genes.update(df_filtered[1].tolist())
            
            # Add syntax fix due to scanpy and var_names_make_unique
            valid_genes.update(gene.lower() for gene in df_filtered[1].tolist())
            valid_genes.update(gene.split('-')[0] for gene in df_filtered[1].tolist() if '-' in gene)
            valid_genes.update(gene.lower().split('-')[0] for gene in df_filtered[1].tolist() if '-' in gene)
            
            # Create a mask for filtering genes
            mask = []
            for gene in sample_adata.var_names:
                # Check different possible variations of the gene name
                if (gene in valid_genes or 
                    gene.lower() in valid_genes or
                    (('-' in gene) and gene.split('-')[0] in valid_genes) or
                    (('-' in gene) and gene.lower().split('-')[0] in valid_genes)):
                    mask.append(True)
                else:
                    mask.append(False)
            
            sample_adata = sample_adata[:, mask].copy()
            
            # Report retention statistics
            after = sample_adata.n_vars
            print(f"Retained {after} features out of {before} ({after/before*100:.2f}%) after removing Peaks.", flush=True)
        else:
            print(f"No 'Peaks' found, keeping all features.", flush=True)
    else:
        print(f"Only 2 rows, keeping all features.", flush=True)
    
    return sample_adata

def fix_and_gzip_10x(base_path: str) -> str:

    # Create a new folder to put gzipped files in
    out_dir = os.path.join(base_path, "adjusted_for_scanpy")
    os.makedirs(out_dir, exist_ok=True)

    # Read the genes path
    p = os.path.join(base_path, "genes.tsv")
    df = pd.read_csv(
        p, sep="\t", header=None,
        compression=None,
        dtype=str
    )

    # trim to 3 columns
    if df.shape[1] > 3:
        df = df.iloc[:, :3]

    # write features.tsv.gz
    out_feats = os.path.join(out_dir, "features.tsv.gz")
    with gzip.open(out_feats, "wt") as f:
        df.to_csv(f, sep="\t", header=False, index=False)

    # gzip barcodes.tsv
    src = os.path.join(base_path, "barcodes.tsv")
    if os.path.exists(src):
        dst = os.path.join(out_dir, "barcodes.tsv.gz")
        with open(src, "rb") as fi, gzip.open(dst, "wb") as fo:
            shutil.copyfileobj(fi, fo)

    # 3) gzip matrix.mtx
    src = os.path.join(base_path, "matrix.mtx")
    if os.path.exists(src):
        dst = os.path.join(out_dir, "matrix.mtx.gz")
        with open(src, "rb") as fi, gzip.open(dst, "wb") as fo:
            shutil.copyfileobj(fi, fo)

    return out_dir

def trim_gzip_to_three_columns(base_path: str):

    # Path to genes
    gz_file    = os.path.join(base_path, "features.tsv.gz")
    backup_file = gz_file + ".backup"
    
    # Skip if already backed up or already ≤3 columns
    if os.path.exists(backup_file):
        return
    
    # Load all columns as strings
    df = pd.read_csv(gz_file, sep="\t", header=None,
                     compression="gzip", dtype=str)
    if df.shape[1] <= 3:
        return
    
    # Backup original, then trim and write
    shutil.move(gz_file, backup_file)
    df.iloc[:, :3].to_csv(
        gz_file,
        sep="\t",
        header=False,
        index=False,
        compression="gzip"
    )

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
            adata.var_names_make_unique()
            # Remove any non-gex data
            adata = filter_peaks(adata, filtered_feature_path)
            adata.obs['sample'] = sample_name
            adata.write(h5ad_file)
            print(f"Created {h5ad_file}")
        except Exception as e:
            try:
                # Read the adata
                adata = sc.read_10x_mtx(sample_dir)
                adata.var_names_make_unique()
                # Remove any non-gex data
                adata = filter_peaks(adata, filtered_feature_path)
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

def get_protein_activity(adata, network, save_name, new_data_dir, logger, num_cores=1):

    logger.info("Running VIPER with only TFs and coTFs")

    output_dir = os.path.join(new_data_dir, "pyviper_h5ad_outputs")
    os.makedirs(output_dir, exist_ok=True)

    output_path = os.path.join(output_dir, save_name)

    if os.path.exists(output_path):
        print(f"Output file {output_path} already exists. Loading existing file...")
        return sc.read_h5ad(output_path)

    print("Running pyviper...")
    vp_data = pyviper.viper(gex_data=adata, # gene expression signature
                                interactome=network, # gene regulatory network
                                enrichment = "area",
                                output_as_anndata=True,
                                njobs=num_cores,
                                verbose=False)
    
    print(vp_data.shape)
    vp_data.write(output_path)
    print(f"Successfully saved protein activity data to {output_path}")

    return vp_data

def concat_prot_act(vp_data_sc, vp_data_sn, new_data_dir, save_name, logger, harmony=False):

    logger.info("Combining protein activity")

    output_dir = os.path.join(new_data_dir, "pyviper_h5ad_outputs")
    os.makedirs(output_dir, exist_ok=True)

    output_path = os.path.join(output_dir, save_name)

    if os.path.exists(output_path):
        print(f"Output file {output_path} already exists. Loading existing file...")
        return sc.read_h5ad(output_path)

    # This is the data for the vast majority of the analysis
    vp_data = vp_data_sc.concatenate(
        vp_data_sn,
        batch_categories=["scRNASeq","snRNASeq"],
        batch_key="technology",
        join='inner',
        fill_value=0
    )

    print(f"Combined vp_data shape: {vp_data.shape}")

    logger.info(">>> >> Adding new sample nomenclature for the paper")
    df = pd.read_excel(os.path.join(new_data_dir, "Single cell dataset info-use-revised-nature-rebuttal.xlsx"),skiprows=1)
    df.rename(columns={"Name":"sample_id_for_paper","sc-RNA-seq ID":"sample_id"},inplace=True)
    df.set_index("sample_id",inplace=True)
    vp_data.obs = vp_data.obs.join( df , on="sample_id" )

    logger.info(">>> >> Adding layer of -log10(cpm) data to [mLog10] layer")
    vp_data.layers['mLog10'] = -1*np.log10(norm.sf( vp_data.X ))

    if harmony:
        # Apply harmony batch correction
        sc.tl.pca(vp_data, svd_solver='arpack', random_state=my_random_seed)
        sc.external.pp.harmony_integrate(vp_data, basis="X_pca" , key='technology', random_state=my_random_seed)

        ## 3 Cluster Solution
        n_pcs = 50
        n_neighbors=9
        resolution = 0.06
        seed_from_acdc = 1

        logger.info(">>> >> Computing Nearest Neighbors with n=%s and total PCs=%s" , n_neighbors , n_pcs)
        
        sc.pp.neighbors(vp_data, n_neighbors=n_neighbors, n_pcs=n_pcs,
                        use_rep="X_pca_harmony",
                        random_state=seed_from_acdc)

        logger.info(">>> >> Cluster Analysis wiht Leiden ...")
        sc.tl.leiden(vp_data, random_state=seed_from_acdc, resolution=resolution, key_added="leiden_pas")
        logger.info(">>> >> Cluster Analysis wiht Leiden | Solution from ACDC | knn=%s , PC=%s , resolution=%s , seed=%s" , n_neighbors , n_pcs , resolution, seed_from_acdc)
        vp_data.obs['leiden_pas'].cat.categories

    # Make sure all sata types are str for saving to h5ad
    for col in vp_data.obs.select_dtypes(include="object").columns:
        vp_data.obs[col] = vp_data.obs[col].astype(str)

    vp_data.write(output_path)
    print(f"Successfully saved protein activity data to {output_path}")

    return vp_data

def human_vp_conversion(new_data_dir, save_name, vp_data_full, logger):

    logger.info("Human gene pathway analysis")

    output_dir = os.path.join(new_data_dir, "pyviper_h5ad_outputs")
    os.makedirs(output_dir, exist_ok=True)

    output_path = os.path.join(output_dir, save_name)

    if os.path.exists(output_path):
        print(f"Output file {output_path} already exists. Loading existing file...")
        return

    logger.info(">>> >> Using VP Data with SIG and SURF for Pathway Analysis")
    vp_data_human = copy.deepcopy(vp_data_full)
    pyviper.pp.translate(vp_data_human, desired_format = "human_symbol")

    MSigDB_H_Pathways_regulon = pyviper.load.msigdb_regulon("h")
    # MSigDB_H_Pathways_regulon.filter_targets(vp_data_human.var_names)
    vp_cancer_hallmarks = pyviper.viper( gex_data=vp_data_human , interactome=MSigDB_H_Pathways_regulon , enrichment="area" , 
                                        min_targets=10,
                                        output_as_anndata=True , njobs=8 , verbose=False )
    
    print("Adding layer of -log10(cpm) data to [mLog10] layer to Beltran and Hieronymus pathway enrichment analysis")
    vp_cancer_hallmarks.layers['mLog10'] = -1*np.log10(norm.sf( vp_cancer_hallmarks.X ))

    # Make sure all sata types are str for saving to h5ad
    for col in vp_data_human.obs.select_dtypes(include="object").columns:
        vp_data_human.obs[col] = vp_data_human.obs[col].astype(str)

    vp_data_human.write(output_path)
    print(f"Successfully saved protein activity data to {output_path}")

    return
