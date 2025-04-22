import pandas as pd
import scanpy as sc
import os
import gc
import pyviper
import matplotlib
import warnings

from processing_functions import *

## These lines are to fix saving PDF with text recognizible by Adobe Illustrator
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

warnings.filterwarnings("ignore", category=UserWarning, module="matplotlib")
warnings.filterwarnings("ignore", category=UserWarning, module="scanpy")
sc.settings.verbosity = 0  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()

my_random_seed = 666
n_cells_to_subsample = 3e3
n_cores = 16

sc.settings.n_jobs = int(n_cores)
sc.settings.set_figure_params(dpi=80,vector_friendly=False)

### --- MAIN --- ###
def main():

    # Set the directories
    experiment_name = "nsd2-paper-experiment"
    new_data_dir = "/shares/vasciaveo_lab/data/nepc_organoid_project/new_data/"
    logging_path = "logs/"

    # set up logging
    print(f"Creating logger file at {logging_path}")
    logger = set_up_logger(logging_path, experiment_name)

    # Read the 10x data and save h5ad files for each sample in the new-data/samples_h5ad_files directory
    create_h5ad_for_samples(new_data_dir)

    # Get each of the h5ad file paths and save to a dict with sample_id as key
    sample_file_dict = get_sample_h5ad_dict(new_data_dir)

    # Load all the h5ad files for the samples and concat
    sc_samples = ["MJ002", "MJ004", "MJ005", "MJ007", "MJ008", "MJ014", "MJ015"]
    sc_adata_name = "adata_scRNASeq.h5ad"
    sc_adata = load_concat_adata(sample_file_dict, sc_samples, logger, new_data_dir, sc_adata_name)

    sn_samples = ["MJ018", "MJ019", "MJ020", "MJ021", "MJ022", "MJ023", "MJ024", "MJ025"]
    sn_adata_name = "adata_snRNASeq.h5ad"
    sn_adata = load_concat_adata(sample_file_dict, sn_samples, logger, new_data_dir, sn_adata_name)


    # Run pyviper to get protein activity data
    print("Inferring protein activity data...")
    network_path = os.path.join(new_data_dir, 'networks/MJ-02-04-05-07-metacells-with-chga-merged.tsv')
    aracne_network = pd.read_csv(network_path, delimiter="\t")
    network_interactome_full = pyviper.Interactome('organoids-network', aracne_network) # convert to class Interactome

    logger.info("Pruning Network with only TFs and coTFs")
    network_interactome_tfs_cotfs_pruned = copy.deepcopy(network_interactome_full)
    network_interactome_full_pruned = network_interactome_full

    tf_list = pyviper.load.TFs(species="mouse")
    cotf_list = pyviper.load.coTFs(species="mouse")

    network_interactome_tfs_cotfs_pruned.filter_regulators( regulators_keep = tf_list + cotf_list + ["Chga","Vim"] , )

    network_interactome_tfs_cotfs_pruned.prune(max_targets=50,eliminate=True) # prune interactome to have exactly 50 targets
    network_interactome_full_pruned.prune(max_targets=50,eliminate=True) # prune interactome to have exactly 50 targets

    logger.info(">>> >> Pruning Network with max_targets=50,eliminate=True and regulators_keep = tf_list + cotf_list ")

    print(network_interactome_tfs_cotfs_pruned.size())
    print(network_interactome_full_pruned.size())

    print("Running pyviper on sc...")
    sc_prot_act_name = "sc_prot_act_pruned.h5ad"
    sc_prot_act = get_protein_activity(sc_adata, network_interactome_tfs_cotfs_pruned, sc_prot_act_name, new_data_dir, logger, num_cores=n_cores)

    print("Running pyviper on sc with full interactome...")
    sc_prot_act_full_name = "sc_prot_act_full_pruned.h5ad"
    sc_prot_act_full = get_protein_activity(sc_adata, network_interactome_full_pruned, sc_prot_act_full_name, new_data_dir, logger, num_cores=n_cores)

    print("Running pyviper on sn...")
    sn_prot_act_name = "sn_prot_act_pruned.h5ad"
    sn_prot_act = get_protein_activity(sn_adata, network_interactome_tfs_cotfs_pruned, sn_prot_act_name, new_data_dir, logger, num_cores=n_cores)

    print("Running pyviper on sn with full interactome...")
    sn_prot_act_full_name = "sn_prot_act_full_pruned.h5ad"
    sn_prot_act_full = get_protein_activity(sn_adata, network_interactome_full_pruned, sn_prot_act_full_name, new_data_dir, logger, num_cores=n_cores)
    
    # Concatenate the protein activity data if not already created
    sc_adata.obs["RFP"] = [ 1 if i > 0 else 0 for i in sc_adata[:,"addgene26001"].X ]
    sc_adata.obs["RFP_int"] = sc_adata.obs["RFP"]
    sc_adata.obs["RFP"] = sc_adata.obs["RFP"].astype("category")

    sc_prot_act.obs = sc_adata.obs.join( sc_prot_act.obs , rsuffix="_pas")
    sn_prot_act.obs = sn_adata.obs.join( sn_prot_act.obs , rsuffix="_pas")
    sc_prot_act_full.obs = sc_adata.obs.join( sc_prot_act_full.obs , rsuffix="_pas")
    sn_prot_act_full.obs = sn_adata.obs.join( sn_prot_act_full.obs , rsuffix="_pas")

    combined_prot_act_name = "prot_act_concatenated.h5ad"
    concat_prot_act(sc_prot_act, sn_prot_act, new_data_dir, combined_prot_act_name, logger, harmony=True) 

    combined_prot_act_full_name = "prot_act_full_concatenated.h5ad"
    prot_act_full_concat = concat_prot_act(sc_prot_act_full, sn_prot_act_full, new_data_dir, combined_prot_act_full_name, logger) 

    human_prot_act_name = "human_prot_act.h5ad"
    human_vp_conversion(new_data_dir, human_prot_act_name, prot_act_full_concat, logger)

if __name__ == '__main__':
    main()
