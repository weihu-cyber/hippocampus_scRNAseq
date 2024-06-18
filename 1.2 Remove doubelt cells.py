import sys
import os
import scrublet as scr
import scipy.io
import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import umap
import numba
import numba.typed
import scanpy

## Set up parameters and variables ##
sample_dir = "/data/whu/home/whr_brain_project_2/data/matrix"
counts_matrix_dir = "/data/whu/home/whr_brain_project_SoupX_2_splited_adjusted/results/base/splited_strainedCounts/" ## Change this based on the path on your system
outdir = "/data/whu/home/whr_brain_project_SoupX_2_splited_adjusted/data/scrublet/" ## Change this based on the path on your system

if not os.path.isdir(outdir):
  os.mkdir(outdir)

plt.rc('font', size=14)
plt.rcParams['pdf.fonttype'] = 42

samples = os.listdir(sample_dir)

for s in samples:
    ## Basic run with scrublet
    counts_matrix = scanpy.read_10x_h5(counts_matrix_dir + s + ".h5")
    barcodes_df=pd.Series(counts_matrix.obs_names, name="Barcode")

    print('Counts matrix shape: {} rows, {} columns'.format(counts_matrix.shape[0], counts_matrix.shape[1]))
    scrub = scr.Scrublet(counts_matrix.X, expected_doublet_rate=0.06, sim_doublet_ratio = 2)
    doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=3,
                                                              min_cells=3,
                                                              min_gene_variability_pctl=85,
                                                              n_prin_comps=30)
    
    ### Plotting and saving
    scrub.plot_histogram();
    plt.savefig(os.path.join(outdir, s + '_doublet_score_histogram.png'))
    print('Running UMAP...')
    scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
    print('Done.')
    scrub.plot_embedding('UMAP', order_points=True);
    plt.savefig(os.path.join(outdir, s + '_UMAP.png'))

    results = pd.Series(scrub.predicted_doublets_, name="scrublet_DropletType")
    scores = pd.Series(scrub.doublet_scores_obs_, name="scrublet_Scores")
    dataframe = pd.concat([barcodes_df, results, scores], axis=1)
    dataframe.scrublet_DropletType = dataframe.scrublet_DropletType.replace(True, "doublet")
    dataframe.scrublet_DropletType = dataframe.scrublet_DropletType.replace(False, "singlet")

    dataframe.to_csv(os.path.join(outdir, s + '_scrublet_results.tsv'), sep = "\t", index = False)


    ### Make summary of singlets and doublets and write to file ###
    summary = pd.DataFrame(dataframe.scrublet_DropletType.value_counts())
    summary.index.name = 'Classification'
    summary.reset_index(inplace=True)
    summary = summary.rename({'scrublet_DropletType': 'Droplet N'}, axis=1)

    summary.to_csv(os.path.join(outdir, s + '_scrublet_summary.tsv'), sep = "\t", index = False)
    
