

import spatialdata
import sopa.io
from matplotlib import pyplot as plt
from tqdm import tqdm
import anndata
import numpy as np
import sopa
from pathlib import Path
import pandas as pd
#from sopa.segmentation.comseg import git
from spatialdata import SpatialData, read_zarr
import json
from sopa import segmentation

### load data
image_key = "image"
points_key = "transcripts" # (ignore this for multiplex imaging)
gene_column = "genes" # (optional) column of sdata[points_key] containing the gene names
sdata = sopa.io.uniform()



###############
# Cellpose segmentation
################
patches = sopa.segmentation.Patches2D(sdata, image_key, patch_width=1200, patch_overlap=50)
patches.write()
channels = ["DAPI"]
method = sopa.segmentation.methods.cellpose_patch(diameter=35, channels=channels, flow_threshold=2, cellprob_threshold=-6)
segmentation = sopa.segmentation.StainingSegmentation(sdata, method, channels, min_area=2500)
# The cellpose boundaries will be temporary saved here. You can choose a different path
cellpose_temp_dir = "tuto.zarr/.sopa_cache/cellpose2"
segmentation.write_patches_cells(cellpose_temp_dir)
cells = sopa.segmentation.StainingSegmentation.read_patches_cells(cellpose_temp_dir)
cells = sopa.segmentation.shapes.solve_conflicts(cells)
shapes_key = "cellpose_boundaries" # name of the key given to the cells in sdata.shapes
sopa.segmentation.StainingSegmentation.add_shapes(sdata, cells, image_key, shapes_key)


###############
# create comseg input
################

### create patch folder with files [transcipts.csv, centroid.csv]




##########################
#### apply comseg to the patches
##########################


MEAN_CELL_DIAMETER = 15  # in micrometer
MAX_CELL_RADIUS = 50  # in micrometer

config =  {
    "dict_scale" : {"x": 1, 'y': 1, "z": 1},
    "mean_cell_diameter" : MEAN_CELL_DIAMETER,
    "alpha" :  0.5,
    "min_rna_per_cell" : 5,
    "gene_column" : "genes", ## harcoded for sopa
    "image_csv_files" : ["transcripts.csv"], ## harcoded for sopa
    "centroid_csv_files" : ["centroid.csv"], ## harcoded for sopa
    "prior_name" : "cell", ## harcoded for sopa
}

path_transcript = "./tuto.zarr/.sopa_cache/comseg_test_centroid_csv"
from comseg import dataset as ds
from comseg import dictionary
import importlib
importlib.reload(ds)
importlib.reload(dictionary)
for patch_index in tqdm(range(4)):




    path_dataset_folder = Path(path_transcript) / str(patch_index)

    dataset = ds.ComSegDataset(
        path_dataset_folder=path_dataset_folder,
        dict_scale=config["dict_scale"],  # scale factor for the coordinates [x, y, z
        mean_cell_diameter = config["mean_cell_diameter"],
        gene_column = config["gene_column"],
        image_csv_files= config["image_csv_files"],
        centroid_csv_files=config["centroid_csv_files"],
        path_cell_centroid=path_dataset_folder,
        )

    dico_proba_edge, count_matrix = dataset.compute_edge_weight()

    Comsegdict = dictionary.ComSegDict(
        dataset=dataset,
        mean_cell_diameter=config["mean_cell_diameter"],
        prior_name="cell",
    )

    Comsegdict.run_all(max_cell_radius=MAX_CELL_RADIUS)

    #########
    ## save result for each patch
    #########
    anndata_comseg, json_dict = Comsegdict.anndata_from_comseg_result(
                    return_polygon = True,
                    alpha = config["alpha"],
                    min_rna_per_cell = config["min_rna_per_cell"])

    anndata_comseg.write_loom(path_dataset_folder / 'segmentation_counts.loom')
    with open(path_dataset_folder / "segmentation_polygons.json", 'w') as f:
        json.dump(json_dict['transcripts'], f)



################
# Aggregate result
#################


from sopa.segmentation.baysor.resolve import resolve
resolve(sdata, path_transcript, gene_column, min_area=10)
shapes_key = "baysor_boundaries"
aggregator = sopa.segmentation.Aggregator(sdata, image_key=image_key, shapes_key=shapes_key)
aggregator.compute_table(gene_column=gene_column, average_intensities=True)

#########
# Plot
##########
import spatialdata_plot
sdata\
    .pl.render_points(size=0.01, color="r")\
    .pl.render_images()\
    .pl.render_shapes(shapes_key, outline=True, fill_alpha=0, outline_color="w")\
    .pl.show("global")
plt.show()