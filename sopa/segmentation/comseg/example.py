

import spatialdata
import sopa.io
from matplotlib import pyplot as plt
from tqdm import tqdm
import anndata
import numpy as np
import sopa
from pathlib import Path
import pandas as pd
from sopa.segmentation.comseg import add_centroids_to_sdata
from spatialdata import SpatialData, read_zarr
import json

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
config_comseg = {}
baysor_temp_dir = "tuto.zarr/.sopa_cache/comseg"
patches = sopa.segmentation.Patches2D(sdata, points_key, patch_width=200, patch_overlap=50)
valid_indices = patches.patchify_transcripts(baysor_temp_dir, config=config_comseg, use_prior=True)


#################
### Compute centroid
#################
sdata_centroid = SpatialData()
sdata_centroid['cellpose_boundaries'] = sdata['cellpose_boundaries']
sdata_centroid['transcripts'] = sdata['transcripts']
sdata_centroid = add_centroids_to_sdata(sdata_centroid,
                           points_key="transcripts",
                           shapes_key='cellpose_boundaries',
                               z_constant=1)


baysor_temp_dir = "tuto.zarr/.sopa_cache/comseg_centroid"
config_comseg ={}
points_key = "centroid"
patches = sopa.segmentation.Patches2D(sdata_centroid, points_key, patch_width=200, patch_overlap=50)
valid_indices = patches.patchify_transcripts(baysor_temp_dir, config=config_comseg, use_prior=True)


##########################
#### apply comseg to the patches
##########################


MEAN_CELL_DIAMETER = 15  # in micrometer
MAX_CELL_RADIUS = 50  # in micrometer

path_transcript = "tuto.zarr/.sopa_cache/comseg"
path_centroid = "tuto.zarr/.sopa_cache/comseg_centroid"
from comseg import dataset as ds
from comseg import dictionary
for patch_index in tqdm(list(range(len(patches.ilocs)))):




    path_dataset_folder = Path(path_transcript) / str(patch_index)
    path_dataset_folder_centroid = Path(path_centroid) / str(patch_index)
    dataset = ds.ComSegDataset(
        path_dataset_folder=path_dataset_folder,
        dict_scale={"x": 1, 'y': 1, "z": 1},
        mean_cell_diameter = MEAN_CELL_DIAMETER,
        gene_column = "genes",
        )

    dico_proba_edge, count_matrix = dataset.compute_edge_weight(  # in micrometer
        images_subset=None,
        n_neighbors=40,
        sampling=True,
        sampling_size=10000
    )

    Comsegdict = dictionary.ComSegDict(
        dataset=dataset,
        mean_cell_diameter=MEAN_CELL_DIAMETER,
        community_detection="with_prior",
        prior_name="cell",
    )

    Comsegdict.run_all(max_cell_radius=MAX_CELL_RADIUS,
                       path_dataset_folder_centroid=path_dataset_folder_centroid,
                       file_extension=".csv")

    #########
    ## save result for each patch
    #########
    anndata_comseg, json_dict = Comsegdict.anndata_from_comseg_result(
                    return_polygon = True,
                    alpha = 0.5,
                    min_rna_per_cell = 5)
    random_string = str(np.random.randint(0, 10000))
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

sdata\
    .pl.render_points(size=0.01, color="r")\
    .pl.render_images()\
    .pl.render_shapes(shapes_key, outline=True, fill_alpha=0, outline_color="w")\
    .pl.show("global")
plt.show()