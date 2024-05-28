
from spatialdata.models import PointsModel
import pandas as pd
from sopa._sdata import to_intrinsic
from spatialdata import SpatialData

gene_column: str | None = None,
area_ratio_threshold: float = 0.25,


def add_centroids_to_sdata(sdata: SpatialData,
                           points_key: str = "transcripts",
                           shapes_key: str = 'cellpose_boundaries',
                           z_constant: float | None = None) -> SpatialData:
    """
    Function to add centroids to the sdata
    Args:
        sdata:
        points_key: name of the RNA coodirdiante in sdata
        shapes_key: name of the cellpose boundaries in sdata
        z_constant: this value corresponds to the z value of the cell centroid.

    Returns: sdata with the centroids

    """
    centroid = sdata[shapes_key].geometry.centroid
    x_centroid = list(centroid.geometry.x)
    y_centroid = list(centroid.geometry.y)
    if z_constant is not None:
        z_centroid = [z_constant] * len(y_centroid)
        coords = pd.DataFrame({"x": x_centroid, "y": y_centroid, "z": z_centroid})

    else:
        if "z" in sdata[points_key].columns:
            z = list(sdata[points_key].z.unique().compute())
            assert len(z)==1, "3D point cloud with 2D segmentation, manually set z_constant"
            z_centroid = [z[0]] * len(y_centroid)
            coords = pd.DataFrame({"x": x_centroid, "y": y_centroid, "z": z_centroid})
        else:
            coords = pd.DataFrame({"x": x_centroid, "y": y_centroid})
    points = PointsModel.parse(coords)
    sdata['centroid'] = points
    sdata['centroid'] = to_intrinsic(sdata, sdata['centroid'], points_key)
    return sdata

