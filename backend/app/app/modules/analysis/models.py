from typing import List, Optional, Any

from pydantic import BaseModel

# Properties to receive via API on creation
from app.modules.channel.models import FilterModel, ScalebarModel, ChannelSettingsModel


class SegmentationSettingsModel(BaseModel):
    algorithm: str
    iterations: int
    kernel_size: int
    mask_color: str
    result_type: str


class AnalysisModel(BaseModel):
    filter: FilterModel
    scalebar: ScalebarModel
    channels: List[ChannelSettingsModel]
    format: Optional[str] = "png"
    settings: SegmentationSettingsModel


class PlotSeriesModel(BaseModel):
    """
    Scatter plot axis model
    """

    label: str
    data: List[float]


class ScatterPlotModel(BaseModel):
    """
    Scatter plot model
    """

    x: PlotSeriesModel
    y: PlotSeriesModel
    z: Optional[PlotSeriesModel]
    heatmap: Optional[PlotSeriesModel]


class BoxPlotModel(BaseModel):
    """
    Box plot model
    """

    series: List[PlotSeriesModel]


class PCAModel(BaseModel):
    """
    PCA plot model
    """

    x: PlotSeriesModel
    y: PlotSeriesModel
    z: Optional[PlotSeriesModel]
    heatmap: Optional[PlotSeriesModel]
    explained_variance_ratio: Optional[List[float]]


class TSNESubmissionModel(BaseModel):
    """
    t-SNE submission model
    """

    dataset_id: int
    acquisition_ids: List[int]
    n_components: int
    perplexity: int
    learning_rate: int
    iterations: int
    theta: float
    init: str
    markers: List[str]


class TSNEModel(BaseModel):
    """
    t-SNE result model
    """

    x: PlotSeriesModel
    y: PlotSeriesModel
    z: Optional[PlotSeriesModel]
    heatmap: Optional[PlotSeriesModel]


class UMAPSubmissionModel(BaseModel):
    """
    UMAP submission model
    """

    dataset_id: int
    acquisition_ids: List[int]
    n_components: int
    n_neighbors: int
    metric: str
    min_dist: float
    markers: List[str]


class UMAPModel(BaseModel):
    """
    UMAP result model
    """

    x: PlotSeriesModel
    y: PlotSeriesModel
    z: Optional[PlotSeriesModel]
    heatmap: Optional[PlotSeriesModel]


class PhenographSubmissionModel(BaseModel):
    """
    PhenoGraph submission model
    """

    dataset_id: int
    acquisition_ids: List[int]
    markers: List[str]
    nearest_neighbors: int
    jaccard: bool
    primary_metric: str
    min_cluster_size: int


class PhenographModel(BaseModel):
    """
    PhenoGraph result model
    """

    x: PlotSeriesModel
    y: PlotSeriesModel
    z: Optional[PlotSeriesModel]
    heatmap: Optional[PlotSeriesModel]


class RegionStatsSubmissionModel(BaseModel):
    """
    Region's stats submission model
    """

    experiment_id: int
    acquisition_id: int
    region_polygon: List[Any]


class RegionChannelStatsModel(BaseModel):
    """
    Region's channel stats model
    """

    metal: str
    min: float
    max: float
    mean: float
