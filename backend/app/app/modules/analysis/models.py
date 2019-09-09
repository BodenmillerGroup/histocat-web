from typing import List, Optional

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
    format: Optional[str] = 'png'
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
    acquisition_id: int
    n_components: int
    markers: List[str]
    heatmap: Optional[str]


class TSNEModel(BaseModel):
    """
    t-SNE result model
    """
    x: PlotSeriesModel
    y: PlotSeriesModel
    z: Optional[PlotSeriesModel]
    heatmap: Optional[PlotSeriesModel]
