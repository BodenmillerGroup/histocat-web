from typing import Any, Optional, Sequence

from pydantic import BaseModel


class PlotSeriesDto(BaseModel):
    """Scatter plot axis model."""
    label: str
    data: Sequence[float]


class ScatterPlotDto(BaseModel):
    """Scatter plot model."""
    x: PlotSeriesDto
    y: PlotSeriesDto
    z: Optional[PlotSeriesDto]
    heatmap: Optional[PlotSeriesDto]


class BoxPlotDto(BaseModel):
    """Box plot model."""
    series: Sequence[PlotSeriesDto]


class PcaDto(BaseModel):
    """PCA plot model."""

    x: PlotSeriesDto
    y: PlotSeriesDto
    z: Optional[PlotSeriesDto]
    heatmap: Optional[PlotSeriesDto]
    explained_variance_ratio: Optional[Sequence[float]]


class TsneSubmissionDto(BaseModel):
    """t-SNE submission model."""

    dataset_id: int
    acquisition_ids: Sequence[int]
    n_components: int
    perplexity: int
    learning_rate: int
    iterations: int
    theta: float
    init: str
    markers: Sequence[str]


class TsneDto(BaseModel):
    """t-SNE result model."""

    x: PlotSeriesDto
    y: PlotSeriesDto
    z: Optional[PlotSeriesDto]
    heatmap: Optional[PlotSeriesDto]


class UmapSubmissionDto(BaseModel):
    """UMAP submission model."""

    dataset_id: int
    acquisition_ids: Sequence[int]
    n_components: int
    n_neighbors: int
    metric: str
    min_dist: float
    markers: Sequence[str]


class UmapDto(BaseModel):
    """UMAP result model."""

    x: PlotSeriesDto
    y: PlotSeriesDto
    z: Optional[PlotSeriesDto]
    heatmap: Optional[PlotSeriesDto]


class PhenographSubmissionDto(BaseModel):
    """PhenoGraph submission model."""

    dataset_id: int
    acquisition_ids: Sequence[int]
    markers: Sequence[str]
    clustering_algo: str
    nearest_neighbors: int
    jaccard: bool
    primary_metric: str
    min_cluster_size: int


class PhenographDto(BaseModel):
    """PhenoGraph result model."""

    x: PlotSeriesDto
    y: PlotSeriesDto
    z: Optional[PlotSeriesDto]
    heatmap: Optional[PlotSeriesDto]


class RegionStatsSubmissionDto(BaseModel):
    """Region's stats submission model."""

    project_id: int
    acquisition_id: int
    region_polygon: Sequence[Any]


class RegionChannelStatsDto(BaseModel):
    """Region's channel stats model."""

    metal: str
    min: float
    max: float
    mean: float
