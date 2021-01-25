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
