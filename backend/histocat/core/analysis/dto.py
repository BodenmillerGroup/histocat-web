from typing import Any, Dict, Optional, Sequence

from pydantic import BaseModel


class PlotSeriesDto(BaseModel):
    """Scatter plot axis model."""

    label: str
    data: Sequence[float]


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


class ClassifyCellsSubmissionDto(BaseModel):
    """Cell classification submission model."""

    dataset_id: int
    result_id: Optional[int]
    channels: Sequence[str]
    thresholds: Dict
    n_estimators: int
    cell_classes: Optional[Any]
    annotations: Optional[Any]


class ClassifyCellsDto(BaseModel):
    """Cell classification result model."""

    cellClasses: Dict
    annotations: Sequence[Any]
