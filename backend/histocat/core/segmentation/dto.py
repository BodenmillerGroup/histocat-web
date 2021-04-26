from typing import List, Optional

from pydantic import BaseModel


class SegmentationPreprocessingSettingsDto(BaseModel):
    """Segmentation preprocessing settings model."""

    threshold: Optional[bool]
    percentile: Optional[float]
    normalize: bool
    kernel_size: int


class SegmentationPostprocessingSettingsDto(BaseModel):
    """Segmentation postprocessing settings model."""

    radius: int
    maxima_threshold: float
    interior_threshold: float
    small_objects_threshold: int
    fill_holes_threshold: int
    interior_model: str
    maxima_model: str
    interior_model_smooth: int
    maxima_model_smooth: int
    pixel_expansion: Optional[int]


class SegmentationSubmissionDto(BaseModel):
    """Segmentation submission model."""

    dataset_name: Optional[str]
    dataset_description: Optional[str]
    model_id: int
    acquisition_ids: List[int]
    channels: List[str]
    nuclei_channels: List[str]
    cytoplasm_channels: List[str]
    preprocessing: SegmentationPreprocessingSettingsDto
    postprocessing: SegmentationPostprocessingSettingsDto
