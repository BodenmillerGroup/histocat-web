from typing import Optional, Sequence

from pydantic import BaseModel


class SegmentationFilterSettingsDto(BaseModel):
    """Segmentation filter settings model."""

    apply: bool
    type: str
    kernel_size: Optional[int]
    sigma: Optional[float]

class SegmentationSubmissionDto(BaseModel):
    """Segmentation submission model."""

    model_id: int
    acquisition_ids: Sequence[int]
    nuclei_channels: Sequence[str]
    cytoplasm_channels: Sequence[str]
    scaling_factor: float
    upper_limit: float
    overlap: float
    filter_settings: SegmentationFilterSettingsDto
