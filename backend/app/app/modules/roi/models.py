from datetime import datetime
from typing import Optional, List

from pydantic import BaseModel

from app.modules.acquisition.models import AcquisitionDatasetModel
from app.modules.roi_point.models import ROIPointModel


# Properties to receive via API on creation
class ROICreateModel(BaseModel):
    panorama_id: int
    roi_type: Optional[str]


# Shared properties
class ROIModel(BaseModel):
    id: int
    panorama_id: int
    roi_type: Optional[str]
    location: Optional[str]
    created_at: datetime

    class Config:
        orm_mode = True


# Full ROI dataset
class ROIDatasetModel(ROIModel):
    acquisitions: List[AcquisitionDatasetModel]
    roi_points: List[ROIPointModel]

    class Config:
        orm_mode = True
