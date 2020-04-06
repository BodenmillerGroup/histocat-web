from typing import Dict, List, Optional

from pydantic import BaseModel

from app.modules.acquisition.dto import AcquisitionDatasetDto


# Properties to receive via API on creation
class ROICreateModel(BaseModel):
    panorama_id: int
    origin_id: int
    meta: Dict[str, Optional[str]]


# Shared properties
class ROIModel(BaseModel):
    id: int
    panorama_id: int
    origin_id: int
    meta: Dict[str, Optional[str]]

    class Config:
        orm_mode = True


# Full ROI dataset
class ROIDatasetModel(ROIModel):
    acquisitions: List[AcquisitionDatasetDto]

    class Config:
        orm_mode = True
