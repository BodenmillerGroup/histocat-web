from typing import List, Optional

from pydantic import BaseModel

# Properties to receive via API on creation
from app.modules.channel.models import FilterModel, ScalebarModel, ChannelSettingsModel


class AnalysisModel(BaseModel):
    filter: FilterModel
    scalebar: ScalebarModel
    channels: List[ChannelSettingsModel]
    format: Optional[str] = 'png'
