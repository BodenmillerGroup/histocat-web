from typing import Optional

import numpy as np
from fastapi import HTTPException
from imctools.io.ometiff.ometiffparser import OmeTiffParser
from sqlalchemy.orm import Session
from starlette.status import HTTP_404_NOT_FOUND

from histocat.core.acquisition import service as acquisition_service
from histocat.core.acquisition.dto import ChannelStackDto
from histocat.core.image import colorize, scale_image


def get_additive_image(db: Session, params: ChannelStackDto):
    additive_image: Optional[np.ndarray] = None

    acquisition = acquisition_service.get_by_id(db, id=params.acquisitionId)
    if not acquisition:
        raise HTTPException(status_code=HTTP_404_NOT_FOUND, detail="Acquisition not found.")
    parser = OmeTiffParser(acquisition.location)
    acq = parser.get_acquisition_data()

    for channel in params.channels:
        data = acq.get_image_by_name(channel.name)
        item = acquisition.channels.get(channel.name)
        levels = (
            (channel.min, channel.max)
            if channel.min is not None and channel.max is not None
            else (item.get("min_intensity"), item.get("max_intensity"))
        )
        data = scale_image(data, levels)

        color = channel.color if channel.color else "#ffffff"
        image = colorize(data, color)

        if additive_image is None:
            additive_image = np.zeros(shape=(data.shape[0], data.shape[1], 4), dtype=np.float32)
        additive_image += image
    return np.clip(additive_image, 0, 1, out=additive_image)
