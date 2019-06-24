from __future__ import annotations

import os

import numpy as np
from imctools.io.txtparser import TxtParser
from sqlalchemy.orm import Session

from app.core.utils import timeit
from app.modules.acquisition import crud as acquisition_crud
from app.modules.acquisition.models import AcquisitionCreateModel
from app.modules.channel import crud as channel_crud
from app.modules.channel.models import ChannelCreateModel
from app.modules.slide import crud as slide_crud
from app.modules.slide.models import SlideCreateModel


class TextLoader:
    @classmethod
    @timeit
    def load(cls, db: Session, uri: str, experiment_id: int):
        txt = TxtParser(uri)
        file_name = os.path.basename(uri)
        imc_acquisition = txt.get_imc_acquisition()
        slide_meta = dict()
        slide_meta["Origin"] = imc_acquisition.origin
        slide_meta["ImageDescription"] = imc_acquisition.image_description
        slide_params = SlideCreateModel(
            experiment_id=experiment_id,
            name=file_name,
            filename=file_name,
            description=slide_meta["ImageDescription"],
            meta=slide_meta,
        )
        slide = slide_crud.create(db, params=slide_params)

        acquisition_meta = dict()
        acquisition_meta["Origin"] = imc_acquisition.origin
        acquisition_meta["ImageDescription"] = imc_acquisition.image_description
        acquisition_meta["MaxX"] = imc_acquisition.shape[0]
        acquisition_meta["MaxY"] = imc_acquisition.shape[1]

        acquisition_params = AcquisitionCreateModel(
            slide_id=slide.id,
            name=file_name,
            width=imc_acquisition.shape[0],
            height=imc_acquisition.shape[1],
            description=acquisition_meta["ImageDescription"],
            meta=acquisition_meta,
        )
        acquisition = acquisition_crud.create(db, params=acquisition_params)

        for i in range(imc_acquisition.n_channels):
            img = imc_acquisition.get_img_by_label(
                imc_acquisition.channel_labels[i]
            )
            channel_meta = dict()
            channel_meta["Label"] = imc_acquisition.channel_labels[i]
            channel_meta["Metal"] = imc_acquisition.channel_metals[i]
            channel_meta["Mass"] = imc_acquisition.channel_mass[i]
            channel_params = ChannelCreateModel(
                acquisition_id=acquisition.id,
                name=channel_meta["Label"],
                metal=channel_meta["Metal"],
                mass=channel_meta["Mass"],
                max_intensity=img.max(),
                min_intensity=img.min(),
                meta=channel_meta,
            )
            channel = channel_crud.create(db, params=channel_params)
            np.save(os.path.join(channel.location, "origin.npy"), img)
