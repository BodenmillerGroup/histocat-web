from __future__ import annotations

import os
from tempfile import NamedTemporaryFile

import h5py
from fastapi import UploadFile
from imctools.io.ometiffparser import OmetiffParser
from sqlalchemy.orm import Session

from app.modules.acquisition import crud as acquisition_crud
from app.modules.acquisition.models import AcquisitionCreateModel
from app.modules.channel import crud as channel_crud
from app.modules.channel.models import ChannelCreateModel
from app.modules.experiment.db import Experiment
from app.modules.slide import crud as slide_crud
from app.modules.slide.models import SlideCreateModel


class OmeTiffLoader:
    @classmethod
    def load(cls, upload_file: UploadFile, db: Session, experiment: Experiment):
        with NamedTemporaryFile() as file:
            file.write(upload_file.file.read())
            file.seek(0)
            ome = OmetiffParser(file)
            file_name = os.path.basename(upload_file.filename)
            slide_meta = ome.meta_dict
            slide_params = SlideCreateModel(
                experiment_id=experiment.id,
                name=slide_meta["image_ID"],
                filename=file_name,
                description=ome.origin,
                meta=slide_meta,
            )
            slide = slide_crud.create(db, params=slide_params)

            imc_acquisition = ome.get_imc_acquisition()
            acquisition_meta = dict()
            acquisition_meta["Origin"] = imc_acquisition.origin
            acquisition_meta["ImageDescription"] = imc_acquisition.image_description
            acquisition_meta["MaxX"] = imc_acquisition.shape[0]
            acquisition_meta["MaxY"] = imc_acquisition.shape[1]

            acquisition_params = AcquisitionCreateModel(
                slide_id=slide.id,
                name=slide_meta["image_ID"],
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
                with h5py.File(os.path.join(channel.location, "origin.h5"), "w") as f:
                    f.create_dataset("image", data=img)
                # cv2.imwrite(os.path.join(channel.location, 'thumbnail.png'), img)
