from __future__ import annotations

import os

import numpy as np
from fastapi import UploadFile
from imctools.io import mcdparser
from sqlalchemy.orm import Session

from app.core.utils import timeit
from app.modules.acquisition import crud as acquisition_crud
from app.modules.acquisition.models import AcquisitionCreateModel
from app.modules.channel import crud as channel_crud
from app.modules.channel.models import ChannelCreateModel
from app.modules.experiment.db import Experiment
from app.modules.slide import crud as slide_crud
from app.modules.slide.models import SlideCreateModel


class McdLoader:
    @classmethod
    @timeit
    def load(cls, upload_file: UploadFile, db: Session, experiment: Experiment):
        with mcdparser.McdParser(
            upload_file.filename, filehandle=upload_file.file
        ) as mcd:
            slide_item = mcd.meta.objects["Slide"]["0"]
            os.path.basename(upload_file.filename)
            slide_params = SlideCreateModel(
                experiment_id=experiment.id,
                name=slide_item.properties["Name"]
                if "Name" in slide_item.properties
                else slide_item.properties["Description"],
                filename=slide_item.properties["Filename"],
                width_um=slide_item.properties["WidthUm"],
                height_um=slide_item.properties["HeightUm"],
                description=slide_item.properties["Description"],
                meta=slide_item.properties,
            )
            slide = slide_crud.create(db, params=slide_params)
            for panorama_item in slide_item.childs["Panorama"].values():
                if "AcquisitionROI" in panorama_item.childs:
                    for acquisition_roi_item in panorama_item.childs[
                        "AcquisitionROI"
                    ].values():
                        for acquisition_item in acquisition_roi_item.childs[
                            "Acquisition"
                        ].values():
                            acquisition_params = AcquisitionCreateModel(
                                slide_id=slide.id,
                                name=acquisition_item.properties["Description"],
                                width=acquisition_item.properties["MaxX"],
                                height=acquisition_item.properties["MaxY"],
                                description=acquisition_item.properties["Description"],
                                meta=acquisition_item.properties,
                            )
                            acquisition = acquisition_crud.create(
                                db, params=acquisition_params
                            )
                            imc_acquisition = mcd.get_imc_acquisition(
                                acquisition_item.properties["ID"]
                            )
                            for i in range(imc_acquisition.n_channels):
                                meta = dict()
                                label = imc_acquisition.channel_labels[i]
                                meta["Label"] = label
                                meta["Metal"] = imc_acquisition.channel_metals[i]
                                meta["Mass"] = imc_acquisition.channel_mass[i]
                                img = imc_acquisition.get_img_by_label(label)
                                channel_params = ChannelCreateModel(
                                    acquisition_id=acquisition.id,
                                    name=label,
                                    metal=meta["Metal"],
                                    mass=meta["Mass"],
                                    max_intensity=img.max(),
                                    min_intensity=img.min(),
                                    meta=meta,
                                )
                                channel = channel_crud.create(db, params=channel_params)
                                np.save(
                                    os.path.join(channel.location, "origin.npy"), img
                                )
                                # cv2.imwrite(os.path.join(channel.location, 'thumbnail.png'), img)
