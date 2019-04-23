from __future__ import annotations

import os
from typing import List

from fastapi import UploadFile
from imctools.io import mcdparser
from sqlalchemy.orm import Session

from app.modules.experiment.db import Experiment
from app.modules.slide import crud as slide_crud
from app.modules.slide.models import SlideCreateModel


class McdLoader:

    @classmethod
    async def load(cls, upload_file: UploadFile, db: Session, experiment: Experiment):
        with mcdparser.McdParser(upload_file.filename, filehandle=upload_file.file) as mcd:
            slide_item = mcd.meta.objects['Slide']['0']
            file_name = os.path.basename(upload_file.filename)
            slide_params = SlideCreateModel(
                experiment_id = experiment.id,
                name = file_name,
                description = upload_file.filename,
                meta=slide_item.properties
            )
            slide = slide_crud.create(db, params=slide_params)
            pass
