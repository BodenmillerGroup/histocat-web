import csv
import os
import logging
import shutil
import re
import zipfile
from collections import namedtuple
from typing import Dict

from sqlalchemy.orm import Session

from imctools.scripts import convertfolder2imcfolder

from app.modules.dataset.db import Dataset
from app.modules.acquisition import crud as acquisition_crud
from app.modules.panorama.db import Panorama
from app.modules.roi.db import ROI
from app.modules.slide.db import Slide

logger = logging.getLogger(__name__)


# 1) input folders
# the folders with the ziped acquisition files for the analysis
# -> If you want to analyse  your own data, put the zipped inuput files (see above `Data requirements`)
#    into a directory and add it to the folder
# Example: if you put it into a folder 'data' which is a subdirectory of the ÌmcSegmentationPipeline folder change this to
#  folders = ['../data']
# Example2: if you put your data into a folder C://Users/dummy/mydata change this to
# folders = ['C://Users/dummy/mydata']

folders = ['../example_data']

# 4) pannel
# pannel:
# This CSV file is specific to the pannel used for the Acquisitions
# It is a comma seperated file that contains metadata about the antibodies and isotopes measured
# This absolutely **needs** to be adapted if you work with your own files1
# Please look at the example!
csv_pannel = '../config/example_pannel.csv'
# Three columns are obligatory
csv_pannel_metal = 'Metal Tag' # Contains the isotope channel measured int he form (Metal)(Mass), e.g. Ir191, Yb171 etc.
# ilastik columm: Bool, either 0 or 1: this selects channels to be used for cellular segmentation
# It is recommended to choose all channels/isotope that gave a signal that could help for cell identification,
# e.g. nuclear, cytoplasmic or membranous signal
# The more channels selected, the slower the pixel classification
csv_pannel_ilastik = 'ilastik'
# full column: Contains the channels that should be quantified/measured in cellprofiler
csv_pannel_full = 'full'

suffix_full = '_full'
suffix_ilastik = '_ilastik'
suffix_ilastik_scale = '_s2'
suffix_mask = '_mask.tiff'
suffix_probablities = '_Probabilities'

failed_images = list()

# Make a list of all the analysis stacks with format:
# (CSV_NAME, SUFFIX, ADDSUM)
# CSV_NAME: name of the column in the CSV to be used
# SUFFIX: suffix of the tiff
# ADDSUM: BOOL, should the sum of all channels be added as the first channel?
list_analysis_stacks =[
    (csv_pannel_ilastik, suffix_ilastik, 1),
    (csv_pannel_full, suffix_full, 0)
]


FileLink = namedtuple('FileLink', ['src', 'dst'])


def _save_meta_csv(items: Dict[int, dict], filename: str):
    """
    Writes the xml data as csv tables

    """
    with open(filename, 'w') as f:
        cols = next(iter(items.values())).keys()
        writer = csv.DictWriter(f, sorted(cols))
        writer.writeheader()
        for row in items.values():
            writer.writerow(row)


def prepare_dataset(db: Session, dataset: Dataset):
    # output folder where the output files will be saved
    folder_base = dataset.location

    # parameters for resizing the images for ilastik
    folder_analysis = os.path.join(folder_base, 'tiffs')
    folder_ilastik = os.path.join(folder_base, 'ilastik')
    folder_ome = os.path.join(folder_base, 'ometiff')
    folder_cp = os.path.join(folder_base, 'cpout')
    folder_histocat = os.path.join(folder_base, 'histocat')
    folder_uncertainty = os.path.join(folder_base, 'uncertainty')

    for folder in [folder_base, folder_analysis, folder_ilastik, folder_ome, folder_cp, folder_histocat, folder_uncertainty]:
        if not os.path.exists(folder):
            os.makedirs(folder)

    input = dataset.meta['input']
    acquisition_ids = input['acquisition_ids']
    metals = input['metals']
    channel_settings = input['channel_settings']

    slide_csv: Dict[int, dict] = dict()
    panorama_csv: Dict[int, dict] = dict()
    roi_csv: Dict[int, dict] = dict()
    acquisition_csv: Dict[int, dict] = dict()
    channel_csv: Dict[int, dict] = dict()

    short_slide_name = None

    for acquisition_id in acquisition_ids:
        acquisition = acquisition_crud.get(db, id=acquisition_id)
        roi: ROI = acquisition.roi
        panorama: Panorama = roi.panorama
        slide: Slide = panorama.slide

        short_slide_name = re.sub('_s\d+$', '', slide.metaname)

        slide_csv[slide.id] = slide.meta
        panorama_csv[panorama.id] = panorama.meta
        # roi_csv[roi.id] = roi.meta
        acquisition_csv[acquisition.id] = acquisition.meta

        links = [
            FileLink(
                src=os.path.join(slide.location, short_slide_name + '_schema.xml'),
                dst=os.path.join(folder_ome, short_slide_name + '_schema.xml')
            ),
            # FileLink(
            #     src=os.path.join(slide.location, short_slide_name + '_Slide_meta.csv'),
            #     dst=os.path.join(folder_ome, short_slide_name + '_Slide_meta.csv')
            # ),
            # FileLink(
            #     src=os.path.join(slide.location, short_slide_name + '_AcquisitionChannel_meta.csv'),
            #     dst=os.path.join(folder_ome, short_slide_name + '_AcquisitionChannel_meta.csv')
            # ),
            # FileLink(
            #     src=os.path.join(slide.location, short_slide_name + '_Acquisition_meta.csv'),
            #     dst=os.path.join(folder_ome, short_slide_name + '_Acquisition_meta.csv')
            # ),
            # FileLink(
            #     src=os.path.join(slide.location, short_slide_name + '_AcquisitionROI_meta.csv'),
            #     dst=os.path.join(folder_ome, short_slide_name + '_AcquisitionROI_meta.csv')
            # ),
            # FileLink(
            #     src=os.path.join(slide.location, short_slide_name + '_Panorama_meta.csv'),
            #     dst=os.path.join(folder_ome, short_slide_name + '_Panorama_meta.csv')
            # ),
            # FileLink(
            #     src=os.path.join(slide.location, short_slide_name + '_ROIPoint_meta.csv'),
            #     dst=os.path.join(folder_ome, short_slide_name + '_ROIPoint_meta.csv')
            # ),
            FileLink(
                src=os.path.join(slide.location, slide.metaname + '_slide.png'),
                dst=os.path.join(folder_ome, slide.metaname + '_slide.png')
            ),
            FileLink(
                src=os.path.join(panorama.location, panorama.metaname + '_pano.png'),
                dst=os.path.join(folder_ome, panorama.metaname + '_pano.png')
            ),
        ]
        for link in links:
            if not os.path.exists(link.dst):
                shutil.copy2(link.src, link.dst)

    _save_meta_csv(slide_csv, os.path.join(folder_ome, short_slide_name + '_Slide_meta.csv'))