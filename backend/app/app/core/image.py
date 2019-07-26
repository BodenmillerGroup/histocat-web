import logging
from typing import Tuple, List

import cv2
import numpy as np
from matplotlib.colors import to_rgb, LinearSegmentedColormap
from skimage import filters

from app.modules.channel.models import FilterModel, ScalebarModel, LegendModel

logger = logging.getLogger(__name__)


def apply_filter(image: np.ndarray, filter: FilterModel):
    if filter.type == 'gaussian':
        sigma = filter.settings.get('sigma')
        sigma = float(sigma) if sigma is not None and sigma != '' else 1.0
        mode = filter.settings.get('mode')
        mode = mode if mode is not None and mode != '' else 'nearest'
        return filters.gaussian(image, sigma=sigma, mode=mode, output=image, preserve_range=True)
    elif filter.type == 'median':
        mode = filter.settings.get('mode')
        mode = mode if mode is not None and mode != '' else 'nearest'
        return filters.median(image, behavior='ndimage', mode=mode, out=image)


def colorize(image: np.ndarray, color: str):
    try:
        channel_color = to_rgb(color)
    except:
        channel_color = to_rgb('#ffffff')
    channel_colormap = LinearSegmentedColormap.from_list(None, [(0, 0, 0), channel_color])
    result = channel_colormap(image)
    return result * 255.0


def scale_image(image: np.ndarray, levels: Tuple[float, float]):
    channel_image = image - levels[0]
    channel_image /= levels[1] - levels[0]
    return np.clip(channel_image, 0, 1, out=channel_image)


def draw_scalebar(image: np.ndarray, scalebar: ScalebarModel):
    width, height, _ = image.shape
    length = 64
    cv2.line(
        image,
        (width - 60, height - 60),
        (width - 60 - length, height - 60),
        (255, 255, 255),
        2,
        cv2.LINE_4
    )
    cv2.line(
        image,
        (width - 60, height - 55),
        (width - 60, height - 65),
        (255, 255, 255),
        2,
        cv2.LINE_4
    )
    cv2.line(
        image,
        (width - 60 - length, height - 55),
        (width - 60 - length, height - 65),
        (255, 255, 255),
        2,
        cv2.LINE_4
    )

    scale_text = length
    if scalebar.settings is not None and 'scale' in scalebar.settings:
        scale = scalebar.settings.get('scale')
        if scale is not None and scale != '':
            scale_text = int(length * float(scale))
    cv2.putText(
        image,
        f'{scale_text} um',
        (width - 60 - length, height - 30),
        cv2.FONT_HERSHEY_PLAIN,
        1,
        (255, 255, 255),
        1,
        cv2.LINE_4
    )
    return image


def draw_legend(image: np.ndarray, legend_labels: List[Tuple[str, str]], legend: LegendModel):
    for i, label in enumerate(legend_labels):
        (label_width, label_height), baseline = cv2.getTextSize(label[0], cv2.FONT_HERSHEY_DUPLEX, legend.fontScale, 1)
        cv2.rectangle(
            image,
            (
                5,
                (label_height + 20) * (i + 1) + 5
            ),
            (
                15 + label_width,
                (label_height + 20) * (i + 1) - label_height - 5
            ),
            (0, 0, 0),
            cv2.FILLED,
            cv2.LINE_AA
        )

        b, g, r = tuple([255 * x for x in to_rgb(label[1])])
        color = (r, g, b)
        cv2.putText(
            image,
            label[0],
            (
                10,
                (label_height + 20) * (i + 1)
            ),
            cv2.FONT_HERSHEY_DUPLEX,
            legend.fontScale,
            color,
            1,
            cv2.LINE_AA
        )
    return image
