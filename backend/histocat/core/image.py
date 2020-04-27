import logging
import sys
from typing import List, Tuple, Sequence

import cv2
import numpy as np
import tifffile
from mahotas import bwperim
from matplotlib import cm
from matplotlib.colors import LinearSegmentedColormap, rgb2hex, to_rgb
from skimage.color import label2rgb

from histocat.modules.acquisition.dto import (
    FilterDto,
    LegendDto,
    MaskSettingsDto,
    ScalebarDto,
)
from histocat.modules.analysis.dto import SegmentationSettingsDto

EPSILON = sys.float_info.epsilon  # Smallest possible difference.

logger = logging.getLogger(__name__)

OTSU_GRAYSCALE = "Otsu Grayscale"
OTSU_HUE_ALGORITHM = "Otsu Hue"
OTSU_SATURATION_ALGORITHM = "Otsu Saturation"
OTSU_LIGHTNESS_ALGORITHM = "Otsu Lightness"


def apply_filter(image: np.ndarray, filter: FilterDto):
    if filter.type == "gaussian":
        sigma = filter.settings.get("sigma")
        sigma = float(sigma) if sigma is not None and sigma != "" else 1.0
        kernel_size = filter.settings.get("kernel_size")
        kernel_size = (kernel_size, kernel_size) if kernel_size is not None and kernel_size != "" else (0, 0)
        return cv2.GaussianBlur(image, kernel_size, sigma)
    elif filter.type == "median":
        kernel_size = filter.settings.get("kernel_size")
        kernel_size = kernel_size if kernel_size is not None and kernel_size != "" else 3
        return cv2.medianBlur(image, kernel_size)


def colorize(image: np.ndarray, color: str):
    try:
        channel_color = to_rgb(color)
    except:
        channel_color = to_rgb("#ffffff")
    channel_colormap = LinearSegmentedColormap.from_list(None, [(0, 0, 0), channel_color])
    result = channel_colormap(image)
    return result * 255.0


def scale_image(image: np.ndarray, levels: Tuple[float, float]):
    channel_image = image - levels[0]
    channel_image /= levels[1] - levels[0]
    return np.clip(channel_image, 0, 1, out=channel_image)


def draw_mask(image: np.ndarray, mask_settings: MaskSettingsDto):
    mask = tifffile.imread(mask_settings.location)
    if mask_settings.colorize:
        return label2rgb(mask, image=image, alpha=0.3, bg_label=0, image_alpha=1, kind="avg")
    else:
        if mask_settings.gated:
            return mask_gated_img(image, mask, mask_settings.cell_ids)
        else:
            return mask_color_img(image, mask)


def mask_gated_img(image: np.ndarray, mask: np.ndarray, cell_ids: Sequence[int], color=(255, 255, 100), alpha=0.3):
    """
    img: cv2 image
    mask: bool or np.where
    color: BGR triplet [_, _, _]. Default: [0, 255, 255] is yellow.
    alpha: float [0, 1].

    Ref: http://www.pyimagesearch.com/2016/03/07/transparent-overlays-with-opencv/
    """
    m = np.isin(mask, cell_ids)
    mask[~m] = 0
    mask = bwperim(mask, n=2)
    # bit_mask = mask == 0
    mask_layer = image.copy()
    mask_layer[mask] = color
    cv2.addWeighted(mask_layer, alpha, image, 1 - alpha, 0, image)
    # cv2.add(image, mask_layer, image)
    return image


def mask_color_img(image: np.ndarray, mask: np.ndarray, color=(255, 255, 100), alpha=0.3):
    """
    img: cv2 image
    mask: bool or np.where
    color: BGR triplet [_, _, _]. Default: [0, 255, 255] is yellow.
    alpha: float [0, 1].

    Ref: http://www.pyimagesearch.com/2016/03/07/transparent-overlays-with-opencv/
    """
    mask = bwperim(mask, n=2)
    # bit_mask = mask == 0
    mask_layer = image.copy()
    mask_layer[mask] = color
    cv2.addWeighted(mask_layer, alpha, image, 1 - alpha, 0, image)
    # cv2.add(image, mask_layer, image)
    return image


def draw_scalebar(image: np.ndarray, scalebar: ScalebarDto):
    height, width, _ = image.shape
    length = 64
    cv2.line(
        image, (width - 60, height - 60), (width - 60 - length, height - 60), (255, 255, 255), 2, cv2.LINE_4,
    )

    scale_text = length
    if scalebar.settings is not None and "scale" in scalebar.settings:
        scale = scalebar.settings.get("scale")
        if scale is not None and scale != "":
            scale_text = int(length * float(scale))
    cv2.putText(
        image,
        f"{scale_text} um",
        (width - 60 - length, height - 30),
        cv2.FONT_HERSHEY_PLAIN,
        1,
        (255, 255, 255),
        1,
        cv2.LINE_4,
    )
    return image


def draw_legend(image: np.ndarray, legend_labels: List[Tuple[str, str, float]], legend: LegendDto):
    for i, label in enumerate(legend_labels):
        text = f"{label[0]} - {int(label[2])}" if legend.showIntensity else f"{label[0]}"
        (label_width, label_height), baseline = cv2.getTextSize(text, cv2.FONT_HERSHEY_DUPLEX, legend.fontScale, 1)
        cv2.rectangle(
            image,
            (5, (label_height + 20) * (i + 1) + 5),
            (15 + label_width, (label_height + 20) * (i + 1) - label_height - 5),
            (0, 0, 0),
            cv2.FILLED,
            cv2.LINE_AA,
        )

        b, g, r = tuple([255 * x for x in to_rgb(label[1])])
        color = (r, g, b)
        cv2.putText(
            image,
            text,
            (10, (label_height + 20) * (i + 1)),
            cv2.FONT_HERSHEY_DUPLEX,
            legend.fontScale,
            color,
            1,
            cv2.LINE_AA,
        )
    return image


def otsu_grayscale(image_rgb):
    image_gray = cv2.cvtColor(image_rgb, cv2.COLOR_BGR2GRAY)
    _, mask = cv2.threshold(image_gray.astype(np.uint8), 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)
    return mask


def _otsu_hls(image_rgb, channel_name: str, flip: bool):
    image_gray = cv2.cvtColor(image_rgb, cv2.COLOR_BGR2HLS)
    hue, lightness, saturation = np.split(image_gray, 3, axis=2)

    hsl = locals()[channel_name]
    hsl = hsl.reshape((hsl.shape[0], hsl.shape[1]))

    _, mask = cv2.threshold(hsl.astype(np.uint8), 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)

    if flip:
        mask = ~mask

    return mask


def otsu_hue(image_rgb):
    return _otsu_hls(image_rgb, channel_name="hue", flip=False)


def otsu_saturation(image_rgb):
    return _otsu_hls(image_rgb, channel_name="saturation", flip=False)


def otsu_lightness(image_rgb):
    return _otsu_hls(image_rgb, channel_name="lightness", flip=False)


def get_mask(image_rgb, settings: SegmentationSettingsDto):
    if settings.algorithm == OTSU_GRAYSCALE:
        return otsu_grayscale(image_rgb)
    elif settings.algorithm == OTSU_HUE_ALGORITHM:
        return otsu_hue(image_rgb)
    elif settings.algorithm == OTSU_SATURATION_ALGORITHM:
        return otsu_saturation(image_rgb)
    elif settings.algorithm == OTSU_LIGHTNESS_ALGORITHM:
        return otsu_lightness(image_rgb)
    return image_rgb


def apply_morphology(mask, settings: SegmentationSettingsDto):
    kernel = np.ones((settings.kernel_size, settings.kernel_size), np.uint8)
    mask = cv2.morphologyEx(mask, cv2.MORPH_OPEN, kernel, iterations=settings.iterations)
    mask = cv2.erode(mask, kernel, iterations=2)
    mask = cv2.dilate(mask, kernel, iterations=5)
    return mask


def get_heatmap_colors(values: np.ndarray, categorical_values: bool):
    keys = values.unique().tolist()
    color_range = np.linspace(0, 1, len(keys), endpoint=False).tolist()
    colors = [cm.tab20b(x) for x in color_range]
    color_dict = dict(zip(keys, colors))
    result = []
    for v in values:
        result.append(rgb2hex(color_dict[v]))
    return result
