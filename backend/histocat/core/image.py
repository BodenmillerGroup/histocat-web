import logging
from typing import Optional, Tuple

import cv2
import numpy as np
import tifffile
from matplotlib import cm
from matplotlib.colors import LinearSegmentedColormap, to_rgb
from skimage import img_as_ubyte, io
from skimage.color import label2rgb
from skimage.segmentation import find_boundaries

from histocat.core.acquisition.dto import FilterDto, MaskSettingsDto, ScalebarDto

logger = logging.getLogger(__name__)


def get_sequential_colors():
    return cm.ScalarMappable(None, "jet")


def get_qualitative_colors():
    return cm.ScalarMappable(None, "Accent")


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
    return result


def scale_image(image: np.ndarray, levels: Tuple[float, float]):
    channel_image = image - levels[0]
    channel_image /= levels[1] - levels[0]
    return np.clip(channel_image, 0, 1, out=channel_image)


def draw_mask(image: np.ndarray, mask_settings: MaskSettingsDto, heatmap_dict: Optional[dict] = None):
    mask = tifffile.imread(mask_settings.location)

    if mask_settings.gated:
        object_numbers = []
        for ann_object_numbers in mask_settings.cells.values():
            object_numbers.extend(ann_object_numbers)
        m = np.isin(mask, object_numbers)
        mask[~m] = 0

    if heatmap_dict:
        if mask_settings.gated:
            heatmap_dict = {k: heatmap_dict[k] for k in np.unique(mask) if k != 0}

        colors = heatmap_dict.values()
        # alpha = 0.3 if mask_settings.gated else 1
        img = label2rgb(
            label=mask,
            image=image,
            colors=colors,
            alpha=mask_settings.opacity,
            bg_label=0,
            image_alpha=1,
            kind="overlay",
        )
        return img
    else:
        boundary = find_boundaries(mask, connectivity=1, mode="inner")
        image[boundary > 0] = mask_settings.opacity
        return image


def draw_overlay(mask_settings: MaskSettingsDto):
    filename = mask_settings.location.replace("mask.tiff", "origin.png")
    overlay = io.imread(filename)
    return overlay


def draw_scalebar(image: np.ndarray, scalebar: ScalebarDto):
    image = img_as_ubyte(image)
    height, width, _ = image.shape
    length = 64
    thickness = 2
    if scalebar.settings is not None and "length" in scalebar.settings:
        new_length = scalebar.settings.get("length")
        if new_length is not None and new_length != "":
            length = int(new_length)
    cv2.line(
        image,
        (width - 60, height - 60),
        (width - 60 - length, height - 60),
        (255, 255, 255),
        thickness,
        cv2.LINE_4,
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


# def draw_legend(image: np.ndarray, legend_labels: List[Tuple[str, str, float]], legend: LegendDto):
#     for i, label in enumerate(legend_labels):
#         text = f"{label[0]} - {int(label[2])}" if legend.showIntensity else f"{label[0]}"
#         (label_width, label_height), baseline = cv2.getTextSize(text, cv2.FONT_HERSHEY_DUPLEX, legend.fontScale, 1)
#         cv2.rectangle(
#             image,
#             (5, (label_height + 20) * (i + 1) + 5),
#             (15 + label_width, (label_height + 20) * (i + 1) - label_height - 5),
#             (0, 0, 0),
#             cv2.FILLED,
#             cv2.LINE_AA,
#         )
#
#         b, g, r = tuple([255 * x for x in to_rgb(label[1])])
#         color = (r, g, b)
#         cv2.putText(
#             image,
#             text,
#             (10, (label_height + 20) * (i + 1)),
#             cv2.FONT_HERSHEY_DUPLEX,
#             legend.fontScale,
#             color,
#             1,
#             cv2.LINE_AA,
#         )
#     return image
