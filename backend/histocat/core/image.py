import logging
import sys
from typing import Optional, Sequence, Tuple

import cv2
import numpy as np
import tifffile
from mahotas import bwperim
from matplotlib import cm
from matplotlib.colors import LinearSegmentedColormap, rgb2hex, to_rgb
from skimage import img_as_ubyte
from skimage.color import label2rgb

from histocat.modules.acquisition.dto import FilterDto, MaskSettingsDto, ScalebarDto

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
    return result


def scale_image(image: np.ndarray, levels: Tuple[float, float]):
    channel_image = image - levels[0]
    channel_image /= levels[1] - levels[0]
    return np.clip(channel_image, 0, 1, out=channel_image)


def normalize_image(image: np.ndarray):
    """Normalize image to [0, 1] range"""
    return (image - np.min(image)) / (np.max(image) - np.min(image))


def draw_mask(image: np.ndarray, mask_settings: MaskSettingsDto, heatmap_dict: Optional[dict] = None):
    mask = tifffile.imread(mask_settings.location)
    if not mask_settings.colorize:
        if mask_settings.gated:
            m = np.isin(mask, mask_settings.cell_ids)
            mask[~m] = 0
            # if heatmap_dict:
            #     mask = replace_with_dict(mask, heatmap_dict)
            # nlabel = len(np.unique(mask))
            # colors = [tuple(map(tuple, np.random.rand(1, 3)))[0] for i in range(0, nlabel)]
            colors = ("darkorange", "darkorange")
            img = label2rgb(
                label=mask, image=image, colors=colors, alpha=0.3, bg_label=0, image_alpha=1, kind="overlay"
            )
            return img
        else:
            # if heatmap_dict:
            #     mask = replace_with_dict(mask, heatmap_dict)
            # nlabel = len(np.unique(mask))
            # colors = [tuple(map(tuple, np.random.rand(1, 3)))[0] for i in range(0, nlabel)]
            colors = ("darkorange", "darkorange")
            img = label2rgb(
                label=mask, image=image, colors=colors, alpha=0.3, bg_label=0, image_alpha=1, kind="overlay"
            )
            return img
    else:
        if mask_settings.gated:
            return mask_gated_img(image, mask, mask_settings.cell_ids)
        else:
            return mask_color_img(image, mask)


def mask_gated_img(image: np.ndarray, mask: np.ndarray, cell_ids: Sequence[int], color=(0, 146, 63, 100), alpha=0.8):
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
    # cv2.addWeighted(mask_layer, alpha, image, 1 - alpha, 0, image)
    # cv2.add(image, mask_layer, image)
    return np.clip(mask_layer + image, 0, 1)


def mask_color_img(image: np.ndarray, mask: np.ndarray, color=(0, 146, 63, 100), alpha=0.3):
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
    # cv2.addWeighted(mask_layer, alpha, image, 1 - alpha, 0, image)
    # cv2.add(image, mask_layer, image)
    return np.clip(mask_layer + image, 0, 1)


def draw_scalebar(image: np.ndarray, scalebar: ScalebarDto):
    image = img_as_ubyte(image)
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


def get_heatmap_colors(values: np.ndarray, categorical_values: bool):
    keys = values.unique().tolist()
    color_range = np.linspace(0, 1, len(keys), endpoint=False).tolist()
    colors = [cm.tab20b(x) for x in color_range]
    color_dict = dict(zip(keys, colors))
    result = []
    for v in values:
        result.append(rgb2hex(color_dict[v]))
    return result


def normalize_embedding(embedding):
    """Normalize embedding layout to meet client assumptions.
    Embedding is an ndarray, shape (n_obs, n)., where n is normally 2.
    """

    # scale isotropically
    min = embedding.min(axis=0)
    max = embedding.max(axis=0)
    scale = np.amax(max - min)
    normalized_layout = (embedding - min) / scale

    # translate to center on both axis
    # translate = 0.5 - ((max - min) / scale / 2)
    # TODO: improve translation
    translate = -1 + ((max - min) / scale / 2)
    normalized_layout = normalized_layout + translate

    normalized_layout = normalized_layout.astype(dtype=np.float32)
    return normalized_layout


def replace_with_dict(ar, dic):
    # Extract out keys and values
    k = np.array(list(dic.keys()))
    v = np.array(list(dic.values()))

    # Get argsort indices
    sidx = k.argsort()

    # Drop the magic bomb with searchsorted to get the corresponding
    # places for a in keys (using sorter since a is not necessarily sorted).
    # Then trace it back to original order with indexing into sidx
    # Finally index into values for desired output.
    return v[sidx[np.searchsorted(k, ar, sorter=sidx)]]
