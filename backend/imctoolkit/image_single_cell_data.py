import numpy as np
import pandas as pd
import tifffile
import xarray as xr

from enum import Enum
from functools import cached_property
from pathlib import Path
from scipy.ndimage import distance_transform_edt
from skimage import measure
from typing import Any, Callable, Optional, Sequence, Union

from imctoolkit import utils
from imctoolkit.multichannel_image import MultichannelImage
from imctoolkit.spatial_single_cell_data import SpatialSingleCellData


class ImageSingleCellData(SpatialSingleCellData):
    """Single-cell data accessor for multi-channel images

    :ivar img: intensity image, as :class:`xarray.DataArray` with dimensions ``(c, y, x)`` and, optionally, channel
        names as coordinate for the ``c`` dimension
    :ivar mask: cell mask, as :class:`numpy.ndarray` of shape ``(y, x)``, where ``0`` indicates background pixels and
        non-zero pixels indicate the cell ID
    :ivar region_properties: list of :class:`RegionProperties` computed by the current instance, see :func:`regionprops`
    """

    class RegionProperties(Enum):
        """Enumeration of regionprops properties supported by this class, see
        https://scikit-image.org/docs/0.17.x/api/skimage.measure.html#skimage.measure.regionprops
        """
        AREA = 'area'
        BBOX = 'bbox'
        BBOX_AREA = 'bbox_area'
        CONVEX_AREA = 'convex_area'
        CONVEX_IMAGE = 'convex_image'
        COORDS = 'coords'
        ECCENTRICITY = 'eccentricity'
        EQUIVALENT_DIAMETER = 'equivalent_diameter'
        EULER_NUMBER = 'euler_number'
        EXTENT = 'extent'
        FILLED_AREA = 'filled_area'
        FILLED_IMAGE = 'filled_image'
        IMAGE = 'image'
        INERTIA_TENSOR = 'inertia_tensor'
        INERTIA_TENSOR_EIGVALS = 'inertia_tensor_eigvals'
        LOCAL_CENTROID = 'local_centroid'
        MAJOR_AXIS_LENGTH = 'major_axis_length'
        MINOR_AXIS_LENGTH = 'minor_axis_length'
        MOMENTS = 'moments'
        MOMENTS_CENTRAL = 'moments_central'
        MOMENTS_HU = 'moments_hu'
        MOMENTS_NORMALIZED = 'moments_normalized'
        ORIENTATION = 'orientation'
        PERIMETER = 'perimeter'
        SLICE = 'slice'
        SOLIDITY = 'solidity'

    DEFAULT_REGION_PROPERTIES = [
        RegionProperties.AREA,
        RegionProperties.ECCENTRICITY,
        RegionProperties.MAJOR_AXIS_LENGTH,
        RegionProperties.MINOR_AXIS_LENGTH,
        RegionProperties.ORIENTATION,
    ]  #: List of :class:`RegionProperties` computed by default, see :func:`__init__`

    _REGIONPROPS_CENTROID_COLUMNS = ['centroid-0', 'centroid-1']

    def __init__(self, img, mask, channel_names: Optional[Sequence[str]] = None,
                 region_properties: Optional[Sequence[RegionProperties]] = None):
        """

        :param img: intensity image, shape: ``(c, y, x)``
        :type img: MultichannelImage or DataArray-like
        :param mask: (path to) cell mask of shape: ``(y, x)``
        :type mask: mask file path or array-like
        :param channel_names: channel names
        :param region_properties: list of :class:`RegionProperties` to compute, defaults to
            :attr:`DEFAULT_REGION_PROPERTIES` when ``None``
        """
        super(ImageSingleCellData, self).__init__()
        if region_properties is None:
            region_properties = self.DEFAULT_REGION_PROPERTIES
        if isinstance(img, MultichannelImage):
            img = img.data
        if not isinstance(img, xr.DataArray):
            img = xr.DataArray(data=img, dims=('c', 'y', 'x'))
        if isinstance(mask, str) or isinstance(mask, Path):
            mask = tifffile.imread(mask).squeeze()
        mask = np.asarray(mask)
        if img.dims != ('c', 'y', 'x'):
            raise ValueError(f'Invalid image dimensions: expected ("c", "y", "x"), got {img.dims}')
        if channel_names is not None:
            img.coords['c'] = channel_names
        if mask.shape != img.shape[1:]:
            raise ValueError(f'Inconsistent mask {mask.shape} and image {img.shape[1:]} shapes')
        self.img = img
        self.mask = mask
        self.region_properties = list(region_properties)
        self._cell_ids = np.unique(mask[mask != 0])

    @property
    def image_width(self) -> int:
        """Image width in pixels"""
        return self.img.sizes['x']

    @property
    def image_height(self) -> int:
        """Image height in pixels"""
        return self.img.sizes['y']

    @property
    def num_channels(self) -> int:
        return self.img.sizes['c']

    @property
    def channel_names(self) -> np.ndarray:
        return self.img.coords['c'].values

    @property
    def num_cells(self) -> int:
        return len(self.cell_ids)

    @property
    def cell_ids(self) -> np.ndarray:
        return self._cell_ids

    @property
    def cell_centroids(self) -> xr.DataArray:
        return self._regionprops_with_centroids.loc[:, self._REGIONPROPS_CENTROID_COLUMNS]

    @cached_property
    def min_intensities(self) -> xr.DataArray:
        """Minimum cell intensities

        :return: DataArray with coordinates ``(cell IDs, channel names)``
        """
        return self.compute_cell_intensities(np.nanmin)

    @property
    def min_intensities_table(self) -> pd.DataFrame:
        """Minimum cell intensities

        :return: DataFrame (index: cell IDs, columns: channel names)
        """
        return utils.to_table(self.min_intensities)

    @cached_property
    def max_intensities(self) -> xr.DataArray:
        """Maximum cell intensities

        :return: DataArray with coordinates ``(cell IDs, channel names)``
        """
        return self.compute_cell_intensities(np.nanmax)

    @property
    def max_intensities_table(self) -> pd.DataFrame:
        """Maximum cell intensities

        :return: DataFrame (index: cell IDs, columns: channel names)
        """
        return utils.to_table(self.max_intensities)

    @cached_property
    def mean_intensities(self) -> xr.DataArray:
        """Mean cell intensities

        :return: DataArray with coordinates ``(cell IDs, channel names)``
        """
        return self.compute_cell_intensities(np.nanmean)

    @property
    def mean_intensities_table(self) -> pd.DataFrame:
        """Mean cell intensities

        :return: DataFrame (index: cell IDs, columns: channel names)
        """
        return utils.to_table(self.mean_intensities)

    @cached_property
    def median_intensities(self) -> xr.DataArray:
        """Median cell intensities

        :return: DataArray with coordinates ``(cell IDs, channel names)``
        """
        return self.compute_cell_intensities(np.nanmedian)

    @property
    def median_intensities_table(self) -> pd.DataFrame:
        """Median cell intensities

        :return: DataFrame (index: cell IDs, columns: channel names)
        """
        return utils.to_table(self.median_intensities)

    @cached_property
    def std_intensities(self) -> xr.DataArray:
        """Standard deviations of cell intensities

        :return: DataArray with coordinates ``(cell IDs, channel names)``
        """
        return self.compute_cell_intensities(np.nanstd)

    @property
    def std_intensities_table(self) -> pd.DataFrame:
        """Standard deviations of cell intensities

        :return: DataFrame (index: cell IDs, columns: channel names)
        """
        return utils.to_table(self.std_intensities)

    @cached_property
    def var_intensities(self) -> xr.DataArray:
        """Variances of cell intensities

        :return: DataArray with coordinates ``(cell IDs, channel names)``
        """
        return self.compute_cell_intensities(np.nanvar)

    @property
    def var_intensities_table(self) -> pd.DataFrame:
        """Variances of cell intensities

        :return: DataFrame (index: cell IDs, columns: channel names)
        """
        return utils.to_table(self.var_intensities)

    @cached_property
    def _regionprops_with_centroids(self) -> xr.DataArray:
        regionprops_properties = ['label', 'centroid'] + [rp.value for rp in self.region_properties]
        regionprops_dict = measure.regionprops_table(self.mask, properties=regionprops_properties)
        df = pd.DataFrame(regionprops_dict, index=regionprops_dict.pop('label'))
        return xr.DataArray(data=df, dims=('cell', 'property'))

    @property
    def regionprops(self) -> xr.DataArray:
        """Region properties

        For a list of computed properties, see :attr:`region_properties`.

        :return: DataArray with coordinates ``(cell_id, property_name)``
        """
        return self._regionprops_with_centroids.drop_sel(property=self._REGIONPROPS_CENTROID_COLUMNS)

    @property
    def regionprops_table(self) -> pd.DataFrame:
        """Region properties

        For a list of computed properties, see :attr:`region_properties`.

        :return: DataFrame (index: cell IDs, columns: regionprops property names)
        """
        return utils.to_table(self.regionprops)

    def compute_cell_intensities(self, aggr: Callable[[np.ndarray], Any]) -> xr.DataArray:
        """Compute cell intensity values

        :param aggr: function for aggregating the pixel values of a cell
        :return: DataArray with coordinates ``(cell IDs, channel names)``
        """
        arr = xr.DataArray(dims=('cell', 'channel'), coords={'cell': self.cell_ids, 'channel': self.channel_names})
        for channel_name in self.channel_names:
            channel_img = self.img.loc[channel_name].values
            arr.loc[:, channel_name] = [aggr(channel_img[self.mask == cell_id]) for cell_id in self.cell_ids]
        return arr

    def compute_cell_border_distances(self) -> xr.DataArray:
        """Compute the pairwise Euclidean distances between cell borders

        :return: symmetric border distance matrix
        """
        # TODO speed up computation, e.g. by only computing distances between pixels belonging to cells
        dist_mat = np.zeros((self.num_cells, self.num_cells))
        cell_masks = [self.mask == cell_id for cell_id in self.cell_ids]
        for i, i_id in enumerate(self.cell_ids[:-1]):
            i_dist = distance_transform_edt(self.mask != i_id)
            dist_mat[i, (i + 1):] = [np.amin(i_dist[cell_masks[j]]) for j in range(i + 1, self.num_cells)]
        dist_mat += dist_mat.transpose()
        return xr.DataArray(data=dist_mat, dims=('cell_i', 'cell_j'),
                            coords={'cell_i': self.cell_ids, 'cell_j': self.cell_ids})

    def to_dataset(self, cell_properties: Union[bool, Sequence[str]] = False,
                   cell_channel_properties: Union[bool, Sequence[str]] = False) -> xr.Dataset:
        if not cell_properties and not cell_channel_properties:
            raise ValueError('At least one of cell_properties, cell_channel_properties must be specified')
        data_vars = {}
        if cell_properties:
            if isinstance(cell_properties, Sequence):
                for cell_property in cell_properties:
                    data_vars[cell_property] = getattr(self, cell_property)
            else:
                data_vars = {
                    **data_vars,
                    'regionprops': self.regionprops,
                }
        if cell_channel_properties:
            if isinstance(cell_channel_properties, Sequence):
                for cell_channel_property in cell_channel_properties:
                    data_vars[cell_channel_property] = getattr(self, cell_channel_property)
            else:
                data_vars = {
                    **data_vars,
                    'min_intensities': self.min_intensities,
                    'max_intensities': self.max_intensities,
                    'mean_intensities': self.mean_intensities,
                    'median_intensities': self.median_intensities,
                    'std_intensities': self.std_intensities,
                    'var_intensities': self.var_intensities,
                }
        return xr.Dataset(data_vars=data_vars)
