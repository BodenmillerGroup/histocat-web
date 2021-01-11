import numpy as np
import pandas as pd
import xarray as xr

from abc import ABC, abstractmethod
from pathlib import Path
from scipy.spatial import distance
from typing import Optional, Sequence, Union

from imctoolkit import utils

try:
    import anndata
except:
    anndata = None

try:
    import fcswrite
except:
    fcswrite = None


class SpatialSingleCellData(ABC):
    """Base class for in situ single-cell data"""

    @property
    @abstractmethod
    def num_channels(self) -> int:
        """Number of channels"""
        pass

    @property
    @abstractmethod
    def channel_names(self) -> np.ndarray:
        """Channel names"""
        pass

    @property
    @abstractmethod
    def num_cells(self) -> int:
        """Number of cells"""
        pass

    @property
    @abstractmethod
    def cell_ids(self) -> np.ndarray:
        """Integer cell IDs"""
        pass

    @property
    @abstractmethod
    def cell_centroids(self) -> xr.DataArray:
        """Cell centroids, shape: ``(cells, dimensions=2)``"""
        pass

    @abstractmethod
    def to_dataset(self, cell_properties: Union[bool, Sequence[str]] = False,
                   cell_channel_properties: Union[bool, Sequence[str]] = False) -> xr.Dataset:
        """Returns an :class:`xarray.Dataset` representation of the current instance

        At least one of :paramref:`cell_properties` or :paramref:`cell_channel_properties` has to be ``True``

        :param cell_properties: list of cell properties (e.g. regionprops) to include; set to ``True`` to include all
        :param cell_channel_properties: list of cell channel properties (e.g. intensity values) to include; set to
            ``True`` to include all
        :return: Dataset with cell ID, channel name and/or property name as coordinates
        """
        pass

    def to_dataframe(self, cell_properties: Union[bool, Sequence[str]] = False,
                     cell_channel_properties: Union[bool, Sequence[str]] = False) -> pd.DataFrame:
        """Returns a :class:`pandas.DataFrame` representation of the current instance

        Column names for cell channel properties (e.g. intensity values) are prefixed by the property name, e.g. a
        ``mean_intensities`` property would be included as ``'mean_intensities_{channel_name}'``. Column names for
        cell properties (e.g. regionprops) are not prefixed.

        :param cell_properties: list of cell properties (e.g. regionprops) to include; set to ``True`` to include all
        :param cell_channel_properties: list of cell channel properties (e.g. intensity values) to include; set to
            ``True`` to include all
        :return: DataFrame (index: cell IDs, columns: see description)
        """
        df = pd.DataFrame(index=pd.Index(self.cell_ids, name='cell'))
        if cell_properties:
            cell_property_dataset = self.to_dataset(cell_properties=cell_properties)
            for da in cell_property_dataset.data_vars.values():
                df = pd.merge(df, utils.to_table(da), left_index=True, right_index=True)
        if cell_channel_properties:
            cell_channel_property_dataset = self.to_dataset(cell_channel_properties=cell_channel_properties)
            for property_name, da in cell_channel_property_dataset.data_vars.items():
                df = pd.merge(df, utils.to_table(da).add_prefix(f'{property_name}_'), left_index=True, right_index=True)
        df.columns.name = 'feature'
        return df

    def to_anndata(self, cell_properties: Union[bool, Sequence[str]] = False,
                   cell_channel_properties: Union[bool, Sequence[str]] = False) -> 'anndata.AnnData':
        """Returns an :class:`anndata.AnnData` representation of the current instance

        :param cell_properties: list of cell properties (e.g. regionprops) to include; set to ``True`` to include all
        :param cell_channel_properties: list of cell channel properties (e.g. intensity values) to include; set to
            ``True`` to include all
        :return: AnnData object, in which cell channel properties (e.g. intensity values) are stored as layers and cell
            properties (e.g. regionprops) are stored as observations
        """
        if anndata is None:
            raise RuntimeError('anndata is not installed')
        obs_data = None
        if cell_properties:
            cell_property_dataset = self.to_dataset(cell_properties=cell_properties)
            obs_data = utils.to_table(xr.concat(cell_property_dataset.data_vars.values(), 'property'))
        layers = None
        if cell_channel_properties:
            cell_channel_property_dataset = self.to_dataset(cell_channel_properties=cell_channel_properties)
            layers = {property_name: da.values for property_name, da in cell_channel_property_dataset.data_vars.items()}
        return anndata.AnnData(
            obs=pd.DataFrame(index=pd.Index(data=self.cell_ids.astype(str), name='cell'), data=obs_data),
            var=pd.DataFrame(index=pd.Index(data=self.channel_names, name='channel')),
            layers=layers,
            shape=(self.num_cells, self.num_channels)
        )

    def compute_cell_centroid_distances(self, metric: str = 'euclidean') -> xr.DataArray:
        """Compute the pairwise distances between cell centroids

        :param metric: the distance metric to use, see
            https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.pdist.html
        :return: symmetric cell centroid distance matrix
        """
        dist_mat = distance.squareform(distance.pdist(self.cell_centroids.values, metric=metric))
        return xr.DataArray(data=dist_mat, dims=('cell_i', 'cell_j'),
                            coords={'cell_i': self.cell_ids, 'cell_j': self.cell_ids})

    def write_csv(self, path: Union[str, Path], **kwargs):
        """Writes a CSV file, see :func:`to_dataframe` for format specifications

        :param path: path to the .csv file to be written
        :param kwargs: other arguments passed to :func:`pandas.DataFrame.to_csv`
        """
        self.to_dataframe().to_csv(path, **kwargs)

    def write_fcs(self, path: Union[str, Path], **kwargs):
        """Writes an FCS file, see :func:`to_dataframe` for format specifications

        Uses :func:`fcswrite.write_fcs` for writing FCS 3.0 files.

        :param path: path to the .fcs file to be written
        :param kwargs: other arguments passed to :func:`fcswrite.write_fcs`
        """
        if fcswrite is None:
            raise RuntimeError('fcswrite is not installed')
        fcswrite.write_fcs(path, self.channel_names, self.to_dataframe().values, **kwargs)
