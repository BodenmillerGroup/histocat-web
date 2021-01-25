import numpy as np
import pandas as pd
import tifffile
import xarray as xr
import xtiff

from imctools.io.mcd.mcdparser import McdParser
from imctools.io.txt.txtparser import TxtParser
from pathlib import Path
from typing import Optional, Sequence, Union
from xml.etree import ElementTree


class MultichannelImage:
    """Functionality for processing multi-channel images

    :ivar data: raw image data, as :class:`xarray.DataArray` with dimensions ``(c, y, x)`` and, optionally, channel
        names as coordinate for the ``c`` dimension
    """

    def __init__(self, data, channel_names: Optional[Sequence[str]] = None):
        """

        :param data: raw image data, shape: ``(c, y, x)``
        :type data: DataArray-like
        :param channel_names: channel names
        """
        if not isinstance(data, xr.DataArray):
            data = xr.DataArray(data=data, dims=('c', 'y', 'x'))
        if data.dims != ('c', 'y', 'x'):
            raise ValueError(f'Invalid image dimensions: expected ("c", "y", "x"), got {data.dims}')
        if channel_names is not None:
            data.coords['c'] = channel_names
        self.data = data

    @property
    def width(self) -> int:
        """Image width in pixels"""
        return self.data.sizes['x']

    @property
    def height(self) -> int:
        """Image height in pixels"""
        return self.data.sizes['y']

    @property
    def num_channels(self) -> int:
        """Number of channels"""
        return self.data.sizes['c']

    @property
    def channel_names(self) -> np.ndarray:
        """Channel names"""
        return self.data.coords['c'].values

    def copy(self) -> 'MultichannelImage':
        """Creates a copy of the current instance

        :return: a deep copy of the current instance
        """
        return MultichannelImage(self.data.copy(deep=True))

    def write_ome_tiff(self, path: Union[str, Path], **kwargs):
        """Writes an OME-TIFF file using :func:`xtiff.to_tiff`

        Note that the written OME-TIFF file will not be identical to any OME-TIFF file from which the current instance
        was potentially created initially.

        :param path: path to the .ome.tiff file to be written
        :param kwargs: other arguments passed to :func:`xtiff.to_tiff`
        """
        xtiff.to_tiff(self.data, path, **kwargs)

    @classmethod
    def read(cls, path: Union[str, Path], **kwargs) -> 'MultichannelImage':
        """Creates a new :class:`MultichannelImage` from the specified file

        The file type is determined based on the file extension.

        :param path: path to the file
        :param kwargs: optional keyword arguments passed to the respective ``read_*`` function
        :return: a new :class:`MultichannelImage` instance
        """
        path_suffix = Path(path).suffix
        if path_suffix == '.txt':
            return cls.read_imc_txt(path, **kwargs)
        elif path_suffix == '.mcd':
            return cls.read_imc_mcd(path, **kwargs)
        elif path_suffix == '.tif' or path_suffix == '.tiff':
            return cls.read_tiff(path, **kwargs)
        raise ValueError(f'Unsupported file extension: {path_suffix}')

    @staticmethod
    def read_imc_txt(path: Union[str, Path], channel_names_attr: str = 'channel_names') -> 'MultichannelImage':
        """Creates a new :class:`MultichannelImage` from the specified Fluidigm(TM) TXT file

        Uses :class:`imctools.io.txt.txtparser.TxtParser` for reading .txt files.

        :param path: path to the .txt file
        :param channel_names_attr: :class:`imctools.data.AcquisitionData` attribute from which the channel names will be
            taken, e.g. ``'channel_labels'``
        :return: a new :class:`MultichannelImage` instance
        """
        path = Path(path)
        parser = TxtParser(path)
        acquisition_data = parser.get_acquisition_data()
        img_data = xr.DataArray(data=acquisition_data.image_data, dims=('c', 'y', 'x'), name=path.name)
        return MultichannelImage(img_data, channel_names=getattr(acquisition_data, channel_names_attr))

    @staticmethod
    def read_imc_mcd(path: Union[str, Path], acquisition_id: int,
                     channel_names_attr: str = 'channel_names') -> 'MultichannelImage':
        """Creates a new :class:`MultichannelImage` from the specified Fluidigm(TM) MCD file

        Uses :class:`imctools.io.txt.mcdparser.McdParser` for reading .mcd files.

        :param path: path to the .mcd file
        :param acquisition_id: acquisition ID to read (unique across slides)
        :param channel_names_attr: :class:`imctools.data.AcquisitionData` attribute from which the channel names will be
            taken, e.g. ``'channel_labels'``
        :return: a new :class:`MultichannelImage` instance
        """
        path = Path(path)
        parser = McdParser(path)
        acquisition_data = parser.get_acquisition_data(acquisition_id)
        img_data = xr.DataArray(data=acquisition_data.image_data, dims=('c', 'y', 'x'), name=path.name)
        return MultichannelImage(img_data, channel_names=getattr(acquisition_data, channel_names_attr))

    @staticmethod
    def read_tiff(path, panel=None, panel_channel_index_col: Optional[str] = None,
                  panel_channel_name_col: str = 'channel_name', channel_names: Optional[Sequence[str]] = None,
                  ome_channel_name_attr: str = 'Name') -> 'MultichannelImage':
        """Creates a new :class:`MultichannelImage` from the specified TIFF/OME-TIFF file

        Uses :class:`tifffile.TiffFile` for reading .tiff files. When reading an OME-TIFF file and :paramref:`panel` is
        not specified, the channel names are taken from the embedded OME-XML.

        :param path: path to the .tiff file
        :type path: tifffile.TiffFile path-argument
        :param panel: panel that maps channel indices to channel names. If specified, overrides channel names
            extracted from OME-XML and acts as a channel selector.
        :type panel: None, pandas.DataFrame or pandas.read_csv path-argument
        :param panel_channel_index_col: channel index column in panel, or None for defaulting to `range(n_channels)`
        :param panel_channel_name_col: channel name column in panel
        :param channel_names: channel names, matching the number of image channels. If specified for OME-TIFFs or
            together with :paramref:`panel`, acts as a channel selector.
        :param ome_channel_name_attr: name of the OME-XML Channel element attribute from which the channel names will
            be taken, see https://www.openmicroscopy.org/Schemas/Documentation/Generated/OME-2016-06/ome.html.
        :return: a new :class:`MultichannelImage` instance
        """
        with tifffile.TiffFile(path) as tiff:
            img_data = xr.DataArray(data=tiff.asarray().squeeze(), dims=('c', 'y', 'x'))
            ome_metadata = tiff.ome_metadata
        if panel is not None:
            if not isinstance(panel, pd.DataFrame):
                panel = pd.read_csv(panel)
            if panel_channel_index_col is not None and panel_channel_index_col not in panel.columns:
                raise ValueError(f'Column {panel_channel_index_col} not found in panel')
            if panel_channel_name_col not in panel.columns:
                raise ValueError(f'Column {panel_channel_name_col} not found in panel')
            if panel_channel_index_col is not None:
                panel = panel.loc[panel[panel_channel_index_col].notna(), :]
                img_data = img_data[panel[panel_channel_index_col].astype(int).values, :, :]
            img_data.coords['c'] = panel[panel_channel_name_col].astype(str).values
        elif ome_metadata is not None:
            element_tree = ElementTree.fromstring(ome_metadata)
            ome_namespaces = {'ome': 'http://www.openmicroscopy.org/Schemas/OME/2016-06'}
            channel_elems = element_tree.findall('./ome:Image/ome:Pixels/ome:Channel', namespaces=ome_namespaces)
            if len(channel_elems) == img_data.sizes['c']:
                channel_elems.sort(key=lambda channel_elem: int(channel_elem.attrib['ID'].rpartition(':')[-1]))
                img_data.coords['c'] = [channel_elem.attrib[ome_channel_name_attr] for channel_elem in channel_elems]
        elif channel_names is not None:
            img_data.coords['c'] = channel_names
        if channel_names is not None:
            img_data = img_data.loc[channel_names, :, :]
        return MultichannelImage(img_data)
