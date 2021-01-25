import numpy as np
import pandas as pd
import xarray as xr

from typing import Collection, Sequence, Union

from imctoolkit.spatial_single_cell_data import SpatialSingleCellData

try:
    import networkx as nx
except:
    nx = None

try:
    import igraph
except:
    igraph = None


class SpatialCellGraph:
    """Spatial cell graph constructed from single-cell data

    :ivar data: single-cell data, as :class:`pandas.DataFrame`, with cell IDs as index and feature names as columns
    :ivar adj_mat: boolean adjacency matrix, as :class:`xarray.DataArray` with shape ``(cell_i, cell_j)`` and
        coordinates ``(cell IDs, cell IDs)``. A cell `j` is a neighbor of cell `i`, iff ``adj_mat[i, j] == True``.
    """

    def __init__(self, data, adj_mat, cell_properties: Union[bool, Sequence[str]] = False,
                 cell_channel_properties: Union[bool, Sequence[str]] = False, _skip_data_preparation: bool = False):
        """

        :param data: single-cell data (rows: cell IDs, columns: feature names)
        :type data: SingleCellData or DataFrame-like
        :param adj_mat: boolean adjacency matrix, shape: ``(cells, features)``
        :type adj_mat: DataArray-like
        :param cell_properties: list of cell properties (e.g. regionprops) to include as node attributes when using a
            SingleCellData object for :paramref:`data`; set to ``True`` to include all
        :param cell_channel_properties: list of cell channel properties (e.g. intensity values) to include as node
            attributes when using a SingleCellData object for :paramref:`data`; set to ``True`` to include all
        """
        if not _skip_data_preparation:
            data, adj_mat = self._prepare_data(data, adj_mat, cell_properties, cell_channel_properties)
        self.data = data
        self.adj_mat = adj_mat

    @staticmethod
    def _prepare_data(data, mat, cell_properties: Union[bool, Sequence[str]],
                      cell_channel_properties: Union[bool, Sequence[str]]):
        if isinstance(data, SpatialSingleCellData):
            data = data.to_dataframe(cell_properties=cell_properties, cell_channel_properties=cell_channel_properties)
        if not isinstance(data, pd.DataFrame):
            data = pd.DataFrame(data)
        if not isinstance(mat, xr.DataArray):
            mat = xr.DataArray(np.asarray(mat), dims=('cell_i', 'cell_j'))
        if not np.all(mat.coords['cell_i'].values == mat.coords['cell_j'].values):
            raise ValueError('Inconsistent row and column indices')
        mat_cells_in_data = np.in1d(mat.coords['cell_i'].values, data.index.values)
        if not np.all(mat_cells_in_data):
            missing_cell_ids = mat.coords['cell_i'][mat_cells_in_data].values.tolist()
            raise ValueError(f'Missing cell data for cell IDs in matrix: {missing_cell_ids}')
        return data, mat

    @property
    def num_cells(self) -> int:
        """Number of cells"""
        return len(self.cell_ids)

    @property
    def cell_ids(self) -> Collection[int]:
        """Cell IDs"""
        return self.data.index.tolist()

    @property
    def num_features(self) -> int:
        """Number of features"""
        return len(self.feature_names)

    @property
    def feature_names(self) -> np.ndarray:
        """Feature names"""
        return self.data.columns.values

    @property
    def is_undirected(self) -> bool:
        """``True``, if :attr:`adj_mat` is symmetric, ``False`` otherwise"""
        return np.allclose(self.adj_mat, self.adj_mat.transpose())

    def to_dataset(self) -> xr.Dataset:
        """Returns an :class:`xarray.Dataset` representation of the current instance

        :return: Dataset with :attr:`data` and :attr:`adj_mat` as members and ``(cell IDs, feature names)`` as
            coordinates
        """
        return xr.Dataset(data_vars={
            'data': xr.DataArray(self.data, dims=('cell', 'feature'),
                                 coords={'cell': self.cell_ids, 'feature': self.feature_names}),
            'adj_mat': self.adj_mat,
        })

    def to_networkx(self, weight_mat=None, create_using=None) -> 'nx.Graph':
        """Returns a :class:`networkx.Graph` representation of the current instance

        :param weight_mat: optional edge weight matrix, will be multiplied by :attr:`adj_mat`
        :param create_using: type of graph to create, defaults to :class:`networkx.Graph` for undirected graphs and
            :class:`networkx.DiGraph` for directed graphs when ``None``
        :type create_using: see :func:`networkx.from_numpy_array`
        :return: Graph or DiGraph with cell IDs as node labels and features as node attributes
        """
        if nx is None:
            raise RuntimeError('networkx is not installed')
        if create_using is None:
            create_using = nx.Graph if self.is_undirected else nx.DiGraph
        adj_mat = self.adj_mat.values
        if weight_mat is not None:
            adj_mat *= weight_mat
        graph: nx.Graph = nx.from_numpy_array(adj_mat, create_using=create_using)
        graph = nx.relabel_nodes(graph, mapping=dict(zip(graph, self.cell_ids)), copy=False)
        if weight_mat is None:
            for n1, n2, edge_attrs in graph.edges(data=True):
                edge_attrs.clear()
        node_attributes = {}
        for cell_id in self.cell_ids:
            node_attributes[cell_id] = self.data.loc[cell_id, :].to_dict()
        nx.set_node_attributes(graph, node_attributes)
        return graph

    def to_igraph(self, mode=None) -> 'igraph.Graph':
        """Returns an :class:`igraph.Graph` representation of the current instance

        :param mode: graph mode, defaults to :attr:`igraph.ADJ_UNDIRECTED` for undirected graphs and
            :attr:`igraph.ADJ_DIRECTED` for directed graphs when ``None``
        :type mode: see :func:`igraph.Graph.Adjacency`
        :return: Graph with cell IDs and features as vertex attributes
        """
        if igraph is None:
            raise RuntimeError('python-igraph is not installed')
        if mode is None:
            mode = igraph.ADJ_UNDIRECTED if self.is_undirected else igraph.ADJ_DIRECTED
        graph = igraph.Graph.Adjacency(self.adj_mat.values.tolist(), mode=mode)
        graph.vs['cell_id'] = self.cell_ids
        for feature_name in self.feature_names:
            graph.vs[feature_name] = self.data.loc[self.cell_ids, feature_name].values.tolist()
        return graph

    @staticmethod
    def load_dataset(dataset: xr.Dataset) -> 'SpatialCellGraph':
        """Creates a new :class:`SpatialCellGraph` from its dataset representation, see

        :param dataset: Dataset of the same format as created by :func:`to_dataset`
        :return: a new :class:`SpatialCellGraph` instance
        """
        data = pd.DataFrame(dataset['data'], index=pd.Index(data=dataset['data'].coords['cell'].values, name='cell'),
                            columns=pd.Index(data=dataset['data'].coords['feature'].values, name='feature'))
        return SpatialCellGraph(data, dataset['adj_mat'])

    @classmethod
    def construct_knn_graph(cls, data, dist_mat, k: int, cell_properties: Union[bool, Sequence[str]] = False,
                            cell_channel_properties: Union[bool, Sequence[str]] = False) -> 'SpatialCellGraph':
        """Constructs a new k-nearest cell neighbor graph

        :param data: single-cell data (rows: cell IDs, columns: feature names)
        :type data: SingleCellData or DataFrame-like
        :param dist_mat: symmetric distance matrix, shape: ``(cells, cells)``
        :type dist_mat: DataArray-like
        :param k: number of nearest neighbors for the graph construction
        :param cell_properties: list of cell properties (e.g. regionprops) to include as node attributes; set to
            ``True`` to include all
        :param cell_channel_properties: list of cell channel properties (e.g. intensity values) to include  as node
            attributes; set to ``True`` to include all
        :return: a directed k-nearest cell neighbor graph
        """
        data, dist_mat = cls._prepare_data(data, dist_mat, cell_properties, cell_channel_properties)
        adj_mat = xr.zeros_like(dist_mat, dtype='bool')
        knn_indices = np.argpartition(dist_mat.values, k + 1, axis=1)[:, :(k + 1)]
        for current_index, current_knn_indices in enumerate(knn_indices):
            adj_mat[current_index, current_knn_indices] = True
        np.fill_diagonal(adj_mat.values, False)
        return SpatialCellGraph(data, adj_mat, _skip_data_preparation=True)

    @classmethod
    def construct_dist_graph(cls, data, dist_mat, dist_thres: float,
                             cell_properties: Union[bool, Sequence[str]] = False,
                             cell_channel_properties: Union[bool, Sequence[str]] = False) -> 'SpatialCellGraph':
        """Constructs a new cell neighborhood graph by distance thresholding

        :param data: single-cell data (rows: cell IDs, columns: feature names)
        :type data: SingleCellData or DataFrame-like
        :param dist_mat: symmetric distance matrix, shape: ``(cells, cells)``
        :type dist_mat: DataArray-like
        :param dist_thres: distance hot_pixel_thres, (strictly) below which cells are considered neighbors
        :param cell_properties: list of cell properties (e.g. regionprops) to include as node attributes; set to
            ``True`` to include all
        :param cell_channel_properties: list of cell channel properties (e.g. intensity values) to include  as node
            attributes; set to ``True`` to include all
        :return: an undirected cell neighborhood graph
        """
        data, dist_mat = cls._prepare_data(data, dist_mat, cell_properties, cell_channel_properties)
        adj_mat = xr.DataArray(dist_mat < dist_thres)
        np.fill_diagonal(adj_mat.values, False)
        return SpatialCellGraph(data, adj_mat, _skip_data_preparation=True)
