export interface IPlotSeries {
  label: string;
  data: number[];
}

export interface IScatterPlotData {
  x: IPlotSeries;
  y: IPlotSeries;
  z?: IPlotSeries;
  heatmap?: IPlotSeries;
}

export interface IChart3DData {
  x: IPlotSeries;
  y: IPlotSeries;
  z?: IPlotSeries;
  heatmap?: IPlotSeries;
}

export interface IChart2DData {
  acquisitionIds: number[];
  cellIds: number[];
  objectNumbers: number[];
  x: IPlotSeries;
  y: IPlotSeries;
  z?: IPlotSeries;
  heatmap?: IPlotSeries;
}

export interface IPCAData {
  cell_ids: any[];
  x: IPlotSeries;
  y: IPlotSeries;
  z?: IPlotSeries;
  heatmap?: IPlotSeries;
  explained_variance_ratio?: number[];
}

export interface ITSNESubmission {
  dataset_id: number;
  acquisition_ids: number[];
  markers: string[];
  n_components: number;
  perplexity: number;
  learning_rate: number;
  iterations: number;
  theta: number;
  init: string;
}

export interface ITSNEData {
  x: IPlotSeries;
  y: IPlotSeries;
  z?: IPlotSeries;
  heatmap?: IPlotSeries;
}

export interface IUMAPSubmission {
  dataset_id: number;
  acquisition_ids: number[];
  markers: string[];
  n_components: number;
  n_neighbors: number;
  min_dist: number;
  metric: string;
}

export interface IUMAPData {
  x: IPlotSeries;
  y: IPlotSeries;
  z?: IPlotSeries;
  heatmap?: IPlotSeries;
}

export interface IPhenoGraphSubmission {
  dataset_id: number;
  acquisition_ids: number[];
  markers: string[];
  clustering_algo: string;
  nearest_neighbors: number;
  jaccard: boolean;
  primary_metric: string;
  min_cluster_size: number;
}

export interface IPhenoGraphData {
  community: IPlotSeries;
}

export interface IRegionChannelData {
  metal: string;
  min: number;
  max: number;
  mean: number;
}

export interface IRegionStatsSubmission {
  project_id: number;
  acquisition_id: number;
  region_polygon: any[];
}
