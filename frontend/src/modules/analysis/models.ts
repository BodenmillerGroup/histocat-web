import { IImageFilter, IImageScalebar } from "@/modules/settings/models";

export type ImageResultType = "origin" | "mask";

export interface IImageSegmentationSettings {
  algorithm: string;
  iterations: number;
  kernel_size: number;
  mask_color: string;
  result_type: ImageResultType;
  suppressBroadcast?: boolean;
}

export interface IImageSegmentationSubmission {
  format?: string;
  filter: IImageFilter;
  scalebar: IImageScalebar;
  channels: {
    name: string;
    color?: string;
    customLabel?: string;
    min?: number;
    max?: number;
  }[];
  settings: IImageSegmentationSettings;
}

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

export interface IPCASubmission {
  dataset_id: number;
  acquisition_ids: number[];
  n_components: number;
  markers: string[];
  heatmapType: string;
  heatmap: string;
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
  experiment_id: number;
  acquisition_id: number;
  region_polygon: any[];
}

export interface ICentroidsSubmission {
  dataset_id: number;
  acquisition_id: number;
}

export interface ICentroidsData {
  cell_ids: any[];
  coords: any[];
}
