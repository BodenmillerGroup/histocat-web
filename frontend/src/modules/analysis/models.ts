import { IImageFilter, IImageScalebar } from "@/modules/settings/models";

export type ImageResultType = "origin" | "mask";

export interface IImageSegmentationSettings {
  algorithm: string;
  iterations: number;
  kernel_size: number;
  mask_color: string;
  result_type: ImageResultType;
}

export interface IImageSegmentationSubmission {
  format?: string;
  filter: IImageFilter;
  scalebar: IImageScalebar;
  channels: Array<{
    id: number;
    color?: string;
    customLabel?: string;
    min?: number;
    max?: number;
  }>;
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
  acquisition_id: number;
  n_components: number;
  markers: string[];
  heatmap: string;
}

export interface IPCAData {
  x: IPlotSeries;
  y: IPlotSeries;
  z?: IPlotSeries;
  heatmap?: IPlotSeries;
  explained_variance_ratio?: number[];
}

export interface ITSNESubmission {
  dataset_id: number;
  acquisition_id: number;
  markers: string[];
  n_components: number;
  perplexity: number;
  learning_rate: number;
  iterations: number;
}

export interface ITSNEData {
  x: IPlotSeries;
  y: IPlotSeries;
  z?: IPlotSeries;
  heatmap?: IPlotSeries;
}

export interface IUMAPSubmission {
  dataset_id: number;
  acquisition_id: number;
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
