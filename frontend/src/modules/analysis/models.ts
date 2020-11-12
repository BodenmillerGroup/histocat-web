export interface IPlotSeries {
  label: string;
  data: number[];
}

export interface IScatterPlotData {
  x: IPlotSeries;
  y: IPlotSeries;
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
