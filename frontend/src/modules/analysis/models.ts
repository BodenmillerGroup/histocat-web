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
