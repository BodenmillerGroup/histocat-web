import { IAnnotation } from "@/modules/annotations/models";

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

export interface IClassifyCellsData {
  gateId: number;
}

export interface IClassifyCellsSubmission {
  dataset_id: number;
  result_id?: number;
  cell_classes: { [name: string]: string };
  annotations: IAnnotation[];
}
