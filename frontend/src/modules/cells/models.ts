export interface ICentroidsSubmission {
  datasetId: number;
}

export interface ICentroidsData {
  acquisitionIds: number[];
  cellIds: string[];
  objectNumbers: number[];
  x: number[];
  y: number[];
  colors: string[];
}

export interface ICell {
  cellId: string;
  objectNumber: number;
  acquisitionId: number;
  xy: [number, number];
  color: string;
  defaultColor: string;
  mappings: { [name: string]: [number, number] };
}

export interface IResult {
  id: number;
  dataset_id: number;
  status: string;
  name: string;
  description: string;
  pipeline: any;
  input: number[];
  output: any;
  location: string;
  created_at: string;
}

export interface IResultUpdate {
  name?: string | null;
  description?: string | null;
}

export interface IAxisData {
  label: string;
  data: number[];
}

export interface IResultDataMapping {
  x: IAxisData;
  y: IAxisData;
}

export interface IResultDataColors {
  type: string;
  name: string;
  data: string[];
}

export interface IRawResultData {
  cellIds: string[];
  markers: string[];
  mappings: { [key: string]: IResultDataMapping };
}

export interface IRawColorsData {
  cellIds: string[];
  colors: IResultDataColors;
}

export interface IPlotSeries {
  label: string;
  data: number[];
}

export interface IRawScatterData {
  cellIds: string[];
  x: IPlotSeries;
  y: IPlotSeries;
}
