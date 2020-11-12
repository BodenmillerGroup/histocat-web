export interface ICentroidsSubmission {
  datasetId: number;
}

export interface ICentroidsData {
  acquisitionIds: number[];
  cellIds: string[];
  objectNumbers: number[];
  x: number[];
  y: number[];
}
