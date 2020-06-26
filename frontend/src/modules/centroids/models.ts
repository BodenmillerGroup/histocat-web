export interface ICentroidsSubmission {
  datasetId: number;
}

export interface ICentroidsData {
  acquisitionIds: number[];
  cellIds: number[];
  objectNumbers: number[];
  x: number[];
  y: number[];
}
