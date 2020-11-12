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
  data: any[];
}

export interface IRawResultData {
  cellIds: string[];
  acquisitionIds: number[];
  objectNumbers: number[];
  x: number[];
  y: number[];
  mappings: { [key: string]: IResultDataMapping };
  colors?: IResultDataColors;
}

export interface ICellData {
  cellId: string;
  acquisitionId: number;
  objectNumber: number;
  x: number;
  y: number;
  mappings?: {
    [key: string]: {
      x: number;
      y: number;
    };
  };
  color?: any;
}

export interface ICellPoint {
  cellId: string;
  acquisitionId: number;
  objectNumber: number;
  x: number;
  y: number;
  color: any;
}

export interface ISelectedCell {
  cellId: string;
  acquisitionId: number;
  objectNumber: number;
}
