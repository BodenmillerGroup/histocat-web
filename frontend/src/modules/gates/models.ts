export interface IGateCreate {
  dataset_id: number;
  name: string;
  description?: string;
  acquisition_ids: number[];
  indices: number[];
  cell_ids: number[];
}

export interface IGateUpdate {
  name?: string | null;
  description?: string | null;
}

export interface IGate {
  id: number;
  dataset_id: number;
  name: string;
  description?: string;
  acquisition_ids: number[];
  indices: number[];
  cell_ids: number[];
  created_at: string;
}
