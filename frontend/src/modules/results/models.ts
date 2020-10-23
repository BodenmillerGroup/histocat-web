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
