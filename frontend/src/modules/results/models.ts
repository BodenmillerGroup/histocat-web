export interface IResult {
  id: number;
  dataset_id: number;
  status: string;
  name: string;
  description: string;
  pipeline: any;
  input: number[];
  location: string;
  created_at: string;
}
