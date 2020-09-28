export interface IResult {
  id: number;
  dataset_id: number;
  parent_id?: number;
  type: string;
  name: string;
  description: string;
  status: string;
  params: any;
  location: string;
  created_at: string;
}
