export interface IModelCreate {
  name: string;
  description: string | null;
}

export interface IModelUpdate {
  name: string;
  description: string | null;
}

export interface IModel {
  id: number;
  group_id: number;
  name: string;
  description: string;
  location: string;
  meta: object;
  created_at: string;
}
