export interface IModelUpdate {
  name: string;
  description: string | null;
}

export interface IModel {
  id: number;
  name: string;
  application: string;
  description: string;
  location: string;
  meta: object;
  created_at: string;
}
