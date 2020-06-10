export interface IGateCreate {
  experiment_id: number;
  name?: string;
  description?: string;
  data?: object;
}

export interface IGate {
  id: number;
  experiment_id: number;
  name: string;
  description: string;
  data: object;
  created_at: string;
}
