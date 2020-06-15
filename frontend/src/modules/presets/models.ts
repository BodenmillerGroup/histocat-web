export interface IPresetCreate {
  experiment_id: number;
  name: string;
  description?: string;
  data?: object;
}

export interface IPreset {
  id: number;
  experiment_id: number;
  name: string;
  description: string;
  data: object;
  created_at: string;
}
