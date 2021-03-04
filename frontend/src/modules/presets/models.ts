export interface IPresetCreate {
  project_id: number;
  name: string;
  description?: string;
  data?: object;
}

export interface IPreset {
  id: number;
  project_id: number;
  name: string;
  description: string;
  data: object;
  created_at: string;
}
