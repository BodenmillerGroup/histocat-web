import { IAnnotation } from "@/modules/annotations/models";

export interface IGateCreate {
  dataset_id: number;
  name: string;
  description?: string;
  cell_classes: { [name: string]: string };
  annotations: IAnnotation[];
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
  cell_classes: { [name: string]: string };
  annotations: IAnnotation[];
  created_at: string;
}
