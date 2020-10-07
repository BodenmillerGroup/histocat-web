export type PipelineStepType = "transformation" | "scale" | "regressOut" | "tsne" | "umap" | "pca";

export interface IPipelineCreate {
  project_id: number;
  name: string;
  description?: string;
  steps: any[];
}

export interface IPipelineUpdate {
  name: string;
  description?: string;
  steps: any[];
}

export interface IPipeline {
  id: number;
  project_id: number;
  name: string;
  description: string;
  steps: any[];
  created_at: string;
}
