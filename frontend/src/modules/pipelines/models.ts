export type PipelineStepType = "markersFilter" | "transformation" | "scale" | "neighbors" | "tsne" | "umap" | "pca" | "leiden" | "louvain";

export interface IPipelineCreate {
  project_id: number;
  name: string;
  description?: string;
  steps: any[];
}

export interface IPipelineUpdate {
  name?: string | null;
  description?: string | null;
}

export interface IPipeline {
  id: number;
  project_id: number;
  name: string;
  description: string;
  steps: any[];
  created_at: string;
}

export interface IProcessPipeline {
  dataset_id: number;
  acquisition_ids: number[];
  steps: any[];
}
