export interface IDatasetTSNEOutput {
  name: string;
  location: string;
  params: any;
}

export interface IDatasetUMAPOutput {
  name: string;
  location: string;
  params: any;
}

export interface IDatasetPhenoGraphOutput {
  name: string;
  location: string;
  params: any;
}

export interface IDatasetUpdate {
  name?: string | null;
  description?: string | null;
}

export interface IDataset {
  id: number;
  project_id: number;
  user_id: number;
  uid: string;
  name: string;
  description: string;
  status: string;
  output?: {
    tsne: { [name: string]: IDatasetTSNEOutput };
    umap: { [name: string]: IDatasetUMAPOutput };
    phenograph: { [name: string]: IDatasetPhenoGraphOutput };
  };
  meta: {
    acquisition_metadata?: {
      location: string;
    };
    cell?: {
      location: string;
    };
    image?: {
      location: string;
    };
    object_relationships?: {
      location: string;
    };
    probability_masks?: {
      [id: number]: {
        location: string;
        slide: {
          id: number;
          origin_id: number;
        };
        panorama: {
          id: number;
          origin_id: number;
        };
        roi: {
          id: number;
          origin_id: number;
        };
        acquisition: {
          id: number;
          origin_id: number;
        };
      };
    };
  };
  location: string;
  created_at: string;
}
