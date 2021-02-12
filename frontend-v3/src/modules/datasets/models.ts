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
  origin: string;
  acquisition_ids: number[];
  channels: string[];
  status: string;
  output?: {
    tsne: { [name: string]: IDatasetTSNEOutput };
    umap: { [name: string]: IDatasetUMAPOutput };
    phenograph: { [name: string]: IDatasetPhenoGraphOutput };
  };
  meta: {
    masks?: {
      [id: number]: {
        location: string;
        slide: {
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
