export interface IExperiment {
  id: number;
  name: string;
  description: string;
  root_directory: string;
  location: string;
  meta: object;
  created_at: string;
}

export interface IExperimentUpdate {
  name?: string;
  description?: string;
}

export interface IExperimentCreate {
  name: string;
  description?: string;
}

export interface IExperimentDataset {
  id: number;
  name: string;
  description: string;
  root_directory: string;
  location: string;
  meta: object;
  created_at: string;
  slides: ISlide[];
}

export interface ISlide {
  id: number;
  experiment_id: number;
  name: string;
  description: string;
  location: string;
  meta: object;
  created_at: string;
  acquisitions: IAcquisition[];
}

export interface IAcquisition {
  id: number;
  slide_id: number;
  name: string;
  description: string;
  location: string;
  width: number;
  height: number;
  meta: object;
  created_at: string;
  channels: IChannel[];
}

export interface IChannel {
  id: number;
  acquisition_id: number;
  name: string;
  metal: string;
  mass: number;
  max_intensity: number;
  min_intensity: number;
  location: string;
  meta: object;
  created_at: string;

  levels?: {
    min: number,
    max: number
  }
}

export interface IChannelStats {
  hist: number[];
  bins: number[];
}
