import { IChannelSettings } from '@/modules/settings/models';

export interface IExperiment {
  id: number;
  user_id: number;
  name: string;
  description: string;
  meta: object;
  tags: string[];
  location: string;
  created_at: string;

  slides?: ISlide[];
}

export interface IExperimentUpdate {
  name?: string;
  description?: string;
  tags?: string[];
}

export interface IExperimentCreate {
  name: string;
  description?: string;
  tags?: string[];
}

export type Status = 'pending' | 'processing' | 'terminated';

export interface IDataset {
  id: number;
  experiment_id: number;
  user_id: number;
  uid: string;
  name: string;
  description: string;
  status: Status;
  meta: {
    input: {
      acquisition_ids: number[];
      metals: string[];
      channel_settings: IChannelSettings[];
    }
  };
  location: string;
  created_at: string;
  updated_at: string;
}

export interface IDatasetCreate {
  experiment_id: number;
  name: string;
  description?: string;
  meta?: {
    input: {
      acquisition_ids: number[];
      metals: string[];
      channel_settings: IChannelSettings[];
    }
  };
}

export interface ISlide {
  id: number;
  experiment_id: number;
  uid: string;
  description: string;
  filename: string;
  slide_type: string;
  width_um: number;
  height_um: number;
  image_end_offset: number;
  image_start_offset: number;
  image_file: string;
  meta: object;
  location: string;
  created_at: string;
  panoramas: IPanorama[];
}

export interface IPanorama {
  id: number;
  slide_id: number;
  description: string;
  slide_x1_pos_um: number;
  slide_y1_pos_um: number;
  slide_x2_pos_um: number;
  slide_y2_pos_um: number;
  slide_x3_pos_um: number;
  slide_y3_pos_um: number;
  slide_x4_pos_um: number;
  slide_y4_pos_um: number;
  image_end_offset: number;
  image_start_offset: number;
  pixel_width: number;
  pixel_height: number;
  image_format: string;
  pixel_scale_coef: number;
  meta: object;
  location: string;
  created_at: string;
  rois: IRoi[];
}

export interface IRoi {
  id: number;
  panorama_id: number;
  roi_type: string;
  location: string;
  created_at: string;
  acquisitions: IAcquisition[];
  roi_points: IRoiPoint[];
}

export interface IRoiPoint {
  id: number;
  roi_id: number;
  order_number: number;
  slide_x_pos_um: number;
  slide_y_pos_um: number;
  panorama_pixel_x_pos: number;
  panorama_pixel_y_pos: number;
  created_at: string;
}

export interface IAcquisition {
  id: number;
  roi_id: number;
  description: string;
  order_number: number;
  ablation_power: number;
  ablation_distance_between_shots_x: number;
  ablation_distance_between_shots_y: number;
  ablation_frequency: number;
  signal_type: string;
  dual_count_start: number;
  data_start_offset: number;
  data_end_offset: number;
  start_timestamp: string;
  end_timestamp: string;
  after_ablation_image_end_offset: number;
  after_ablation_image_start_offset: number;
  before_ablation_image_end_offset: number;
  before_ablation_image_start_offset: number;
  roi_start_x_pos_um: number;
  roi_start_y_pos_um: number;
  roi_end_x_pos_um: number;
  roi_end_y_pos_um: number;
  movement_type: string;
  segment_data_format: string;
  value_bytes: number;
  max_y: number;
  max_x: number;
  plume_start: number;
  plume_end: number;
  template: string;
  meta: object;
  location: string;
  created_at: string;
  channels: IChannel[];
}

export interface IChannel {
  id: number;
  acquisition_id: number;
  metal: string;
  label: string;
  mass: number;
  max_intensity: number;
  min_intensity: number;
  meta: object;
  location: string;
  created_at: string;
}

export interface IChannelStats {
  hist: number[];
  bins: number[];
}
