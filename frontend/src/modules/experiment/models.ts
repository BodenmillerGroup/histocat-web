import { IImageFilter, IImageScalebar, IMaskSettings } from "@/modules/settings/models";

export interface IExperimentCreate {
  group_id: number;
  name: string;
  description?: string;
  tags?: string[];
}

export interface IExperimentUpdate {
  name?: string;
  description?: string;
  tags?: string[];
}

export interface IExperiment {
  id: number;
  group_id: number;
  user_id: number;
  name: string;
  description: string;
  meta: object;
  tags: string[];
  location: string;
  created_at: string;
}

export interface IExperimentData extends IExperiment {
  slides?: ISlide[];
}

export interface ISlide {
  id: number;
  experiment_id: number;
  origin_id: number;
  name: string;
  width_um: number;
  height_um: number;
  has_slide_image: boolean;
  meta: object;
  created_at: string;

  panoramas: IPanorama[];
  acquisitions: IAcquisition[];
}

export interface IPanorama {
  id: number;
  slide_id: number;
  origin_id: number;
  image_type: string;
  description: string;
  start_position_x: number;
  start_position_y: number;
  width: number;
  height: number;
  rotation_angle: number;
  meta: object;
}

export interface IAcquisition {
  id: number;
  slide_id: number;
  origin_id: number;
  description: string;
  max_x: number;
  max_y: number;
  signal_type: string;
  segment_data_format: string;
  ablation_frequency: number;
  ablation_power: number;
  start_timestamp: string;
  end_timestamp: string;
  movement_type: string;
  ablation_distance_between_shots_x: number;
  ablation_distance_between_shots_y: number;
  template: string;
  roi_start_x_pos_um: number;
  roi_start_y_pos_um: number;
  roi_end_x_pos_um: number;
  roi_end_y_pos_um: number;
  has_before_ablation_image: boolean;
  has_after_ablation_image: boolean;
  is_valid: boolean;
  meta: object;

  channels: { [name: string]: IChannel };
}

export interface IChannel {
  acquisition_id: number;
  origin_id: number;
  order_number: number;
  name: string;
  label: string;
  mass: number;
  max_intensity: number;
  min_intensity: number;
}

export interface IChannelStats {
  bins: number[];
}

export interface IChannelStack {
  acquisitionId: number | null;
  datasetId?: number;
  format?: string;
  filter: IImageFilter;
  scalebar: IImageScalebar;
  mask?: IMaskSettings;
  channels: {
    name: string;
    color?: string;
    min?: number;
    max?: number;
  }[];
}

export type ExportFormat = "tiff" | "png";
