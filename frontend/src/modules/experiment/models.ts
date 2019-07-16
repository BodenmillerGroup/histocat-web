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
  input: {
    acquisition_ids: number[];
    metals: string[];
    channel_settings: IChannelSettings[];
  }
  meta: object;
  location: string;
  created_at: string;
  updated_at: string;
}

export interface IDatasetCreate {
  experiment_id: number;
  name: string;
  description?: string;
  input: {
    acquisition_ids: number[];
    metals: string[];
    channel_settings: IChannelSettings[];
  }
  meta?: object;
}

export interface ISlide {
  id: number;
  experiment_id: number;
  metaname: string;
  original_id: number;
  original_metadata: string;
  meta: {
    UID: string;
    Description: string;
    FileName: string;
    SlideType: string;
    WidthUm: string;
    HeightUm: string;
    ImageEndOffset: string;
    ImageStartOffset: string;
    ImageFile: string;
  };
  location: string;
  created_at: string;
  panoramas: IPanorama[];
}

export interface IPanorama {
  id: number;
  slide_id: number;
  metaname: string;
  original_id: number;
  meta: {
    Description: string;
    SlideX1PosUm: string;
    SlideY1PosUm: string;
    SlideX2PosUm: string;
    SlideY2PosUm: string;
    SlideX3PosUm: string;
    SlideY3PosUm: string;
    SlideX4PosUm: string;
    SlideY4PosUm: string;
    ImageEndOffset: string;
    ImageStartOffset: string;
    PixelWidth: string;
    PixelHeight: string;
    ImageFormat: string;
    PixelScaleCoef: string;
  };
  location: string;
  created_at: string;
  rois: IRoi[];
}

export interface IRoi {
  id: number;
  panorama_id: number;
  metaname: string;
  original_id: number;
  meta: {
    ROIType: string;
  };
  location: string;
  created_at: string;
  acquisitions: IAcquisition[];
  roi_points: IRoiPoint[];
}

export interface IRoiPoint {
  id: number;
  roi_id: number;
  metaname: string;
  original_id: number;
  meta: {
    OrderNumber: string;
    SlideXPosUm: string;
    SlideYPosUm: string;
    PanoramaPixelXPos: string;
    PanoramaPixelYPos: string;
  };
  created_at: string;
}

export interface IAcquisition {
  id: number;
  roi_id: number;
  metaname: string;
  original_id: number;
  meta: {
    Description: string;
    OrderNumber: string;
    AblationPower: string;
    AblationDistanceBetweenShotsX: string;
    AblationDistanceBetweenShotsY: string;
    AblationFrequency: string;
    SignalType: string;
    DualCountStart: string;
    DataStartOffset: string;
    DataEndOffset: string;
    StartTimestamp: string;
    EndTimestamp: string;
    AfterAblationImageEndOffset: string;
    AfterAblationImageStartOffset: string;
    BeforeAblationImageEndOffset: string;
    BeforeAblationImageStartOffset: string;
    ROIStartXPosUm: string;
    ROIStartYPosUm: string;
    ROIEndXPosUm: string;
    ROIEndYPosUm: string;
    MovementType: string;
    SegmentDataFormat: string;
    ValueBytes: string;
    MaxY: string;
    MaxX: string;
    PlumeStart: string;
    PlumeEnd: string;
    Template: string;
  };
  location: string;
  created_at: string;
  channels: IChannel[];
}

export interface IChannel {
  id: number;
  acquisition_id: number;
  metaname: string;
  original_id: number;
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
  edges: number[];
}

export interface IChannelStack {
  channels: Array<{
    id: number;
    color?: string;
    min?: number;
    max?: number;
  }>;
}
