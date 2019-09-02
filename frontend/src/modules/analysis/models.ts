import { IImageFilter, IImageScalebar } from '@/modules/settings/models';

export type ImageResultType = 'origin' | 'mask';

export interface IImageSegmentationSettings {
  algorithm: string;
  iterations: number;
  kernel_size: number;
  mask_color: string;
  result_type: ImageResultType
}

export interface IImageSegmentationSubmission {
  format?: string;
  filter: IImageFilter;
  scalebar: IImageScalebar;
  channels: Array<{
    id: number;
    color?: string;
    customLabel?: string;
    min?: number;
    max?: number;
  }>;
  settings: IImageSegmentationSettings;
}

export interface IScatterPlotData {
  x: number[];
  y: number[];
}
