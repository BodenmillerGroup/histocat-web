import { IImageFilter, IImageScalebar } from '@/modules/settings/models';

export interface IAnalysisSubmission {
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
}
