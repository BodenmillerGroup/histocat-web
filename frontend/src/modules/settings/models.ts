export interface IChannelSettings {
  id: number;
  levels?: {
    min: number,
    max: number
  }
}

export interface IImageFilter {
  apply: boolean;
  type: string;
  settings?: any;
}

export interface IImageLegend {
  apply: boolean;
  settings?: any;
}

export interface IImageScalebar {
  apply: boolean;
  settings?: any;
}
