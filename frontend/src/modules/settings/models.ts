export interface IChannelSettings {
  id: number;
  customLabel: string;
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
  fontScale: number;
  showIntensity: boolean;
}

export interface IImageScalebar {
  apply: boolean;
  settings?: any;
}

export interface IMaskSettings {
  apply: boolean;
  location?: string;
  settings?: any;
}
