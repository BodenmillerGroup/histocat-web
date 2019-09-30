export interface IChannelSettings {
  id: number;
  customLabel: string;
  levels?: {
    min: number;
    max: number;
  };
  suppressBroadcast?: boolean;
}

export interface IImageFilter {
  apply: boolean;
  type: string;
  settings?: any;
  suppressBroadcast?: boolean;
}

export interface IImageLegend {
  apply: boolean;
  fontScale: number;
  showIntensity: boolean;
  suppressBroadcast?: boolean;
}

export interface IImageScalebar {
  apply: boolean;
  settings?: any;
  suppressBroadcast?: boolean;
}

export interface IMaskSettings {
  apply: boolean;
  location?: string;
  settings?: any;
  suppressBroadcast?: boolean;
}
