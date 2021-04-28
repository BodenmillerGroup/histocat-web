export interface IChannelSettings {
  color: string;
  levels?: {
    min: number;
    max: number;
  };
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
  mode: "raw" | "mask" | "origin";
  opacity: number;
  gated?: boolean;
  cells: {
    [color: string]: number[];
  };
  resultId?: number;
  location?: string;
  settings?: any;
  colorsType?: string;
  colorsName?: string;
}
