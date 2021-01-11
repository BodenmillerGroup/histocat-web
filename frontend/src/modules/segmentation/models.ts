export interface ISegmentationFilterSettings {
  apply: boolean;
  type: string;
  kernel_size?: number;
  sigma?: number;
}

export interface ISegmentationSubmission {
  model_id: number;
  acquisition_ids: readonly number[];
  nuclei_channels: string[];
  cytoplasm_channels: string[];
  scaling_factor: number;
  upper_limit: number;
  overlap: number;
  filter_settings: ISegmentationFilterSettings;
}
