export interface ISegmentationPreprocessingSettings {
  channels_normalization?: string;
  threshold?: boolean;
  percentile?: number;
  normalize: boolean;
  kernel_size: number;
}

export interface ISegmentationPostprocessingSettings {
  radius: number;
  maxima_threshold: number;
  interior_threshold: number;
  small_objects_threshold: number;
  fill_holes_threshold: number;
  interior_model: string;
  maxima_model: string;
  interior_model_smooth: number;
  maxima_model_smooth: number;
  pixel_expansion: number | null;
}

export interface ISegmentationSubmission {
  dataset_name: string | null;
  dataset_description: string | null;
  model_id: number;
  compartment: string;
  acquisition_ids: readonly number[];
  channels: string[];
  nuclei_channels: string[];
  cytoplasm_channels: string[];
  preprocessing: ISegmentationPreprocessingSettings;
  postprocessing: ISegmentationPostprocessingSettings;
}
