import { IChannelSettings } from '@/modules/settings/models';

export type Status = 'pending' | 'processing' | 'terminated';

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

export interface IDataset {
  id: number;
  experiment_id: number;
  user_id: number;
  uid: string;
  name: string;
  description: string;
  status: Status;
  artifacts: {
    acquisition_metadata?: {
      location: string;
    }
    cell?: {
      location: string;
    }
    image?: {
      location: string;
    }
    object_relationships?: {
      location: string;
    }
    probability_masks?: {
      [id: number]: {
        location: string;
        slide: {
          id: number;
          origin_id: number;
        }
        panorama: {
          id: number;
          origin_id: number;
        }
        roi: {
          id: number;
          origin_id: number;
        }
        acquisition: {
          id: number;
          origin_id: number;
        }
      }
    }
  }
  errors?: object;
  meta?: object;
  location: string;
  created_at: string;
  updated_at: string;
}
