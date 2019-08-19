import { IChannelSettings } from '@/modules/settings/models';

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
