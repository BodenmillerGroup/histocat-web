import { IChannel, IExperiment } from './models';

export interface ExperimentsState {
  experiments: IExperiment[];
  tags: string[];
  selectedExperimentId?: number;
  selectedAcquisitionId?: number;
  channels: IChannel[];
  selectedMetals: string[];
  metalColorMap: { [metal: string]: string }
}
