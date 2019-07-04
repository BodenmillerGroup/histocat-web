import { IAcquisition, IChannel, IExperiment } from './models';

export interface ExperimentsState {
  experiments: IExperiment[];
  tags: string[];
  selectedExperimentId?: number;
  channels: IChannel[];
  selectedAcquisition?: IAcquisition;
  selectedMetals: string[];
  metalColorMap: { [metal: string]: string }
}
