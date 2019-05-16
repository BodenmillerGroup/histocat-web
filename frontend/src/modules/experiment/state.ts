import { IAcquisition, IChannel, IExperiment, IExperimentDataset } from './models';

export interface ExperimentsState {
  experiments: IExperiment[];
  selectedExperimentId?: number;
  dataset?: IExperimentDataset;
  selectedMeta?: object;
  channels: IChannel[];
  selectedAcquisition?: IAcquisition;
  selectedMetals: string[];
  metalColorMap: { [metal: string]: string }
}
