import { IExperiment, IExperimentDataset } from './models';

export interface ExperimentsState {
    experiments: IExperiment[];
    activeExperimentId?: number;
    dataset?: IExperimentDataset;
    activeMeta?: object;
}
