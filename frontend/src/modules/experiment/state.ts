import { IExperiment } from './models';

export interface ExperimentsState {
    experiments: IExperiment[];
    activeExperimentId?: number;
}
