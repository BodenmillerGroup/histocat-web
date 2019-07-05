import { RootState } from '@/store';
import { getStoreAccessors } from 'typesafe-vuex';
import { ExperimentsState } from '.';
import { IExperiment } from './models';

export const getters = {
  experiments: (state: ExperimentsState) => state.experiments,
  datasets: (state: ExperimentsState) => state.datasets,
  adminOneExperiment: (state: ExperimentsState) => (id: number) => {
    return state.experiments.find((item) => item.id === id);
  },
  selectedExperimentId: (state: ExperimentsState) => {
    return state.selectedExperimentId;
  },
  selectedExperiment: (state: ExperimentsState) => {
    return state.experiments.find((item) => item.id === state.selectedExperimentId);
  },
  selectedAcquisition: (state: ExperimentsState, self) => {
    const experiment: IExperiment = self.selectedExperiment;
    if (experiment && experiment.slides) {
      for (const slide of experiment.slides) {
        for (const panorama of slide.panoramas) {
          for (const roi of panorama.rois) {
            const acquisition = roi.acquisitions.find((item) => item.id === state.selectedAcquisitionId);
            if (acquisition) {
              return acquisition;
            }
          }
        }
      }
    }
    return undefined;
  },
  selectedMetals: (state: ExperimentsState) => {
    return state.selectedMetals;
  },
  selectedChannels: (state: ExperimentsState, self) => {
    const acquisition = self.selectedAcquisition;
    if (acquisition) {
      const channels = acquisition.channels.filter((channel) => {
        if (state.selectedMetals.includes(channel.metal)) {
          return channel;
        }
      });
      return channels;
    }
  },
  tags: (state: ExperimentsState) => {
    return state.tags;
  },
};

const { read } = getStoreAccessors<ExperimentsState, RootState>('');

export const readAdminOneExperiment = read(getters.adminOneExperiment);
export const readExperiments = read(getters.experiments);
export const readSelectedExperiment = read(getters.selectedExperiment);
export const readSelectedExperimentId = read(getters.selectedExperimentId);
export const readTags = read(getters.tags);
export const readSelectedAcquisition = read(getters.selectedAcquisition);
export const readSelectedMetals = read(getters.selectedMetals);
export const readSelectedChannels = read(getters.selectedChannels);
export const readDatasets = read(getters.datasets);
