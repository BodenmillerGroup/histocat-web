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
  activeExperiment: (state: ExperimentsState) => {
    return state.experiments.find((item) => item.id === state.activeExperimentId);
  },
  activeAcquisition: (state: ExperimentsState, self) => {
    const experiment: IExperiment = self.activeExperiment;
    if (experiment && experiment.slides) {
      for (const slide of experiment.slides) {
        for (const panorama of slide.panoramas) {
          for (const roi of panorama.rois) {
            const acquisition = roi.acquisitions.find((item) => item.id === state.activeAcquisitionId);
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
    const acquisition = self.activeAcquisition;
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
  selectedAcquisitionIds: (state: ExperimentsState) => {
    return state.selectedAcquisitionIds;
  },
  activeExperimentId: (state: ExperimentsState) => {
    return state.activeExperimentId;
  },
  activeSlideId: (state: ExperimentsState) => {
    return state.activeSlideId;
  },
  activePanoramaId: (state: ExperimentsState) => {
    return state.activePanoramaId;
  },
};

const { read } = getStoreAccessors<ExperimentsState, RootState>('');

export const readAdminOneExperiment = read(getters.adminOneExperiment);
export const readExperiments = read(getters.experiments);
export const readActiveExperiment = read(getters.activeExperiment);
export const readTags = read(getters.tags);
export const readActiveAcquisition = read(getters.activeAcquisition);
export const readSelectedMetals = read(getters.selectedMetals);
export const readSelectedChannels = read(getters.selectedChannels);
export const readDatasets = read(getters.datasets);
export const readSelectedAcquisitionIds = read(getters.selectedAcquisitionIds);
export const readActiveExperimentId = read(getters.activeExperimentId);
export const readActiveSlideId = read(getters.activeSlideId);
export const readActivePanoramaId = read(getters.activePanoramaId);
