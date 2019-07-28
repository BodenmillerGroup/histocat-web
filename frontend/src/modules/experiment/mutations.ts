import { Mutations } from 'vuex-smart-module';
import { ExperimentState } from '.';
import { IDataset, IExperiment, IShare } from './models';


export class ExperimentMutations extends Mutations<ExperimentState> {
  setExperiments(experiments: IExperiment[]) {
    this.state.experiments = experiments;
  }

  setDatasets(datasets: IDataset[]) {
    this.state.datasets = datasets;
  }

  setShares(shares: IShare[]) {
    this.state.shares = shares;
  }

  setTags(tags: string[]) {
    this.state.tags = tags;
  }

  setExperiment(experiment: IExperiment) {
    const items = this.state.experiments.filter((item) => item.id !== experiment.id);
    items.push(experiment);
    this.state.experiments = items;
  }

  setDataset(dataset: IDataset) {
    const items = this.state.datasets.filter((item) => item.id !== dataset.id);
    items.push(dataset);
    this.state.datasets = items;
  }

  deleteExperiment(id: number) {
    this.state.experiments = this.state.experiments.filter((item) => item.id !== id);
  }

  deleteDataset(id: number) {
    this.state.datasets = this.state.datasets.filter((item) => item.id !== id);
  }

  setSelectedMetals(metals: string[]) {
    this.state.selectedMetals = metals;
  }

  setSelectedAcquisitionIds(ids: number[]) {
    this.state.selectedAcquisitionIds = ids;
  }

  setActiveExperimentId(id?: number) {
    this.state.activeExperimentId = id;
  }

  setActiveSlideId(id?: number) {
    this.state.activeSlideId = id;
  }

  setActivePanoramaId(id?: number) {
    this.state.activePanoramaId = id;
  }

  setActiveAcquisitionId(id?: number) {
    this.state.activeAcquisitionId = id;
  }

  setActiveWorkspaceNode(node?: { id: number, type: string }) {
    this.state.activeWorkspaceNode = node;
    if (node) {
      if (node.type === 'slide') {
        this.setActiveSlideId(node.id);
      } else if (node.type === 'panorama') {
        this.setActivePanoramaId(node.id);
      } else if (node.type === 'acquisition') {
        this.setActiveAcquisitionId(node.id);
      }
    }
  }

  setChannelStackImage(base64Image: string | ArrayBuffer | null) {
    this.state.channelStackImage = base64Image;
  }

  resetExperiment() {
    this.setActiveExperimentId(undefined);
    this.setActiveAcquisitionId(undefined);
    this.setSelectedMetals([]);
  }
}
