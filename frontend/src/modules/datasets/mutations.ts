import { Mutations } from 'vuex-smart-module';
import { DatasetState } from '.';
import { IDataset } from './models';


export class DatasetMutations extends Mutations<DatasetState> {
  setDatasets(datasets: IDataset[]) {
    this.state.datasets = datasets;
  }

  setDataset(dataset: IDataset) {
    const items = this.state.datasets.filter((item) => item.id !== dataset.id);
    items.push(dataset);
    this.state.datasets = items;
  }

  deleteDataset(id: number) {
    this.state.datasets = this.state.datasets.filter((item) => item.id !== id);
  }
}
