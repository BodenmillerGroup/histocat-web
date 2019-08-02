import { Getters } from 'vuex-smart-module';
import { DatasetState } from '.';

export class DatasetGetters extends Getters<DatasetState> {
  get datasets() {
    return this.state.datasets;
  }
}
