import { Actions, Context } from "vuex-smart-module";
import { Store } from "vuex";
import { analysisModule } from "@/modules/analysis";
import { datasetModule } from "@/modules/datasets";
import { experimentModule } from "@/modules/experiment";
import {selectionModule} from "@/modules/selection";

export class RootActions extends Actions {
  analysis?: Context<typeof analysisModule>;
  dataset?: Context<typeof datasetModule>;
  experiment?: Context<typeof experimentModule>;
  selection?: Context<typeof selectionModule>;

  // Called after the module is initialized
  $init(store: Store<any>): void {
    this.analysis = analysisModule.context(store);
    this.dataset = datasetModule.context(store);
    this.experiment = experimentModule.context(store);
    this.selection = selectionModule.context(store);
  }

  // Reset global store
  reset() {
    this.analysis?.mutations.reset();
    this.dataset?.mutations.reset();
    this.experiment?.mutations.reset();
    this.selection?.mutations.reset();
  }
}
