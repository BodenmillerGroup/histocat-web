import { Actions, Context } from "vuex-smart-module";
import { Store } from "vuex";
import { analysisModule } from "@/modules/analysis";
import { datasetModule } from "@/modules/datasets";
import { experimentModule } from "@/modules/experiment";
import { selectionModule } from "@/modules/selection";
import { centroidsModule } from "@/modules/centroids";
import { presetModule } from "@/modules/presets";
import { gateModule } from "@/modules/gates";

export class RootActions extends Actions {
  analysis?: Context<typeof analysisModule>;
  dataset?: Context<typeof datasetModule>;
  experiment?: Context<typeof experimentModule>;
  selection?: Context<typeof selectionModule>;
  centroids?: Context<typeof centroidsModule>;
  presets?: Context<typeof presetModule>;
  gates?: Context<typeof gateModule>;

  // Called after the module is initialized
  $init(store: Store<any>): void {
    this.analysis = analysisModule.context(store);
    this.dataset = datasetModule.context(store);
    this.experiment = experimentModule.context(store);
    this.selection = selectionModule.context(store);
    this.centroids = centroidsModule.context(store);
    this.presets = presetModule.context(store);
    this.gates = gateModule.context(store);
  }

  // Reset global store
  reset() {
    this.analysis?.mutations.reset();
    this.dataset?.mutations.reset();
    this.experiment?.mutations.reset();
    this.selection?.mutations.reset();
    this.centroids?.mutations.reset();
    this.presets?.mutations.reset();
    this.gates?.mutations.reset();
  }
}
