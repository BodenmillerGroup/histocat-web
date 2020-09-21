import { Actions, Context } from "vuex-smart-module";
import { Store } from "vuex";
import { analysisModule } from "@/modules/analysis";
import { datasetsModule } from "@/modules/datasets";
import { projectsModule } from "@/modules/projects";
import { selectionModule } from "@/modules/selection";
import { centroidsModule } from "@/modules/centroids";
import { presetsModule } from "@/modules/presets";
import { gatesModule } from "@/modules/gates";
import { groupModule } from "@/modules/group";
import { memberModule } from "@/modules/member";
import { resultsModule } from "@/modules/results";

export class RootActions extends Actions {
  group?: Context<typeof groupModule>;
  member?: Context<typeof memberModule>;
  analysis?: Context<typeof analysisModule>;
  datasets?: Context<typeof datasetsModule>;
  results?: Context<typeof resultsModule>;
  projects?: Context<typeof projectsModule>;
  selection?: Context<typeof selectionModule>;
  centroids?: Context<typeof centroidsModule>;
  presets?: Context<typeof presetsModule>;
  gates?: Context<typeof gatesModule>;

  // Called after the module is initialized
  $init(store: Store<any>): void {
    this.group = groupModule.context(store);
    this.member = memberModule.context(store);
    this.analysis = analysisModule.context(store);
    this.datasets = datasetsModule.context(store);
    this.results = resultsModule.context(store);
    this.projects = projectsModule.context(store);
    this.selection = selectionModule.context(store);
    this.centroids = centroidsModule.context(store);
    this.presets = presetsModule.context(store);
    this.gates = gatesModule.context(store);
  }

  // Reset project store
  resetProject() {
    this.member?.mutations.reset();
    this.analysis?.mutations.reset();
    this.datasets?.mutations.reset();
    this.results?.mutations.reset();
    this.projects?.mutations.reset();
    this.selection?.mutations.reset();
    this.centroids?.mutations.reset();
    this.presets?.mutations.reset();
    this.gates?.mutations.reset();
  }

  // Reset global store
  reset() {
    this.resetProject();
    this.group?.mutations.reset();
  }
}
