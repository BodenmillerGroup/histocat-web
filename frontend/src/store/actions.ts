import { Actions, Context } from "vuex-smart-module";
import { Store } from "vuex";
import { analysisModule } from "@/modules/analysis";
import { datasetsModule } from "@/modules/datasets";
import { projectsModule } from "@/modules/projects";
import { presetsModule } from "@/modules/presets";
import { gatesModule } from "@/modules/gates";
import { groupModule } from "@/modules/group";
import { memberModule } from "@/modules/member";
import { pipelinesModule } from "@/modules/pipelines";
import { modelsModule } from "@/modules/models";
import { segmentationModule } from "@/modules/segmentation";
import { annotationsModule } from "@/modules/annotations";
import { cellsModule } from "@/modules/cells";
import { uiModule } from "@/modules/ui";

export class RootActions extends Actions {
  ui?: Context<typeof uiModule>;
  group?: Context<typeof groupModule>;
  member?: Context<typeof memberModule>;
  analysis?: Context<typeof analysisModule>;
  datasets?: Context<typeof datasetsModule>;
  projects?: Context<typeof projectsModule>;
  presets?: Context<typeof presetsModule>;
  gates?: Context<typeof gatesModule>;
  pipelines?: Context<typeof pipelinesModule>;
  models?: Context<typeof modelsModule>;
  segmentation?: Context<typeof segmentationModule>;
  annotations?: Context<typeof annotationsModule>;
  cells?: Context<typeof cellsModule>;

  // Called after the module is initialized
  $init(store: Store<any>): void {
    this.ui = uiModule.context(store);
    this.group = groupModule.context(store);
    this.member = memberModule.context(store);
    this.analysis = analysisModule.context(store);
    this.datasets = datasetsModule.context(store);
    this.projects = projectsModule.context(store);
    this.presets = presetsModule.context(store);
    this.gates = gatesModule.context(store);
    this.pipelines = pipelinesModule.context(store);
    this.models = modelsModule.context(store);
    this.segmentation = segmentationModule.context(store);
    this.annotations = annotationsModule.context(store);
    this.cells = cellsModule.context(store);
  }

  // Reset project store
  resetProject() {
    this.member?.mutations.reset();
    this.analysis?.mutations.reset();
    this.datasets?.mutations.reset();
    this.projects?.mutations.reset();
    this.presets?.mutations.reset();
    this.gates?.mutations.reset();
    this.pipelines?.mutations.reset();
    this.segmentation?.mutations.reset();
    this.annotations?.mutations.reset();
    this.cells?.mutations.reset();
  }

  // Reset global store
  reset() {
    this.resetProject();
    this.models?.mutations.reset();
    this.group?.mutations.reset();
    this.ui?.mutations.reset();
  }
}
