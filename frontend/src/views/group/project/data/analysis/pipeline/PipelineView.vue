<template>
  <div>
    <v-toolbar dense flat>
      <v-menu offset-y>
        <template v-slot:activator="{ on }">
          <v-btn v-on="on" small elevation="1">
            <v-icon left small>mdi-pipe</v-icon>
            Add step
          </v-btn>
        </template>
        <v-list dense>
          <v-list-item @click="addStep('transformation')">
            <v-list-item-title>Transformation</v-list-item-title>
          </v-list-item>
          <v-list-item @click="addStep('tsne')">
            <v-list-item-title>t-SNE</v-list-item-title>
          </v-list-item>
          <v-list-item @click="addStep('umap')">
            <v-list-item-title>UMAP</v-list-item-title>
          </v-list-item>
        </v-list>
      </v-menu>
      <v-tooltip bottom>
        <template v-slot:activator="{ on }">
          <v-btn small elevation="1" v-on="on" @click="clearPipeline" class="ml-2" :disabled="steps.length === 0">
            Clear pipeline
          </v-btn>
        </template>
        <span>Clear pipeline</span>
      </v-tooltip>
    </v-toolbar>
    <v-timeline dense class="mx-2">
      <v-slide-x-reverse-transition group hide-on-leave>
        <v-timeline-item v-for="step in steps" :key="step.id" :color="step.color" small fill-dot>
          <TransformationStepEditor v-if="step.type === `transformation`" :step="step" :deleteStep="deleteStep" />
          <TsneStepEditor v-else-if="step.type === `tsne`" :step="step" :deleteStep="deleteStep" />
        </v-timeline-item>
      </v-slide-x-reverse-transition>
    </v-timeline>
  </div>
</template>

<script lang="ts">
import { analysisModule } from "@/modules/analysis";
import { datasetsModule } from "@/modules/datasets";
import { projectsModule } from "@/modules/projects";
import { settingsModule } from "@/modules/settings";
import { Component, Vue } from "vue-property-decorator";
import { mainModule } from "@/modules/main";
import { pipelinesModule } from "@/modules/pipelines";
import { PipelineStepType } from "@/modules/pipelines/models";
import TransformationStepEditor from "@/views/group/project/data/analysis/pipeline/steps/TransformationStepEditor.vue";
import TsneStepEditor from "@/views/group/project/data/analysis/pipeline/steps/TsneStepEditor.vue";

@Component({
  components: {TsneStepEditor, TransformationStepEditor}
})
export default class PipelineView extends Vue {
  readonly mainContext = mainModule.context(this.$store);
  readonly projectsContext = projectsModule.context(this.$store);
  readonly analysisContext = analysisModule.context(this.$store);
  readonly settingsContext = settingsModule.context(this.$store);
  readonly datasetsContext = datasetsModule.context(this.$store);
  readonly pipelinesContext = pipelinesModule.context(this.$store);

  steps: any[] = [];

  get dataset() {
    return this.datasetsContext.getters.activeDataset;
  }

  addStep(type: PipelineStepType) {
    switch (type) {
      case "transformation": {
        const step = {
          id: "transformation",
          type: "transformation",
          color: "blue",
          icon: "mdi-information",
          cofactor: 5,
        };
        this.steps.push(step);
        break;
      }
      case "tsne": {
        const step = {
          id: "tsne",
          type: "tsne",
          color: "blue",
          icon: "mdi-information",
        };
        this.steps.push(step);
        break;
      }
    }
  }

  deleteStep(step) {
    this.steps = this.steps.filter((item) => item !== step);
  }

  clearPipeline() {
    this.steps = [];
  }
}
</script>
