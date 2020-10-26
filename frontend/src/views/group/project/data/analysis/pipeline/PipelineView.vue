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
          <v-subheader>Preprocessing</v-subheader>
          <v-list-item @click="addStep('markersFilter')">
            <v-list-item-title>Markers Filter</v-list-item-title>
          </v-list-item>
          <v-list-item @click="addStep('transformation')">
            <v-list-item-title>Transformation</v-list-item-title>
          </v-list-item>
          <v-list-item @click="addStep('scale')">
            <v-list-item-title>Scale</v-list-item-title>
          </v-list-item>
          <v-list-item @click="addStep('neighbors')">
            <v-list-item-title>Neighbors</v-list-item-title>
          </v-list-item>
          <v-subheader>Embeddings</v-subheader>
          <v-list-item @click="addStep('pca')">
            <v-list-item-title>PCA</v-list-item-title>
          </v-list-item>
          <v-list-item @click="addStep('tsne')">
            <v-list-item-title>t-SNE</v-list-item-title>
          </v-list-item>
          <v-list-item @click="addStep('umap')">
            <v-list-item-title>UMAP</v-list-item-title>
          </v-list-item>
          <v-subheader>Clustering</v-subheader>
          <v-list-item @click="addStep('leiden')">
            <v-list-item-title>Leiden</v-list-item-title>
          </v-list-item>
          <v-list-item @click="addStep('louvain')">
            <v-list-item-title>Louvain</v-list-item-title>
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
      <v-tooltip bottom>
        <template v-slot:activator="{ on }">
          <v-btn small elevation="1" v-on="on" @click="printPipeline" class="ml-2" :disabled="steps.length === 0">
            Print
          </v-btn>
        </template>
        <span>Print pipeline</span>
      </v-tooltip>
      <v-dialog v-model="dialog" scrollable max-width="600px">
        <template v-slot:activator="{ on, attrs }">
          <v-btn small elevation="1" color="primary" v-bind="attrs" v-on="on" class="ml-2" :disabled="steps.length === 0">Process</v-btn>
        </template>
        <v-card>
          <v-card-title>Select Acquisitions</v-card-title>
          <v-divider />
          <v-card-text style="height: 300px">
            <AcquisitionsSelector />
          </v-card-text>
          <v-divider />
          <v-card-actions>
            <v-spacer />
            <v-btn text @click="dialog = false">Cancel</v-btn>
            <v-btn color="primary" text @click="processPipeline">Process</v-btn>
          </v-card-actions>
        </v-card>
      </v-dialog>
    </v-toolbar>
    <v-timeline dense class="mr-2">
      <v-slide-x-reverse-transition group hide-on-leave>
        <v-timeline-item v-for="step in steps" :key="step.id" color="blue" :icon="icons[step.type]" fill-dot>
          <MarkersFilterStepEditor v-if="step.type === `markersFilter`" :step="step" :deleteStep="deleteStep" />
          <TransformationStepEditor v-else-if="step.type === `transformation`" :step="step" :deleteStep="deleteStep" />
          <ScaleStepEditor v-else-if="step.type === `scale`" :step="step" :deleteStep="deleteStep" />
          <NeighborsStepEditor v-else-if="step.type === `neighbors`" :step="step" :deleteStep="deleteStep" />
          <PcaStepEditor v-else-if="step.type === `pca`" :step="step" :deleteStep="deleteStep" />
          <TsneStepEditor v-else-if="step.type === `tsne`" :step="step" :deleteStep="deleteStep" />
          <UmapStepEditor v-else-if="step.type === `umap`" :step="step" :deleteStep="deleteStep" />
          <LeidenStepEditor v-else-if="step.type === `leiden`" :step="step" :deleteStep="deleteStep" />
          <LouvainStepEditor v-else-if="step.type === `louvain`" :step="step" :deleteStep="deleteStep" />
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
import { v4 as uuidv4 } from "uuid";
import TransformationStepEditor from "@/views/group/project/data/analysis/pipeline/steps/TransformationStepEditor.vue";
import TsneStepEditor from "@/views/group/project/data/analysis/pipeline/steps/TsneStepEditor.vue";
import UmapStepEditor from "@/views/group/project/data/analysis/pipeline/steps/UmapStepEditor.vue";
import PcaStepEditor from "@/views/group/project/data/analysis/pipeline/steps/PcaStepEditor.vue";
import ScaleStepEditor from "@/views/group/project/data/analysis/pipeline/steps/ScaleStepEditor.vue";
import MarkersFilterStepEditor from "@/views/group/project/data/analysis/pipeline/steps/MarkersFilterStepEditor.vue";
import NeighborsStepEditor from "@/views/group/project/data/analysis/pipeline/steps/NeighborsStepEditor.vue";
import AcquisitionsSelector from "@/views/group/project/data/analysis/pipeline/AcquisitionsSelector.vue";
import LeidenStepEditor from "@/views/group/project/data/analysis/pipeline/steps/LeidenStepEditor.vue";
import LouvainStepEditor from "@/views/group/project/data/analysis/pipeline/steps/LouvainStepEditor.vue";

@Component({
  components: {
    LouvainStepEditor,
    LeidenStepEditor,
    AcquisitionsSelector,
    NeighborsStepEditor,
    MarkersFilterStepEditor,
    ScaleStepEditor,
    PcaStepEditor,
    UmapStepEditor,
    TsneStepEditor,
    TransformationStepEditor,
  },
})
export default class PipelineView extends Vue {
  readonly mainContext = mainModule.context(this.$store);
  readonly projectsContext = projectsModule.context(this.$store);
  readonly analysisContext = analysisModule.context(this.$store);
  readonly settingsContext = settingsModule.context(this.$store);
  readonly datasetsContext = datasetsModule.context(this.$store);
  readonly pipelinesContext = pipelinesModule.context(this.$store);

  readonly icons = {
    markersFilter: "mdi-filter-outline",
    transformation: "mdi-label-variant-outline",
    scale: "mdi-label-variant-outline",
    neighbors: "mdi-graph-outline",
    pca: "mdi-chart-scatter-plot",
    tsne: "mdi-chart-scatter-plot",
    umap: "mdi-chart-scatter-plot",
    leiden: "mdi-dots-hexagon",
    louvain: "mdi-dots-hexagon",
  };

  dialog = false;

  get steps() {
    return this.pipelinesContext.getters.steps;
  }

  get dataset() {
    return this.datasetsContext.getters.activeDataset;
  }

  addMarkersFilterStep() {
    return {
      id: uuidv4(),
      type: "markersFilter",
      markers: [],
    };
  }

  addTransformationStep() {
    return {
      id: uuidv4(),
      type: "transformation",
      mode: "arcsinh",
      cofactor: 5,
    };
  }

  addScaleStep() {
    return {
      id: uuidv4(),
      type: "scale",
      zeroCenter: true,
      maxValue: undefined,
    };
  }

  addNeighborsStep() {
    return {
      id: uuidv4(),
      type: "neighbors",
      nNeighbors: 15,
      knn: true,
      method: "umap",
      metric: "euclidean",
    };
  }

  addPcaStep() {
    return {
      id: uuidv4(),
      type: "pca",
      svdSolver: "arpack",
    };
  }

  addTsneStep() {
    return {
      id: uuidv4(),
      type: "tsne",
      perplexity: 30,
      earlyExaggeration: 12,
      learningRate: 1000,
    };
  }

  addUmapStep() {
    return {
      id: uuidv4(),
      type: "umap",
      minDist: 0.5,
      spread: 1.0,
    };
  }

  addLeidenStep() {
    return {
      id: uuidv4(),
      type: "leiden",
      resolution: 1.0,
      directed: true,
      useWeights: true,
    };
  }

  addLouvainStep() {
    return {
      id: uuidv4(),
      type: "louvain",
      resolution: 1.0,
      flavor: "vtraag",
      directed: true,
      useWeights: false,
    };
  }

  addStep(type: PipelineStepType) {
    switch (type) {
      case "markersFilter": {
        this.steps.push(this.addMarkersFilterStep());
        break;
      }
      case "transformation": {
        this.steps.push(this.addTransformationStep());
        break;
      }
      case "scale": {
        this.steps.push(this.addScaleStep());
        break;
      }
      case "neighbors": {
        this.steps.push(this.addNeighborsStep());
        break;
      }
      case "pca": {
        this.steps.push(this.addPcaStep());
        break;
      }
      case "tsne": {
        this.steps.push(this.addTsneStep());
        break;
      }
      case "umap": {
        this.steps.push(this.addUmapStep());
        break;
      }
      case "leiden": {
        this.steps.push(this.addLeidenStep());
        break;
      }
      case "louvain": {
        this.steps.push(this.addLouvainStep());
        break;
      }
    }
  }

  deleteStep(step) {
    this.pipelinesContext.mutations.setSteps(this.steps.filter((item) => item !== step));
  }

  clearPipeline() {
    this.pipelinesContext.mutations.setSteps([]);
  }

  printPipeline() {
    console.log(this.steps);
  }

  async processPipeline() {
    this.dialog = false;
    await this.pipelinesContext.actions.processPipeline();
  }
}
</script>
