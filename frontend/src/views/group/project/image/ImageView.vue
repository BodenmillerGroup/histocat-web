<template>
  <div>
    <v-toolbar dense flat>
      <v-menu offset-y>
        <template v-slot:activator="{ on }">
          <v-btn v-on="on" small elevation="1" :disabled="!activeAcquisition">
            <v-icon left small>mdi-download</v-icon>
            Export
          </v-btn>
        </template>
        <v-list dense>
          <v-list-item @click="exportImage('tiff')">
            <v-list-item-title>Export TIFF</v-list-item-title>
          </v-list-item>
          <v-list-item @click="exportImage('png')">
            <v-list-item-title>Export PNG</v-list-item-title>
          </v-list-item>
          <v-list-item @click="exportImage('ome-tiff')">
            <v-list-item-title>Export OME-TIFF</v-list-item-title>
          </v-list-item>
        </v-list>
      </v-menu>
      <v-btn-toggle v-model="mouseMode" dense mandatory class="ml-4">
        <v-btn value="panZoom" small>
          <v-icon>mdi-arrow-top-left</v-icon>
        </v-btn>
        <v-btn value="lasso" small>
          <v-icon>mdi-lasso</v-icon>
        </v-btn>
      </v-btn-toggle>
      <v-switch
        v-model="showMask"
        label="Mask"
        hide-details
        inset
        class="ml-4"
        dense
        :disabled="!activeAcquisition || !hasMask"
      />
      <span class="opacityContainer ml-4">
        <v-slider
          label="Opacity"
          :value="maskOpacity"
          @end="maskOpacityHandler"
          hide-details
          min="0"
          max="1"
          step="0.01"
          dense
          :disabled="!activeAcquisition || !hasMask"
        />
      </span>
      <v-spacer />
      <v-switch v-model="regionsEnabled" label="Region" hide-details inset dense />
    </v-toolbar>
    <ImageViewer ref="imageViewer" />
  </div>
</template>

<script lang="ts">
import { datasetsModule } from "@/modules/datasets";
import { projectsModule } from "@/modules/projects";
import { ExportFormat } from "@/modules/projects/models";
import { Component, Vue } from "vue-property-decorator";
import ImageViewer from "@/components/ImageViewer.vue";
import { analysisModule } from "@/modules/analysis";
import { uiModule } from "@/modules/ui";

@Component({
  components: { ImageViewer },
})
export default class ImageView extends Vue {
  readonly uiContext = uiModule.context(this.$store);
  readonly projectsContext = projectsModule.context(this.$store);
  readonly datasetsContext = datasetsModule.context(this.$store);
  readonly analysisContext = analysisModule.context(this.$store);

  get showMask() {
    return this.uiContext.getters.showMask;
  }

  set showMask(value: boolean) {
    this.uiContext.mutations.setShowMask(value);
    this.projectsContext.actions.getChannelStackImage();
  }

  get maskOpacity() {
    return this.uiContext.getters.maskOpacity;
  }

  maskOpacityHandler(value: number) {
    this.uiContext.mutations.setMaskOpacity(value);
    this.projectsContext.actions.getChannelStackImage();
  }

  get mouseMode() {
    return this.uiContext.getters.mouseMode;
  }

  set mouseMode(value: "panZoom" | "lasso" | "rotate") {
    this.uiContext.mutations.setMouseMode(value);
  }

  get regionsEnabled() {
    return this.analysisContext.getters.regionsEnabled;
  }

  set regionsEnabled(value: boolean) {
    this.analysisContext.mutations.setRegionsEnabled(value);
  }

  get activeProjectId() {
    return this.projectsContext.getters.activeProjectId;
  }

  get activeAcquisition() {
    return this.projectsContext.getters.activeAcquisition;
  }

  get activeDataset() {
    return this.datasetsContext.getters.activeDataset;
  }

  get hasMask() {
    let hasMask = false;
    if (this.activeAcquisition && this.activeDataset && this.activeDataset.meta.masks) {
      hasMask = !!this.activeDataset.meta.masks[this.activeAcquisition.id];
    }
    return hasMask;
  }

  exportImage(format: ExportFormat) {
    this.projectsContext.actions.exportChannelStackImage(format);
  }

  refresh() {
    (this.$refs.imageViewer as any).refresh();
  }
}
</script>

<style scoped>
.opacityContainer {
  width: 200px;
}
</style>
