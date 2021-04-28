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
      <v-menu offset-y open-on-hover>
        <template v-slot:activator="{ on, attrs }">
          <v-btn v-bind="attrs" v-on="on" small elevation="1" :disabled="!activeAcquisition || !hasMask" class="ml-2">
            <v-icon left small>mdi-cog-outline</v-icon>
            Mode
          </v-btn>
        </template>
        <v-list dense>
          <v-list-item-group v-model="maskMode" color="primary">
            <v-list-item value="raw">
              <v-list-item-title>Show Raw Image</v-list-item-title>
            </v-list-item>
            <v-list-item value="mask">
              <v-list-item-title>Show Mask Overlay</v-list-item-title>
            </v-list-item>
            <v-list-item value="origin" :disabled="!activeDataset || activeDataset.origin !== 'DeepCell'">
              <v-list-item-title>Show Mask Source</v-list-item-title>
            </v-list-item>
          </v-list-item-group>
        </v-list>
        <v-list dense subheader>
          <v-subheader>Mask Opacity</v-subheader>
          <v-list-item>
            <v-list-item-content>
              <v-slider :value="maskOpacity" @end="maskOpacityHandler" hide-details min="0" max="1" step="0.01" />
            </v-list-item-content>
          </v-list-item>
        </v-list>
      </v-menu>
      <v-btn-toggle v-model="mouseMode" dense mandatory class="ml-8">
        <v-btn value="panZoom" small>
          <v-icon>mdi-arrow-top-left</v-icon>
        </v-btn>
        <v-btn value="lasso" small>
          <v-icon>mdi-lasso</v-icon>
        </v-btn>
      </v-btn-toggle>
      <v-switch v-model="regionsEnabled" label="Region statistics" hide-details inset class="ml-8" dense />
    </v-toolbar>
    <ImageViewer />
  </div>
</template>

<script lang="ts">
import { datasetsModule } from "@/modules/datasets";
import { projectsModule } from "@/modules/projects";
import { ExportFormat } from "@/modules/projects/models";
import { Component, Vue } from "vue-property-decorator";
import ImageViewer from "@/components/ImageViewer.vue";
import { mainModule } from "@/modules/main";
import { analysisModule } from "@/modules/analysis";

@Component({
  components: { ImageViewer },
})
export default class BlendTab extends Vue {
  readonly mainContext = mainModule.context(this.$store);
  readonly projectsContext = projectsModule.context(this.$store);
  readonly datasetsContext = datasetsModule.context(this.$store);
  readonly analysisContext = analysisModule.context(this.$store);

  get maskMode() {
    return this.mainContext.getters.maskMode;
  }

  set maskMode(value: "raw" | "mask" | "origin") {
    this.mainContext.mutations.setMaskMode(value);
    this.projectsContext.actions.getChannelStackImage();
  }

  get maskOpacity() {
    return this.mainContext.getters.maskOpacity;
  }

  maskOpacityHandler(value: number) {
    this.mainContext.mutations.setMaskOpacity(value);
    this.projectsContext.actions.getChannelStackImage();
  }

  get mouseMode() {
    return this.mainContext.getters.mouseMode;
  }

  set mouseMode(value: "panZoom" | "lasso" | "rotate") {
    this.mainContext.mutations.setMouseMode(value);
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
}
</script>
