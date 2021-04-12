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
          <v-list-item-group v-model="mode" color="primary">
            <v-list-item value="raw">
              <v-list-item-title>Raw</v-list-item-title>
            </v-list-item>
            <v-list-item value="mask">
              <v-list-item-title>Mask</v-list-item-title>
            </v-list-item>
            <v-list-item value="origin" :disabled="!activeDataset || activeDataset.origin !== 'DeepCell'">
              <v-list-item-title>Mask Origin</v-list-item-title>
            </v-list-item>
          </v-list-item-group>
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
import { settingsModule } from "@/modules/settings";
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
  readonly settingsContext = settingsModule.context(this.$store);
  readonly datasetsContext = datasetsModule.context(this.$store);
  readonly analysisContext = analysisModule.context(this.$store);

  get mode() {
    return this.settingsContext.getters.maskSettings.mode;
  }

  set mode(value: "raw" | "mask" | "origin") {
    this.settingsContext.mutations.setMaskSettings({
      ...this.settingsContext.getters.maskSettings,
      mode: value,
    });
    this.projectsContext.actions.getChannelStackImage();
  }

  get mouseMode() {
    return this.settingsContext.getters.mouseMode;
  }

  set mouseMode(value: "panZoom" | "lasso" | "rotate") {
    this.settingsContext.mutations.setMouseMode(value);
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
