<template>
  <div>
    <v-toolbar dense flat>
      <v-menu offset-y>
        <template v-slot:activator="{ on }">
          <v-btn v-on="on" small elevation="1">
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
        </v-list>
      </v-menu>
      <v-radio-group label="Mode:" v-model="mode" row dense mandatory hide-details :disabled="!hasMask" class="ml-8">
        <v-radio label="Raw" value="raw" />
        <v-radio label="Mask" value="mask" />
        <v-radio label="Mask Origin" value="origin" :disabled="!activeDataset || activeDataset.origin !== 'DeepCell'" />
      </v-radio-group>
      <v-switch v-model="regionsEnabled" label="Region statistics" hide-details inset class="ml-8" dense />
      <v-btn-toggle v-model="mouseMode" dense mandatory class="ml-8">
        <v-btn value="panZoom" small>
          <v-icon>mdi-arrow-top-left</v-icon>
        </v-btn>
        <v-btn value="lasso" small>
          <v-icon>mdi-lasso</v-icon>
        </v-btn>
      </v-btn-toggle>
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
    this.settingsContext.actions.setMaskSettings({
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
