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
      <v-switch v-model="applyMask" label="Mask overlay" hide-details inset class="ml-8" :disabled="!hasMask" dense />
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

@Component({
  components: { ImageViewer },
})
export default class BlendTab extends Vue {
  readonly mainContext = mainModule.context(this.$store);
  readonly projectsContext = projectsModule.context(this.$store);
  readonly settingsContext = settingsModule.context(this.$store);
  readonly datasetsContext = datasetsModule.context(this.$store);

  get applyMask() {
    return this.settingsContext.getters.maskSettings.apply;
  }

  set applyMask(value: boolean) {
    this.settingsContext.actions.setMaskSettings({
      ...this.settingsContext.getters.maskSettings,
      apply: value,
    });
    this.projectsContext.actions.getChannelStackImage();
  }

  get activeProjectId() {
    return this.projectsContext.getters.activeProjectId;
  }

  get activeAcquisition() {
    return this.projectsContext.getters.activeAcquisition;
  }

  get dataset() {
    return this.datasetsContext.getters.activeDataset;
  }

  get hasMask() {
    let hasMask = false;
    if (this.activeAcquisition && this.dataset && this.dataset.meta.probability_masks) {
      hasMask = !!this.dataset.meta.probability_masks[this.activeAcquisition.id];
    }
    return hasMask;
  }

  exportImage(format: ExportFormat) {
    this.projectsContext.actions.exportChannelStackImage(format);
  }
}
</script>
