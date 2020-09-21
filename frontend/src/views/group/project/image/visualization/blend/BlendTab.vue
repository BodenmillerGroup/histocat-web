<template>
  <v-row no-gutters>
    <v-col>
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
        <v-switch v-model="regionsEnabled" label="Enable regions" hide-details inset class="ml-8" dense />
        <v-tooltip bottom>
          <template v-slot:activator="{ on }">
            <v-btn small elevation="1" v-on="on" @click="deleteRegions" class="ml-2" :disabled="!selectedRegion">
              Delete region(s)
            </v-btn>
          </template>
          <span>Delete selected regions</span>
        </v-tooltip>
        <v-tooltip bottom>
          <template v-slot:activator="{ on }">
            <v-btn small elevation="1" v-on="on" @click="calculateRegionStats" class="ml-2" :disabled="!selectedRegion">
              Region stats
            </v-btn>
          </template>
          <span>Calculate region's statistics</span>
        </v-tooltip>
        <!--        <v-tooltip bottom>-->
        <!--          <template v-slot:activator="{ on }">-->
        <!--            <v-btn-->
        <!--              small-->
        <!--              elevation="1"-->
        <!--              v-on="on"-->
        <!--              @click="getColorizedMaskImage"-->
        <!--              class="ml-2"-->
        <!--              :loading="colorizeMaskInProgress"-->
        <!--              :disabled="colorizeMaskInProgress"-->
        <!--            >-->
        <!--              Colorize mask-->
        <!--            </v-btn>-->
        <!--          </template>-->
        <!--          <span>Request calculation of colorized cell mask</span>-->
        <!--        </v-tooltip>-->
        <v-switch v-model="applyMask" label="Mask overlay" hide-details inset class="ml-8" :disabled="!hasMask" dense />
      </v-toolbar>
      <v-row no-gutters>
        <v-col>
          <BlendView ref="blendView" class="blend-view" />
        </v-col>
        <IntensityView />
      </v-row>
    </v-col>
  </v-row>
</template>

<script lang="ts">
import { analysisModule } from "@/modules/analysis";
import { IRegionStatsSubmission } from "@/modules/analysis/models";
import { datasetsModule } from "@/modules/datasets";
import { projectsModule } from "@/modules/projects";
import { ExportFormat } from "@/modules/projects/models";
import { settingsModule } from "@/modules/settings";
import BlendView from "@/views/group/project/image/visualization/blend/BlendView.vue";
import IntensityView from "@/views/group/project/image/visualization/blend/IntensityView.vue";
import Polygon from "ol/geom/Polygon";
import { Component, Vue } from "vue-property-decorator";

@Component({
  components: { IntensityView, BlendView },
})
export default class BlendTab extends Vue {
  readonly projectsContext = projectsModule.context(this.$store);
  readonly analysisContext = analysisModule.context(this.$store);
  readonly settingsContext = settingsModule.context(this.$store);
  readonly datasetContext = datasetsModule.context(this.$store);

  get colorizeMaskInProgress() {
    return this.projectsContext.getters.colorizeMaskInProgress;
  }

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

  get regionsEnabled() {
    return this.analysisContext.getters.regionsEnabled;
  }

  set regionsEnabled(value: boolean) {
    this.analysisContext.mutations.setRegionsEnabled(value);
  }

  get selectedRegion() {
    return this.analysisContext.getters.selectedRegion;
  }

  get activeProjectId() {
    return this.projectsContext.getters.activeProjectId;
  }

  get activeAcquisition() {
    return this.projectsContext.getters.activeAcquisition;
  }

  get dataset() {
    return this.datasetContext.getters.activeDataset;
  }

  get hasMask() {
    let hasMask = false;
    if (this.activeAcquisition && this.dataset && this.dataset.meta.probability_masks) {
      hasMask = !!this.dataset.meta.probability_masks[this.activeAcquisition.id];
    }
    return hasMask;
  }

  calculateRegionStats() {
    if (!this.activeProjectId || !this.activeAcquisition || !this.selectedRegion) {
      return;
    }
    const polygon = this.selectedRegion.getGeometry()! as Polygon;
    const coords = polygon.getCoordinates()[0];
    const height = this.activeAcquisition.max_y;
    // TODO: Y axis flip
    const convertedCoords = coords.map((point) => [point[0], Math.abs(point[1] - height)]);
    const params: IRegionStatsSubmission = {
      project_id: this.activeProjectId,
      acquisition_id: this.activeAcquisition.id,
      region_polygon: convertedCoords,
    };
    return this.analysisContext.actions.calculateRegionStats(params);
  }

  exportImage(format: ExportFormat) {
    this.projectsContext.actions.exportChannelStackImage(format);
  }

  getColorizedMaskImage() {
    this.projectsContext.actions.getColorizedMaskImage();
  }

  deleteRegions() {
    (this.$refs.blendView as any).deleteRegions();
  }
}
</script>

<style scoped>
.blend-view {
  height: calc(100vh - 174px);
}
</style>
