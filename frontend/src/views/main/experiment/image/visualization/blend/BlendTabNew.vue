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
          <ImageViewer ref="imageViewer" :canvas-width="canvasWidth" :canvas-height="canvasHeight" />
        </v-col>
        <IntensityView />
      </v-row>
    </v-col>
  </v-row>
</template>

<script lang="ts">
import { analysisModule } from "@/modules/analysis";
import { IRegionStatsSubmission } from "@/modules/analysis/models";
import { datasetModule } from "@/modules/datasets";
import { experimentModule } from "@/modules/experiment";
import { ExportFormat } from "@/modules/experiment/models";
import { settingsModule } from "@/modules/settings";
import IntensityView from "@/views/main/experiment/image/visualization/blend/IntensityView.vue";
import Polygon from "ol/geom/Polygon";
import { Component, Vue } from "vue-property-decorator";
import ImageViewer from "@/components/ImageViewer.vue";
import { responsiveModule } from "@/modules/responsive";
import { mainModule } from "@/modules/main";
import { BroadcastManager } from "@/utils/BroadcastManager";
import { SET_MASK_SETTINGS } from "@/modules/settings/events";

@Component({
  components: { ImageViewer, IntensityView },
})
export default class BlendTabNew extends Vue {
  readonly mainContext = mainModule.context(this.$store);
  readonly experimentContext = experimentModule.context(this.$store);
  readonly analysisContext = analysisModule.context(this.$store);
  readonly settingsContext = settingsModule.context(this.$store);
  readonly datasetContext = datasetModule.context(this.$store);
  readonly responsiveContext = responsiveModule.context(this.$store);

  get showOptions() {
    return this.mainContext.getters.showOptions;
  }

  get canvasWidth() {
    return this.showOptions
      ? this.responsiveContext.getters.responsive.width! - 840
      : this.responsiveContext.getters.responsive.width! - 460;
  }

  get canvasHeight() {
    return this.responsiveContext.getters.responsive.height! - 174;
  }

  get colorizeMaskInProgress() {
    return this.experimentContext.getters.colorizeMaskInProgress;
  }

  get applyMask() {
    return this.settingsContext.getters.maskSettings.apply;
  }

  set applyMask(value: boolean) {
    BroadcastManager.publish(SET_MASK_SETTINGS, {
      ...this.settingsContext.getters.maskSettings,
      apply: value,
    });
    this.experimentContext.actions.getChannelStackImage();
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

  get activeExperimentId() {
    return this.experimentContext.getters.activeExperimentId;
  }

  get activeAcquisition() {
    return this.experimentContext.getters.activeAcquisition;
  }

  get dataset() {
    return this.datasetContext.getters.activeDataset;
  }

  get hasMask() {
    let hasMask = false;
    if (this.activeAcquisition && this.dataset && this.dataset.input && this.dataset.input.probability_masks) {
      hasMask = !!this.dataset.input.probability_masks[this.activeAcquisition.id];
    }
    return hasMask;
  }

  calculateRegionStats() {
    if (!this.activeExperimentId || !this.activeAcquisition || !this.selectedRegion) {
      return;
    }
    const polygon = this.selectedRegion.getGeometry()! as Polygon;
    const coords = polygon.getCoordinates()[0];
    const height = this.activeAcquisition.max_y;
    // TODO: Y axis flip
    const convertedCoords = coords.map((point) => [point[0], Math.abs(point[1] - height)]);
    const params: IRegionStatsSubmission = {
      experiment_id: this.activeExperimentId,
      acquisition_id: this.activeAcquisition.id,
      region_polygon: convertedCoords,
    };
    return this.analysisContext.actions.calculateRegionStats(params);
  }

  exportImage(format: ExportFormat) {
    this.experimentContext.actions.exportChannelStackImage(format);
  }

  getColorizedMaskImage() {
    this.experimentContext.actions.getColorizedMaskImage();
  }

  deleteRegions() {
    (this.$refs.blendView as any).deleteRegions();
  }
}
</script>
