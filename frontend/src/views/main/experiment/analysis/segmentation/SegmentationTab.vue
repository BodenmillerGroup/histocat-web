<template>
  <v-flex column>
    <v-toolbar dense flat>
      <v-menu offset-y>
        <template v-slot:activator="{ on }">
          <v-btn
            v-on="on"
            small
            color="primary lighten-2"
            elevation="1"
          >
            <v-icon left small>mdi-download</v-icon>
            Export image
          </v-btn>
        </template>
        <v-list dense>
          <v-list-item
            @click="exportImage('tiff')"
          >
            <v-list-item-title>Export TIFF</v-list-item-title>
          </v-list-item>
          <v-list-item
            @click="exportImage('png')"
          >
            <v-list-item-title>Export PNG</v-list-item-title>
          </v-list-item>
        </v-list>
      </v-menu>
      <v-switch
        v-model="showSettings"
        label="Show settings"
        hide-details
        class="ml-2"
      ></v-switch>
    </v-toolbar>
    <v-layout row>
      <v-flex :class="mainClass">
        <keep-alive>
          <SegmentationView
            class="segmentation-view"
            ref="segmentationView"
          ></SegmentationView>
        </keep-alive>
      </v-flex>
      <v-flex
        v-if="showSettings"
        md4
      >
        <SegmentationSettingsView/>
      </v-flex>
    </v-layout>
  </v-flex>
</template>

<script lang="ts">
  import { analysisModule } from '@/modules/analysis';
  import { ExportFormat } from '@/modules/experiment/models';
  import { settingsModule } from '@/modules/settings';
  import SegmentationSettingsView from '@/views/main/experiment/analysis/segmentation/SegmentationSettingsView.vue';
  import SegmentationView from '@/views/main/experiment/analysis/segmentation/SegmentationView.vue';
  import { Component, Vue } from 'vue-property-decorator';

  @Component({
    components: { SegmentationSettingsView, SegmentationView },
  })
  export default class SegmentationTab extends Vue {
    readonly analysisContext = analysisModule.context(this.$store);
    readonly settingsContext = settingsModule.context(this.$store);

    showSettings = true;

    get mainClass() {
      if (this.showSettings) {
        return 'md8';
      }
      return 'md12';
    }

    exportImage(format: ExportFormat) {
      const settings = this.settingsContext.getters.segmentationSettings;
      this.analysisContext.actions.exportSegmentationImage({
        format: format,
        settings: settings,
      });
    }
  }
</script>

<style scoped>
  .segmentation-view {
    height: calc(100vh - 212px);
  }
</style>
