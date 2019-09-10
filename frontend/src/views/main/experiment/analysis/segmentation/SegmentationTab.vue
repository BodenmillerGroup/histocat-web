<template>
  <v-col>
    <v-toolbar dense flat>
      <v-menu offset-y>
        <template v-slot:activator="{ on }">
          <v-btn
            v-on="on"
            small
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
    </v-toolbar>
    <v-row no-gutters>
      <v-col :cols="columns">
        <keep-alive>
          <SegmentationView
            class="segmentation-view"
            ref="segmentationView"
          ></SegmentationView>
        </keep-alive>
      </v-col>
      <v-col
        v-if="showOptions"
        cols="3"
      >
        <SegmentationSettingsView/>
      </v-col>
    </v-row>
  </v-col>
</template>

<script lang="ts">
  import { analysisModule } from '@/modules/analysis';
  import { ExportFormat } from '@/modules/experiment/models';
  import { mainModule } from '@/modules/main';
  import { settingsModule } from '@/modules/settings';
  import SegmentationSettingsView from '@/views/main/experiment/analysis/segmentation/SegmentationSettingsView.vue';
  import SegmentationView from '@/views/main/experiment/analysis/segmentation/SegmentationView.vue';
  import { Component, Vue } from 'vue-property-decorator';

  @Component({
    components: { SegmentationSettingsView, SegmentationView },
  })
  export default class SegmentationTab extends Vue {
    readonly mainContext = mainModule.context(this.$store);
    readonly analysisContext = analysisModule.context(this.$store);
    readonly settingsContext = settingsModule.context(this.$store);

    get showOptions() {
      return this.mainContext.getters.showOptions;
    }

    get columns() {
      return this.showOptions ? 9 : 12;
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
