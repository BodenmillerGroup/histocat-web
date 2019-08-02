<template>
  <v-flex column>
    <v-toolbar dense flat>
      <v-menu offset-y>
        <template v-slot:activator="{ on }">
          <v-btn
            text
            small
            color="primary"
            v-on="on"
          >
            <v-icon left>mdi-cloud-download-outline</v-icon>
            Export image
          </v-btn>
        </template>
        <v-list dense>
          <v-list-item
            @click="download('tiff')"
          >
            <v-list-item-title>Export TIFF</v-list-item-title>
          </v-list-item>
          <v-list-item
            @click="download('png')"
          >
            <v-list-item-title>Export PNG</v-list-item-title>
          </v-list-item>
        </v-list>
      </v-menu>
      <v-divider
        class="mr-4"
        inset
        vertical
      ></v-divider>
      <v-switch
        v-model="showSettings"
        label="Show settings"
        hide-details
      ></v-switch>
    </v-toolbar>
    <v-layout row>
      <v-flex :class="mainClass">
        <SegmentationView class="segmentation-view" ref="segmentationView"></SegmentationView>
      </v-flex>
      <v-flex v-if="showSettings" md4>
        <SegmentationSettingsView></SegmentationSettingsView>
      </v-flex>
    </v-layout>
  </v-flex>
</template>

<script lang="ts">
  import { experimentModule } from '@/modules/experiment';
  import { ExportTypes } from '@/modules/experiment/models';
  import { mainModule } from '@/modules/main';
  import SegmentationSettingsView from '@/views/main/experiment/analysis/segmentation/SegmentationSettingsView.vue';
  import SegmentationView from '@/views/main/experiment/analysis/segmentation/SegmentationView.vue';
  import { Component, Vue } from 'vue-property-decorator';

  @Component({
    components: { SegmentationSettingsView, SegmentationView },
  })
  export default class SegmentationTab extends Vue {
    readonly mainContext = mainModule.context(this.$store);
    readonly experimentContext = experimentModule.context(this.$store);

    showSettings = true;

    get mainClass() {
      if (this.showSettings) {
        return 'md8';
      }
      return 'md12';
    }

    download(type: ExportTypes) {
      this.experimentContext.actions.exportChannelStackImage(type);
    }
  }
</script>

<style scoped>
  .segmentation-view {
    height: calc(100vh - 162px);
  }
</style>
