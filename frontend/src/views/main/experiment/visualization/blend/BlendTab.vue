<template>
  <v-row no-gutters>
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
        <v-switch
          v-model="applyMask"
          label="Mask overlay"
          hide-details
          class="ml-2"
        ></v-switch>
      </v-toolbar>
      <v-row no-gutters>
        <v-col>
          <BlendView class="blend-view"/>
        </v-col>
        <IntensityView/>
      </v-row>
    </v-col>
  </v-row>
</template>

<script lang="ts">
  import { experimentModule } from '@/modules/experiment';
  import { ExportFormat } from '@/modules/experiment/models';
  import { settingsModule } from '@/modules/settings';
  import BlendView from '@/views/main/experiment/visualization/blend/BlendView.vue';
  import IntensityView from '@/views/main/experiment/visualization/blend/IntensityView.vue';
  import { Component, Vue } from 'vue-property-decorator';

  @Component({
    components: { IntensityView, BlendView },
  })
  export default class BlendTab extends Vue {
    readonly experimentContext = experimentModule.context(this.$store);
    readonly settingsContext = settingsModule.context(this.$store);

    get colorizeMaskInProgress() {
      return this.experimentContext.getters.colorizeMaskInProgress;
    }

    get applyMask() {
      return this.settingsContext.getters.maskSettings.apply;
    }

    set applyMask(value: boolean) {
      this.settingsContext.mutations.setMaskSettings({
        ...this.settingsContext.getters.maskSettings,
        apply: value,
      });
      this.experimentContext.actions.getChannelStackImage();
    }

    exportImage(format: ExportFormat) {
      this.experimentContext.actions.exportChannelStackImage(format);
    }

    getColorizedMaskImage() {
      this.experimentContext.actions.getColorizedMaskImage();
    }
  }
</script>

<style scoped>
  .blend-view {
    height: calc(100vh - 204px);
  }
</style>
