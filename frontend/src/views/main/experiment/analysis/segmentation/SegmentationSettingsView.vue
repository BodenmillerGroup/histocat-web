<template>
  <v-card tile>
    <v-card-title>Segmentation Settings</v-card-title>
    <v-card-text>
      <v-select
        class="input-row"
        :items="resultTypes"
        v-model="resultType"
        label="Result image"
        hint="Type of result image."
        persistent-hint
      ></v-select>
      <v-select
        class="input-row"
        :items="algorithms"
        v-model="algorithm"
        label="Algorithm"
        hint="The masking algorithm to use for image pre-labelling."
        persistent-hint
      ></v-select>
      <v-text-field
        class="input-row"
        type="number"
        min="0"
        step="1"
        label="Iterations"
        v-model="iterations"
        :rules="[required]"
        hint="Morphology iterations."
        persistent-hint
      ></v-text-field>
      <v-text-field
        class="input-row"
        type="number"
        min="0"
        step="1"
        label="Kernel size"
        v-model="kernelSize"
        :rules="[required]"
        hide-details
      ></v-text-field>
      <v-color-picker
        v-model="color"
        hide-canvas
        hide-inputs
      ></v-color-picker>
    </v-card-text>
    <v-card-actions>
      <v-btn
        @click="process"
        color="primary"
      >
        Process
      </v-btn>
      <v-btn
        @click="produceContours"
        color="primary"
      >
        Produce contours
      </v-btn>
    </v-card-actions>
  </v-card>
</template>

<script lang="ts">
  import { analysisModule } from '@/modules/analysis';
  import { ImageResultType } from '@/modules/analysis/models';
  import { experimentModule } from '@/modules/experiment';
  import { settingsModule } from '@/modules/settings';
  import { required } from '@/utils';
  import { Component, Vue } from 'vue-property-decorator';

  @Component
  export default class SegmentationSettingsView extends Vue {
    readonly experimentContext = experimentModule.context(this.$store);
    readonly analysisContext = analysisModule.context(this.$store);
    readonly settingsContext = settingsModule.context(this.$store);

    readonly resultTypes: ImageResultType[] = ['origin', 'mask'];
    readonly algorithms = ['Otsu Grayscale', 'Otsu Hue', 'Otsu Saturation', 'Otsu Lightness'];
    resultType: ImageResultType = 'origin';
    algorithm = 'Otsu Hue';
    iterations = 1;
    kernelSize = 3;
    color = '#00AAFF40';

    readonly required = required;

    get selectedChannels() {
      return this.experimentContext.getters.selectedChannels;
    }

    process() {
      const settings = {
        result_type: this.resultType,
        algorithm: this.algorithm,
        iterations: this.iterations,
        kernel_size: this.kernelSize,
        mask_color: this.color,
      };
      this.settingsContext.mutations.setSegmentationSettings(settings);
      this.analysisContext.actions.getSegmentationImage(settings);
    }

    produceContours() {
      const settings = {
        result_type: this.resultType,
        algorithm: this.algorithm,
        iterations: this.iterations,
        kernel_size: this.kernelSize,
        mask_color: this.color,
      };
      this.settingsContext.mutations.setSegmentationSettings(settings);
      this.analysisContext.actions.produceSegmentationContours(settings);
    }
  }
</script>

<style scoped>
  .input-row {
    margin-bottom: 32px;
  }
</style>
