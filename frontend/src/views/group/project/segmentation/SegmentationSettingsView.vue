<template>
  <div class="px-2 pt-1">
    <h4>Settings</h4>
    <v-form v-model="valid" ref="form">
      <v-tabs v-model="tab">
        <v-tab>General</v-tab>
        <v-tab>Pre</v-tab>
        <v-tab>Post</v-tab>
        <v-tab-item>
          <v-select
            label="Model"
            v-model="modelId"
            :items="models"
            item-value="id"
            item-text="name"
            :rules="modelRules"
          />
          <v-radio-group v-model="compartment" label="Compartment" dense>
            <v-radio label="Whole-Cell" value="whole-cell" />
            <v-radio label="Nuclear" value="nuclear" />
          </v-radio-group>
          <v-text-field label="Dataset Name" hint="Resulting dataset name" v-model="datasetName" />
          <v-text-field label="Dataset Description" hint="Resulting dataset description" v-model="datasetDescription" />
        </v-tab-item>
        <v-tab-item>
          <v-radio-group v-model="channelsNormalization" label="Per-channel normalization" dense>
            <v-radio label="None" value="none" />
            <v-radio label="Min-Max" value="minmax" />
            <v-radio label="Z-score" value="zscore" />
          </v-radio-group>
          <v-switch v-model="threshold" label="Threshold" inset class="ml-1" dense />
          <v-text-field
            label="Percentile"
            v-model.number="percentile"
            type="number"
            min="0"
            :rules="percentileRules"
            :disabled="!threshold"
          />
          <v-switch v-model="normalize" label="Normalize" inset class="ml-1" dense />
          <v-text-field
            label="Kernel Size"
            v-model.number="kernelSize"
            type="number"
            min="0"
            step="1"
            :rules="kernelSizeRules"
            :disabled="!normalize"
          />
        </v-tab-item>
        <v-tab-item>
          <v-text-field
            label="Radius"
            hint="Radius of disk used to search for maxima"
            v-model.number="radius"
            type="number"
            min="0"
            step="1"
            :rules="radiusRules"
          />
          <v-text-field
            label="Maxima Threshold"
            hint="Threshold for the maxima prediction"
            v-model.number="maximaThreshold"
            type="number"
            min="0"
            :rules="maximaThresholdRules"
          />
          <v-text-field
            label="Interior Threshold"
            hint="Threshold for the interior prediction"
            v-model.number="interiorThreshold"
            type="number"
            min="0"
            :rules="interiorThresholdRules"
          />
          <v-text-field
            label="Small Objects Threshold"
            hint="Removes objects smaller than this size"
            v-model.number="smallObjectsThreshold"
            type="number"
            min="0"
            step="1"
            :rules="smallObjectsThresholdRules"
          />
          <v-text-field
            label="Fill Holes Threshold"
            hint="Maximum size for holes within segmented objects to be filled"
            v-model.number="fillHolesThreshold"
            type="number"
            min="0"
            step="1"
            :rules="fillHolesThresholdRules"
          />
          <v-select
            :items="validModelNames"
            v-model="interiorModel"
            label="Interior Model"
            hint="Semantic head to use to predict interior of each object"
            :rules="interiorModelRules"
          />
          <v-select
            :items="validModelNames"
            v-model="maximaModel"
            label="Maxima Model"
            hint="Semantic head to use to predict maxima of each object"
            :rules="maximaModelRules"
          />
          <v-text-field
            label="Interior Model Smooth"
            hint="Smoothing factor to apply to interior model predictions"
            v-model.number="interiorModelSmooth"
            type="number"
            min="0"
            step="1"
            :rules="interiorModelSmoothRules"
          />
          <v-text-field
            label="Maxima Model Smooth"
            hint="Smoothing factor to apply to maxima model predictions"
            v-model.number="maximaModelSmooth"
            type="number"
            min="0"
            step="1"
            :rules="maximaModelSmoothRules"
          />
          <v-text-field
            label="Pixel Expansion"
            hint="Optional number of pixels to expand segmentation labels"
            v-model.number="pixelExpansion"
            type="number"
            min="0"
            step="1"
            clearable
            :rules="pixelExpansionRules"
          />
        </v-tab-item>
      </v-tabs>
    </v-form>
    <v-btn @click="submit" color="primary" block>Submit</v-btn>
  </div>
</template>

<script lang="ts">
import { projectsModule } from "@/modules/projects";
import { Component, Vue } from "vue-property-decorator";
import { segmentationModule } from "@/modules/segmentation";
import { modelsModule } from "@/modules/models";
import { required, positiveNumber, percentNumber, nonNegativeNumber } from "@/utils/validators";

@Component
export default class SegmentationSettingsView extends Vue {
  readonly projectsContext = projectsModule.context(this.$store);
  readonly segmentationContext = segmentationModule.context(this.$store);
  readonly modelsContext = modelsModule.context(this.$store);

  readonly modelRules = [required];
  readonly percentileRules = [percentNumber];
  readonly kernelSizeRules = [positiveNumber];
  readonly radiusRules = [positiveNumber];
  readonly maximaThresholdRules = [nonNegativeNumber];
  readonly interiorThresholdRules = [nonNegativeNumber];
  readonly smallObjectsThresholdRules = [nonNegativeNumber];
  readonly fillHolesThresholdRules = [nonNegativeNumber];
  readonly interiorModelRules = [required];
  readonly maximaModelRules = [required];
  readonly interiorModelSmoothRules = [nonNegativeNumber];
  readonly maximaModelSmoothRules = [nonNegativeNumber];
  readonly pixelExpansionRules = [(value) => value === null || value >= 0 || "Should be null or non-negative number"];

  readonly validModelNames = ["inner-distance", "outer-distance", "fgbg-fg", "pixelwise-interior"];

  valid = false;
  tab = 0;

  datasetName: string | null = null;
  datasetDescription: string | null = null;
  modelId: number | null = null;
  compartment = "whole-cell";

  channelsNormalization = "none";
  threshold = true;
  percentile = 99.9;
  normalize = true;
  kernelSize = 128;

  radius = 2;
  maximaThreshold = 0.1;
  interiorThreshold = 0.3;
  smallObjectsThreshold = 15;
  fillHolesThreshold = 15;
  interiorModel = "pixelwise-interior";
  maximaModel = "inner-distance";
  interiorModelSmooth = 2;
  maximaModelSmooth = 0;
  pixelExpansion: number | null = null;

  get selectedAcquisitionIds() {
    return this.segmentationContext.getters.selectedAcquisitionIds;
  }

  get projectData() {
    return this.projectsContext.getters.projectData!;
  }

  get models() {
    return this.modelsContext.getters.models;
  }

  async submit() {
    if ((this.$refs.form as any).validate()) {
      const acquisitionIds = this.segmentationContext.getters.selectedAcquisitionIds;
      const channels = this.segmentationContext.getters.channels;
      const nucleiChannels = this.segmentationContext.getters.nucleiChannels;
      const cytoplasmChannels = this.segmentationContext.getters.cytoplasmChannels;

      if (nucleiChannels.length === 0 || cytoplasmChannels.length === 0) {
        self.alert("Select at least one nuclear and one cytoplasm channel!");
        return;
      }

      const preprocessingParams = {
        channels_normalization: this.channelsNormalization,
        threshold: this.threshold,
        percentile: this.percentile,
        normalize: this.normalize,
        kernel_size: this.kernelSize,
      };

      const postprocessingParams = {
        radius: this.radius,
        maxima_threshold: this.maximaThreshold,
        interior_threshold: this.interiorThreshold,
        small_objects_threshold: this.smallObjectsThreshold,
        fill_holes_threshold: this.fillHolesThreshold,
        interior_model: this.interiorModel,
        maxima_model: this.maximaModel,
        interior_model_smooth: this.interiorModelSmooth,
        maxima_model_smooth: this.maximaModelSmooth,
        pixel_expansion: this.pixelExpansion,
      };

      await this.segmentationContext.actions.processSegmentation({
        dataset_name: this.datasetName,
        dataset_description: this.datasetDescription,
        model_id: this.modelId!,
        compartment: this.compartment,
        acquisition_ids: acquisitionIds,
        channels: channels,
        nuclei_channels: nucleiChannels,
        cytoplasm_channels: cytoplasmChannels,
        preprocessing: preprocessingParams,
        postprocessing: postprocessingParams,
      });
    }
  }

  async mounted() {
    await this.modelsContext.actions.getModels(+this.$router.currentRoute.params.groupId);
  }
}
</script>
