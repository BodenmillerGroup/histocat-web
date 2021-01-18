<template>
  <v-card tile class="ma-1">
    <v-card-title>Settings</v-card-title>
    <v-card-text>
      <v-form v-model="valid" ref="form">
        <v-select
          label="Model"
          v-model="modelId"
          :items="models"
          item-value="id"
          item-text="name"
          :rules="modelRules"
        />
        <v-text-field
          label="Scaling factor"
          hint="This parameter is commonly used in the IMCSegmentation pipeline"
          v-model.number="scalingFactor"
          type="number"
          min="0"
          :rules="scalingFactorRules"
        />
        <v-text-field
          label="Upper limit [%]"
          hint="Histogram rescaled from [0-upper_limit]"
          v-model.number="upperLimit"
          type="number"
          min="0"
          max="100"
          :rules="upperLimitRules"
        />
        <v-text-field
          label="Overlap"
          hint="Overlap between tiles"
          v-model.number="overlap"
          type="number"
          min="0"
          :rules="overlapRules"
        />
        <v-switch v-model="applyFilter" label="Apply Filter" inset class="mt-2" />
        <v-select :items="filterTypes" v-model="filterType" label="Filter Type" :disabled="!applyFilter" />
        <v-select
          v-if="filterType === 'median'"
          :items="kernelSizes"
          v-model.number="kernelSize"
          label="Kernel Size"
          :disabled="!applyFilter"
        />
        <v-text-field
          v-if="filterType === 'gaussian'"
          type="number"
          min="0"
          step="0.1"
          label="Sigma"
          v-model.number="sigma"
          :rules="sigmaRules"
          :disabled="!applyFilter"
        />
      </v-form>
    </v-card-text>
    <v-card-actions>
      <v-btn @click="submit" color="primary" block>Submit</v-btn>
    </v-card-actions>
  </v-card>
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
  readonly scalingFactorRules = [positiveNumber];
  readonly upperLimitRules = [percentNumber];
  readonly overlapRules = [nonNegativeNumber];
  readonly sigmaRules = [required];

  readonly filterTypes = ["gaussian", "median"];
  readonly kernelSizes = [3, 5, 7, 9];

  search = "";
  valid = false;

  modelId: number | null = null;
  scalingFactor = 1;
  upperLimit = 98;
  overlap = 0.0;

  applyFilter = true;
  filterType = "gaussian";
  kernelSize = 3;
  sigma = 1.0;

  readonly headers = [
    {
      text: "ID",
      sortable: true,
      filterable: false,
      value: "id",
      align: "end",
      width: 70,
    },
    {
      text: "Slide ID",
      sortable: true,
      filterable: false,
      value: "slide_id",
      align: "end",
      width: 100,
    },
    {
      text: "Descriptions",
      sortable: true,
      value: "description",
    },
  ];

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
      const nucleiChannels = this.segmentationContext.getters.nucleiChannels;
      const cytoplasmChannels = this.segmentationContext.getters.cytoplasmChannels;

      await this.segmentationContext.actions.processSegmentation({
        model_id: this.modelId!,
        acquisition_ids: acquisitionIds,
        nuclei_channels: nucleiChannels,
        cytoplasm_channels: cytoplasmChannels,
        scaling_factor: this.scalingFactor,
        upper_limit: this.upperLimit,
        overlap: this.overlap,
        filter_settings: {
          apply: this.applyFilter,
          type: this.filterType,
          kernel_size: this.kernelSize,
          sigma: this.sigma,
        },
      });
    }
  }

  async mounted() {
    await this.modelsContext.actions.getModels(+this.$router.currentRoute.params.groupId);
  }
}
</script>
