<template>
  <v-expansion-panel>
    <v-expansion-panel-header>Filter</v-expansion-panel-header>
    <v-expansion-panel-content class="ma-0 pa-0">
      <v-switch v-model="apply" label="Apply Filter" inset hide-details class="ma-0 pa-0"></v-switch>
      <v-select :items="filterTypes" v-model="filterType" label="Filter Type" hide-details></v-select>
      <v-select :items="kernelSizes" v-model.number="kernelSize" label="Kernel Size" hide-details></v-select>
      <v-text-field
        v-if="filterType === 'gaussian'"
        type="number"
        min="0"
        step="0.1"
        label="Sigma"
        v-model.number="sigma"
        :rules="[required]"
        hide-details
      ></v-text-field>
    </v-expansion-panel-content>
  </v-expansion-panel>
</template>

<script lang="ts">
import { experimentModule } from "@/modules/experiment";
import { settingsModule } from "@/modules/settings";
import { Component, Vue } from "vue-property-decorator";
import { required } from "@/utils/validators";

@Component
export default class FilterSettingsView extends Vue {
  readonly settingsContext = settingsModule.context(this.$store);
  readonly experimentContext = experimentModule.context(this.$store);

  filterTypes = ["gaussian", "median"];

  readonly required = required;

  get kernelSizes() {
    return this.filterType === "gaussian" ? [0, 1, 3, 5, 7, 9] : [3, 5, 7, 9];
  }

  get apply() {
    return this.settingsContext.getters.filter.apply;
  }

  set apply(value: boolean) {
    this.settingsContext.mutations.setFilter({
      ...this.settingsContext.getters.filter,
      apply: value,
    });
    this.experimentContext.actions.getChannelStackImage();
  }

  get filterType() {
    return this.settingsContext.getters.filter.type;
  }

  set filterType(value: string) {
    this.settingsContext.mutations.setFilter({
      ...this.settingsContext.getters.filter,
      type: value,
    });
    if (this.apply) {
      this.experimentContext.actions.getChannelStackImage();
    }
  }

  get sigma() {
    return this.settingsContext.getters.filter.settings.sigma
      ? this.settingsContext.getters.filter.settings.sigma
      : 1.0;
  }

  set sigma(value: number) {
    this.settingsContext.mutations.setFilter({
      ...this.settingsContext.getters.filter,
      settings: {
        sigma: value,
      },
    });
    if (this.apply) {
      this.experimentContext.actions.getChannelStackImage();
    }
  }

  get kernelSize() {
    return this.settingsContext.getters.filter.settings.kernel_size
      ? this.settingsContext.getters.filter.settings.kernel_size
      : 1;
  }

  set kernelSize(value: string) {
    this.settingsContext.mutations.setFilter({
      ...this.settingsContext.getters.filter,
      settings: {
        kernel_size: value,
      },
    });
    if (this.apply) {
      this.experimentContext.actions.getChannelStackImage();
    }
  }
}
</script>
