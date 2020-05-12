<template>
  <v-expansion-panel>
    <v-expansion-panel-header>Legend</v-expansion-panel-header>
    <v-expansion-panel-content class="ma-0 pa-0">
      <v-switch v-model="apply" label="Show Legend" hide-details inset class="ma-0 pa-0"></v-switch>
      <v-text-field
        type="number"
        label="Font Scale"
        v-model.number="legendFontScale"
        :rules="[required]"
        min="0"
        step="0.05"
        hide-details
      ></v-text-field>
      <v-switch v-model="showIntensity" label="Show Intensity" hide-details inset></v-switch>
    </v-expansion-panel-content>
  </v-expansion-panel>
</template>

<script lang="ts">
import { experimentModule } from "@/modules/experiment";
import { settingsModule } from "@/modules/settings";
import { Component, Vue } from "vue-property-decorator";
import { required } from "@/utils/validators";

@Component
export default class LegendSettingsView extends Vue {
  readonly settingsContext = settingsModule.context(this.$store);
  readonly experimentContext = experimentModule.context(this.$store);

  readonly required = required;

  get apply() {
    return this.settingsContext.getters.legend.apply;
  }

  set apply(value: boolean) {
    this.settingsContext.mutations.setLegend({
      ...this.settingsContext.getters.legend,
      apply: value,
    });
    this.experimentContext.actions.getChannelStackImage();
  }

  get legendFontScale() {
    return this.settingsContext.getters.legend.fontScale;
  }

  set legendFontScale(value: number) {
    this.settingsContext.mutations.setLegend({
      ...this.settingsContext.getters.legend,
      fontScale: value,
    });
    if (this.apply) {
      this.experimentContext.actions.getChannelStackImage();
    }
  }

  get showIntensity() {
    return this.settingsContext.getters.legend.showIntensity;
  }

  set showIntensity(value: boolean) {
    this.settingsContext.mutations.setLegend({
      ...this.settingsContext.getters.legend,
      showIntensity: value,
    });
    this.experimentContext.actions.getChannelStackImage();
  }
}
</script>
