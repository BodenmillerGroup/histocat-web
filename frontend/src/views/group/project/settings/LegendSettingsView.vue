<template>
  <v-expansion-panel>
    <v-expansion-panel-header>Legend</v-expansion-panel-header>
    <v-expansion-panel-content class="ma-0 pa-0">
      <v-switch v-model="apply" label="Show Legend" dense hide-details inset class="ma-0 pa-0" />
      <v-text-field
        type="number"
        label="Font Size (pt)"
        v-model.number="legendFontScale"
        :rules="[required]"
        min="0"
        step="1"
        hide-details
      />
      <v-switch v-model="showIntensity" label="Show Intensity" dense hide-details inset />
    </v-expansion-panel-content>
  </v-expansion-panel>
</template>

<script lang="ts">
import { settingsModule } from "@/modules/settings";
import { Component, Vue } from "vue-property-decorator";
import { required } from "@/utils/validators";

@Component
export default class LegendSettingsView extends Vue {
  readonly settingsContext = settingsModule.context(this.$store);

  readonly required = required;

  get apply() {
    return this.settingsContext.getters.legend.apply;
  }

  set apply(value: boolean) {
    this.settingsContext.mutations.setLegend({
      ...this.settingsContext.getters.legend,
      apply: value,
    });
  }

  get legendFontScale() {
    return this.settingsContext.getters.legend.fontScale;
  }

  set legendFontScale(value: number) {
    this.settingsContext.mutations.setLegend({
      ...this.settingsContext.getters.legend,
      fontScale: value,
    });
  }

  get showIntensity() {
    return this.settingsContext.getters.legend.showIntensity;
  }

  set showIntensity(value: boolean) {
    this.settingsContext.mutations.setLegend({
      ...this.settingsContext.getters.legend,
      showIntensity: value,
    });
  }
}
</script>
