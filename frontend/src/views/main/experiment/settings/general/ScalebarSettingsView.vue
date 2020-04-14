<template>
  <v-expansion-panel>
    <v-expansion-panel-header>Scalebar</v-expansion-panel-header>
    <v-expansion-panel-content class="ma-0 pa-0">
      <v-switch v-model="apply" label="Show Scalebar" hide-details inset class="ma-0 pa-0"></v-switch>
      <v-text-field
        type="number"
        :rules="[required]"
        min="0"
        step="0.1"
        label="Scale"
        v-model.number="scale"
        persistent-hint
        hint="1px to μm"
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
export default class ScalebarSettingsView extends Vue {
  readonly settingsContext = settingsModule.context(this.$store);
  readonly experimentContext = experimentModule.context(this.$store);

  readonly required = required;

  get apply() {
    return this.settingsContext.getters.scalebar.apply;
  }

  set apply(value: boolean) {
    this.settingsContext.mutations.setScalebar({
      ...this.settingsContext.getters.scalebar,
      apply: value,
    });
    this.experimentContext.actions.getChannelStackImage();
  }

  get scale() {
    return this.settingsContext.getters.scalebar.settings.scale
      ? this.settingsContext.getters.scalebar.settings.scale
      : 1.0;
  }

  set scale(value: number) {
    this.settingsContext.mutations.setScalebar({
      ...this.settingsContext.getters.scalebar,
      settings: {
        scale: value,
      },
    });
    if (this.apply) {
      this.experimentContext.actions.getChannelStackImage();
    }
  }
}
</script>
