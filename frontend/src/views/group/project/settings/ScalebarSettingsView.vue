<template>
  <v-expansion-panel>
    <v-expansion-panel-header>Scalebar</v-expansion-panel-header>
    <v-expansion-panel-content class="ma-0 pa-0">
      <v-switch v-model="apply" label="Show Scalebar" dense hide-details inset class="ma-0 pa-0" />
      <v-text-field
        type="number"
        :rules="[required]"
        min="0"
        step="0.1"
        label="Scale"
        v-model.number="scale"
        hint="1px to μm"
      />
      <v-text-field
        type="number"
        :rules="[required]"
        min="1"
        max="300"
        step="1"
        label="Length"
        v-model.number="length"
        hint="Scalebar length"
      />
    </v-expansion-panel-content>
  </v-expansion-panel>
</template>

<script lang="ts">
import { projectsModule } from "@/modules/projects";
import { settingsModule } from "@/modules/settings";
import { Component, Vue } from "vue-property-decorator";
import { required } from "@/utils/validators";

@Component
export default class ScalebarSettingsView extends Vue {
  readonly settingsContext = settingsModule.context(this.$store);
  readonly projectsContext = projectsModule.context(this.$store);

  readonly required = required;

  get apply() {
    return this.settingsContext.getters.scalebar.apply;
  }

  set apply(value: boolean) {
    this.settingsContext.mutations.setScalebar({
      ...this.settingsContext.getters.scalebar,
      apply: value,
    });
    this.projectsContext.actions.getChannelStackImage();
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
      this.projectsContext.actions.getChannelStackImage();
    }
  }

  get length() {
    return this.settingsContext.getters.scalebar.settings.length
      ? this.settingsContext.getters.scalebar.settings.length
      : 64;
  }

  set length(value: number) {
    this.settingsContext.mutations.setScalebar({
      ...this.settingsContext.getters.scalebar,
      settings: {
        length: value,
      },
    });
    if (this.apply) {
      this.projectsContext.actions.getChannelStackImage();
    }
  }
}
</script>
