<template>
  <v-expansion-panels
    v-model="panel"
    multiple
    accordion
  >
    <v-expansion-panel>
      <v-expansion-panel-header>
        <v-switch
          v-model="legendApply"
          label="Show Legend"
          hide-details
        ></v-switch>
      </v-expansion-panel-header>
    </v-expansion-panel>
    <v-expansion-panel>
      <v-expansion-panel-header>
        <v-switch
          v-model="scalebarApply"
          label="Show Scalebar"
          hide-details
        ></v-switch>
      </v-expansion-panel-header>
      <v-expansion-panel-content>
        <v-text-field
          type="number"
          min="0"
          label="Scale"
          v-model="scale"
          persistent-hint
          hint="1px to μm"
        ></v-text-field>
      </v-expansion-panel-content>
    </v-expansion-panel>
    <v-expansion-panel>
      <v-expansion-panel-header>
        <v-switch
          v-model="filterApply"
          label="Apply Filter"
          hide-details
        ></v-switch>
      </v-expansion-panel-header>
      <v-expansion-panel-content>
        <v-select
          :items="filterTypes"
          v-model="filterType"
          label="Filter Type"
          hide-details
        ></v-select>
        <v-text-field
          type="number"
          label="Sigma"
          v-model="sigma"
          hide-details
        ></v-text-field>
      </v-expansion-panel-content>
    </v-expansion-panel>
  </v-expansion-panels>
</template>

<script lang="ts">
  import { experimentModule } from '@/modules/experiment';
  import { settingsModule } from '@/modules/settings';
  import ChannelHistogramView from '@/views/main/experiment/ChannelHistogramView.vue';
  import { Component, Vue } from 'vue-property-decorator';

  @Component({
    components: { ChannelHistogramView },
  })
  export default class GeneralSettingsView extends Vue {
    readonly settingsContext = settingsModule.context(this.$store);
    readonly experimentContext = experimentModule.context(this.$store);

    panel = [0];
    filterTypes = ['gaussian'];

    get legendApply() {
      return this.settingsContext.getters.legend.apply;
    }

    set legendApply(value: boolean) {
      this.settingsContext.mutations.setLegend({
        ...this.settingsContext.getters.legend,
        apply: value,
      });
      this.experimentContext.actions.getChannelStackImage();
    }

    get scalebarApply() {
      return this.settingsContext.getters.scalebar.apply;
    }

    set scalebarApply(value: boolean) {
      this.settingsContext.mutations.setScalebar({
        ...this.settingsContext.getters.scalebar,
        apply: value,
      });
      this.experimentContext.actions.getChannelStackImage();
    }

    get filterApply() {
      return this.settingsContext.getters.filter.apply;
    }

    set filterApply(value: boolean) {
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
      this.experimentContext.actions.getChannelStackImage();
    }

    get sigma() {
      return this.settingsContext.getters.filter.settings.sigma ? this.settingsContext.getters.filter.settings.sigma : 1.0;
    }

    set sigma(value: number) {
      this.settingsContext.mutations.setFilter({
        ...this.settingsContext.getters.filter,
        settings: {
          sigma: value,
        },
      });
      this.experimentContext.actions.getChannelStackImage();
    }

    get scale() {
      return this.settingsContext.getters.scalebar.settings.scale ? this.settingsContext.getters.scalebar.settings.scale : 1;
    }

    set scale(value: number) {
      this.settingsContext.mutations.setScalebar({
        ...this.settingsContext.getters.scalebar,
        settings: {
          scale: value,
        },
      });
      this.experimentContext.actions.getChannelStackImage();
    }
  }
</script>

<style scoped>
  .expansion-panel-header {

  }
</style>