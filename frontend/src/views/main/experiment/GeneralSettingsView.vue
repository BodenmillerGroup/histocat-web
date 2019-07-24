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
        ></v-switch>
      </v-expansion-panel-header>
    </v-expansion-panel>
    <v-expansion-panel>
      <v-expansion-panel-header>
        <v-switch
          v-model="filterApply"
          label="Apply Filter"
        ></v-switch>
      </v-expansion-panel-header>
      <v-expansion-panel-content>
        <v-select
          :items="filterTypes"
          v-model="filterType"
          label="Filter Type"
        ></v-select>
        <v-text-field
          type="number"
          label="Sigma"
          v-model="sigma"
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
      return this.settingsContext.getters.filter.settings.sigma ? this.settingsContext.getters.filter.settings.sigma : 1;
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
  }
</script>
