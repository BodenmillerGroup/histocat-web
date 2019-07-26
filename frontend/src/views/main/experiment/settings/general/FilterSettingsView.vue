<template>
  <v-expansion-panel>
    <v-expansion-panel-header>
      <v-switch
        v-model="apply"
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
        v-if="filterType === 'gaussian'"
        type="number"
        min="0"
        step="0.1"
        label="Sigma"
        v-model="sigma"
        :rules="[required]"
        hide-details
      ></v-text-field>
      <v-select
        :items="modes"
        v-model="mode"
        label="Mode"
        hide-details
      ></v-select>
    </v-expansion-panel-content>
  </v-expansion-panel>
</template>

<script lang="ts">
  import { experimentModule } from '@/modules/experiment';
  import { settingsModule } from '@/modules/settings';
  import { Component, Vue } from 'vue-property-decorator';

  @Component({
    components: {},
  })
  export default class FilterSettingsView extends Vue {
    readonly settingsContext = settingsModule.context(this.$store);
    readonly experimentContext = experimentModule.context(this.$store);

    filterTypes = ['gaussian', 'median'];
    modes = ['reflect', 'constant', 'nearest', 'mirror', 'wrap'];

    required = value => !!value || 'Required';

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
      return this.settingsContext.getters.filter.settings.sigma ? this.settingsContext.getters.filter.settings.sigma : 1.0;
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

    get mode() {
      return this.settingsContext.getters.filter.settings.mode ? this.settingsContext.getters.filter.settings.mode : 'nearest';
    }

    set mode(value: string) {
      this.settingsContext.mutations.setFilter({
        ...this.settingsContext.getters.filter,
        settings: {
          mode: value,
        },
      });
      if (this.apply) {
        this.experimentContext.actions.getChannelStackImage();
      }
    }
  }
</script>
