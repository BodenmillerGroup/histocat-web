<template>
  <v-card tile>
    <v-card-title>Settings</v-card-title>
    <v-card-text>
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
    </v-card-text>
  </v-card>
</template>

<script lang="ts">
  import { experimentModule } from '@/modules/experiment';
  import ChannelSettingsView from '@/views/main/experiment/settings/channel/ChannelSettingsView.vue';
  import GeneralSettingsView from '@/views/main/experiment/settings/general/GeneralSettingsView.vue';
  import { Component, Vue } from 'vue-property-decorator';

  @Component({
    components: { GeneralSettingsView, ChannelSettingsView },
  })
  export default class SegmentationSettingsView extends Vue {
    readonly experimentContext = experimentModule.context(this.$store);

    tabSettings = 0;

    get selectedChannels() {
      return this.experimentContext.getters.selectedChannels;
    }
  }
</script>

<style scoped>
  .channel-settings-view {
    height: calc(50vh - 92px);
  }
</style>
