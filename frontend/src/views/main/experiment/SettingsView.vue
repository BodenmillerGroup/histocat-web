<template>
  <v-card tile>
    <v-card-title class="card-title"><h4>Settings</h4></v-card-title>
    <v-tabs v-model="tabSettings">
      <v-tab>Channels</v-tab>
      <v-tab>General</v-tab>
      <v-tab-item class="scroll-y scroll-view">
        <ChannelSettingsView v-for="channel in selectedChannels" :key="channel.id" :channel="channel"/>
      </v-tab-item>
      <v-tab-item lazy class="scroll-y scroll-view">
        <GeneralSettingsView/>
      </v-tab-item>
    </v-tabs>
  </v-card>
</template>

<script lang="ts">
  import { experimentModule } from '@/modules/experiment';
  import ChannelSettingsView from '@/views/main/experiment/ChannelSettingsView.vue';
  import GeneralSettingsView from '@/views/main/experiment/GeneralSettingsView.vue';
  import { Component, Vue } from 'vue-property-decorator';

  @Component({
    components: { GeneralSettingsView, ChannelSettingsView },
  })
  export default class SettingsView extends Vue {
    readonly experimentContext = experimentModule.context(this.$store);

    tabSettings = 0;

    get selectedChannels() {
      return this.experimentContext.getters.selectedChannels;
    }
  }
</script>

<style scoped>
  .card-title {
    padding-bottom: 4px;
  }

  .scroll-view {
    height: calc(50vh - 121px);
  }
</style>