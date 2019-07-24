<template>
  <v-card tile id="settings-container">
    <v-tabs v-model="tabSettings">
      <v-tab>Channels</v-tab>
      <v-tab>General</v-tab>
      <v-tab-item class="overflow-y-auto scroll-view">
        <v-expansion-panels>
          <ChannelSettingsView v-for="channel in selectedChannels" :key="channel.id" :channel="channel"/>
        </v-expansion-panels>
      </v-tab-item>
      <v-tab-item class="overflow-y-auto scroll-view">
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
  .scroll-view {
    height: calc(50vh - 92px);
  }
</style>
