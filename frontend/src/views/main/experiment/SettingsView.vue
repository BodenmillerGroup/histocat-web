<template>
  <v-card tile class="scroll-y">
    <v-card-title><h4>Settings</h4></v-card-title>
    <v-divider></v-divider>
    <ChannelSettingsView v-for="channel in selectedChannels" :key="channel.id" :channel="channel"/>
  </v-card>
</template>

<script lang="ts">
  import { Component, Vue, Watch } from 'vue-property-decorator';
  import { readSelectedAcquisition, readSelectedChannels } from '@/modules/experiment/getters';
  import { IChannel } from '@/modules/experiment/models';
  import ChannelSettingsView from '@/views/main/experiment/ChannelSettingsView.vue';

  @Component({
    components: { ChannelSettingsView },
  })
  export default class SettingsView extends Vue {

    get selectedChannels() {
      return readSelectedChannels(this.$store);
    }

    @Watch('selectedChannels')
    onSelectedChannelsChanged(channels: IChannel[]) {
      const acquisition = readSelectedAcquisition(this.$store);
      if (!acquisition) {
        return;
      }
    }
  }
</script>
