<template>
  <v-card tile>
    <v-card-title><h4>Settings</h4></v-card-title>
    <v-divider></v-divider>
    <v-expansion-panel
      v-model="panel"
      expand
      class="scroll-y local-height"
    >
      <v-expansion-panel-content v-for="channel in selectedChannels" :key="channel.id">
        <template v-slot:header>
          <div>{{channel.name}}</div>
        </template>
        <ChannelSettingsView/>
      </v-expansion-panel-content>
    </v-expansion-panel>
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

    panel = [];

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

<style scoped>
  .local-height {
    max-height: 35vh;
  }
</style>
