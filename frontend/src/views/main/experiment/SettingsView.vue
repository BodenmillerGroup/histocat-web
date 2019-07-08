<template>
  <v-card tile class="scroll-y">
    <v-card-title><h4>Settings</h4></v-card-title>
    <v-divider></v-divider>
    <ChannelSettingsView v-for="channel in selectedChannels" :key="channel.id" :channel="channel"/>
  </v-card>
</template>

<script lang="ts">
  import { experimentModule } from '@/modules/experiment';
  import ChannelSettingsView from '@/views/main/experiment/ChannelSettingsView.vue';
  import { Component, Vue } from 'vue-property-decorator';

  @Component({
    components: { ChannelSettingsView },
  })
  export default class SettingsView extends Vue {
    readonly experimentContext = experimentModule.context(this.$store);

    get selectedChannels() {
      return this.experimentContext.getters.selectedChannels;
    }
  }
</script>
