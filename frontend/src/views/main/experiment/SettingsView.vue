<template>
  <v-card tile class="scroll-y">
    <v-card-title><h4>Settings</h4></v-card-title>
    <v-divider></v-divider>
    <ChannelSettingsView v-for="channel in selectedChannels" :key="channel.id" :channel="channel"/>
  </v-card>
</template>

<script lang="ts">
  import { Component, Vue } from 'vue-property-decorator';
  import { readSelectedChannels } from '@/modules/experiment/getters';
  import ChannelSettingsView from '@/views/main/experiment/ChannelSettingsView.vue';

  @Component({
    components: { ChannelSettingsView },
  })
  export default class SettingsView extends Vue {

    get selectedChannels() {
      return readSelectedChannels(this.$store);
    }
  }
</script>
