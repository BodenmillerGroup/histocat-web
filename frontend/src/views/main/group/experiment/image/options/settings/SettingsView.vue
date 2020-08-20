<template>
  <v-card tile id="settings-container">
    <v-tabs v-model="tabSettings">
      <v-tab>Channels</v-tab>
      <v-tab>General</v-tab>
      <v-tab>Presets</v-tab>
      <v-tab-item class="overflow-y-auto channel-settings-view">
        <BrushableHistogram v-for="channel in selectedChannels" :key="channel.id" :channel="channel" />
      </v-tab-item>
      <v-tab-item class="overflow-y-auto channel-settings-view">
        <GeneralSettingsView />
      </v-tab-item>
      <v-tab-item class="overflow-y-auto channel-settings-view">
        <PresetsView />
      </v-tab-item>
    </v-tabs>
  </v-card>
</template>

<script lang="ts">
import { experimentModule } from "@/modules/experiment";
import GeneralSettingsView from "@/views/main/group/experiment/image/options/settings/general/GeneralSettingsView.vue";
import { Component, Vue } from "vue-property-decorator";
import BrushableHistogram from "@/components/BrushableHistogram.vue";
import PresetsView from "@/views/main/group/experiment/image/options/settings/PresetsView.vue";

@Component({
  components: { PresetsView, BrushableHistogram, GeneralSettingsView },
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
.channel-settings-view {
  height: calc(50vh - 82px);
}
</style>
