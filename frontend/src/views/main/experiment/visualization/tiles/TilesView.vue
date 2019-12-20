<template>
  <v-container fluid class="overflow-y-auto tiles-view">
    <v-row>
      <v-col v-for="item in items" :key="item.id" class="d-flex child-flex" :cols="cols">
        <v-card flat>
          <v-card-title class="subtitle-2 pa-0">{{ item.caption }}</v-card-title>
          <v-img :src="`${item.url}`" aspect-ratio="1" class="grey lighten-2" eager>
            <template v-slot:placeholder>
              <v-row class="fill-height ma-0" align="center" justify="center">
                <v-progress-circular indeterminate color="grey lighten-5"></v-progress-circular>
              </v-row>
            </template>
          </v-img>
        </v-card>
      </v-col>
    </v-row>
  </v-container>
</template>

<script lang="ts">
import { apiUrl } from "@/env";
import { experimentModule } from "@/modules/experiment";
import { settingsModule } from "@/modules/settings";
import { Component, Vue } from "vue-property-decorator";

@Component
export default class TilesView extends Vue {
  readonly experimentContext = experimentModule.context(this.$store);
  readonly settingsContext = settingsModule.context(this.$store);

  get items() {
    return this.experimentContext.getters.selectedChannels.map(channel => {
      let color = this.metalColorMap.get(channel.metal);
      if (color) {
        color = color.replace("#", "");
      }
      const channelSettings = this.settingsContext.getters.getChannelSettings(channel.id);
      const min = channelSettings && channelSettings.levels ? channelSettings.levels.min : "";
      const max = channelSettings && channelSettings.levels ? channelSettings.levels.max : "";
      return {
        id: channel.id,
        url: `${apiUrl}/channels/${channel.id}/image?color=${color}&min=${min}&max=${max}`,
        caption: channelSettings && channelSettings.customLabel ? channelSettings.customLabel : channel.label
      };
    });
  }

  get metalColorMap() {
    return this.settingsContext.getters.metalColorMap;
  }

  get cols() {
    const l = this.items.length;
    if (l === 1) {
      return 12;
    } else if (l < 5) {
      return 6;
    }
    return 4;
  }
}
</script>

<style scoped>
.tiles-view {
  height: calc(100vh - 154px);
}
</style>
