<template>
  <v-container fluid class="overflow-y-auto tiles-view">
    <v-row>
      <v-col v-for="item in items" :key="item.name" class="d-flex child-flex" :cols="cols">
        <v-card flat>
          <v-card-title class="text-subtitle-2 pa-0">{{ item.caption }}</v-card-title>
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
import { projectsModule } from "@/modules/projects";
import { settingsModule } from "@/modules/settings";
import { Component, Vue } from "vue-property-decorator";

@Component
export default class TilesView extends Vue {
  readonly projectsContext = projectsModule.context(this.$store);
  readonly settingsContext = settingsModule.context(this.$store);

  get channelsSettings() {
    return this.settingsContext.getters.channelsSettings;
  }

  get items() {
    const acquisitionId = this.projectsContext.getters.activeAcquisitionId;
    if (!acquisitionId) {
      return [];
    }
    return this.projectsContext.getters.selectedChannels.map((channel) => {
      let url = `${apiUrl}/acquisitions/${acquisitionId}/${channel.name}/image?`;
      const channelSettings = this.channelsSettings[channel.name];
      if (channelSettings) {
        const color = channelSettings.color.replace("#", "");
        url += `color=${color}`;
      }
      if (channelSettings && channelSettings.levels) {
        url += `&min=${channelSettings.levels.min}&max=${channelSettings.levels.max}`;
      }
      return {
        name: channel.name,
        url: url,
        caption: channel.customLabel,
      };
    });
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
  height: calc(100vh - 84px);
}
</style>
