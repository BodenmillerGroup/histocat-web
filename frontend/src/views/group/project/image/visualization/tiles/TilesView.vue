<template>
  <div v-intersect="onIntersect">
    <v-toolbar dense flat>
      <v-menu offset-y>
        <template v-slot:activator="{ on }">
          <v-btn v-on="on" small elevation="1">
            <v-icon left small>mdi-download</v-icon>
            Export
          </v-btn>
        </template>
        <v-list dense>
          <v-list-item @click="exportImage">
            <v-list-item-title>Export PNG</v-list-item-title>
          </v-list-item>
        </v-list>
      </v-menu>
    </v-toolbar>
    <v-container ref="tilesContainer" fluid>
      <v-row>
        <v-col v-for="item in items" :key="item.name" :cols="cols">
          <span class="text-subtitle-2">{{ item.caption }}</span>
          <v-img :src="`${item.url}`" aspect-ratio="1" eager>
            <template v-slot:placeholder>
              <v-row class="fill-height ma-0" align="center" justify="center">
                <v-progress-circular indeterminate color="grey lighten-5" />
              </v-row>
            </template>
          </v-img>
        </v-col>
      </v-row>
    </v-container>
  </div>
</template>

<script lang="ts">
import { apiUrl } from "@/env";
import { projectsModule } from "@/modules/projects";
import { settingsModule } from "@/modules/settings";
import { Component, Vue } from "vue-property-decorator";
import html2canvas from "html2canvas";
import { saveAs } from "file-saver";

@Component
export default class TilesView extends Vue {
  readonly projectsContext = projectsModule.context(this.$store);
  readonly settingsContext = settingsModule.context(this.$store);

  isVisible = false;

  get channelsSettings() {
    return this.settingsContext.getters.channelsSettings;
  }

  get items() {
    const acquisitionId = this.projectsContext.getters.activeAcquisitionId;
    if (!acquisitionId || !this.isVisible) {
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

  onIntersect(entries, observer, isIntersecting) {
    this.isVisible = isIntersecting;
  }

  async exportImage(e) {
    const canvas = await html2canvas(this.$refs.tilesContainer as HTMLElement, {
      useCORS: true,
      backgroundColor: null,
    });
    canvas.toBlob(
      (blob) => {
        if (blob) {
          saveAs(blob);
        }
      },
      "image/png",
      1.0
    );
  }
}
</script>
