<template>
  <v-card flat>
    <v-container grid-list-sm fluid>
      <v-layout row wrap>
        <v-flex
          v-for="url in urls"
          :key="url"
          xs4
          flex-grow-1
        >
          <v-card flat tile flex-grow-1>
            <v-img
              :src="`${url}`"
              aspect-ratio="1"
              class="grey lighten-2"
            >
              <template v-slot:placeholder>
                <v-layout
                  fill-height
                  align-center
                  justify-center
                  ma-0
                >
                  <v-progress-circular indeterminate color="grey lighten-5"></v-progress-circular>
                </v-layout>
              </template>
            </v-img>
          </v-card>
        </v-flex>
      </v-layout>
    </v-container>
  </v-card>
</template>

<script lang="ts">
  import { apiUrl } from '@/env';
  import { experimentModule } from '@/modules/experiment';
  import { settingsModule } from '@/modules/settings';
  import 'ol/ol.css';
  import { Component, Vue } from 'vue-property-decorator';

  @Component
  export default class TilesView extends Vue {
    readonly experimentContext = experimentModule.context(this.$store);
    readonly settingsContext = settingsModule.context(this.$store);

    get urls() {
      return this.experimentContext.getters.selectedChannels.map((channel) => {
        let color = this.metalColorMap.get(channel.metal);
        if (color) {
          color = color.replace('#', '');
        }
        const channelSettings = this.settingsContext.getters.channelSettings(channel.id);
        const min = channelSettings && channelSettings.levels ? channelSettings.levels.min : '';
        const max = channelSettings && channelSettings.levels ? channelSettings.levels.max : '';
        return `${apiUrl}/api/v1/channels/${channel.id}/image?color=${color}&min=${min}&max=${max}`;
      });
    }

    get metalColorMap() {
      return this.settingsContext.getters.metalColorMap;
    }
  }
</script>
