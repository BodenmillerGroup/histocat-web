<template>
  <v-tabs v-model="tab">
    <v-tab>Blend</v-tab>
    <v-tab>Regions</v-tab>
    <v-tab>Tiles</v-tab>
    <v-tab-item>
      <BlendTab />
    </v-tab-item>
    <v-tab-item>
      <RegionsTab />
    </v-tab-item>
    <v-tab-item>
      <TilesView />
    </v-tab-item>
  </v-tabs>
</template>

<script lang="ts">
import RegionsTab from "@/views/group/project/image/visualization/regions/RegionsTab.vue";
import TilesView from "@/views/group/project/image/visualization/tiles/TilesView.vue";
import { Component, Vue } from "vue-property-decorator";
import BlendTab from "@/views/group/project/image/visualization/blend/BlendTab.vue";

@Component({
  components: {
    BlendTab,
    RegionsTab,
    TilesView,
  },
})
export default class VisualizationView extends Vue {
  tab = 0;
}
</script>
