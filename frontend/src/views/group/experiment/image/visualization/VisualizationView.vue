<template>
  <v-tabs v-model="tab">
    <v-tab>Blend New</v-tab>
    <v-tab>Blend</v-tab>
    <v-tab>Tiles</v-tab>
    <v-tab-item>
      <BlendTabNew />
    </v-tab-item>
    <v-tab-item>
      <BlendTab />
    </v-tab-item>
    <v-tab-item>
      <TilesView />
    </v-tab-item>
  </v-tabs>
</template>

<script lang="ts">
import BlendTab from "@/views/group/experiment/image/visualization/blend/BlendTab.vue";
import TilesView from "@/views/group/experiment/image//visualization/tiles/TilesView.vue";
import { Component, Vue } from "vue-property-decorator";
import BlendTabNew from "@/views/group/experiment/image/visualization/blend/BlendTabNew.vue";

@Component({
  components: {
    BlendTabNew,
    BlendTab,
    TilesView,
  },
})
export default class VisualizationView extends Vue {
  tab = 0;
}
</script>
