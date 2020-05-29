<template>
  <div :class="layoutClass">
    <v-tabs v-model="mainTab">
      <v-tab>Blend</v-tab>
      <v-tab>Blend New</v-tab>
      <v-tab>Tiles</v-tab>
      <v-tab-item>
        <BlendTab />
      </v-tab-item>
      <v-tab-item>
        <BlendTabNew />
      </v-tab-item>
      <v-tab-item>
        <TilesView />
      </v-tab-item>
    </v-tabs>
    <div v-show="showOptions">
      <v-row no-gutters>
        <v-col>
          <v-tabs v-model="secondaryTab">
            <v-tab>Channels</v-tab>
            <v-tab>Region</v-tab>
            <v-tab-item>
              <ChannelsView />
            </v-tab-item>
            <v-tab-item>
              <RegionView />
            </v-tab-item>
          </v-tabs>
        </v-col>
      </v-row>
      <v-row dense>
        <v-col>
          <SettingsView />
        </v-col>
      </v-row>
    </div>
  </div>
</template>

<script lang="ts">
import { mainModule } from "@/modules/main";
import ChannelsView from "@/views/main/experiment/image/ChannelsView.vue";
import RegionView from "@/views/main/experiment/image/RegionView.vue";
import SettingsView from "@/views/main/experiment/image/settings/SettingsView.vue";
import BlendTab from "@/views/main/experiment/image/visualization/blend/BlendTab.vue";
import TilesView from "@/views/main/experiment/image//visualization/tiles/TilesView.vue";
import { Component, Vue } from "vue-property-decorator";
import BlendTabNew from "@/views/main/experiment/image/visualization/blend/BlendTabNew.vue";

@Component({
  components: {
    BlendTabNew,
    RegionView,
    SettingsView,
    ChannelsView,
    BlendTab,
    TilesView,
  },
})
export default class VisualizationView extends Vue {
  readonly mainContext = mainModule.context(this.$store);

  mainTab = 0;
  secondaryTab = 0;

  get showOptions() {
    return this.mainContext.getters.showOptions;
  }

  get layoutClass() {
    return this.showOptions ? "layout-options" : null;
  }
}
</script>

<style scoped>
.layout-options {
  display: grid;
  grid-template-columns: 1fr 380px;
  grid-template-rows: auto;
}
</style>
