<template>
  <v-row no-gutters>
    <v-col :cols="columns">
      <v-tabs v-model="mainTab">
        <v-tab>Blend</v-tab>
        <v-tab>Tiles</v-tab>
        <v-tab-item>
          <BlendTab />
        </v-tab-item>
        <v-tab-item>
          <TilesView />
        </v-tab-item>
      </v-tabs>
    </v-col>
    <v-col v-show="showOptions" cols="3">
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
    </v-col>
  </v-row>
</template>

<script lang="ts">
import { mainModule } from "@/modules/main";
import ChannelsView from "@/views/main/experiment/ChannelsView.vue";
import RegionView from "@/views/main/experiment/RegionView.vue";
import SettingsView from "@/views/main/experiment/settings/SettingsView.vue";
import BlendTab from "@/views/main/experiment/visualization/blend/BlendTab.vue";
import TilesView from "@/views/main/experiment/visualization/tiles/TilesView.vue";
import { Component, Vue } from "vue-property-decorator";

@Component({
  components: {
    RegionView,
    SettingsView,
    ChannelsView,
    BlendTab,
    TilesView
  }
})
export default class VisualizationView extends Vue {
  readonly mainContext = mainModule.context(this.$store);

  mainTab = 0;
  secondaryTab = 0;

  get showOptions() {
    return this.mainContext.getters.showOptions;
  }

  get columns() {
    if (this.showOptions) {
      return 9;
    }
    return 12;
  }
}
</script>
