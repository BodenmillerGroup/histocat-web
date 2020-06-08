<template>
  <div>
    <v-tabs v-model="tab">
      <v-tab>Channels</v-tab>
      <v-tab>Region</v-tab>
      <v-tab-item>
        <ChannelsView />
      </v-tab-item>
      <v-tab-item>
        <RegionView />
      </v-tab-item>
    </v-tabs>
    <SettingsView class="mt-1" />
  </div>
</template>

<script lang="ts">
import ChannelsView from "@/views/main/experiment/image/options/ChannelsView.vue";
import RegionView from "@/views/main/experiment/image/options/RegionView.vue";
import SettingsView from "@/views/main/experiment/image/options/settings/SettingsView.vue";
import { Component, Vue } from "vue-property-decorator";

@Component({
  components: {
    RegionView,
    SettingsView,
    ChannelsView,
  },
})
export default class OptionsView extends Vue {
  tab = 0;
}
</script>
