<template>
  <div class="settings-view">
    <v-expansion-panels v-model="panel" multiple accordion>
      <LegendSettingsView></LegendSettingsView>
      <ScalebarSettingsView></ScalebarSettingsView>
      <FilterSettingsView></FilterSettingsView>
    </v-expansion-panels>
  </div>
</template>

<script lang="ts">
import FilterSettingsView from "@/views/group/project/settings/FilterSettingsView.vue";
import LegendSettingsView from "@/views/group/project/settings/LegendSettingsView.vue";
import ScalebarSettingsView from "@/views/group/project/settings/ScalebarSettingsView.vue";
import { Component, Vue } from "vue-property-decorator";

@Component({
  components: { FilterSettingsView, ScalebarSettingsView, LegendSettingsView },
})
export default class GeneralSettingsView extends Vue {
  panel = [];
}
</script>

<style scoped>
.settings-view {
  width: 100%;
  height: 100%;
  overflow-y: auto;
}
</style>
