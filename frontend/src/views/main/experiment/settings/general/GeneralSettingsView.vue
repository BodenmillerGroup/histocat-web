<template>
  <v-expansion-panels v-model="panel" multiple accordion>
    <LegendSettingsView></LegendSettingsView>
    <ScalebarSettingsView></ScalebarSettingsView>
    <FilterSettingsView></FilterSettingsView>
  </v-expansion-panels>
</template>

<script lang="ts">
import FilterSettingsView from "@/views/main/experiment/settings/general/FilterSettingsView.vue";
import LegendSettingsView from "@/views/main/experiment/settings/general/LegendSettingsView.vue";
import ScalebarSettingsView from "@/views/main/experiment/settings/general/ScalebarSettingsView.vue";
import { Component, Vue } from "vue-property-decorator";

@Component({
  components: { FilterSettingsView, ScalebarSettingsView, LegendSettingsView },
})
export default class GeneralSettingsView extends Vue {
  panel = [];
}
</script>
