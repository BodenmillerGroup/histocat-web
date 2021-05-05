<template>
  <v-data-table
    :headers="headers"
    :items="items"
    :search="search"
    hide-default-footer
    dense
    disable-pagination
    no-data-text="Please select a region"
    class="root"
  >
    <template v-slot:top>
      <v-text-field v-model="search" label="Search" clearable single-line dense>
        <template v-slot:append>
          <v-icon>mdi-magnify</v-icon>
        </template>
      </v-text-field>
    </template>
  </v-data-table>
</template>

<script lang="ts">
import { analysisModule } from "@/modules/analysis";
import { settingsModule } from "@/modules/settings";
import { Component, Vue } from "vue-property-decorator";

@Component
export default class RegionsView extends Vue {
  readonly settingsModule = settingsModule.context(this.$store);
  readonly analysisContext = analysisModule.context(this.$store);

  search = "";

  readonly headers = [
    {
      text: "Metal",
      sortable: true,
      value: "metal",
      align: "start",
    },
    {
      text: "Min",
      sortable: true,
      value: "min",
      align: "end",
    },
    {
      text: "Max",
      sortable: true,
      value: "max",
      align: "end",
    },
    {
      text: "Mean",
      sortable: true,
      value: "mean",
      align: "end",
    },
  ];

  get items() {
    return this.analysisContext.getters.selectedRegionStats;
  }
}
</script>

<style scoped>
.root {
  width: 100%;
  height: 100%;
  overflow-y: auto;
  padding: 8px;
  margin-right: auto;
  margin-left: auto;
}
</style>
