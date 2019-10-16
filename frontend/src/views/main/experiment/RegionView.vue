<template>
  <v-card tile>
    <v-card-title>
      <v-text-field v-model="search" append-icon="mdi-magnify" label="Search" single-line hide-details clearable />
    </v-card-title>
    <v-data-table
      :headers="headers"
      :items="items"
      :search="search"
      hide-default-footer
      class="overflow-y-auto scroll-view"
      dense
      disable-pagination
      no-data-text="Please select a region"
    >
    </v-data-table>
  </v-card>
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
      align: "start"
    },
    {
      text: "Min",
      sortable: true,
      value: "min",
      align: "end"
    },
    {
      text: "Max",
      sortable: true,
      value: "max",
      align: "end"
    },
    {
      text: "Mean",
      sortable: true,
      value: "mean",
      align: "end"
    }
  ];

  get items() {
    return this.analysisContext.getters.selectedRegionStats;
  }
}
</script>

<style scoped>
table.v-table tbody td,
table.v-table tbody th {
  height: 35px;
}

.scroll-view {
  height: calc(50vh - 120px);
}
</style>

<style>
.channels-table table {
  table-layout: fixed;
}
</style>
