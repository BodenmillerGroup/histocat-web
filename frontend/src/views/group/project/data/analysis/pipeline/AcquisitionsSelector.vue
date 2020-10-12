<template>
  <v-data-table
    :headers="headers"
    :items="items"
    :search="search"
    v-model="selected"
    show-select
    hide-default-footer
    dense
    disable-pagination
    no-data-text="No available single-cell data"
  >
    <template v-slot:top>
      <v-text-field v-model="search" label="Search" clearable single-line>
        <template v-slot:append>
          <v-icon>mdi-magnify</v-icon>
        </template>
      </v-text-field>
    </template>
  </v-data-table>
</template>

<script lang="ts">
import { projectsModule } from "@/modules/projects";
import { IAcquisition } from "@/modules/projects/models";
import { equals } from "rambda";
import { Component, Vue } from "vue-property-decorator";
import { datasetsModule } from "@/modules/datasets";
import {pipelinesModule} from "@/modules/pipelines";

@Component
export default class AcquisitionsSelector extends Vue {
  readonly projectsContext = projectsModule.context(this.$store);
  readonly datasetsContext = datasetsModule.context(this.$store);
  readonly pipelinesContext = pipelinesModule.context(this.$store);

  search = "";

  readonly headers = [
    {
      text: "ID",
      sortable: true,
      filterable: false,
      value: "id",
      align: "end",
      width: 70,
    },
    {
      text: "Slide ID",
      sortable: true,
      filterable: false,
      value: "slide_id",
      align: "end",
      width: 100,
    },
    {
      text: "Descriptions",
      sortable: true,
      value: "description",
    },
  ];

  get selectedAcquisitionIds() {
    return this.pipelinesContext.getters.selectedAcquisitionIds;
  }

  get projectData() {
    return this.projectsContext.getters.projectData!;
  }

  get activeDataset() {
    return this.datasetsContext.getters.activeDataset;
  }

  get items() {
    const acquisitions: IAcquisition[] = [];
    if (this.projectData.slides) {
      for (const s of this.projectData.slides) {
        for (const a of s.acquisitions) {
          let hasMask = false;
          if (this.activeDataset && this.activeDataset.meta.probability_masks) {
            hasMask = !!this.activeDataset.meta.probability_masks[a.id];
          }
          if (hasMask) {
            acquisitions.push(a);
          }
        }
      }
    }
    return acquisitions;
  }

  get selected() {
    return this.items.filter((item) => {
      if (this.selectedAcquisitionIds.includes(item.id)) {
        return item;
      }
    }) as any;
  }

  set selected(items: IAcquisition[]) {
    const selectedIds = items.map((item) => item.id);
    if (!equals(this.selectedAcquisitionIds, selectedIds)) {
      this.pipelinesContext.mutations.setSelectedAcquisitionIds(selectedIds);
    }
  }
}
</script>
