<template>
  <v-card tile class="ma-1">
    <v-card-title>Acquisitions</v-card-title>
    <v-card-text>
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
    </v-card-text>
  </v-card>
</template>

<script lang="ts">
import { projectsModule } from "@/modules/projects";
import { IAcquisition } from "@/modules/projects/models";
import { isEqual } from "lodash-es";
import { Component, Vue } from "vue-property-decorator";
import { segmentationModule } from "@/modules/segmentation";

@Component
export default class AcquisitionsView extends Vue {
  readonly projectsContext = projectsModule.context(this.$store);
  readonly segmentationContext = segmentationModule.context(this.$store);

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
    return this.segmentationContext.getters.selectedAcquisitionIds;
  }

  get projectData() {
    return this.projectsContext.getters.projectData!;
  }

  get items() {
    const acquisitions: IAcquisition[] = [];
    if (this.projectData.slides) {
      for (const s of this.projectData.slides) {
        for (const a of s.acquisitions) {
          acquisitions.push(a);
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
    if (!isEqual(this.selectedAcquisitionIds, selectedIds)) {
      this.segmentationContext.mutations.setSelectedAcquisitionIds(selectedIds);
    }
  }
}
</script>
