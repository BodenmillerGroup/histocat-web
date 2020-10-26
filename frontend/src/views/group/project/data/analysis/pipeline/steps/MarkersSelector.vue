<template>
  <v-card tile>
    <v-card-title>
      Markers
      <v-spacer />
      <v-text-field v-model="search" label="Search" clearable single-line hide-details>
        <template v-slot:append>
          <v-icon>mdi-magnify</v-icon>
        </template>
      </v-text-field>
    </v-card-title>
    <v-data-table
      :headers="headers"
      :items="channels"
      :search="search"
      v-model="selected"
      item-key="name"
      show-select
      hide-default-footer
      class="overflow-y-auto scroll-view"
      dense
      disable-pagination
      no-data-text="Please first select an acquisition"
    />
  </v-card>
</template>

<script lang="ts">
import { projectsModule } from "@/modules/projects";
import { IChannel } from "@/modules/projects/models";
import { settingsModule } from "@/modules/settings";
import { equals } from "rambda";
import { Component, Prop, Vue } from "vue-property-decorator";
import { datasetsModule } from "@/modules/datasets";

@Component
export default class MarkersSelector extends Vue {
  readonly settingsContext = settingsModule.context(this.$store);
  readonly datasetsContext = datasetsModule.context(this.$store);
  readonly projectsContext = projectsModule.context(this.$store);

  @Prop(Object) step;

  search = "";

  readonly headers = [
    {
      text: "Name",
      sortable: true,
      value: "name",
    },
    {
      text: "Label",
      sortable: true,
      value: "customLabel",
    },
    {
      text: "Mass",
      sortable: true,
      filterable: false,
      value: "mass",
      align: "end",
    },
  ];

  get activeAcquisitionId() {
    return this.projectsContext.getters.activeAcquisitionId;
  }

  get availableMarkers() {
    return this.datasetsContext.getters.channels;
  }

  get channels() {
    const acquisition = this.projectsContext.getters.activeAcquisition;
    const availableMarkers = this.availableMarkers;
    return acquisition
      ? Object.values(acquisition.channels)
          .filter((v) => availableMarkers.includes(v.name))
          .sort((a, b) => a.mass - b.mass)
      : [];
  }

  get selected() {
    return this.channels.filter((channel) => {
      if (this.step.markers.includes(channel.name)) {
        return channel;
      }
    }) as any;
  }

  set selected(items: IChannel[]) {
    const selectedMarkers = items.map((item) => item.name);
    if (!equals(this.step.markers, selectedMarkers)) {
      this.step.markers = selectedMarkers;
    }
  }
}
</script>

<style scoped>
.scroll-view {
  height: 300px;
}
</style>
