<template>
  <v-data-table
    :headers="headers"
    :items="items"
    :search="search"
    v-model="selected"
    item-key="name"
    show-select
    hide-default-footer
    dense
    disable-pagination
    no-data-text="Please first select an acquisition"
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
import { IChannel } from "@/modules/projects/models";
import { isEqual } from "lodash-es";
import { Component, Vue } from "vue-property-decorator";
import { datasetsModule } from "@/modules/datasets";

@Component
export default class ChannelSelector extends Vue {
  readonly datasetsContext = datasetsModule.context(this.$store);
  readonly projectsContext = projectsModule.context(this.$store);

  private selectedChannels: string[] = [];

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

  get channels() {
    return this.datasetsContext.getters.channels;
  }

  get items() {
    const acquisition = this.projectsContext.getters.activeAcquisition;
    const channels = this.channels;
    return acquisition
      ? Object.values(acquisition.channels)
          .filter((v) => channels.includes(v.name))
          .sort((a, b) => a.mass - b.mass)
      : [];
  }

  get selected() {
    return this.items.filter((channel) => {
      if (this.selectedChannels.includes(channel.name)) {
        return channel;
      }
    }) as any;
  }

  set selected(items: IChannel[]) {
    const selectedMarkers = items.map((item) => item.name);
    if (!isEqual(this.selectedChannels, selectedMarkers)) {
      this.selectedChannels = selectedMarkers;
    }
  }

  public getChannels() {
    return Object.freeze(this.selectedChannels);
  }
}
</script>
