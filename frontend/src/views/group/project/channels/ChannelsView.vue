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
    no-data-text="Please first select an acquisition"
    class="root"
  >
    <template v-slot:top>
      <v-text-field v-model="search" label="Search" clearable single-line dense>
        <template v-slot:append>
          <v-icon>mdi-magnify</v-icon>
        </template>
      </v-text-field>
    </template>
    <template v-slot:item.customLabel="props">
      <v-edit-dialog :return-value.sync="props.item.customLabel" @save="save(props.item)">
        {{ props.item.customLabel }}
        <template v-slot:input>
          <v-text-field
            v-model="props.item.customLabel"
            :rules="[max25chars]"
            label="Edit"
            single-line
            counter
          ></v-text-field>
        </template>
      </v-edit-dialog>
    </template>
  </v-data-table>
</template>

<script lang="ts">
import { projectsModule } from "@/modules/projects";
import { IChannel } from "@/modules/projects/models";
import { isEqual } from "lodash-es";
import { Component, Vue } from "vue-property-decorator";

@Component
export default class ChannelsView extends Vue {
  readonly projectsContext = projectsModule.context(this.$store);

  search = "";

  readonly headers = [
    {
      text: "Name",
      sortable: true,
      value: "name",
      align: "start",
    },
    {
      text: "Label",
      sortable: true,
      value: "customLabel",
      align: "start",
    },
  ];

  max25chars = (v) => v.length <= 25 || "Input too long!";

  get activeAcquisitionId() {
    return this.projectsContext.getters.activeAcquisitionId;
  }

  get channels() {
    const acquisition = this.projectsContext.getters.activeAcquisition;
    return acquisition ? Object.values(acquisition.channels).sort((a, b) => a.mass - b.mass) : [];
  }

  get items() {
    if (!this.activeAcquisitionId) {
      return [];
    }
    return this.channels.map((channel) => {
      return {
        id: channel.name,
        customLabel: channel.customLabel,
        name: channel.name,
        mass: channel.mass,
      };
    });
  }

  get selectedMetals() {
    return this.projectsContext.getters.selectedMetals;
  }

  get selected() {
    return this.items.filter((channel) => {
      if (this.selectedMetals.includes(channel.name)) {
        return channel;
      }
    }) as any;
  }

  set selected(items: IChannel[]) {
    const selectedMetals = items.map((item) => item.name);
    if (!isEqual(this.selectedMetals, selectedMetals)) {
      this.projectsContext.mutations.setSelectedMetals(selectedMetals);
      this.projectsContext.actions.getChannelStackImage();
    }
  }

  save(item) {
    if (!this.activeAcquisitionId) {
      return;
    }
    this.projectsContext.actions.updateChannel({
      acquisitionId: this.activeAcquisitionId,
      data: {
        name: item.name,
        customLabel: item.customLabel,
      },
    });
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
