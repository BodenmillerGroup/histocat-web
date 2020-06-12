<template>
  <v-card tile>
    <v-toolbar dense flat>
      <v-text-field
        v-model="search"
        label="Search"
        single-line
        hide-details
        clearable
        dense
      >
        <template v-slot:append-outer>
          <v-icon dense>mdi-magnify</v-icon>
        </template>
      </v-text-field>
    </v-toolbar>
    <v-data-table
      :headers="headers"
      :items="items"
      :search="search"
      v-model="selected"
      show-select
      hide-default-footer
      class="overflow-y-auto scroll-view"
      dense
      disable-pagination
      no-data-text="Please first select an acquisition"
    >
      <template v-slot:item.label="props">
        <v-edit-dialog :return-value.sync="props.item.label" @save="save">
          {{ props.item.label }}
          <template v-slot:input>
            <v-text-field
              v-model="props.item.label"
              :rules="[max25chars]"
              label="Edit"
              single-line
              counter
            ></v-text-field>
          </template>
        </v-edit-dialog>
      </template>
    </v-data-table>
  </v-card>
</template>

<script lang="ts">
import { experimentModule } from "@/modules/experiment";
import { IChannel } from "@/modules/experiment/models";
import { settingsModule } from "@/modules/settings";
import { equals } from "rambda";
import { Component, Vue } from "vue-property-decorator";
import { IChannelSettings } from "@/modules/settings/models";
import { BroadcastManager } from "@/utils/BroadcastManager";
import { SET_SHARED_CHANNEL_SETTINGS } from "@/modules/settings/events";
import { SET_SELECTED_METALS } from "@/modules/experiment/events";

@Component
export default class ChannelsView extends Vue {
  readonly settingsModule = settingsModule.context(this.$store);
  readonly experimentContext = experimentModule.context(this.$store);

  search = "";

  readonly headers = [
    {
      text: "Name",
      sortable: true,
      value: "name",
      align: "start",
      width: "30%",
    },
    {
      text: "Label",
      sortable: true,
      value: "label",
      align: "start",
      width: "50%",
    },
  ];

  max25chars = (v) => v.length <= 25 || "Input too long!";

  get activeAcquisitionId() {
    return this.experimentContext.getters.activeAcquisitionId;
  }

  get channels() {
    const acquisition = this.experimentContext.getters.activeAcquisition;
    return acquisition ? Object.values(acquisition.channels).sort((a, b) => a.mass - b.mass) : [];
  }

  get items() {
    if (!this.activeAcquisitionId) {
      return [];
    }
    return this.channels.map((channel) => {
      const settings = this.settingsModule.getters.getChannelSettings(this.activeAcquisitionId!, channel.name);
      return {
        id: channel.name,
        label: settings && settings.customLabel ? settings.customLabel : channel.label,
        name: channel.name,
        mass: channel.mass,
      };
    });
  }

  get selectedMetals() {
    return this.experimentContext.getters.selectedMetals;
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
    if (!equals(this.selectedMetals, selectedMetals)) {
      BroadcastManager.publish(SET_SELECTED_METALS, selectedMetals);
      this.experimentContext.actions.getChannelStackImage();
    }
  }

  save() {
    if (!this.activeAcquisitionId) {
      return;
    }
    let allSettings: IChannelSettings[] = [];
    this.items.forEach((item) => {
      const settings = this.settingsModule.getters.getChannelSettings(this.activeAcquisitionId!, item.name);
      if (!settings) {
        allSettings.push({
          acquisitionId: this.activeAcquisitionId!,
          name: item.name,
          customLabel: item.label,
        });
      } else {
        if (settings.customLabel !== item.label) {
          allSettings.push({
            ...settings,
            customLabel: item.label,
          });
        }
      }
    });
    if (allSettings.length > 0) {
      BroadcastManager.publish(SET_SHARED_CHANNEL_SETTINGS, allSettings);
      this.experimentContext.actions.getChannelStackImage();
    }
  }
}
</script>

<style scoped>
table.v-table tbody td,
table.v-table tbody th {
  height: 35px;
}

.scroll-view {
  height: calc(50vh - 94px);
}
</style>

<style>
.channels-table table {
  table-layout: fixed;
}
</style>
