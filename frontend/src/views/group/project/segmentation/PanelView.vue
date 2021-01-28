<template>
  <v-card tile class="ma-1">
    <v-card-title>Channels</v-card-title>
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
        no-data-text="Please first select an acquisition"
      >
        <template v-slot:top>
          <v-text-field v-model="search" label="Search" clearable single-line>
            <template v-slot:append>
              <v-icon>mdi-magnify</v-icon>
            </template>
          </v-text-field>
        </template>
        <template v-slot:item.type="props">
          <v-radio-group v-model="props.item.type" row dense mandatory hide-details class="ma-0">
            <v-radio label="" value="none" @click="typeChanged" />
            <v-radio label="nuclear" value="nuclear" @click="typeChanged" />
            <v-radio label="cytoplasm" value="cytoplasm" @click="typeChanged" />
          </v-radio-group>
        </template>
      </v-data-table>
    </v-card-text>
  </v-card>
</template>

<script lang="ts">
import { projectsModule } from "@/modules/projects";
import { IChannel } from "@/modules/projects/models";
import { settingsModule } from "@/modules/settings";
import { intersectionBy } from "lodash-es";
import { Component, Vue } from "vue-property-decorator";
import { segmentationModule } from "@/modules/segmentation";
import { channelTypeEnum } from "@/utils/enums";
import { channelTypeToString } from "@/utils/converters";

@Component
export default class PanelView extends Vue {
  readonly settingsContext = settingsModule.context(this.$store);
  readonly projectsContext = projectsModule.context(this.$store);
  readonly segmentationContext = segmentationModule.context(this.$store);

  search = "";

  readonly channelTypes = channelTypeEnum;
  readonly channelTypeToString = channelTypeToString;

  readonly headers = [
    {
      text: "Name",
      sortable: true,
      value: "name",
      align: "start",
      width: 100,
    },
    {
      text: "Label",
      sortable: true,
      value: "label",
      align: "start",
    },
    {
      text: "Mass",
      sortable: true,
      value: "mass",
      align: "end",
      width: 100,
    },
    {
      text: "Type",
      sortable: true,
      value: "type",
      align: "start",
      width: 350,
    },
  ];

  get selectedAcquisitionIds() {
    return this.segmentationContext.getters.selectedAcquisitionIds;
  }

  get channels() {
    const acquisitionChannels: any[] = [];
    for (const slide of this.projectsContext.getters.projectData?.slides!) {
      for (const acquisition of slide.acquisitions) {
        if (this.selectedAcquisitionIds.includes(acquisition.id)) {
          acquisitionChannels.push(Array.from(Object.values(acquisition.channels)));
        }
      }
    }
    const commonChannels = intersectionBy<IChannel>(...acquisitionChannels, "name");
    return commonChannels.sort((a, b) => a.mass - b.mass);
  }

  get items() {
    return this.channels.map((channel) => {
      return {
        id: channel.name,
        label: channel.customLabel,
        name: channel.name,
        mass: channel.mass,
        type: "none",
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
    // const selectedMetals = items.map((item) => item.name);
    // if (!isEqual(this.selectedMetals, selectedMetals)) {
    //   this.projectsContext.actions.setSelectedMetals(selectedMetals);
    //   this.projectsContext.actions.getChannelStackImage();
    // }
  }

  typeChanged() {
    const nucleiChannels = this.items.filter((item) => item.type === "nuclear").map((item) => item.name);
    const cytoplasmChannels = this.items.filter((item) => item.type === "cytoplasm").map((item) => item.name);
    this.segmentationContext.mutations.setNucleiChannels(nucleiChannels);
    this.segmentationContext.mutations.setCytoplasmChannels(cytoplasmChannels);
  }
}
</script>
