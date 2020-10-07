<template>
  <v-card tile>
    <v-card-title>
      Variables
      <v-spacer />
      <v-text-field v-model="search" label="Search" clearable single-line hide-details>
        <template v-slot:append>
          <v-icon>mdi-magnify</v-icon>
        </template>
      </v-text-field>
    </v-card-title>
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
    />
  </v-card>
</template>

<script lang="ts">
import { projectsModule } from "@/modules/projects";
import { IChannel } from "@/modules/projects/models";
import { settingsModule } from "@/modules/settings";
import { equals } from "rambda";
import { Component, Prop, Vue } from "vue-property-decorator";

@Component
export default class VariablesSelector extends Vue {
  readonly settingsContext = settingsModule.context(this.$store);
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

  get variables() {
    return this.step.variables;
  }

  get selected() {
    return this.items.filter((channel) => {
      if (this.variables.includes(channel.name)) {
        return channel;
      }
    }) as any;
  }

  set selected(items: IChannel[]) {
    const selectedVariables = items.map((item) => item.name);
    if (!equals(this.variables, selectedVariables)) {
      this.step.variables = selectedVariables;
    }
  }
}
</script>

<style scoped>
.scroll-view {
  height: 300px;
}
</style>
