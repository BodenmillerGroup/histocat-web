<template>
  <v-card tile>
    <v-toolbar flat dense color="grey lighten-4">
      <v-btn @click="createGate" color="primary" elevation="1" small>
        Create gate
      </v-btn>
      <v-spacer></v-spacer>
      <v-tooltip bottom>
        <template v-slot:activator="{ on }">
          <v-btn icon small v-on="on" @click="refreshGates">
            <v-icon small>mdi-refresh</v-icon>
          </v-btn>
        </template>
        <span>Refresh gates</span>
      </v-tooltip>
    </v-toolbar>
    <v-list dense two-line class="overflow-y-auto scroll-view pa-0">
      <v-list-item-group v-model="selected" color="primary">
        <v-list-item v-for="item in items" :key="item.id">
          <v-list-item-content>
            <v-list-item-title>{{ item.name }}</v-list-item-title>
            <v-list-item-subtitle>{{ item.createdAt }}</v-list-item-subtitle>
          </v-list-item-content>

          <v-list-item-action>
            <v-tooltip bottom>
              <template v-slot:activator="{ on }">
                <v-btn icon v-on="on" download color="primary lighten-3" @click="applyGate($event, item.id)">
                  <v-icon>mdi-refresh-circle</v-icon>
                </v-btn>
              </template>
              <span>Apply gate</span>
            </v-tooltip>
          </v-list-item-action>
          <v-list-item-action>
            <v-tooltip bottom>
              <template v-slot:activator="{ on }">
                <v-btn icon v-on="on" color="secondary lighten-3" @click.stop="deleteGate($event, item.id)">
                  <v-icon>mdi-delete</v-icon>
                </v-btn>
              </template>
              <span>Delete gate</span>
            </v-tooltip>
          </v-list-item-action>
        </v-list-item>
      </v-list-item-group>
    </v-list>
  </v-card>
</template>

<script lang="ts">
import { experimentModule } from "@/modules/experiment";
import { Component, Vue, Watch } from "vue-property-decorator";
import {gateModule} from "@/modules/gates";
import {datasetModule} from "@/modules/datasets";

@Component
export default class GatesView extends Vue {
  readonly experimentContext = experimentModule.context(this.$store);
  readonly gateContext = gateModule.context(this.$store);
  readonly datasetContext = datasetModule.context(this.$store);

  selected?: number | null = null;

  @Watch("selected")
  gateChanged(index: number | null) {
    if (index !== null && index !== undefined) {
      const preset = this.gates[index];
      // BroadcastManager.publish(SET_ACTIVE_DATASET, dataset);
      // this.centroidsContext.actions.getCentroids({ datasetId: dataset.id });
    } else {
      // BroadcastManager.publish(SET_ACTIVE_DATASET, undefined);
    }
  }

  get gates() {
    return this.gateContext.getters.gates;
  }

  get items() {
    return this.gates.map((gate) => {
      return Object.assign({}, gate, {
        createdAt: new Date(gate.created_at).toUTCString(),
      });
    });
  }

  async applyGate(event, id: number) {
    await this.gateContext.actions.applyGate(id);
    await this.experimentContext.actions.getChannelStackImage();
  }

  async deleteGate(event, id: number) {
    if (self.confirm("Do you really want to delete gate?")) {
      await this.gateContext.actions.deleteGate(id);
    }
  }

  async createGate() {
    const name = self.prompt("Please enter gate name:");
    if (name) {
      await this.gateContext.actions.createGate(name);
    }
  }

  async refreshGates() {
    const dataset = this.datasetContext.getters.activeDataset;
    if (dataset) {
      await this.gateContext.actions.getGates(dataset.id);
    }
  }

  mounted() {
    this.refreshGates();
  }
}
</script>

<style scoped>
.scroll-view {
  height: calc(100vh - 132px);
}
</style>
