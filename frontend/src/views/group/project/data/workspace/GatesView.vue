<template>
  <v-card tile>
    <v-toolbar flat dense color="grey lighten-4">
      <v-btn @click="createGate" color="primary" elevation="1" small>Save gate</v-btn>
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
            <v-list-item-subtitle v-if="item.description">{{ item.description }}</v-list-item-subtitle>
            <v-list-item-subtitle class="font-weight-light">{{ item.createdAt }}</v-list-item-subtitle>
          </v-list-item-content>

          <v-list-item-action>
            <v-row>
              <v-tooltip bottom>
                <template v-slot:activator="{ on }">
                  <v-btn icon small v-on="on" download color="primary lighten-2" @click="applyGate(item.id)">
                    <v-icon small>mdi-refresh-circle</v-icon>
                  </v-btn>
                </template>
                <span>Apply gate</span>
              </v-tooltip>
              <v-dialog v-model="dialog" scrollable max-width="600px">
                <template v-slot:activator="{ on, attrs }">
                  <v-btn
                    icon
                    small
                    v-bind="attrs"
                    v-on="on"
                    color="primary lighten-2"
                    @click="
                      name = item.name;
                      description = item.description;
                    "
                  >
                    <v-icon small>mdi-pencil-outline</v-icon>
                  </v-btn>
                </template>
                <v-card>
                  <v-card-title>Edit Gate</v-card-title>
                  <v-card-text>
                    <v-text-field label="Name" v-model="name" />
                    <v-text-field label="Description" v-model="description" />
                  </v-card-text>
                  <v-card-actions>
                    <v-spacer />
                    <v-btn text @click="dialog = false">Cancel</v-btn>
                    <v-btn color="primary" text @click="updateGate(item.id)">Update</v-btn>
                  </v-card-actions>
                </v-card>
              </v-dialog>
              <v-tooltip bottom>
                <template v-slot:activator="{ on }">
                  <v-btn icon small v-on="on" color="secondary lighten-2" @click.stop="deleteGate(item.id)">
                    <v-icon small>mdi-delete-outline</v-icon>
                  </v-btn>
                </template>
                <span>Delete gate</span>
              </v-tooltip>
            </v-row>
          </v-list-item-action>
        </v-list-item>
      </v-list-item-group>
    </v-list>
  </v-card>
</template>

<script lang="ts">
import { projectsModule } from "@/modules/projects";
import { Component, Vue, Watch } from "vue-property-decorator";
import { gatesModule } from "@/modules/gates";
import { datasetsModule } from "@/modules/datasets";

@Component
export default class GatesView extends Vue {
  readonly projectsContext = projectsModule.context(this.$store);
  readonly gateContext = gatesModule.context(this.$store);
  readonly datasetContext = datasetsModule.context(this.$store);

  dialog = false;
  name: string | null = null;
  description: string | null = null;

  selected?: number | null = null;

  @Watch("selected")
  gateChanged(index: number | null) {
    if (index !== null && index !== undefined) {
      const gate = this.gates[index];
      this.gateContext.actions.setActiveGateId(gate.id);
    } else {
      this.gateContext.actions.setActiveGateId(null);
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

  async applyGate(id: number) {
    await this.gateContext.actions.applyGate(id);
    await this.projectsContext.actions.getChannelStackImage();
  }

  async deleteGate(id: number) {
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

  async updateGate(gateId: number) {
    this.dialog = false;
    await this.gateContext.actions.updateGate({
      gateId: gateId,
      data: { name: this.name, description: this.description },
    });
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
