<template>
  <v-card tile>
    <v-toolbar flat dense color="grey lighten-4">
      <v-btn @click="createGate" color="primary" elevation="1" small>Save gate</v-btn>
    </v-toolbar>
    <v-list dense class="overflow-y-auto scroll-view pa-0">
      <v-list-item-group v-model="selected" color="primary">
        <v-list-item v-for="item in items" :key="item.id">
          <v-list-item-content>
            <v-list-item-title>{{ item.name }}</v-list-item-title>
            <v-list-item-subtitle v-if="item.description">{{ item.description }}</v-list-item-subtitle>
          </v-list-item-content>

          <v-list-item-action>
            <v-row>
              <v-tooltip bottom>
                <template v-slot:activator="{ on }">
                  <v-btn icon small v-on="on" download color="primary lighten-2" @click="applyGate(item.id)">
                    <v-icon small>mdi-refresh-circle</v-icon>
                  </v-btn>
                </template>
                <span>Load gate</span>
              </v-tooltip>
              <v-tooltip bottom>
                <template v-slot:activator="{ on }">
                  <v-btn
                    icon
                    small
                    v-on="on"
                    color="primary lighten-2"
                    @click.stop="
                      activeId = item.id;
                      name = item.name;
                      description = item.description;
                      dialog = true;
                    "
                  >
                    <v-icon small>mdi-pencil-outline</v-icon>
                  </v-btn>
                </template>
                <span>Edit gate</span>
              </v-tooltip>
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
    <v-dialog v-model="dialog" scrollable max-width="600px">
      <v-card>
        <v-card-title>Edit Gate</v-card-title>
        <v-card-text>
          <v-text-field label="Name" v-model="name" />
          <v-text-field label="Description" v-model="description" />
        </v-card-text>
        <v-card-actions>
          <v-spacer />
          <v-btn text @click="dialog = false">Cancel</v-btn>
          <v-btn color="primary" text @click="updateGate()">Update</v-btn>
        </v-card-actions>
      </v-card>
    </v-dialog>
  </v-card>
</template>

<script lang="ts">
import { Component, Vue, Watch } from "vue-property-decorator";
import { gatesModule } from "@/modules/gates";

@Component
export default class GatesView extends Vue {
  readonly gateContext = gatesModule.context(this.$store);

  dialog = false;
  activeId: number | null = null;
  name: string | null = null;
  description: string | null = null;

  selected?: number | null = null;

  @Watch("selected")
  gateChanged(index: number | null) {
    if (index !== null && index !== undefined) {
      const gate = this.gates[index];
      this.gateContext.mutations.setActiveGateId(gate.id);
    } else {
      this.gateContext.mutations.setActiveGateId(null);
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
    await this.gateContext.actions.loadGate(id);
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

  async updateGate() {
    this.dialog = false;
    if (this.activeId) {
      await this.gateContext.actions.updateGate({
        gateId: this.activeId,
        data: { name: this.name, description: this.description },
      });
    }
  }
}
</script>

<style scoped>
.scroll-view {
  height: calc(33vh - 100px);
}
</style>
