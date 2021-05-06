<template>
  <div class="presets-view">
    <v-toolbar flat dense color="grey lighten-4">
      <v-btn @click="createPreset" color="primary" elevation="1" small>Create preset</v-btn>
      <v-spacer></v-spacer>
      <v-tooltip bottom>
        <template v-slot:activator="{ on }">
          <v-btn icon small v-on="on" @click="refreshPresets">
            <v-icon small>mdi-refresh</v-icon>
          </v-btn>
        </template>
        <span>Refresh presets</span>
      </v-tooltip>
    </v-toolbar>
    <v-list dense class="pa-0">
      <v-list-item-group v-model="selected" color="primary">
        <v-list-item v-for="item in items" :key="item.id">
          <v-list-item-content>
            <v-list-item-title>{{ item.name }}</v-list-item-title>
            <v-list-item-subtitle>{{ item.createdAt }}</v-list-item-subtitle>
          </v-list-item-content>

          <v-list-item-action>
            <v-tooltip bottom>
              <template v-slot:activator="{ on }">
                <v-btn icon v-on="on" download color="primary lighten-3" @click="applyPreset($event, item.id)">
                  <v-icon>mdi-refresh-circle</v-icon>
                </v-btn>
              </template>
              <span>Apply preset</span>
            </v-tooltip>
          </v-list-item-action>
          <v-list-item-action>
            <v-tooltip bottom>
              <template v-slot:activator="{ on }">
                <v-btn icon v-on="on" color="secondary lighten-3" @click.stop="deletePreset($event, item.id)">
                  <v-icon>mdi-delete</v-icon>
                </v-btn>
              </template>
              <span>Delete preset</span>
            </v-tooltip>
          </v-list-item-action>
        </v-list-item>
      </v-list-item-group>
    </v-list>
  </div>
</template>

<script lang="ts">
import { projectsModule } from "@/modules/projects";
import { Component, Vue, Watch } from "vue-property-decorator";
import { presetsModule } from "@/modules/presets";

@Component
export default class PresetsView extends Vue {
  readonly projectsContext = projectsModule.context(this.$store);
  readonly presetContext = presetsModule.context(this.$store);

  selected?: number | null = null;

  @Watch("selected")
  presetChanged(index: number | null) {
    if (index !== null && index !== undefined) {
      const preset = this.presets[index];
      // BroadcastManager.publish(SET_ACTIVE_DATASET, dataset);
      // this.centroidsContext.actions.getCentroids({ datasetId: dataset.id });
    } else {
      // BroadcastManager.publish(SET_ACTIVE_DATASET, undefined);
    }
  }

  get presets() {
    return this.presetContext.getters.presets;
  }

  get items() {
    return this.presets.map((preset) => {
      return Object.assign({}, preset, {
        createdAt: new Date(preset.created_at).toUTCString(),
      });
    });
  }

  async applyPreset(event, id: number) {
    await this.presetContext.actions.applyPreset(id);
  }

  async deletePreset(event, id: number) {
    if (self.confirm("Do you really want to delete preset?")) {
      await this.presetContext.actions.deletePreset(id);
    }
  }

  async createPreset() {
    const name = self.prompt("Please enter preset name:");
    if (name) {
      await this.presetContext.actions.createPreset(name);
    }
  }

  async refreshPresets() {
    const projectId = this.projectsContext.getters.activeProjectId;
    if (projectId) {
      await this.presetContext.actions.getPresets(projectId);
    }
  }

  mounted() {
    this.refreshPresets();
  }
}
</script>

<style scoped>
.presets-view {
  width: 100%;
  height: 100%;
  overflow-y: auto;
}
</style>
