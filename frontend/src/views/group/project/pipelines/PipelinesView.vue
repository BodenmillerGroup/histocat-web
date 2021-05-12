<template>
  <div class="pipelines-view">
    <v-toolbar flat dense color="grey lighten-4">
      <v-btn @click="savePipeline" color="primary" elevation="1" small>Save pipeline</v-btn>
      <v-spacer></v-spacer>
      <v-tooltip bottom>
        <template v-slot:activator="{ on }">
          <v-btn icon small v-on="on" @click="refreshPipelines">
            <v-icon small>mdi-refresh</v-icon>
          </v-btn>
        </template>
        <span>Refresh pipelines</span>
      </v-tooltip>
    </v-toolbar>
    <v-list dense two-line class="pa-0">
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
                <v-btn icon small v-on="on" download color="primary lighten-2" @click="loadPipeline(item.id)">
                  <v-icon small>mdi-refresh-circle</v-icon>
                </v-btn>
              </template>
              <span>Load pipeline</span>
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
              <span>Edit pipeline</span>
            </v-tooltip>
            <v-tooltip bottom>
              <template v-slot:activator="{ on }">
                <v-btn icon small v-on="on" color="secondary lighten-2" @click.stop="deletePipeline(item.id)">
                  <v-icon small>mdi-delete-outline</v-icon>
                </v-btn>
              </template>
              <span>Delete pipeline</span>
            </v-tooltip>
          </v-row>
        </v-list-item-action>
      </v-list-item>
    </v-list>
    <v-dialog v-model="dialog" scrollable max-width="600px">
      <v-card>
        <v-card-title>Edit Pipeline</v-card-title>
        <v-card-text>
          <v-text-field label="Name" v-model="name" />
          <v-text-field label="Description" v-model="description" />
        </v-card-text>
        <v-card-actions>
          <v-spacer />
          <v-btn text @click="dialog = false">Cancel</v-btn>
          <v-btn color="primary" text @click="updatePipeline()">Update</v-btn>
        </v-card-actions>
      </v-card>
    </v-dialog>
  </div>
</template>

<script lang="ts">
import { projectsModule } from "@/modules/projects";
import { Component, Vue } from "vue-property-decorator";
import { pipelinesModule } from "@/modules/pipelines";

@Component
export default class PipelinesView extends Vue {
  readonly projectsContext = projectsModule.context(this.$store);
  readonly pipelinesContext = pipelinesModule.context(this.$store);

  dialog = false;
  activeId: number | null = null;
  name: string | null = null;
  description: string | null = null;

  get activeProjectId() {
    return this.projectsContext.getters.activeProjectId;
  }

  get pipelines() {
    return this.pipelinesContext.getters.pipelines;
  }

  get items() {
    return this.pipelines.map((pipeline) => {
      return Object.assign({}, pipeline, {
        createdAt: new Date(pipeline.created_at).toUTCString(),
      });
    });
  }

  async loadPipeline(id: number) {
    await this.pipelinesContext.actions.loadPipeline(id);
  }

  async deletePipeline(id: number) {
    if (self.confirm("Do you really want to delete pipeline?")) {
      await this.pipelinesContext.actions.deletePipeline(id);
    }
  }

  async savePipeline() {
    const name = self.prompt("Please enter pipeline name:");
    if (name) {
      await this.pipelinesContext.actions.createPipeline(name);
    }
  }

  async refreshPipelines() {
    if (this.activeProjectId) {
      await this.pipelinesContext.actions.getPipelines(this.activeProjectId);
    }
  }

  async updatePipeline() {
    this.dialog = false;
    if (this.activeId) {
      await this.pipelinesContext.actions.updatePipeline({
        pipelineId: this.activeId,
        data: { name: this.name, description: this.description },
      });
    }
  }

  mounted() {
    this.refreshPipelines();
  }
}
</script>

<style scoped>
.pipelines-view {
  width: 100%;
  height: 100%;
  overflow-y: auto;
}
</style>
