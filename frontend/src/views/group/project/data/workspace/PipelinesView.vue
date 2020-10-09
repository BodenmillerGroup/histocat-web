<template>
  <v-card tile>
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
                <v-btn icon v-on="on" download color="primary lighten-3" @click="loadPipeline(item.id)">
                  <v-icon>mdi-refresh-circle</v-icon>
                </v-btn>
              </template>
              <span>Load pipeline</span>
            </v-tooltip>
          </v-list-item-action>
          <v-list-item-action>
            <v-tooltip bottom>
              <template v-slot:activator="{ on }">
                <v-btn icon v-on="on" color="secondary lighten-3" @click.stop="deletePipeline(item.id)">
                  <v-icon>mdi-delete</v-icon>
                </v-btn>
              </template>
              <span>Delete pipeline</span>
            </v-tooltip>
          </v-list-item-action>
        </v-list-item>
      </v-list-item-group>
    </v-list>
  </v-card>
</template>

<script lang="ts">
import { projectsModule } from "@/modules/projects";
import { Component, Vue, Watch } from "vue-property-decorator";
import { datasetsModule } from "@/modules/datasets";
import { pipelinesModule } from "@/modules/pipelines";

@Component
export default class PipelinesView extends Vue {
  readonly projectsContext = projectsModule.context(this.$store);
  readonly pipelinesContext = pipelinesModule.context(this.$store);
  readonly datasetContext = datasetsModule.context(this.$store);

  selected?: number | null = null;

  get activeProjectId() {
    return this.projectsContext.getters.activeProjectId;
  }

  @Watch("selected")
  pipelineChanged(index: number | null) {
    if (index !== null && index !== undefined) {
      const pipeline = this.pipelines[index];
      this.pipelinesContext.mutations.setActivePipelineId(pipeline.id);
    } else {
      this.pipelinesContext.mutations.setActivePipelineId(null);
    }
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
    // await this.pipelinesContext.actions.(id);
    // await this.projectsContext.actions.getChannelStackImage();
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

  mounted() {
    this.refreshPipelines();
  }
}
</script>

<style scoped>
.scroll-view {
  height: calc(100vh - 132px);
}
</style>
