<template>
  <v-card tile>
    <v-toolbar flat dense color="grey lighten-4">
      <v-spacer></v-spacer>
      <v-tooltip bottom>
        <template v-slot:activator="{ on }">
          <v-btn icon small v-on="on" @click="refreshDatasets">
            <v-icon small>mdi-refresh</v-icon>
          </v-btn>
        </template>
        <span>Refresh datasets</span>
      </v-tooltip>
    </v-toolbar>
    <v-list dense two-line class="overflow-y-auto scroll-view pa-0">
      <v-list-item-group v-model="selected" color="primary">
        <v-list-item v-for="item in items" :key="item.uid">
          <v-list-item-icon>
            <v-tooltip bottom>
              <template v-slot:activator="{ on }">
                <v-icon v-on="on">{{ item.icon }}</v-icon>
              </template>
              <span>Status: {{ item.status }}</span>
            </v-tooltip>
          </v-list-item-icon>

          <v-list-item-content>
            <v-list-item-title>{{ item.name }}</v-list-item-title>
            <v-list-item-subtitle v-if="item.description">{{ item.description }}</v-list-item-subtitle>
            <v-list-item-subtitle class="font-weight-light">{{ item.createdAt }}</v-list-item-subtitle>
          </v-list-item-content>

          <v-list-item-action>
            <v-row>
              <v-tooltip bottom>
                <template v-slot:activator="{ on }">
                  <v-btn
                    icon
                    small
                    v-on="on"
                    download
                    color="primary lighten-2"
                    @click.stop=""
                    :href="`${apiUrl}/datasets/${item.id}/download`"
                  >
                    <v-icon small>mdi-download-outline</v-icon>
                  </v-btn>
                </template>
                <span>Download dataset</span>
              </v-tooltip>
              <v-dialog v-model="dialog" scrollable max-width="600px">
                <template v-slot:activator="{ on, attrs }">
                  <v-btn icon small v-bind="attrs" v-on="on" color="primary lighten-2" @click="name = item.name; description = item.description">
                    <v-icon small>mdi-pencil-outline</v-icon>
                  </v-btn>
                </template>
                <v-card>
                  <v-card-title>Edit Dataset</v-card-title>
                  <v-card-text>
                    <v-text-field label="Name" v-model="name" />
                    <v-text-field label="Description" v-model="description" />
                  </v-card-text>
                  <v-card-actions>
                    <v-spacer />
                    <v-btn text @click="dialog = false">Cancel</v-btn>
                    <v-btn color="primary" text @click="updateDataset(item.id)">Update</v-btn>
                  </v-card-actions>
                </v-card>
              </v-dialog>
              <v-tooltip bottom>
                <template v-slot:activator="{ on }">
                  <v-btn icon small v-on="on" color="secondary lighten-2" @click.stop="deleteDataset(item.id)">
                    <v-icon small>mdi-delete-outline</v-icon>
                  </v-btn>
                </template>
                <span>Delete dataset</span>
              </v-tooltip>
            </v-row>
          </v-list-item-action>
        </v-list-item>
      </v-list-item-group>
    </v-list>
  </v-card>
</template>

<script lang="ts">
import { apiUrl } from "@/env";
import { datasetsModule } from "@/modules/datasets";
import { projectsModule } from "@/modules/projects";
import { Component, Vue, Watch } from "vue-property-decorator";
import { centroidsModule } from "@/modules/centroids";
import { resultsModule } from "@/modules/results";

@Component
export default class DatasetsView extends Vue {
  readonly projectsContext = projectsModule.context(this.$store);
  readonly datasetsContext = datasetsModule.context(this.$store);
  readonly resultContext = resultsModule.context(this.$store);
  readonly centroidsContext = centroidsModule.context(this.$store);

  readonly apiUrl = apiUrl;
  readonly icons = {
    pending: "mdi-progress-clock",
    ready: "mdi-check-circle-outline",
  };

  dialog = false;
  name: string | null = null;
  description: string | null = null;

  selected?: number | null = null;

  @Watch("selected")
  datasetChanged(index?: number | null) {
    if (index !== null && index !== undefined) {
      const dataset = this.datasets[index];
      this.datasetsContext.actions.setActiveDatasetId(dataset.id);
      Promise.all([
        this.centroidsContext.actions.getCentroids({ datasetId: dataset.id }),
        this.resultContext.actions.getDatasetResults(dataset.id),
      ]);
    } else {
      this.datasetsContext.actions.setActiveDatasetId(null);
    }
  }

  get datasets() {
    return this.datasetsContext.getters.datasets;
  }

  get items() {
    return this.datasets.map((dataset) => {
      return Object.assign({}, dataset, {
        icon: this.icons[dataset.status],
        createdAt: new Date(dataset.created_at).toUTCString(),
      });
    });
  }

  async deleteDataset(id: number) {
    if (self.confirm("Do you really want to delete dataset?")) {
      await this.datasetsContext.actions.deleteDataset(id);
    }
  }

  async refreshDatasets() {
    const projectId = this.projectsContext.getters.activeProjectId;
    if (projectId) {
      await this.datasetsContext.actions.getProjectDatasets(projectId);
    }
  }

  async updateDataset(datasetId: number) {
    this.dialog = false;
    await this.datasetsContext.actions.updateDataset({
      datasetId: datasetId,
      data: { name: this.name, description: this.description },
    });
  }

  async mounted() {
    await this.refreshDatasets();
  }
}
</script>

<style scoped>
.scroll-view {
  height: calc(100vh - 132px);
}
</style>
