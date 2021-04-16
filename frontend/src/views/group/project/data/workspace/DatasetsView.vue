<template>
  <v-card tile>
    <v-toolbar flat dense color="grey lighten-4">
      <UploadButton label="Upload dataset" :upload="upload" />
      <v-spacer />
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
            <v-menu :close-on-content-click="false" :nudge-width="200" offset-x>
              <template v-slot:activator="{ on }">
                <v-btn v-on="on" small icon color="grey" @click.stop="">
                  <v-icon>{{ item.icon }}</v-icon>
                </v-btn>
              </template>
              <v-card tile flat class="card overflow-y-auto">
                <TreeView :data="item.meta" :options="{ modifiable: false, rootObjectKey: 'meta' }" />
              </v-card>
            </v-menu>
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
                <span>Edit dataset</span>
              </v-tooltip>
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
    <v-dialog v-model="dialog" scrollable max-width="600px">
      <v-card>
        <v-card-title>Edit Dataset</v-card-title>
        <v-card-text>
          <v-text-field label="Name" v-model="name" />
          <v-text-field label="Description" v-model="description" />
        </v-card-text>
        <v-card-actions>
          <v-spacer />
          <v-btn text @click="dialog = false">Cancel</v-btn>
          <v-btn color="primary" text @click="updateDataset()">Update</v-btn>
        </v-card-actions>
      </v-card>
    </v-dialog>
  </v-card>
</template>

<script lang="ts">
import { apiUrl } from "@/env";
import { datasetsModule } from "@/modules/datasets";
import { projectsModule } from "@/modules/projects";
import { Component, Vue, Watch } from "vue-property-decorator";
import UploadButton from "@/components/UploadButton.vue";
import TreeView from "@/components/vue-json-tree-view/TreeView.vue";
import { cellsModule } from "@/modules/cells";

@Component({
  components: { TreeView, UploadButton },
})
export default class DatasetsView extends Vue {
  readonly projectsContext = projectsModule.context(this.$store);
  readonly datasetsContext = datasetsModule.context(this.$store);
  readonly cellsContext = cellsModule.context(this.$store);

  readonly apiUrl = apiUrl;
  readonly icons = {
    pending: "mdi-progress-clock",
    ready: "mdi-check-circle-outline",
  };

  dialog = false;
  activeId: number | null = null;
  name: string | null = null;
  description: string | null = null;

  selected?: number | null = null;

  @Watch("selected")
  datasetChanged(index?: number | null) {
    this.cellsContext.mutations.reset();
    if (index !== null && index !== undefined) {
      const dataset = this.datasets[index];
      if (dataset.status === "ready") {
        this.datasetsContext.mutations.setActiveDatasetId(dataset.id);
        Promise.all([
          this.cellsContext.actions.initializeCells({ datasetId: dataset.id }),
          this.cellsContext.actions.getDatasetResults(dataset.id),
          this.projectsContext.actions.getChannelStackImage(),
        ]);
      }
    } else {
      this.cellsContext.mutations.setActiveResultId(null);
      this.datasetsContext.mutations.setActiveDatasetId(null);
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

  async updateDataset() {
    this.dialog = false;
    if (this.activeId) {
      await this.datasetsContext.actions.updateDataset({
        datasetId: this.activeId,
        data: { name: this.name, description: this.description },
      });
    }
  }

  async upload(data: FormData) {
    await this.datasetsContext.actions.uploadDataset({ id: this.projectsContext.getters.activeProjectId!, data: data });
  }

  async mounted() {
    await this.refreshDatasets();
  }
}
</script>

<style scoped>
.scroll-view {
  height: calc(50vh - 100px);
}
.card {
  height: 40vh;
}
</style>
