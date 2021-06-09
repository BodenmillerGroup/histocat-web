<template>
  <div class="datasets-view">
    <v-toolbar flat dense color="grey lighten-4">
      <v-menu offset-y>
        <template v-slot:activator="{ on }">
          <v-btn v-on="on" color="primary" elevation="1" small>
            <v-icon left small>mdi-cloud-upload</v-icon>
            Upload dataset
          </v-btn>
        </template>
        <v-list dense>
          <v-list-item @click="openUploadDialog('steinbock')">
            <v-list-item-title>steinbock</v-list-item-title>
          </v-list-item>
          <v-list-item @click="openUploadDialog('ImcSegmentationPipelineV1')">
            <v-list-item-title>ImcSegmentationPipelineV1</v-list-item-title>
          </v-list-item>
          <v-list-item @click="openUploadDialog('ImcSegmentationPipelineV2')">
            <v-list-item-title>ImcSegmentationPipelineV2</v-list-item-title>
          </v-list-item>
          <v-list-item @click="openUploadDialog('masks')">
            <v-list-item-title>Masks</v-list-item-title>
          </v-list-item>
        </v-list>
      </v-menu>
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
    <v-list dense class="pa-0">
      <v-list-item-group v-model="selected" color="primary">
        <v-list-item v-for="item in items" :key="item.id">
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
                    target="_blank"
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
                      dialogEditDataset = true;
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
    <v-dialog v-model="dialogEditDataset" scrollable max-width="600px">
      <v-card>
        <v-card-title>Edit Dataset</v-card-title>
        <v-card-text>
          <v-text-field label="Name" v-model="name" />
          <v-text-field label="Description" v-model="description" />
        </v-card-text>
        <v-card-actions>
          <v-spacer />
          <v-btn text @click="dialogEditDataset = false">Cancel</v-btn>
          <v-btn color="primary" text @click="updateDataset()">Update</v-btn>
        </v-card-actions>
      </v-card>
    </v-dialog>
    <v-dialog v-model="dialogUploadDataset" scrollable max-width="600px">
      <v-card>
        <v-card-title>Upload Dataset</v-card-title>
        <v-card-text>
          <v-form v-model="valid" ref="form">
            <v-text-field
              v-if="type === 'steinbock'"
              label="Masks folder name"
              v-model="masksFolder"
              :rules="requiredRule"
            />
            <v-text-field
              v-if="type === 'steinbock'"
              label="Regionprops folder name"
              v-model="regionpropsFolder"
              :rules="requiredRule"
            />
            <v-text-field
              v-if="type === 'steinbock'"
              label="Intensities folder name"
              v-model="intensitiesFolder"
              :rules="requiredRule"
            />
            <v-file-input v-model="file" label="File" show-size accept=".zip" :rules="requiredRule" />
          </v-form>
        </v-card-text>
        <v-card-actions>
          <v-spacer />
          <v-btn text @click="dialogUploadDataset = false">Cancel</v-btn>
          <v-btn color="primary" text @click="uploadDataset()" :disabled="!valid">Upload</v-btn>
        </v-card-actions>
      </v-card>
    </v-dialog>
  </div>
</template>

<script lang="ts">
import { apiUrl } from "@/env";
import { datasetsModule } from "@/modules/datasets";
import { projectsModule } from "@/modules/projects";
import { Component, Vue } from "vue-property-decorator";
import TreeView from "@/components/vue-json-tree-view/TreeView.vue";
import { cellsModule } from "@/modules/cells";
import { gatesModule } from "@/modules/gates";
import { annotationsModule } from "@/modules/annotations";
import { uiModule } from "@/modules/ui";
import { required } from "@/utils/validators";

@Component({
  components: { TreeView },
})
export default class DatasetsView extends Vue {
  readonly uiContext = uiModule.context(this.$store);
  readonly projectsContext = projectsModule.context(this.$store);
  readonly datasetsContext = datasetsModule.context(this.$store);
  readonly cellsContext = cellsModule.context(this.$store);
  readonly gatesContext = gatesModule.context(this.$store);
  readonly annotationsContext = annotationsModule.context(this.$store);

  readonly apiUrl = apiUrl;
  readonly icons = {
    pending: "mdi-progress-clock",
    ready: "mdi-check-circle-outline",
  };

  readonly requiredRule = [required];

  dialogUploadDataset = false;
  valid = true;
  file: File | null = null;
  type: string | null = null;
  masksFolder = "masks";
  intensitiesFolder = "object_intensities";
  regionpropsFolder = "object_regionprops";

  dialogEditDataset = false;
  activeId: number | null = null;
  name: string | null = null;
  description: string | null = null;

  get selected() {
    return this.datasetsContext.getters.activeDataset
      ? this.datasets.indexOf(this.datasetsContext.getters.activeDataset)
      : null;
  }

  set selected(index: number | null | undefined) {
    this.cellsContext.mutations.reset();
    if (index !== null && index !== undefined) {
      const dataset = this.datasets[index];
      if (dataset.status === "ready") {
        this.datasetsContext.mutations.setActiveDatasetId(dataset.id);
        this.annotationsContext.mutations.reset();
        Promise.all([
          this.cellsContext.actions.initializeCells({ datasetId: dataset.id }),
          this.cellsContext.actions.getDatasetResults(dataset.id),
          this.gatesContext.actions.getGates(),
        ]).then(() => {
          if (this.uiContext.getters.showMask) {
            this.projectsContext.actions.getChannelStackImage();
          }
        });
      }
    } else {
      this.cellsContext.mutations.setActiveResultId(null);
      this.datasetsContext.mutations.setActiveDatasetId(null);
      this.gatesContext.mutations.reset();
      this.annotationsContext.mutations.reset();
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
    this.dialogEditDataset = false;
    if (this.activeId) {
      await this.datasetsContext.actions.updateDataset({
        datasetId: this.activeId,
        data: { name: this.name, description: this.description },
      });
    }
  }

  openUploadDialog(type: string) {
    this.file = null;
    this.type = type;
    this.dialogUploadDataset = true;
  }

  async uploadDataset() {
    this.dialogUploadDataset = false;
    if ((this.$refs.form as any).validate() && this.type && this.file) {
      const formData = new FormData();
      formData.append("type", this.type);
      formData.append("masks_folder", this.masksFolder);
      formData.append("regionprops_folder", this.regionpropsFolder);
      formData.append("intensities_folder", this.intensitiesFolder);
      formData.append("file", this.file);
      await this.datasetsContext.actions.uploadDataset({
        projectId: this.projectsContext.getters.activeProjectId!,
        formData: formData,
      });
    }
  }

  async mounted() {
    await this.refreshDatasets();
  }
}
</script>

<style scoped>
.datasets-view {
  width: 100%;
  height: 100%;
  overflow-y: auto;
}
.card {
  height: 40vh;
}
</style>
