<template>
  <v-banner v-if="!activeDatasetId" icon="mdi-alert-circle-outline">Please select dataset</v-banner>
  <div v-else class="results-view">
    <v-toolbar flat dense color="grey lighten-4">
      <v-select
        :items="heatmaps"
        v-model="heatmap"
        label="Color"
        return-object
        hide-details
        solo
        flat
        clearable
        dense
        item-value="value"
        item-text="label"
      />
      <v-spacer />
      <v-tooltip bottom>
        <template v-slot:activator="{ on }">
          <v-btn icon small v-on="on" @click="refreshResults">
            <v-icon small>mdi-refresh</v-icon>
          </v-btn>
        </template>
        <span>Refresh results</span>
      </v-tooltip>
    </v-toolbar>
    <v-list dense two-line class="pa-0">
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
                  <v-btn
                    icon
                    small
                    v-on="on"
                    download
                    color="primary lighten-2"
                    @click.stop=""
                    :href="`${apiUrl}/results/${item.id}/download`"
                  >
                    <v-icon small>mdi-download-outline</v-icon>
                  </v-btn>
                </template>
                <span>Download result</span>
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
                <span>Edit result</span>
              </v-tooltip>
              <v-tooltip bottom>
                <template v-slot:activator="{ on }">
                  <v-btn icon small v-on="on" color="secondary lighten-2" @click.stop="deleteResult(item.id)">
                    <v-icon small>mdi-delete-outline</v-icon>
                  </v-btn>
                </template>
                <span>Delete result</span>
              </v-tooltip>
            </v-row>
          </v-list-item-action>
        </v-list-item>
      </v-list-item-group>
    </v-list>
    <v-dialog v-model="dialog" scrollable max-width="600px">
      <v-card>
        <v-card-title>Edit Result</v-card-title>
        <v-card-text>
          <v-text-field label="Name" v-model="name" />
          <v-text-field label="Description" v-model="description" />
        </v-card-text>
        <v-card-actions>
          <v-spacer />
          <v-btn text @click="dialog = false">Cancel</v-btn>
          <v-btn color="primary" text @click="updateResult()">Update</v-btn>
        </v-card-actions>
      </v-card>
    </v-dialog>
  </div>
</template>

<script lang="ts">
import { projectsModule } from "@/modules/projects";
import { Component, Vue, Watch } from "vue-property-decorator";
import { datasetsModule } from "@/modules/datasets";
import { apiUrl } from "@/env";
import { cellsModule } from "@/modules/cells";

@Component
export default class ResultsView extends Vue {
  readonly projectsContext = projectsModule.context(this.$store);
  readonly datasetContext = datasetsModule.context(this.$store);
  readonly cellsContext = cellsModule.context(this.$store);

  readonly apiUrl = apiUrl;

  dialog = false;
  activeId: number | null = null;
  name: string | null = null;
  description: string | null = null;

  selected?: number | null = null;

  @Watch("selected")
  resultChanged(index: number | null | undefined) {
    if (index !== null && index !== undefined) {
      const result = this.results[index];
      this.cellsContext.mutations.setActiveResultId(result.id);
      this.getResultData(result.id);
    } else {
      this.cellsContext.mutations.setActiveResultId(null);
      this.cellsContext.mutations.setMarkers([]);
    }
  }

  get channels() {
    const acquisition = this.projectsContext.getters.activeAcquisition;
    return acquisition ? acquisition.channels : null;
  }

  get activeDataset() {
    return this.datasetContext.getters.activeDataset;
  }

  get markers() {
    return this.cellsContext.getters.markers;
  }

  get heatmaps() {
    const activeDataset = this.activeDataset;

    if (!activeDataset) {
      return [];
    }
    const channelItems =
      this.markers.length > 0
        ? this.markers.map((item) => {
            return {
              type: "marker",
              value: item,
              label: this.channels && this.channels[item] ? this.channels[item].customLabel : item,
            };
          })
        : this.datasetContext.getters.channels.map((item) => {
            return {
              type: "marker",
              value: item,
              label: this.channels && this.channels[item] ? this.channels[item].customLabel : item,
            };
          });
    const clusteringItems: any[] = [];
    if (this.cellsContext.getters.activeResult?.output.leiden) {
      clusteringItems.push({
        type: "clustering",
        value: "leiden",
        label: "clustering [leiden]",
      });
    }
    if (this.cellsContext.getters.activeResult?.output.louvain) {
      clusteringItems.push({
        type: "clustering",
        value: "louvain",
        label: "clustering [louvain]",
      });
    }
    const annotationsItems = [
      {
        type: "annotation",
        value: "annotation",
        label: "annotation",
      },
    ];
    return channelItems.concat(clusteringItems, annotationsItems);
  }

  get heatmap() {
    return this.cellsContext.getters.heatmap;
  }

  set heatmap(value) {
    if (value === undefined) {
      value = null;
    }
    this.cellsContext.mutations.setHeatmap(value);
    this.projectsContext.actions.getChannelStackImage();
  }

  @Watch("heatmap")
  async heatmapChanged(value: { type: string; label: string; value: string } | null | undefined) {
    if (value && value.type === "annotation") {
      await this.projectsContext.actions.getAnnotationData();
    } else {
      await this.cellsContext.actions.getColorsData();
    }
  }

  get activeDatasetId() {
    return this.datasetContext.getters.activeDatasetId;
  }

  @Watch("activeDatasetId")
  async activeDatasetIdChanged(value: number | null) {
    this.selected = null;
  }

  get results() {
    return this.cellsContext.getters.results;
  }

  get items() {
    return this.results.map((gate) => {
      return Object.assign({}, gate, {
        createdAt: new Date(gate.created_at).toUTCString(),
      });
    });
  }

  async getResultData(id: number) {
    await this.cellsContext.actions.getResultData(id);
  }

  async deleteResult(id: number) {
    if (self.confirm("Do you really want to delete result?")) {
      await this.cellsContext.actions.deleteResult(id);
    }
  }

  async refreshResults() {
    if (this.activeDatasetId) {
      await this.cellsContext.actions.getDatasetResults(this.activeDatasetId);
    }
  }

  async updateResult() {
    this.dialog = false;
    if (this.activeId) {
      await this.cellsContext.actions.updateResult({
        resultId: this.activeId,
        data: { name: this.name, description: this.description },
      });
    }
  }

  mounted() {
    this.refreshResults();
  }
}
</script>

<style scoped>
.results-view {
  width: 100%;
  height: 100%;
  overflow-y: auto;
}
</style>
