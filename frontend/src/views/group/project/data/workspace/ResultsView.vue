<template>
  <v-banner v-if="!activeDatasetId" icon="mdi-alert-circle-outline">Please select dataset</v-banner>
  <v-card v-else tile>
    <v-toolbar flat dense color="grey lighten-4">
      <v-select
        :items="heatmaps"
        v-model="heatmap"
        label="Heatmap"
        item-text="label"
        return-object
        hide-details
        solo
        flat
        clearable
        dense
      />
      <v-spacer></v-spacer>
      <v-tooltip bottom>
        <template v-slot:activator="{ on }">
          <v-btn icon small v-on="on" @click="refreshResults">
            <v-icon small>mdi-refresh</v-icon>
          </v-btn>
        </template>
        <span>Refresh results</span>
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
                <v-btn icon small v-on="on" download color="primary lighten-2" @click="loadResult(item.id)">
                  <v-icon small>mdi-refresh-circle</v-icon>
                </v-btn>
              </template>
              <span>Load result</span>
            </v-tooltip>
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
              <v-dialog v-model="dialog" scrollable max-width="600px">
                <template v-slot:activator="{ on, attrs }">
                  <v-btn icon small v-bind="attrs" v-on="on" color="primary lighten-2" @click="name = item.name; description = item.description">
                    <v-icon small>mdi-pencil-outline</v-icon>
                  </v-btn>
                </template>
                <v-card>
                  <v-card-title>Edit Result</v-card-title>
                  <v-card-text>
                    <v-text-field label="Name" v-model="name" />
                    <v-text-field label="Description" v-model="description" />
                  </v-card-text>
                  <v-card-actions>
                    <v-spacer />
                    <v-btn text @click="dialog = false">Cancel</v-btn>
                    <v-btn color="primary" text @click="updateResult(item.id)">Update</v-btn>
                  </v-card-actions>
                </v-card>
              </v-dialog>
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
  </v-card>
</template>

<script lang="ts">
import { projectsModule } from "@/modules/projects";
import { Component, Vue, Watch } from "vue-property-decorator";
import { datasetsModule } from "@/modules/datasets";
import { resultsModule } from "@/modules/results";
import { apiUrl } from "@/env";

@Component
export default class ResultsView extends Vue {
  readonly projectsContext = projectsModule.context(this.$store);
  readonly resultsContext = resultsModule.context(this.$store);
  readonly datasetContext = datasetsModule.context(this.$store);

  readonly apiUrl = apiUrl;

  dialog = false;
  name: string | null = null;
  description: string | null = null;

  selected?: number | null = null;

  @Watch("selected")
  resultChanged(index: number | null | undefined) {
    if (index !== null && index !== undefined) {
      const result = this.results[index];
      this.resultsContext.mutations.setActiveResultId(result.id);
    } else {
      this.resultsContext.mutations.setActiveResultId(null);
    }
  }

  get heatmaps() {
    return this.datasetContext.getters.heatmaps;
  }

  get heatmap() {
    return this.resultsContext.getters.heatmap;
  }

  set heatmap(value) {
    if (value === undefined) {
      value = null;
    }
    this.resultsContext.mutations.setHeatmap(value);
  }

  get activeDatasetId() {
    return this.datasetContext.getters.activeDatasetId;
  }

  get results() {
    return this.resultsContext.getters.results;
  }

  get items() {
    return this.results.map((gate) => {
      return Object.assign({}, gate, {
        createdAt: new Date(gate.created_at).toUTCString(),
      });
    });
  }

  async loadResult(id: number) {
    await this.resultsContext.actions.loadResult(id);
  }

  async deleteResult(id: number) {
    if (self.confirm("Do you really want to delete result?")) {
      await this.resultsContext.actions.deleteResult(id);
    }
  }

  async refreshResults() {
    if (this.activeDatasetId) {
      await this.resultsContext.actions.getDatasetResults(this.activeDatasetId);
    }
  }

  async updateResult(resultId: number) {
    this.dialog = false;
    await this.resultsContext.actions.updateResult({
      resultId: resultId,
      data: { name: this.name, description: this.description },
    });
  }

  mounted() {
    this.refreshResults();
  }
}
</script>

<style scoped>
.scroll-view {
  height: calc(100vh - 132px);
}
</style>
