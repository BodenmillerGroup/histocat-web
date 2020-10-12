<template>
  <v-banner v-if="!activeDatasetId" icon="mdi-alert-circle-outline">Please select dataset</v-banner>
  <v-card v-else tile>
    <v-toolbar flat dense color="grey lighten-4">
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
            <v-list-item-subtitle>{{ item.createdAt }}</v-list-item-subtitle>
          </v-list-item-content>

          <v-list-item-action>
            <v-tooltip bottom>
              <template v-slot:activator="{ on }">
                <v-btn icon v-on="on" download color="primary lighten-3" @click="loadResult(item.id)">
                  <v-icon>mdi-refresh-circle</v-icon>
                </v-btn>
              </template>
              <span>Load result</span>
            </v-tooltip>
          </v-list-item-action>
          <v-list-item-action>
            <v-tooltip bottom>
              <template v-slot:activator="{ on }">
                <v-btn icon v-on="on" color="secondary lighten-3" @click.stop="deleteResult(item.id)">
                  <v-icon>mdi-delete</v-icon>
                </v-btn>
              </template>
              <span>Delete result</span>
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
import {resultsModule} from "@/modules/results";

@Component
export default class ResultsView extends Vue {
  readonly projectsContext = projectsModule.context(this.$store);
  readonly resultsContext = resultsModule.context(this.$store);
  readonly datasetContext = datasetsModule.context(this.$store);

  selected?: number | null = null;

  @Watch("selected")
  resultChanged(index: number | null) {
    if (index !== null && index !== undefined) {
      const result = this.results[index];
      this.resultsContext.mutations.setActiveResultId(result.id);
    } else {
      this.resultsContext.mutations.setActiveResultId(null);
    }
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
