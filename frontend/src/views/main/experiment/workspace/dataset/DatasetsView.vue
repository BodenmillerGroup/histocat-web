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
    <v-list
        dense
        two-line
        class="overflow-y-auto scroll-view pa-0"
      >
        <v-list-item-group
          v-model="selected"
          color="primary"
        >
          <v-list-item
            v-for="item in items"
            :key="item.uid"
          >
            <v-list-item-icon>
              <v-tooltip bottom>
                <template v-slot:activator="{ on }">
                  <v-icon v-on="on">{{item.icon}}</v-icon>
                </template>
                <span>Status: {{item.status}}</span>
              </v-tooltip>
            </v-list-item-icon>

            <v-list-item-content>
              <v-list-item-title>Dataset {{item.id}}</v-list-item-title>
              <v-list-item-subtitle>{{item.createdAt}}</v-list-item-subtitle>
            </v-list-item-content>

            <v-list-item-action>
              <v-tooltip bottom>
                <template v-slot:activator="{ on }">
                  <v-btn
                    icon
                    v-on="on"
                    download
                    color="primary lighten-3"
                    @click.stop=""
                    :href="`${apiUrl}/api/v1/datasets/${item.id}/download`"
                  >
                    <v-icon>mdi-download</v-icon>
                  </v-btn>
                </template>
                <span>Download dataset</span>
              </v-tooltip>
            </v-list-item-action>
            <v-list-item-action>
              <v-tooltip bottom>
                <template v-slot:activator="{ on }">
                  <v-btn
                    icon
                    v-on="on"
                    color="secondary lighten-3"
                    @click.stop="deleteDataset($event, item.id)"
                  >
                    <v-icon>mdi-delete</v-icon>
                  </v-btn>
                </template>
                <span>Delete dataset</span>
              </v-tooltip>
            </v-list-item-action>
          </v-list-item>
        </v-list-item-group>
      </v-list>
  </v-card>
</template>

<script lang="ts">
  import { apiUrl } from '@/env';
  import { datasetModule } from '@/modules/datasets';
  import { experimentModule } from '@/modules/experiment';
  import { Component, Vue, Watch } from 'vue-property-decorator';

  @Component
  export default class DatasetsView extends Vue {
    readonly experimentContext = experimentModule.context(this.$store);
    readonly datasetContext = datasetModule.context(this.$store);

    readonly apiUrl = apiUrl;
    readonly icons = {
      pending: 'mdi-progress-clock',
      ready: 'mdi-check-circle-outline',
    };

    selected?: number | null = null;

    @Watch('selected')
    datasetChanged(index: number | null) {
      if (index !== null && index !== undefined) {
        const dataset = this.datasets[index];
        this.datasetContext.mutations.setActiveDataset(dataset);
      } else {
        this.datasetContext.mutations.setActiveDataset(undefined);
      }
    }

    get datasets() {
      return this.datasetContext.getters.datasets;
    }

    get items() {
      return this.datasets.map((dataset) => {
        return Object.assign({}, dataset, {
          icon: this.icons[dataset.status],
          createdAt: new Date(dataset.created_at).toUTCString(),
        });
      });
    }

    async deleteDataset(event, id: number) {
      if (self.confirm('Do you really want to delete dataset?')) {
        await this.datasetContext.actions.deleteDataset(id);
      }
    }

    async refreshDatasets() {
      const experimentId = this.experimentContext.getters.activeExperimentId;
      if (experimentId) {
        await this.datasetContext.actions.getExperimentDatasets(experimentId);
      }
    }

    async mounted() {
      this.refreshDatasets();
    }

    beforeDestroy() {
      this.datasetContext.mutations.reset();
    }
  }
</script>

<style scoped>
  .scroll-view {
    height: calc(100vh - 148px);
  }
</style>
