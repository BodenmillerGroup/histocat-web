<template>
  <v-card tile>
    <v-toolbar flat dense color="grey lighten-4">
      <v-tooltip bottom>
        <template v-slot:activator="{ on }">
          <v-btn icon v-on="on" @click="refreshStatus">
            <v-icon>mdi-refresh</v-icon>
          </v-btn>
        </template>
        <span>Refresh datasets status</span>
      </v-tooltip>
      <CreateDatasetDialog/>
      <UploadDatasetDialog class="ml-2"/>
    </v-toolbar>
    <v-card-title class="card-title">
      <v-text-field
        v-model="search"
        append-icon="mdi-magnify"
        label="Search"
        single-line
        hide-details
        solo-inverted
        flat
      />
    </v-card-title>
    <v-card-text>
      <v-list class="overflow-y-auto scroll-view">
        <v-list-group
          v-for="item in items"
          :key="item.uid"
          :prepend-icon="item.icon"
          no-action
        >
          <template v-slot:activator>
            <v-list-item>
              <v-list-item-content>
                <v-list-item-title>{{ item.name }}</v-list-item-title>
                <v-list-item-subtitle>{{ item.status }}</v-list-item-subtitle>
              </v-list-item-content>
            </v-list-item>
          </template>

          <v-card tile flat>
            <v-img
              src="https://cdn.vuetifyjs.com/images/cards/desert.jpg"
              aspect-ratio="2.75"
            ></v-img>

            <v-card-title primary-title>
              <div>
                <h3>{{ item.name }}</h3>
                <div>{{ item.description }}</div>
                <span class="caption"><v-icon small>mdi-calendar-outline</v-icon> {{ item.createdAt }}</span>
              </div>
            </v-card-title>

            <v-card-actions>
              <v-btn text download :href="`${apiUrl}/api/v1/datasets/${item.id}/download`">Download</v-btn>
              <v-btn text color="error" @click="deleteDataset($event, item.id)">Delete</v-btn>
            </v-card-actions>
          </v-card>
        </v-list-group>
      </v-list>
    </v-card-text>
  </v-card>
</template>

<script lang="ts">
  import InfoCard from '@/components/InfoCard.vue';
  import { apiUrl } from '@/env';
  import { datasetModule } from '@/modules/datasets';
  import { experimentModule } from '@/modules/experiment';
  import UploadDatasetDialog from '@/views/main/experiment/workspace/dataset/UploadDatasetDialog.vue';
  import { Component, Vue } from 'vue-property-decorator';
  import CreateDatasetDialog from '@/views/main/experiment/workspace/dataset/CreateDatasetDialog.vue';

  @Component({
    components: { UploadDatasetDialog, InfoCard, CreateDatasetDialog, },
  })
  export default class DatasetsView extends Vue {
    readonly experimentContext = experimentModule.context(this.$store);
    readonly datasetContext = datasetModule.context(this.$store);

    search = '';
    readonly apiUrl = apiUrl;

    readonly icons = {
      pending: 'mdi-progress-upload',
    };

    get datasets() {
      return this.datasetContext.getters.datasets;
    }

    get items() {
      if (this.datasets) {
        return this.datasets
          .filter((dataset) => dataset.name.includes(this.search))
          .map((dataset) => {
            return Object.assign({}, dataset, {
              icon: this.icons[dataset.status],
              createdAt: new Date(dataset.created_at).toUTCString(),
            });
          });
      }
    }

    async downloadDataset(event, id: number, filename: string) {
      await this.datasetContext.actions.downloadDataset({ datasetId: id, filename: filename });
    }

    async deleteDataset(event, id: number) {
      if (self.confirm('Do you really want to delete dataset?')) {
        await this.datasetContext.actions.deleteDataset(id);
      }
    }

    async refreshStatus() {
      const experimentId = this.experimentContext.getters.activeExperimentId;
      if (experimentId) {
        await this.datasetContext.actions.getExperimentDatasets(experimentId);
      }
    }

    async mounted() {
      this.refreshStatus();
    }

    beforeDestroy() {
      this.datasetContext.mutations.setDatasets([]);
    }
  }
</script>

<style scoped>
  .card-title {
    padding-bottom: 0;
  }

  .scroll-view {
    height: calc(100vh - 220px);
  }
</style>

<style>
  .v-treeview-node__content {
    max-height: 24px;
  }

  .v-treeview-node__label {
    font-size: 10pt;
  }

  .v-treeview-node__root {
    min-height: 24px;
    font-size: 10pt;
  }

  .v-text-field.v-text-field--solo .v-input__control {
    min-height: 28px;
  }
</style>
