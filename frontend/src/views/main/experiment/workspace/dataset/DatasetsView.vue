<template>
  <v-card tile>
    <v-toolbar flat dense color="grey lighten-4">
      <!--      <CreateDatasetDialog/>-->
      <!--      <UploadDatasetDialog class="ml-2"/>-->
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
<!--    <v-toolbar dense flat>-->
<!--      <v-text-field-->
<!--        v-model="search"-->
<!--        append-icon="mdi-magnify"-->
<!--        label="Search"-->
<!--        single-line-->
<!--        hide-details-->
<!--        clearable-->
<!--        flat-->
<!--      />-->
<!--    </v-toolbar>-->
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
                <v-list-item-title>Dataset {{ item.id }}</v-list-item-title>
                <v-list-item-subtitle>
                  <span class="x-small"><v-icon small>mdi-calendar-outline</v-icon> {{ item.createdAt }}</span>
                </v-list-item-subtitle>
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
                <div>{{ item.location }}</div>
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
  import CreateDatasetDialog from '@/views/main/experiment/workspace/dataset/CreateDatasetDialog.vue';
  import UploadDatasetDialog from '@/views/main/experiment/workspace/dataset/UploadDatasetDialog.vue';
  import { Component, Vue } from 'vue-property-decorator';

  @Component({
    components: { UploadDatasetDialog, InfoCard, CreateDatasetDialog },
  })
  export default class DatasetsView extends Vue {
    readonly experimentContext = experimentModule.context(this.$store);
    readonly datasetContext = datasetModule.context(this.$store);

    readonly apiUrl = apiUrl;

    readonly icons = {
      pending: 'mdi-progress-upload',
    };

    get datasets() {
      return this.datasetContext.getters.datasets;
    }

    get items() {
      if (this.datasets) {
        return this.datasets.map((dataset) => {
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
      this.datasetContext.mutations.setDatasets([]);
    }
  }
</script>

<style scoped>
  .card-title {
    padding-bottom: 0;
  }

  .scroll-view {
    height: calc(100vh - 232px);
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
