<template>
  <v-card tile>
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
      <v-list>
        <v-list-group
          v-for="item in items"
          :key="item.uid"
          :prepend-icon="item.icon"
          no-action
        >
          <template v-slot:activator>
            <v-list-tile>
              <v-list-tile-content>
                <v-list-tile-title>{{ item.name }}</v-list-tile-title>
                <v-list-tile-sub-title>{{ item.status }}</v-list-tile-sub-title>
              </v-list-tile-content>
            </v-list-tile>
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
              <v-btn flat color="red" @click="deleteDataset($event, item.id)">Delete</v-btn>
            </v-card-actions>
          </v-card>
        </v-list-group>
      </v-list>
    </v-card-text>
  </v-card>
</template>

<script lang="ts">
  import InfoCard from '@/components/InfoCard.vue';
  import { dispatchDeleteDataset } from '@/modules/experiment/actions';
  import { readDatasets } from '@/modules/experiment/getters';
  import { Component, Vue } from 'vue-property-decorator';

  @Component({
    components: { InfoCard },
  })
  export default class DatasetsView extends Vue {

    search = '';

    icons = {
      pending: 'mdi-progress-upload',
    };

    get datasets() {
      return readDatasets(this.$store);
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

    async deleteDataset(event, id: number) {
      const res = await this.$confirm('Do you really want to delete dataset?', { title: 'Warning' });
      if (res) {
        await dispatchDeleteDataset(this.$store, id);
      }
    }
  }
</script>

<style scoped>
  .card-title {
    padding-bottom: 0;
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
