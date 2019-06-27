<template>
  <v-card tile>
    <v-card-title class="card-title">
      <h4>Channels</h4>
      <v-spacer/>
      <v-text-field
        v-model="search"
        append-icon="mdi-magnify"
        label="Search"
        single-line
        hide-details
        clearable
        solo-inverted
        flat
      />
    </v-card-title>
    <v-data-table
      :headers="headers"
      :items="channels"
      :search="search"
      v-model="selected"
      item-key="label"
      select-all
      disable-initial-sort
      hide-actions
      class="scroll-y channels-table"
    >
      <template slot="items" slot-scope="props">
        <td>
          <v-checkbox
            v-model="props.selected"
            primary
            hide-details
          ></v-checkbox>
        </td>
        <td>{{ props.item.label }}</td>
        <td>{{ props.item.metal }}</td>
        <td>{{ props.item.mass }}</td>
      </template>
    </v-data-table>
  </v-card>
</template>

<script lang="ts">
  import { Component, Vue, Watch } from 'vue-property-decorator';
  import { readChannels } from '@/modules/experiment/getters';
  import { IChannel } from '@/modules/experiment/models';
  import { commitSetSelectedMetals } from '@/modules/experiment/mutations';

  @Component
  export default class ChannelsView extends Vue {

    search = '';
    selected = [];

    headers = [
      {
        text: 'Name',
        sortable: true,
        value: 'name',
        align: 'left',
        width: '50%',
      },
      {
        text: 'Metal',
        sortable: true,
        value: 'metal',
        align: 'left',
        width: '20%',
      },
      {
        text: 'Mass',
        sortable: true,
        value: 'mass',
        align: 'left',
        width: '20%',
      },
    ];

    get channels() {
      return readChannels(this.$store);
    }

    @Watch('selected')
    onSelectedChanged(items: IChannel[]) {
      const selectedMetals = items.map(item => item.metal);
      commitSetSelectedMetals(this.$store, selectedMetals);
    }
  }
</script>

<style scoped>
  .card-title {
    padding-bottom: 4px;
  }

  table.v-table tbody td, table.v-table tbody th {
    height: 35px;
  }

  .channels-table {
    height: calc(50vh - 100px);
  }
</style>

<style>
  .channels-table table {
    table-layout: fixed;
  }
</style>
