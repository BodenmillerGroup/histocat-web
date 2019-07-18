<template>
  <v-card tile>
    <v-card-title>
      Channels
      <v-spacer/>
      <v-text-field
        v-model="search"
        append-icon="mdi-magnify"
        label="Search"
        single-line
        hide-details
        clearable
      />
    </v-card-title>
    <v-data-table
      :headers="headers"
      :items="channels"
      :search="search"
      v-model="selected"
      show-select
      hide-default-footer
      class="scroll-y scroll-view"
      dense
      disable-pagination
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
  import { experimentModule } from '@/modules/experiment';
  import { IChannel } from '@/modules/experiment/models';
  import { Component, Vue, Watch } from 'vue-property-decorator';

  @Component
  export default class ChannelsView extends Vue {
    readonly experimentContext = experimentModule.context(this.$store);

    search = '';
    selected = [];

    headers = [
      {
        text: 'Label',
        sortable: true,
        value: 'label',
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
      const acquisition = this.experimentContext.getters.activeAcquisition;
      return acquisition && acquisition.channels;
    }

    @Watch('selected')
    onSelectedChanged(items: IChannel[]) {
      const selectedMetals = items.map(item => item.metal);
      this.experimentContext.mutations.setSelectedMetals(selectedMetals);
    }
  }
</script>

<style scoped>
  table.v-table tbody td, table.v-table tbody th {
    height: 35px;
  }

  .scroll-view {
    height: calc(50vh - 100px);
  }
</style>

<style>
  .channels-table table {
    table-layout: fixed;
  }
</style>
