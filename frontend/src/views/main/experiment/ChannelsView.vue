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
      no-data-text="Please first select an acquisition"
    >

    </v-data-table>
  </v-card>
</template>

<script lang="ts">
  import { experimentModule } from '@/modules/experiment';
  import { IChannel } from '@/modules/experiment/models';
  import * as R from 'ramda';
  import { Component, Vue } from 'vue-property-decorator';

  @Component
  export default class ChannelsView extends Vue {
    readonly experimentContext = experimentModule.context(this.$store);

    search = '';

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
      return acquisition && acquisition.channels ? acquisition.channels : [];
    }

    get selectedMetals() {
      return this.experimentContext.getters.selectedMetals;
    }

    get selected() {
      return this.channels.filter((channel) => {
        if (this.selectedMetals.includes(channel.metal)) {
          return channel;
        }
      });
    }

    set selected(items: IChannel[]) {
      const selectedMetals = items.map(item => item.metal);
      if (!R.equals(this.selectedMetals, selectedMetals)) {
        this.experimentContext.mutations.setSelectedMetals(selectedMetals);
      }
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
