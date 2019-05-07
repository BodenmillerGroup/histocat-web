<template>
  <v-card tile>
    <v-card-title>
      <h4>Channels</h4>
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
      item-key="name"
      select-all
      disable-initial-sort
      hide-actions
    >
      <template slot="items" slot-scope="props">
        <td>
          <v-checkbox
            v-model="props.selected"
            primary
            hide-details
          ></v-checkbox>
        </td>
        <td>{{ props.item.name }}</td>
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
      },
      {
        text: 'Metal',
        sortable: true,
        value: 'metal',
        align: 'left',
      },
      {
        text: 'Mass',
        sortable: true,
        value: 'mass',
        align: 'left',
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
