<template>
  <v-card tile>
    <v-card-title>
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
      :items="items"
      :search="search"
      v-model="selected"
      show-select
      hide-default-footer
      class="overflow-y-auto scroll-view"
      dense
      disable-pagination
      no-data-text="Please first select an acquisition"
    >
      <template v-slot:item.label="props">
        <v-edit-dialog
          :return-value.sync="props.item.label"
          @save="save"
        > {{ props.item.label }}
          <template v-slot:input>
            <v-text-field
              v-model="props.item.label"
              :rules="[max25chars]"
              label="Edit"
              single-line
              counter
            ></v-text-field>
          </template>
        </v-edit-dialog>
      </template>
    </v-data-table>
  </v-card>
</template>

<script lang="ts">
  import { experimentModule } from '@/modules/experiment';
  import { IChannel } from '@/modules/experiment/models';
  import { settingsModule } from '@/modules/settings';
  import * as R from 'ramda';
  import { Component, Vue } from 'vue-property-decorator';

  @Component
  export default class ChannelsView extends Vue {
    readonly settingsModule = settingsModule.context(this.$store);
    readonly experimentContext = experimentModule.context(this.$store);

    search = '';

    readonly headers = [
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

    max25chars = v => v.length <= 25 || 'Input too long!';

    get channels() {
      const acquisition = this.experimentContext.getters.activeAcquisition;
      return acquisition && acquisition.channels ? acquisition.channels : [];
    }

    get items() {
      return this.channels.map((channel) => {
        const settings = this.settingsModule.getters.channelSettings(channel.id);
        return {
          id: channel.id,
          label: settings && settings.customLabel ? settings.customLabel : channel.label,
          metal: channel.metal,
          mass: channel.mass,
        };
      });
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

    save() {
      this.items.forEach((item) => {
        const settings = this.settingsModule.getters.channelSettings(item.id);
        if (!settings) {
          this.settingsModule.mutations.setChannelSettings({
            id: item.id,
            customLabel: item.label,
          });
        } else {
          if (settings.customLabel !== item.label) {
            this.settingsModule.mutations.setChannelSettings({
              ...settings,
              customLabel: item.label,
            });
          }
        }
      });
      if (this.items.length > 0) {
        this.experimentContext.actions.getChannelStackImage();
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
