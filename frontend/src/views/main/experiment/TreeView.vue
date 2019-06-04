<template>
  <v-card tile>
    <v-card-title class="card-title">
      <h4>Workspace</h4>
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
    <v-card-text>
      <v-treeview
        :items="items"
        :search="search"
        :filter="filter"
        :active.sync="active"
        activatable
        open-on-click
        transition
        return-object
        class="scroll-y workspace-tree"
      >
        <template v-slot:prepend="{ item }">
          <v-icon>
            {{ icons[item.type] }}
          </v-icon>
        </template>
      </v-treeview>
    </v-card-text>
  </v-card>
</template>

<script lang="ts">
  import { Component, Prop, Vue, Watch } from 'vue-property-decorator';
  import { IAcquisition, IExperimentDataset } from '@/modules/experiment/models';
  import {
    commitSetChannels,
    commitSetSelectedAcquisition,
    commitSetSelectedMeta,
  } from '@/modules/experiment/mutations';

  @Component
  export default class TreeView extends Vue {

    @Prop(Object) dataset?: IExperimentDataset;

    search = null;
    open = [];
    active = [];

    icons = {
      slide: 'mdi-folder-outline',
      acquisition: 'mdi-film',
    };

    @Watch('selected')
    onActiveChanged(item: object) {
      if (!item) {
        return;
      }
      if (item.hasOwnProperty('meta')) {
        commitSetSelectedMeta(this.$store, item['meta']);
      }
      if (item.hasOwnProperty('channels')) {
        commitSetSelectedAcquisition(this.$store, item as IAcquisition);
        commitSetChannels(this.$store, item['channels']);
      }
    }

    get filter() {
      return (item, search, textKey) => item[textKey].indexOf(search) > -1;
    }

    get items() {
      if (this.dataset) {
        return this.dataset.slides.map((slide) => {
          const acquisitions = slide.acquisitions.map((acquisition) => {
            return Object.assign({}, acquisition, { type: 'acquisition' });
          });
          return Object.assign({}, slide, { type: 'slide', children: acquisitions });
        });
      }
    }

    get selected() {
      if (!this.active.length) {
        return undefined;
      }
      return this.active[0];
    }
  }
</script>

<style scoped>
  .card-title {
    padding-bottom: 0px;
  }

  .workspace-tree {
    max-height: 37vh;
  }
</style>

<style>
  .v-treeview-node__label {
    font-size: 1em;
  }

  .v-treeview-node__root {
    min-height: 24px;
  }

  .v-text-field.v-text-field--solo .v-input__control {
    min-height: 28px;
  }
</style>
