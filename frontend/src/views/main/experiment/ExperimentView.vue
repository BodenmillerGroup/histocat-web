<template>
  <LoadingView v-if="!dataset" text="Loading..."/>
  <v-container v-else fluid grid-list-md pa-2>
    <v-layout row wrap>
      <div v-if="showWorkspace" max-width="200px">
        <v-layout column wrap>
          <v-flex d-flex>
            <TreeView :dataset="dataset"/>
          </v-flex>
          <v-flex d-flex>
            <InfoView/>
          </v-flex>
        </v-layout>
      </div>
      <v-flex grow>
        <v-card flat tile>
          <v-toolbar card dense>
            <v-toolbar-side-icon></v-toolbar-side-icon>
            <v-toolbar-title>Blend</v-toolbar-title>
            <v-spacer/>
            <v-btn-toggle v-model="toggleMultiple" multiple>
              <v-btn flat value='showWorkspace'>
                <v-icon>mdi-file-tree</v-icon>
                <span>Workspace</span>
              </v-btn>
              <v-btn flat value='showChannels'>
                <v-icon>mdi-format-list-checkbox</v-icon>
                <span>Channels</span>
              </v-btn>
            </v-btn-toggle>
          </v-toolbar>
          <v-tabs v-model="active">
            <v-tab ripple>
              Blend
            </v-tab>
            <v-tab ripple>
              Tiles
            </v-tab>
            <v-tab-item>
              <BlendView/>
            </v-tab-item>
            <v-tab-item>
              <TilesView/>
            </v-tab-item>
          </v-tabs>
        </v-card>
      </v-flex>
      <div v-if="showChannels" max-width="200px">
        <ChannelsView/>
      </div>
    </v-layout>
  </v-container>
</template>

<script lang="ts">
  import { Component, Vue } from 'vue-property-decorator';
  import { readActiveExperiment, readExperimentDataset } from '@/modules/experiment/getters';
  import { dispatchGetExperimentDataset, dispatchSetActiveExperimentId } from '@/modules/experiment/actions';
  import LoadingView from '@/components/LoadingView.vue';
  import TreeView from '@/views/main/experiment/TreeView.vue';
  import InfoView from '@/views/main/experiment/InfoView.vue';
  import ChannelsView from '@/views/main/experiment/ChannelsView.vue';
  import BlendView from '@/views/main/experiment/BlendView.vue';
  import TilesView from '@/views/main/experiment/TilesView.vue';

  @Component({
    components: { ChannelsView, TilesView, BlendView, InfoView, TreeView, LoadingView },
  })
  export default class ExperimentView extends Vue {

    active = undefined;
    toggleMultiple = ['showWorkspace', 'showChannels'];

    get experiment() {
      return readActiveExperiment(this.$store);
    }

    get dataset() {
      return readExperimentDataset(this.$store);
    }

    get showWorkspace() {
      return this.toggleMultiple.includes('showWorkspace');
    }

    get showChannels() {
      return this.toggleMultiple.includes('showChannels');
    }

    async mounted() {
      const experimentId = parseInt(this.$router.currentRoute.params.id, 10);
      await dispatchSetActiveExperimentId(this.$store, { id: experimentId });
      await dispatchGetExperimentDataset(this.$store, { id: experimentId });
    }
  }
</script>
