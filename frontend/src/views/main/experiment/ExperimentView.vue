<template>
  <LoadingView v-if="!dataset" text="Loading..."/>
  <v-container v-else fluid grid-list-md pa-2>
    <v-layout row>
      <v-flex v-if="showWorkspace" md3>
        <v-flex>
          <TreeView :dataset="dataset"/>
        </v-flex>
        <v-flex>
          <InfoView/>
        </v-flex>
      </v-flex>
      <v-flex grow>
        <v-toolbar card dense>
          <v-toolbar-side-icon></v-toolbar-side-icon>
          <v-toolbar-title>Image View</v-toolbar-title>
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
        <BlendView/>
      </v-flex>
      <v-flex v-if="showChannels" md3>
        <v-flex>
          <ChannelsView/>
        </v-flex>
        <v-flex>
          <SettingsView/>
        </v-flex>
      </v-flex>
    </v-layout>
  </v-container>
</template>

<script lang="ts">
  import { Component, Vue } from 'vue-property-decorator';
  import { readExperimentDataset } from '@/modules/experiment/getters';
  import { dispatchGetExperimentDataset } from '@/modules/experiment/actions';
  import LoadingView from '@/components/LoadingView.vue';
  import TreeView from '@/views/main/experiment/TreeView.vue';
  import InfoView from '@/views/main/experiment/InfoView.vue';
  import ChannelsView from '@/views/main/experiment/ChannelsView.vue';
  import SettingsView from '@/views/main/experiment/SettingsView.vue';
  import BlendView from '@/views/main/experiment/BlendView.vue';
  import TilesView from '@/views/main/experiment/TilesView.vue';
  import { commitSetSelectedExperimentId } from '@/modules/experiment/mutations';

  @Component({
    components: { ChannelsView, TilesView, BlendView, InfoView, TreeView, LoadingView, SettingsView },
  })
  export default class ExperimentView extends Vue {

    toggleMultiple = ['showWorkspace', 'showChannels'];

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
      commitSetSelectedExperimentId(this.$store, { id: experimentId });
      await dispatchGetExperimentDataset(this.$store, { id: experimentId });
    }
  }
</script>
