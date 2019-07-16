<template>
  <LoadingView v-if="!experimentData" text="Loading..."/>
  <v-container v-else fluid grid-list-md pa-2>
    <v-layout row>
      <v-flex v-if="showWorkspace" md3>
        <v-tabs v-model="tabWorkspace">
          <v-tab>Workspace</v-tab>
          <v-tab>Datasets</v-tab>
          <v-tab-item>
            <WorkspaceView :experiment="experimentData" class="tree-view"/>
          </v-tab-item>
          <v-tab-item>
            <DatasetsView class="tree-view"/>
          </v-tab-item>
        </v-tabs>
      </v-flex>
      <v-flex :class="viewerClass">
        <v-toolbar card dense>
          <v-toolbar-side-icon></v-toolbar-side-icon>
          <v-menu>
            <template v-slot:activator="{ on }">
              <v-toolbar-title v-on="on">
                <span>Tools</span>
                <v-icon>mdi-menu-down</v-icon>
              </v-toolbar-title>
            </template>
            <v-list>
              <CreateDatasetDialog></CreateDatasetDialog>
            </v-list>
          </v-menu>
          <v-spacer/>
          <v-btn-toggle v-model="toggleMultiple" multiple>
            <v-btn small flat value='showWorkspace'>
              <v-icon>mdi-file-tree</v-icon>
              <span>Workspace</span>
            </v-btn>
            <v-btn small flat value='showChannels'>
              <v-icon>mdi-format-list-checkbox</v-icon>
              <span>Channels</span>
            </v-btn>
          </v-btn-toggle>
        </v-toolbar>
        <v-flex>
          <v-tabs v-model="tabImageView">
            <v-tab>Blend</v-tab>
            <v-tab>Tiles</v-tab>
            <v-tab-item>
              <BlendView class="image-view"/>
            </v-tab-item>
            <v-tab-item>
              <TilesView class="image-view"/>
            </v-tab-item>
          </v-tabs>
        </v-flex>
      </v-flex>
      <v-flex v-if="showChannels" md3>
        <v-flex>
          <ChannelsView class="channels-view"/>
        </v-flex>
        <v-flex>
          <SettingsView class="settings-view"/>
        </v-flex>
      </v-flex>
    </v-layout>
  </v-container>
</template>

<script lang="ts">
  import LoadingView from '@/components/LoadingView.vue';
  import { experimentModule } from '@/modules/experiment';
  import BlendView from '@/views/main/experiment/BlendView.vue';
  import ChannelsView from '@/views/main/experiment/ChannelsView.vue';
  import CreateDatasetDialog from '@/views/main/experiment/CreateDatasetDialog.vue';
  import DatasetsView from '@/views/main/experiment/DatasetsView.vue';
  import SettingsView from '@/views/main/experiment/SettingsView.vue';
  import TilesView from '@/views/main/experiment/TilesView.vue';
  import WorkspaceView from '@/views/main/experiment/WorkspaceView.vue';
  import { Component, Vue } from 'vue-property-decorator';

  @Component({
    components: {
      TilesView,
      DatasetsView,
      CreateDatasetDialog,
      ChannelsView,
      BlendView,
      WorkspaceView,
      LoadingView,
      SettingsView,
    },
  })
  export default class ExperimentView extends Vue {
    readonly experimentContext = experimentModule.context(this.$store);

    tabWorkspace = 0;
    tabImageView = 0;

    toggleMultiple = ['showWorkspace', 'showChannels'];

    get experiment() {
      return this.experimentContext.getters.activeExperiment;
    }

    get experimentData() {
      return this.experiment && this.experiment.slides ? this.experiment : undefined;
    }

    get showWorkspace() {
      return this.toggleMultiple.includes('showWorkspace');
    }

    get showChannels() {
      return this.toggleMultiple.includes('showChannels');
    }

    get viewerClass() {
      if (this.showWorkspace && this.showChannels) {
        return 'md6';
      }
      if (this.showWorkspace || this.showChannels) {
        return 'md9';
      }
      return 'md12';
    }

    async mounted() {
      const experimentId = parseInt(this.$router.currentRoute.params.id, 10);
      this.experimentContext.mutations.setActiveExperimentId(experimentId);
      await this.experimentContext.actions.getExperimentData(experimentId);
    }

    beforeDestroy() {
      this.experimentContext.mutations.setActiveExperimentId(undefined);
    }
  }
</script>

<style scoped>
  .tree-view {
    height: calc(100vh - 115px);
  }

  .image-view {
    height: calc(100vh - 163px);
  }

  .channels-view {
    height: calc(50vh - 50px);
  }

  .settings-view {
    height: calc(50vh - 60px);
  }
</style>