<template>
  <LoadingView v-if="!experimentData" text="Loading..."/>
  <v-container v-else fluid grid-list-md pa-1>
    <v-layout row>
      <v-flex v-if="showWorkspace" md3>
        <v-tabs v-model="tabWorkspace">
          <v-tab>Workspace</v-tab>
          <v-tab>Datasets</v-tab>
          <v-tab-item>
            <WorkspaceView :experiment="experimentData"/>
          </v-tab-item>
          <v-tab-item>
            <DatasetsView/>
          </v-tab-item>
        </v-tabs>
      </v-flex>
      <ImageView/>
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
  import LoadingView from '@/components/LoadingView.vue';
  import { experimentModule } from '@/modules/experiment';
  import { mainModule } from '@/modules/main';
  import ChannelsView from '@/views/main/experiment/ChannelsView.vue';
  import DatasetsView from '@/views/main/experiment/dataset/DatasetsView.vue';
  import ImageView from '@/views/main/experiment/image/ImageView.vue';
  import SettingsView from '@/views/main/experiment/settings/SettingsView.vue';
  import WorkspaceView from '@/views/main/experiment/WorkspaceView.vue';
  import { Component, Vue } from 'vue-property-decorator';

  @Component({
    components: {
      ImageView,
      DatasetsView,
      ChannelsView,
      WorkspaceView,
      LoadingView,
      SettingsView,
    },
  })
  export default class ExperimentView extends Vue {
    readonly mainContext = mainModule.context(this.$store);
    readonly experimentContext = experimentModule.context(this.$store);

    tabWorkspace = 0;

    get experiment() {
      return this.experimentContext.getters.activeExperiment;
    }

    get experimentData() {
      return this.experiment && this.experiment.slides ? this.experiment : undefined;
    }

    get showWorkspace() {
      return this.mainContext.getters.showWorkspace;
    }

    get showChannels() {
      return this.mainContext.getters.showChannels;
    }

    async mounted() {
      const experimentId = parseInt(this.$router.currentRoute.params.id, 10);
      this.experimentContext.mutations.setActiveExperimentId(experimentId);
      await this.experimentContext.actions.getExperimentData(experimentId);
    }

    beforeDestroy() {
      this.experimentContext.mutations.resetExperiment();
    }
  }
</script>
