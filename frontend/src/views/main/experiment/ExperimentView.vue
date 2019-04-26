<template>
  <LoadingView v-if="!experiment" text="Loading..."/>
  <v-container v-else fluid grid-list-md pa-2>
    <v-layout row wrap>
      <v-flex d-flex md3>
        <v-layout column wrap>
          <v-flex d-flex>
            <TreeView/>
          </v-flex>
          <v-flex d-flex>
            <InfoView/>
          </v-flex>
        </v-layout>
      </v-flex>
      <v-flex d-flex md6>
        <v-tabs v-model="active">
          <v-tab ripple>
            Blend
          </v-tab>
          <v-tab ripple>
            Tiles
          </v-tab>
          <v-spacer></v-spacer>
          <v-btn-toggle v-model="toggle_multiple" multiple>
            <v-btn flat>
              <v-icon>format_bold</v-icon>
            </v-btn>
            <v-btn flat>
              <v-icon>format_italic</v-icon>
            </v-btn>
            <v-btn flat>
              <v-icon>format_underlined</v-icon>
            </v-btn>
            <v-btn flat>
              <v-icon>format_color_fill</v-icon>
            </v-btn>
          </v-btn-toggle>
          <v-tab-item>
            <BlendView/>
          </v-tab-item>
          <v-tab-item>
            <TilesView/>
          </v-tab-item>
        </v-tabs>
      </v-flex>
      <v-flex d-flex md3>
        <ChannelsView/>
      </v-flex>
    </v-layout>
  </v-container>
</template>

<script lang="ts">
  import { Component, Vue } from 'vue-property-decorator';
  import { readActiveExperiment } from '@/modules/experiment/getters';
  import { dispatchGetExperiment, dispatchSetActiveExperimentId } from '@/modules/experiment/actions';
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

    toggle_multiple = [];

    get experiment() {
      return readActiveExperiment(this.$store);
    }

    async mounted() {
      const experimentId = parseInt(this.$router.currentRoute.params.id, 10);
      await dispatchSetActiveExperimentId(this.$store, { id: experimentId });
      await dispatchGetExperiment(this.$store, { id: experimentId });
    }
  }
</script>
