<template>
  <v-layout row>
    <v-flex :class="mainClass">
      <v-tabs v-model="tabImageView">
        <v-tab>Blend</v-tab>
        <v-tab>Tiles</v-tab>
        <v-spacer></v-spacer>
        <v-btn-toggle v-model="toggleChannels">
          <v-tooltip bottom>
            <template v-slot:activator="{ on }">
              <v-btn tile v-on="on" value='show'>
                <v-icon>mdi-format-list-checkbox</v-icon>
              </v-btn>
            </template>
            <span v-if="!showChannels">Show channels</span>
            <span v-else>Hide channels</span>
          </v-tooltip>
        </v-btn-toggle>
        <v-tab-item>
          <BlendTab/>
        </v-tab-item>
        <v-tab-item>
          <TilesView/>
        </v-tab-item>
      </v-tabs>
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
</template>

<script lang="ts">
  import { mainModule } from '@/modules/main';
  import ChannelsView from '@/views/main/experiment/ChannelsView.vue';
  import BlendTab from '@/views/main/experiment/image/blend/BlendTab.vue';
  import TilesView from '@/views/main/experiment/image/tiles/TilesView.vue';
  import SettingsView from '@/views/main/experiment/settings/SettingsView.vue';
  import WorkflowTab from '@/views/main/experiment/workflow/WorkflowTab.vue';
  import { Component, Vue, Watch } from 'vue-property-decorator';

  @Component({
    components: {
      SettingsView,
      ChannelsView,
      WorkflowTab,
      BlendTab,
      TilesView,
    },
  })
  export default class ImageView extends Vue {
    readonly mainContext = mainModule.context(this.$store);

    toggleChannels = 'show';
    tabImageView = 0;

    get showChannels() {
      return this.mainContext.getters.showChannels;
    }

    get showWorkspace() {
      return this.mainContext.getters.showWorkspace;
    }

    get mainClass() {
      if (this.showChannels) {
        return 'md9';
      }
      return 'md12';
    }

    @Watch('toggleChannels')
    onToggleChannels(value: string) {
      this.mainContext.mutations.setShowChannels(value === 'show');
    }
  }
</script>
