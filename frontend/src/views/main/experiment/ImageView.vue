<template>
  <v-flex :class="viewerClass">
    <v-tabs v-model="tabImageView">
      <v-tab>Blend</v-tab>
      <v-tab>Tiles</v-tab>
      <v-tab>Workflow</v-tab>
      <v-spacer/>
      <v-btn-toggle v-model="toggleUI" multiple>
        <v-tooltip bottom>
          <template v-slot:activator="{ on }">
            <v-btn icon v-on="on" value='workspace'>
              <v-icon>mdi-file-tree</v-icon>
            </v-btn>
          </template>
          <span>Show workspace</span>
        </v-tooltip>
        <v-tooltip bottom>
          <template v-slot:activator="{ on }">
            <v-btn icon v-on="on" value='channels'>
              <v-icon>mdi-format-list-checkbox</v-icon>
            </v-btn>
          </template>
          <span>Show channels</span>
        </v-tooltip>
      </v-btn-toggle>
      <v-tab-item>
        <BlendTab/>
      </v-tab-item>
      <v-tab-item>
        <TilesView/>
      </v-tab-item>
      <v-tab-item>
        <WorkflowTab/>
      </v-tab-item>
    </v-tabs>
  </v-flex>
</template>

<script lang="ts">
  import { mainModule } from '@/modules/main';
  import BlendTab from '@/views/main/experiment/BlendTab.vue';
  import TilesView from '@/views/main/experiment/TilesView.vue';
  import WorkflowTab from '@/views/main/experiment/workflow/WorkflowTab.vue';
  import { Component, Vue, Watch } from 'vue-property-decorator';

  @Component({
    components: {
      WorkflowTab,
      BlendTab,
      TilesView,
    },
  })
  export default class ImageView extends Vue {
    readonly mainContext = mainModule.context(this.$store);

    tabImageView = 0;
    toggleUI = ['workspace', 'channels'];

    @Watch('toggleUI')
    onToggleMultiple(items: string[]) {
      const showWorkspace = items.includes('workspace');
      const showChannels = items.includes('channels');
      this.mainContext.mutations.setShowWorkspace(showWorkspace);
      this.mainContext.mutations.setShowChannels(showChannels);
    }

    get viewerClass() {
      const showWorkspace = this.mainContext.getters.showWorkspace;
      const showChannels = this.mainContext.getters.showChannels;
      if (showWorkspace && showChannels) {
        return 'md6';
      }
      if (showWorkspace || showChannels) {
        return 'md9';
      }
      return 'md12';
    }
  }
</script>
