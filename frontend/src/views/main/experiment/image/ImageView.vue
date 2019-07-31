<template>
  <v-flex :class="viewerClass">
    <v-tabs v-model="tabImageView">
      <v-tab>Blend</v-tab>
      <v-tab>Tiles</v-tab>
      <v-tab>Workflow</v-tab>
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
  import BlendTab from '@/views/main/experiment/image/blend/BlendTab.vue';
  import TilesView from '@/views/main/experiment/image/tiles/TilesView.vue';
  import WorkflowTab from '@/views/main/experiment/workflow/WorkflowTab.vue';
  import { Component, Vue } from 'vue-property-decorator';

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
