<template>
  <v-layout row>
    <v-flex :class="mainClass">
      <v-tabs v-model="tab">
        <v-tab>Blend</v-tab>
        <v-tab>Tiles</v-tab>
        <v-tab-item>
          <BlendTab/>
        </v-tab-item>
        <v-tab-item>
          <TilesView/>
        </v-tab-item>
      </v-tabs>
    </v-flex>
    <v-flex v-show="showOptions" md3>
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
  import { Component, Vue } from 'vue-property-decorator';

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

    tab = 0;

    get showOptions() {
      return this.mainContext.getters.showOptions;
    }

    get mainClass() {
      if (this.showOptions) {
        return 'md9';
      }
      return 'md12';
    }
  }
</script>
