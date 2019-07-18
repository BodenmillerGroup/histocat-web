<template>
  <v-flex :class="viewerClass">
    <v-toolbar flat dense>
      <CreateDatasetDialog/>
      <v-menu offset-y>
        <template v-slot:activator="{ on }">
          <v-btn
            small
            v-on="on"
          >
            Download
          </v-btn>
        </template>
        <v-list dense>
          <v-list-item
            @click="download('tiff')"
          >
            <v-list-item-title>TIFF</v-list-item-title>
          </v-list-item>
          <v-list-item
            @click="download('png')"
          >
            <v-list-item-title>PNG</v-list-item-title>
          </v-list-item>
        </v-list>
      </v-menu>
      <v-spacer/>
      <v-btn-toggle v-model="toggleUI" multiple>
        <v-btn text tile small value='workspace'>
          <v-icon left>mdi-file-tree</v-icon> Workspace
        </v-btn>
        <v-btn text tile small value='channels'>
          <v-icon left>mdi-format-list-checkbox</v-icon> Channels
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
</template>

<script lang="ts">
  import { experimentModule } from '@/modules/experiment';
  import { mainModule } from '@/modules/main';
  import BlendView from '@/views/main/experiment/BlendView.vue';
  import CreateDatasetDialog from '@/views/main/experiment/CreateDatasetDialog.vue';
  import TilesView from '@/views/main/experiment/TilesView.vue';
  import { Component, Vue, Watch } from 'vue-property-decorator';

  @Component({
    components: {
      TilesView,
      CreateDatasetDialog,
      BlendView,
    },
  })
  export default class ImageView extends Vue {
    readonly mainContext = mainModule.context(this.$store);
    readonly experimentContext = experimentModule.context(this.$store);

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

    download(type: 'png' | 'tiff') {
      this.experimentContext.actions.exportChannelStackImage(type);
    }
  }
</script>

<style scoped>
  .image-view {
    height: calc(100vh - 168px);
  }
</style>