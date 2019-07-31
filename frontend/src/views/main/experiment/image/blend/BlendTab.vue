<template>
  <v-flex column>
    <v-toolbar dense flat>
      <v-menu offset-y>
        <template v-slot:activator="{ on }">
          <v-btn
            text
            small
            v-on="on"
          >
            <v-icon left>mdi-cloud-download-outline</v-icon>
            Export image
          </v-btn>
        </template>
        <v-list dense>
          <v-list-item
            @click="download('tiff')"
          >
            <v-list-item-title>Export TIFF</v-list-item-title>
          </v-list-item>
          <v-list-item
            @click="download('png')"
          >
            <v-list-item-title>Export PNG</v-list-item-title>
          </v-list-item>
        </v-list>
      </v-menu>
    </v-toolbar>
    <v-layout>
      <v-flex pa-0>
        <BlendView class="blend-view"/>
      </v-flex>
      <IntensityView/>
    </v-layout>
  </v-flex>
</template>

<script lang="ts">
  import { experimentModule } from '@/modules/experiment';
  import { ExportTypes } from '@/modules/experiment/models';
  import { mainModule } from '@/modules/main';
  import BlendView from '@/views/main/experiment/image/blend/BlendView.vue';
  import IntensityView from '@/views/main/experiment/image/blend/IntensityView.vue';
  import { Component, Vue } from 'vue-property-decorator';

  @Component({
    components: { IntensityView, BlendView },
  })
  export default class BlendTab extends Vue {
    readonly mainContext = mainModule.context(this.$store);
    readonly experimentContext = experimentModule.context(this.$store);

    get showWorkspace() {
      return this.mainContext.getters.showWorkspace;
    }

    get showChannels() {
      return this.mainContext.getters.showChannels;
    }

    download(type: ExportTypes) {
      this.experimentContext.actions.exportChannelStackImage(type);
    }
  }
</script>

<style scoped>
  .blend-view {
    height: calc(100vh - 162px);
  }
</style>
