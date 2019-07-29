<template>
  <v-flex column>
    <v-toolbar dense flat>
      <v-menu offset-y>
        <template v-slot:activator="{ on }">
          <v-btn
            icon
            v-on="on"
          >
            <v-tooltip bottom>
              <template v-slot:activator="{ on }">
                <v-btn icon v-on="on">
                  <v-icon>mdi-cloud-download-outline</v-icon>
                </v-btn>
              </template>
              <span>Export image</span>
            </v-tooltip>
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
    </v-toolbar>
    <WorkflowView class="workflow-view"/>
  </v-flex>
</template>

<script lang="ts">
  import { experimentModule } from '@/modules/experiment';
  import { mainModule } from '@/modules/main';
  import WorkflowView from '@/views/main/experiment/workflow/WorkflowView.vue';
  import { Component, Vue } from 'vue-property-decorator';

  @Component({
    components: { WorkflowView },
  })
  export default class WorkflowTab extends Vue {
    readonly mainContext = mainModule.context(this.$store);
    readonly experimentContext = experimentModule.context(this.$store);

    get showWorkspace() {
      return this.mainContext.getters.showWorkspace;
    }

    get showChannels() {
      return this.mainContext.getters.showChannels;
    }
  }
</script>

<style scoped>
  .workflow-view {
    height: calc(100vh - 162px);
  }
</style>
