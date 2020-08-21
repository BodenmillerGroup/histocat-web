<template>
  <v-tabs v-model="tab">
    <v-tab>Datasets</v-tab>
    <v-tab>Gates</v-tab>
    <v-tab-item>
      <DatasetsView />
    </v-tab-item>
    <v-tab-item>
      <GatesView />
    </v-tab-item>
  </v-tabs>
</template>

<script lang="ts">
import { Component, Vue } from "vue-property-decorator";
import DatasetsView from "@/views/group/experiment/data/workspace/DatasetsView.vue";
import GatesView from "@/views/group/experiment/data/workspace/GatesView.vue";

@Component({
  components: { GatesView, DatasetsView },
})
export default class DataWorkspaceView extends Vue {
  tab = 0;
}
</script>
