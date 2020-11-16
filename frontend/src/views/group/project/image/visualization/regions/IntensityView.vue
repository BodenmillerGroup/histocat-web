<template>
  <v-col id="intensity-view-container" class="shrink">
    <IntensityBar
      v-for="channel in selectedChannels"
      :key="channel.id"
      :channel="channel"
      class="mb-1 mx-1"
    ></IntensityBar>
  </v-col>
</template>

<script lang="ts">
import IntensityBar from "@/views/group/project/image/visualization/regions/IntensityBar.vue";
import { projectsModule } from "@/modules/projects";
import { Component, Vue } from "vue-property-decorator";

@Component({
  components: { IntensityBar },
})
export default class IntensityView extends Vue {
  readonly projectsContext = projectsModule.context(this.$store);

  get selectedChannels() {
    return this.projectsContext.getters.selectedChannels;
  }
}
</script>
