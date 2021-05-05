<template>
  <div class="root">
    <div>
      <BrushableHistogram
        v-for="channel in selectedChannels"
        :key="channel.acquisition_id + channel.name + containerWidth"
        :channel="channel"
        :containerWidth="containerWidth"
      />
    </div>
  </div>
</template>

<script lang="ts">
import { projectsModule } from "@/modules/projects";
import { Component, Vue } from "vue-property-decorator";
import BrushableHistogram from "@/components/BrushableHistogram.vue";

@Component({
  components: { BrushableHistogram },
})
export default class HistogramView extends Vue {
  readonly projectsContext = projectsModule.context(this.$store);

  containerWidth = 340;

  get selectedChannels() {
    return this.projectsContext.getters.selectedChannels;
  }

  refresh(containerWidth: number) {
    this.containerWidth = containerWidth;
  }
}
</script>

<style scoped>
.root {
  width: 100%;
  height: 100%;
  overflow-y: auto;
}
</style>
