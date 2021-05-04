<template>
  <div id="settings-container">
    <div class="overflow-y-auto channel-settings-view">
      <BrushableHistogram v-for="channel in selectedChannels" :key="channel.name" :channel="channel" />
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

  get selectedChannels() {
    return this.projectsContext.getters.selectedChannels;
  }
}
</script>

<style scoped>
.channel-settings-view {
  height: calc(50vh - 82px);
}
</style>
