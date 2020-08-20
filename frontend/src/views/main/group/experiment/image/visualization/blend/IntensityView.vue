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
import IntensityBar from "@/views/main/group/experiment/image/visualization/blend/IntensityBar.vue";
import { experimentModule } from "@/modules/experiment";
import { Component, Vue } from "vue-property-decorator";

@Component({
  components: { IntensityBar },
})
export default class IntensityView extends Vue {
  readonly experimentContext = experimentModule.context(this.$store);

  get selectedChannels() {
    return this.experimentContext.getters.selectedChannels;
  }
}
</script>
