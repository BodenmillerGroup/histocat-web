<template>
  <div class="segmentation-view">
    <AcquisitionsView />
    <PanelView />
    <SegmentationSettingsView />
  </div>
</template>

<script lang="ts">
import { Component, Vue } from "vue-property-decorator";
import AcquisitionsView from "@/views/group/project/segmentation/AcquisitionsView.vue";
import PanelView from "@/views/group/project/segmentation/PanelView.vue";
import SegmentationSettingsView from "@/views/group/project/segmentation/SegmentationSettingsView.vue";

@Component({
  components: {
    SegmentationSettingsView,
    PanelView,
    AcquisitionsView,
  },
})
export default class SegmentationView extends Vue {}
</script>

<style scoped>
.segmentation-view {
  display: grid;
  grid-template-columns: 30% 1fr 250px;
  grid-template-rows: auto;
  overflow-y: auto;
  width: 100%;
  height: 100%;
}
</style>
