<template>
  <div class="layout-full py-0">
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
.layout-full {
  display: grid;
  grid-template-columns: 30% 1fr 20%;
  grid-template-rows: auto;
}
</style>
