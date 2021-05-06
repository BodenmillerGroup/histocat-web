<template>
  <div class="gating-view">
    <GatesView />
    <AnnotationsView />
    <CellClassesView />
  </div>
</template>

<script lang="ts">
import { Component, Vue } from "vue-property-decorator";
import GatesView from "@/views/group/project/gating/GatesView.vue";
import CellClassesView from "@/views/group/project/gating/CellClassesView.vue";
import AnnotationsView from "@/views/group/project/gating/AnnotationsView.vue";

@Component({
  components: { AnnotationsView, CellClassesView, GatesView },
})
export default class GatingView extends Vue {}
</script>

<style scoped>
.gating-view {
  display: grid;
  grid-template-columns: auto;
  grid-template-rows: 33.33% 33.33% 33.33%;
  width: 100%;
  height: 100%;
}
</style>
