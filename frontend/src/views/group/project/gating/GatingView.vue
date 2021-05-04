<template>
  <div>
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
