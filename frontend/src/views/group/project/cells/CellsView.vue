<template>
  <div class="cells-view">
    <DatasetsView />
    <ResultsView />
  </div>
</template>

<script lang="ts">
import { Component, Vue } from "vue-property-decorator";
import DatasetsView from "@/views/group/project/cells/DatasetsView.vue";
import ResultsView from "@/views/group/project/cells/ResultsView.vue";

@Component({
  components: { ResultsView, DatasetsView },
})
export default class CellsView extends Vue {}
</script>

<style scoped>
.cells-view {
  display: grid;
  grid-template-columns: auto;
  grid-template-rows: 50% 50%;
  width: 100%;
  height: 100%;
}
</style>
