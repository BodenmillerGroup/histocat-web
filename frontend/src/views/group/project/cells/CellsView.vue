<template>
  <div>
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
