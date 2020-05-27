<template>
  <DatasetsView />
</template>

<script lang="ts">
import DatasetsView from "@/views/main/experiment/data/workspace/dataset/DatasetsView.vue";
import { Component, Prop, Vue } from "vue-property-decorator";
import { IExperiment } from "@/modules/experiment/models";

@Component({
  components: {
    DatasetsView,
  },
})
export default class WorkspaceView extends Vue {
  @Prop(Object) readonly experiment!: IExperiment;
}
</script>
