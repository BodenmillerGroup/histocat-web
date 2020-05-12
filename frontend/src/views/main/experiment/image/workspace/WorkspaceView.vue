<template>
  <SlidesTreeView :experiment="experiment" />
</template>

<script lang="ts">
import { IExperiment } from "@/modules/experiment/models";
import SlidesTreeView from "@/views/main/experiment/image/workspace/slide/SlidesTreeView.vue";
import { Component, Prop, Vue } from "vue-property-decorator";

@Component({
  components: {
    SlidesTreeView,
  },
})
export default class WorkspaceView extends Vue {
  @Prop(Object) readonly experiment!: IExperiment;
}
</script>
