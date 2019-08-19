<template>
  <v-flex>
    <!-- You need to add “cwl-workflow” class to the SVG root for cwl-svg rendering -->
    <svg
      id="cwl-workflow-svg"
      class="cwl-workflow"
    ></svg>
  </v-flex>
</template>

<script lang="ts">
  import { experimentModule } from '@/modules/experiment';
  import { settingsModule } from '@/modules/settings';
  import {
    DeletionPlugin,
    SelectionPlugin,
    SVGArrangePlugin,
    SVGEdgeHoverPlugin,
    SVGNodeMovePlugin,
    SVGPortDragPlugin,
    SVGValidatePlugin,
    Workflow,
    ZoomPlugin,
  } from 'cwl-svg';

  import 'cwl-svg/src/assets/styles/theme.scss';
  import 'cwl-svg/src/plugins/port-drag/theme.scss';
  import 'cwl-svg/src/plugins/selection/theme.scss';
  import { WorkflowFactory } from 'cwlts/models';
  import { Component, Vue } from 'vue-property-decorator';
  import sample from './rna-seq-alignment.json';

  @Component
  export default class WorkflowView extends Vue {
    readonly experimentContext = experimentModule.context(this.$store);
    readonly settingsContext = settingsModule.context(this.$store);

    mounted() {
      this.showWorkflow();
    }

    private showWorkflow() {
      const wf = WorkflowFactory.from(sample);
      const svgRoot = document.getElementById('cwl-workflow-svg') as any;
      const workflow = new Workflow({
        model: wf,
        svgRoot: svgRoot,
        editingEnabled: true,
        plugins: [
          new SVGArrangePlugin(),
          new SVGEdgeHoverPlugin(),
          new SVGNodeMovePlugin(),
          new SVGPortDragPlugin(),
          new SelectionPlugin(),
          new ZoomPlugin(),
          new DeletionPlugin(),
          new SVGValidatePlugin(),
        ],
      });

      // workflow.getPlugin(SVGArrangePlugin).arrange();
      window['wf'] = workflow;
    }
  }
</script>

<style scoped>
  #cwl-workflow-svg {
    height: 100%;
  }
</style>
