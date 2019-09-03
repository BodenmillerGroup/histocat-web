<template>
  <v-flex column>
    <v-toolbar dense flat>

    </v-toolbar>
    <WorkflowView class="workflow-view"/>
  </v-flex>
</template>

<script lang="ts">
  import WorkflowView from '@/views/main/experiment/workflow/WorkflowView.vue';
  import { Component, Vue } from 'vue-property-decorator';

  @Component({
    components: { WorkflowView },
  })
  export default class WorkflowTab extends Vue {
  }
</script>

<style scoped>
  .workflow-view {
    height: calc(100vh - 162px);
  }
</style>
