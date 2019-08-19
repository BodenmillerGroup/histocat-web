<template>
  <v-flex
    column
    shrink
    id="intensity-view-container"
  >
    <IntensityBar
      v-for="channel in selectedChannels"
      :key="channel.id"
      :channel="channel"
      class="mb-1"
    ></IntensityBar>
  </v-flex>
</template>

<script lang="ts">
  import IntensityBar from '@/views/main/experiment/image/blend/IntensityBar.vue';
  import { experimentModule } from '@/modules/experiment';
  import { Component, Vue } from 'vue-property-decorator';

  @Component({
    components: { IntensityBar },
  })
  export default class IntensityView extends Vue {
    readonly experimentContext = experimentModule.context(this.$store);

    get selectedChannels() {
      return this.experimentContext.getters.selectedChannels;
    }
  }
</script>
