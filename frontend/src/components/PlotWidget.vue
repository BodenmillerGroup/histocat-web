<template>
  <v-card tile class="widget-container">
    <v-card-title class="widget-title">
      <v-icon left>mdi-drag</v-icon>
      <span class="subtitle-1 font-weight-light">{{ title }}</span>
    </v-card-title>
    <v-divider />
    <v-img :src="src" aspect-ratio="1" class="grey lighten-2" eager>
      <template v-slot:placeholder>
        <v-row class="fill-height ma-0" align="center" justify="center">
          <v-progress-circular indeterminate color="grey lighten-5" />
        </v-row>
      </template>
    </v-img>
  </v-card>
</template>

<script lang="ts">
import {Component, Prop, Vue} from "vue-property-decorator";

@Component
export default class PlotWidget extends Vue {
  @Prop(String) title: string = "";
  @Prop(String) src: string = "";
}
</script>

<style scoped>
.widget-container {
  width: 100%;
  height: 100%;
}
.plot {
  width: 100%;
  height: calc(100% - 36px);
}
</style>
