<template>
  <v-img :src="src" aspect-ratio="1" class="grey lighten-2" eager>
    <template v-slot:placeholder>
      <v-row class="fill-height ma-0" align="center" justify="center">
        <v-progress-circular indeterminate color="grey lighten-5" />
      </v-row>
    </template>
  </v-img>
</template>

<script lang="ts">
import { Component, Prop, Vue } from "vue-property-decorator";

@Component
export default class PlotView extends Vue {
  @Prop(String) title: string = "";
  @Prop(String) src: string = "";
}
</script>
