<template>
  <v-card flat>
    <v-card-title class="widget-title">
      <v-icon left>mdi-drag</v-icon>
      <span class="subtitle-1 font-weight-light">{{ title }}</span>
    </v-card-title>
    <v-divider />
    <v-img :src="src" contain max-height="550px">
      <template v-slot:placeholder>
        <v-row class="fill-height ma-0" align="center" justify="center">
          <v-progress-circular indeterminate color="grey lighten-5" />
        </v-row>
      </template>
    </v-img>
  </v-card>
</template>

<script lang="ts">
import { Component, Prop, Vue } from "vue-property-decorator";

@Component
export default class PlotWidget extends Vue {
  @Prop(String) readonly title!: string;
  @Prop(String) readonly src!: string;
}
</script>
