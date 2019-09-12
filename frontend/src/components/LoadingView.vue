<template>
  <v-container fluid>
    <v-row>
      <v-col align="center" justify="center">
        <v-progress-circular indeterminate color="primary" />
        <h5>{{ text }}</h5>
      </v-col>
    </v-row>
  </v-container>
</template>

<script lang="ts">
import { Component, Prop, Vue } from "vue-property-decorator";

@Component
export default class LoadingView extends Vue {
  @Prop(String) text: string | undefined;
}
</script>
