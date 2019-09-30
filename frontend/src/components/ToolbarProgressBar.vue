<template>
  <v-progress-linear
    :active="processing"
    :color="color"
    :indeterminate="indeterminate"
    :value="progress"
    class="ma-0"
    height="4"
    style="top: -2px;"
  />
</template>

<script lang="ts">
import { Component, Prop, Vue } from "vue-property-decorator";

@Component
export default class ToolbarProgressBar extends Vue {
  @Prop(Boolean) processing?: boolean;
  @Prop(Number) progress?: boolean;
  @Prop(Boolean) indeterminate?: boolean;
  @Prop(String) color?: string;
}
</script>
