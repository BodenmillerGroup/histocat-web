<template>
  <div class="right-sidebar">
    <div v-for="item in items" :key="item">
      {{ item }}
    </div>
  </div>
</template>

<script lang="ts">
import { Component, Prop, Vue } from "vue-property-decorator";

@Component
export default class RightSidebar extends Vue {
  get items() {
    return new Array(100).fill(0).map((v, i) => i)
  }
}
</script>

<style scoped>
.right-sidebar {
  /* x y blur spread color */
  border-left: 1px solid rgb(153, 153, 153);
  display: flex;
  flex-direction: column;
  position: relative;
  overflow-y: inherit;
  height: inherit;
  width: inherit;
}
</style>
