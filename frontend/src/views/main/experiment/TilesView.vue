<template>
  <v-card flat tile>
    <v-toolbar card dense>
      <v-toolbar-side-icon></v-toolbar-side-icon>
      <v-toolbar-title>Tiles</v-toolbar-title>
      <v-spacer/>
      <v-btn icon>
        <v-icon>search</v-icon>
      </v-btn>
      <v-btn icon>
        <v-icon>favorite</v-icon>
      </v-btn>
      <v-btn icon>
        <v-icon>more_vert</v-icon>
      </v-btn>
    </v-toolbar>
  </v-card>
</template>

<script lang="ts">
  import { Component, Vue } from 'vue-property-decorator';

  @Component
  export default class TilesView extends Vue {
  }
</script>
