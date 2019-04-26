<template>
  <v-card tile>
    <v-card-title><h4>Info</h4></v-card-title>
    <v-divider></v-divider>
    <v-list dense>
      <v-list-tile>
        <v-list-tile-content>Calories:</v-list-tile-content>
        <v-list-tile-content class="align-end">x</v-list-tile-content>
      </v-list-tile>
      <v-list-tile>
        <v-list-tile-content>Fat:</v-list-tile-content>
        <v-list-tile-content class="align-end">xx</v-list-tile-content>
      </v-list-tile>
    </v-list>
  </v-card>
</template>

<script lang="ts">
  import { Component, Vue } from 'vue-property-decorator';

  @Component
  export default class InfoView extends Vue {
  }
</script>
