<template>
  <v-card tile class="card scroll-y">
    <v-card-title><h4>Info</h4></v-card-title>
    <v-divider></v-divider>
    <v-list dense>
      <v-list-tile v-for="item in items" :key="item.name" class="info-list-tile">
        <v-list-tile-content>{{item.name}}</v-list-tile-content>
        <v-list-tile-content class="align-end info-list-tile-content">{{item.value}}</v-list-tile-content>
      </v-list-tile>
    </v-list>
  </v-card>
</template>

<script lang="ts">
  import { Component, Prop, Vue } from 'vue-property-decorator';

  @Component
  export default class InfoCard extends Vue {
    @Prop(Object) node;

    get items() {
      if (this.node.item.meta) {
        const items = Object.entries(this.node.item.meta);
        return items.map((item) => {
          return {
            name: item[0],
            value: item[1],
          };
        });
      }
    }
  }
</script>

<style scoped>
  .card {
    height: 400px;
  }

  .info-list-tile {
    height: 25px;
  }

  .info-list-tile-content {
    margin-left: 10px;
  }
</style>