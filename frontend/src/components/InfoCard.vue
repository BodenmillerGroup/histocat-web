<template>
  <v-card tile flat class="card">
    <v-card-title>Info</v-card-title>
    <v-divider></v-divider>
    <v-list dense class="scroll-y list">
      <v-list-tile v-for="item in items" :key="item.name" class="list-tile">
        <v-list-tile-content>{{item.name}}</v-list-tile-content>
        <v-list-tile-content class="align-end list-tile-content">{{item.value}}</v-list-tile-content>
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
    height: 40vh;
  }

  .list {
    height: calc(100% - 55px);
  }

  .list-tile {
    height: 25px;
  }

  .list-tile-content {
    margin-left: 10px;
  }
</style>