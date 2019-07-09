<template>
  <v-card tile flat class="card scroll-y">
    <v-card-title v-if="imageUrl">
      <v-img
        :src="imageUrl"
        class="grey lighten-2"
      >
        <template v-slot:placeholder>
          <v-layout
            fill-height
            align-center
            justify-center
            ma-0
          >
            <v-progress-circular indeterminate color="grey lighten-5"></v-progress-circular>
          </v-layout>
        </template>
      </v-img>
    </v-card-title>
    <v-card-title v-else>
      Info
    </v-card-title>
    <v-divider></v-divider>
    <v-list dense class="list">
      <v-list-tile v-for="item in items" :key="item.name" class="list-tile">
        <v-list-tile-content>{{item.name}}</v-list-tile-content>
        <v-list-tile-content class="align-end list-tile-content">{{item.value}}</v-list-tile-content>
      </v-list-tile>
    </v-list>
  </v-card>
</template>

<script lang="ts">
  import { apiUrl } from '@/env';
  import { Component, Prop, Vue } from 'vue-property-decorator';

  @Component
  export default class InfoCard extends Vue {
    @Prop(Object) node;

    get imageUrl() {
      switch (this.node.item.type) {
        case 'slide':
          return `${apiUrl}/api/v1/slides/${this.node.item.id}/image`;
        case 'panorama':
          return `${apiUrl}/api/v1/panoramas/${this.node.item.id}/image`;
        default:
          return undefined;
      }
    }

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

  .list-tile {
    height: 25px;
  }

  .list-tile-content {
    margin-left: 10px;
  }
</style>