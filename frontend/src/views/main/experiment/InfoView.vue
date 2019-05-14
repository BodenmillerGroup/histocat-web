<template>
  <v-card tile>
    <v-card-title><h4>Info</h4></v-card-title>
    <v-divider></v-divider>
    <v-list dense class="scroll-y local-height">
      <v-list-tile v-for="item in items" :key="item.name">
        <v-list-tile-content>{{item.name}}</v-list-tile-content>
        <v-list-tile-content class="align-end">{{item.value}}</v-list-tile-content>
      </v-list-tile>
    </v-list>
  </v-card>
</template>

<script lang="ts">
  import { Component, Vue } from 'vue-property-decorator';
  import { readSelectedMeta } from '@/modules/experiment/getters';

  @Component
  export default class InfoView extends Vue {

    get meta() {
      return readSelectedMeta(this.$store);
    }

    get items() {
      if (this.meta) {
        const items = Object.entries(this.meta);
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
  .local-height {
    max-height: 46vh;
  }
</style>
