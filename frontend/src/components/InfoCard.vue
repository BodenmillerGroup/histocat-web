<template>
  <v-card tile flat class="card overflow-y-auto">
    <v-card-title>Info</v-card-title>
    <v-divider />
    <v-card-title v-if="imageUrl">
      <v-img :src="imageUrl" class="grey lighten-2" width="380px">
        <template v-slot:placeholder>
          <v-row class="fill-height ma-0" align="center" justify="center">
            <v-progress-circular indeterminate color="grey lighten-5" />
          </v-row>
        </template>
      </v-img>
    </v-card-title>
    <v-list dense>
      <template v-for="item in items">
        <v-list-item dense inactive :key="item.name">
          <v-list-item-content>
            <v-list-item-title v-text="item.name" />
          </v-list-item-content>
          <v-list-item-action>
            <v-list-item-action-text v-text="item.value" />
          </v-list-item-action>
        </v-list-item>
      </template>
    </v-list>
  </v-card>
</template>

<script lang="ts">
import { apiUrl } from "@/env";
import { Component, Prop, Vue } from "vue-property-decorator";

@Component
export default class InfoCard extends Vue {
  @Prop(Object) readonly node;

  get imageUrl() {
    switch (this.node.item.type) {
      case "slide":
        return `${apiUrl}/slides/${this.node.item.id}/image`;
      case "panorama":
        return `${apiUrl}/panoramas/${this.node.item.id}/image`;
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
    } else {
      return [];
    }
  }
}
</script>

<style scoped>
.card {
  height: 40vh;
}
</style>
