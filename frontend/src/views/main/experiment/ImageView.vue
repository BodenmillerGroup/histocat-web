<template>
  <v-tabs v-model="active">
    <v-tab ripple>
      Blend
    </v-tab>
    <v-tab ripple>
      Tiles
    </v-tab>
    <v-spacer></v-spacer>
    <v-btn icon>
      <v-icon>mdi-magnify</v-icon>
    </v-btn>
    <v-tab-item>
      <BlendView/>
    </v-tab-item>
    <v-tab-item>
      <TilesView/>
    </v-tab-item>
  </v-tabs>
</template>

<script lang="ts">
  import { Component, Vue } from 'vue-property-decorator';
  import BlendView from '@/views/main/experiment/BlendView.vue';
  import TilesView from '@/views/main/experiment/TilesView.vue';

  @Component({
    components: { TilesView, BlendView },
  })
  export default class ImageView extends Vue {
  }
</script>
