<template>
  <v-card tile>
    <v-card-text>
      <v-btn-toggle v-model="color">
        <v-btn color="red" small depressed>R</v-btn>
        <v-btn color="green" small depressed>G</v-btn>
        <v-btn color="blue" small depressed>B</v-btn>
        <v-btn color="cyan" small depressed>C</v-btn>
        <v-btn color="#FF00FF" small depressed>M</v-btn>
        <v-btn color="yellow" small depressed>Y</v-btn>
      </v-btn-toggle>
    </v-card-text>
    <v-card-text>
      <ChannelHistogramView></ChannelHistogramView>
    </v-card-text>
  </v-card>
</template>

<script lang="ts">
  import { Component, Vue } from 'vue-property-decorator';
  import ChannelHistogramView from '@/components/ChannelHistogramView.vue';
  @Component({
    components: { ChannelHistogramView },
  })
  export default class ChannelSettingsView extends Vue {
    color = '';
  }
</script>
