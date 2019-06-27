<template>
  <v-flex column pa-3>
    <v-layout justify-space-between>
      <b>{{ channel.label }}</b>
      <v-btn-toggle v-model="colorIndex">
        <v-btn color="red" small depressed>R</v-btn>
        <v-btn color="green" small depressed>G</v-btn>
        <v-btn color="blue" small depressed>B</v-btn>
        <v-btn color="yellow" small depressed>Y</v-btn>
        <v-btn color="cyan" small depressed>C</v-btn>
        <v-btn color="#FF00FF" small depressed>M</v-btn>
      </v-btn-toggle>
    </v-layout>
    <ChannelHistogramView :channel="channel"></ChannelHistogramView>
  </v-flex>
</template>

<script lang="ts">
  import { Component, Prop, Vue, Watch } from 'vue-property-decorator';
  import ChannelHistogramView from '@/views/main/experiment/ChannelHistogramView.vue';
  import { IChannel } from '@/modules/experiment/models';
  import { readMetalColorMap } from '@/modules/experiment/getters';
  import { commitSetMetalColor } from '@/modules/experiment/mutations';
  import { convertColorToIndex, convertIndexToColor } from '@/utils';

  @Component({
    components: { ChannelHistogramView },
  })
  export default class ChannelSettingsView extends Vue {

    @Prop(Object) channel!: IChannel;

    get metalColor() {
      const colorMap = readMetalColorMap(this.$store);
      return colorMap[this.channel.metal];
    }

    colorIndex = this.channel ? convertColorToIndex(this.metalColor) : -1;

    @Watch('colorIndex')
    onColorIndexChanged(colorIndex: number) {
      commitSetMetalColor(this.$store, { metal: this.channel.metal, color: convertIndexToColor(colorIndex) });
    }
  }
</script>
