<template>
  <v-flex column pa-3>
    <v-layout justify-space-between>
      <b>{{ channel.label }}</b>
      <v-btn-toggle v-model="colorIndex">
        <v-btn color="red" x-small>R</v-btn>
        <v-btn color="green" x-small>G</v-btn>
        <v-btn color="blue" x-small>B</v-btn>
        <v-btn color="yellow" x-small>Y</v-btn>
        <v-btn color="cyan" x-small>C</v-btn>
        <v-btn color="#FF00FF" x-small>M</v-btn>
      </v-btn-toggle>
    </v-layout>
    <ChannelHistogramView :channel="channel"></ChannelHistogramView>
  </v-flex>
</template>

<script lang="ts">
  import { IChannel } from '@/modules/experiment/models';
  import { settingsModule } from '@/modules/settings';
  import { convertColorToIndex, convertIndexToColor } from '@/utils';
  import ChannelHistogramView from '@/views/main/experiment/ChannelHistogramView.vue';
  import { Component, Prop, Vue, Watch } from 'vue-property-decorator';

  @Component({
    components: { ChannelHistogramView },
  })
  export default class ChannelSettingsView extends Vue {
    readonly settingsContext = settingsModule.context(this.$store);

    @Prop(Object) channel!: IChannel;

    get metalColor() {
      const colorMap = this.settingsContext.getters.metalColorMap;
      const color = colorMap.get(this.channel.metal);
      return color ? color : '';
    }

    colorIndex = this.channel ? convertColorToIndex(this.metalColor) : -1;

    @Watch('colorIndex')
    onColorIndexChanged(colorIndex: number) {
      this.settingsContext.mutations.setMetalColor({
        metal: this.channel.metal,
        color: convertIndexToColor(colorIndex),
      });
    }
  }
</script>
