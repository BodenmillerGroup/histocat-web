<template>
  <v-expansion-panel>
    <v-expansion-panel-header v-slot="{ open }">
      <v-layout column>
        <v-layout justify-space-between>
          <b>{{ channel.label }}</b>
          <v-btn-toggle v-model="colorIndex">
            <v-btn color="red" x-small @click.stop="">R</v-btn>
            <v-btn color="green" x-small @click.stop="">G</v-btn>
            <v-btn color="blue" x-small @click.stop="">B</v-btn>
            <v-btn color="yellow" x-small @click.stop="">Y</v-btn>
            <v-btn color="cyan" x-small @click.stop="">C</v-btn>
            <v-btn color="#FF00FF" x-small @click.stop="">M</v-btn>
          </v-btn-toggle>
        </v-layout>
        <v-range-slider
          :value="levels"
          :max="channel.max_intensity"
          :min="channel.min_intensity"
          :step="1"
          thumb-label
          :thumb-size="24"
          @click.stop=""
          @end="submitLimit"
          hide-details
        ></v-range-slider>
      </v-layout>
    </v-expansion-panel-header>
    <v-expansion-panel-content>
      <ChannelHistogramView :channel="channel"></ChannelHistogramView>
    </v-expansion-panel-content>
  </v-expansion-panel>
</template>

<script lang="ts">
  import { experimentModule } from '@/modules/experiment';
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
    readonly experimentContext = experimentModule.context(this.$store);

    @Prop(Object) channel!: IChannel;

    colorIndex = this.channel ? convertColorToIndex(this.metalColor) : -1;

    get levels() {
      const settings = this.settingsContext.getters.channelSettings(this.channel.id);
      if (settings && settings.levels) {
        return [settings.levels.min, settings.levels.max];
      } else {
        return [this.channel.min_intensity, this.channel.max_intensity];
      }
    }

    submitLimit(range: number[]) {
      this.settingsContext.mutations.setChannelSettings({
        id: this.channel.id,
        levels: { min: Math.round(range[0]), max: Math.round(range[1]) },
      });
      this.experimentContext.actions.getChannelStackImage();
    }

    get metalColor() {
      const colorMap = this.settingsContext.getters.metalColorMap;
      const color = colorMap.get(this.channel.metal);
      return color ? color : '';
    }

    @Watch('colorIndex')
    onColorIndexChanged(colorIndex: number) {
      this.settingsContext.mutations.setMetalColor({
        metal: this.channel.metal,
        color: convertIndexToColor(colorIndex),
      });
    }
  }
</script>
