<template>
  <v-card
    tile
    width="60"
    elevation="1"
  >
    <v-card-title class="caption">{{caption}}</v-card-title>
    <v-card-text>
      <svg
        ref="svg"
        class="svg"
        shape-rendering="optimizeSpeed"
      ></svg>
    </v-card-text>
  </v-card>
</template>

<script lang="ts">
  import { IChannel } from '@/modules/experiment/models';
  import { settingsModule } from '@/modules/settings';
  import * as d3 from 'd3';
  import { Component, Prop, Vue, Watch } from 'vue-property-decorator';

  @Component
  export default class IntensityBar extends Vue {
    readonly settingsContext = settingsModule.context(this.$store);

    @Prop(Object) readonly channel!: IChannel;

    get caption() {
      const settings = this.settingsContext.getters.getChannelSettings(this.channel.id);
      return settings && settings.levels ?
        settings.levels.max :
        this.channel.max_intensity.toFixed(0);
    }

    get color() {
      const color = this.settingsContext.getters.metalColorMap.get(this.channel.metal);
      return color ? color : '#ffffff';
    }

    mounted() {
      this.draw();
    }

    @Watch('color')
    private draw() {
      // Find svg container
      const svg = d3.select(this.$refs.svg as any);
      if (svg.empty()) {
        return;
      }

      svg.selectAll('*').remove();

      // Append a defs (for definition) element to your SVG
      const defs = svg.append('defs');

      // Append a linearGradient element to the defs and give it a unique id
      const linearGradient = defs
        .append('linearGradient')
        .attr('id', `linear-gradient-${this.channel.id}`);

      // Vertical gradient
      linearGradient
        .attr('x1', '0%')
        .attr('y1', '0%')
        .attr('x2', '0%')
        .attr('y2', '100%');

      // Set the color for the start (0%)
      linearGradient.append('stop')
        .attr('offset', '0%')
        .attr('stop-color', this.color);

      // Set the color for the end (100%)
      linearGradient.append('stop')
        .attr('offset', '100%')
        .attr('stop-color', 'black');

      // Draw the rectangle and fill with gradient
      svg.append('rect')
        .attr('width', 20)
        .attr('height', 100)
        .style('fill', `url(#linear-gradient-${this.channel.id})`);
    }
  }
</script>

<style scoped>
  .svg {
    height: 100px;
  }
</style>
