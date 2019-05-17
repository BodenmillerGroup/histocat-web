<template>
  <svg ref="svg" width="270" height="100"></svg>
</template>

<script lang="ts">
  import { Component, Prop, Vue } from 'vue-property-decorator';
  import * as d3 from 'd3';
  import { IChannel } from '@/modules/experiment/models';
  import { dispatchGetChannelStats } from '@/modules/experiment/actions';

  @Component
  export default class ChannelHistogramView extends Vue {

    @Prop(Object) channel!: IChannel;

    async mounted() {
      const stats = await dispatchGetChannelStats(this.$store, { id: this.channel.id });
      if (!stats) {
        return;
      }

      const data = stats.hist.map((hist, index) => {
        return [stats.bins[index], hist];
      });
      const xMax = Math.max(...stats.bins);
      const yMax = Math.max(...stats.hist);

      const svg = d3.select(this.$refs.svg as any);
      const margin = { top: 0, right: 10, bottom: 20, left: 10 };
      const width = +svg.attr('width') - margin.left - margin.right;
      const height = +svg.attr('height') - margin.top - margin.bottom;
      const g = svg.append('g').attr('transform', 'translate(' + margin.left + ',' + margin.top + ')');

      const x = d3.scaleLinear()
        .domain([0, xMax])
        .range([0, width]);

      const y = d3.scaleLinear()
        .domain([0, yMax])
        .range([height, 0]);

      const brush = d3.brushX()
        .extent([[0, 0], [width, height]])
        .on('start brush', brushed);

      const dot = g.append('g')
        .attr('fill-opacity', 0.2)
        .selectAll('circle')
        .data(data)
        .enter().append('circle')
        .attr('transform', (d) => {
          return 'translate(' + x(d[0]) + ',' + y(d[1]) + ')';
        })
        .attr('r', 1);

      g.append('g')
        .call(brush)
        .call(brush.move, [0, width].map(x))
        .selectAll('.overlay')
        .each((d: any) => {
          d.type = 'selection';
        }); // Treat overlay interaction as move.

      g.append('g')
        .attr('transform', 'translate(0,' + height + ')')
        .call(d3.axisBottom(x));

      function brushed() {
        const extent = d3.event.selection.map(x.invert, x);
        dot.classed('selected', (d) => {
          return extent[0] <= d[0] && d[0] <= extent[1];
        });
      }
    }
  }
</script>
