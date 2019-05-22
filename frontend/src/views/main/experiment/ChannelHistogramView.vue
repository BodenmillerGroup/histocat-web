<template>
  <div ref="svg" height="100" width="300"></div>
</template>

<script lang="ts">
  import { Component, Prop, Vue } from 'vue-property-decorator';
  import * as d3 from 'd3';
  import { IChannel } from '@/modules/experiment/models';
  import { dispatchGetChannelStats } from '@/modules/experiment/actions';
  import { commitSetChannelLevels } from '@/modules/experiment/mutations';

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
      const yMax = Math.max(...stats.hist.slice(1));

      // find svg container
      const div = d3.select(this.$refs.svg as any);
      if (div.empty()) {
        return;
      }

      // set the dimensions and margins of the graph
      const margin = { top: 10, right: 10, bottom: 20, left: 50 };
      const width = +div.attr('width') - margin.left - margin.right;
      const height = +div.attr('height') - margin.top - margin.bottom;

      // append the svg object to the body of the page
      const svg = div
        .append('svg')
        .attr('width', width + margin.left + margin.right)
        .attr('height', height + margin.top + margin.bottom)
        .append('g')
        .attr('transform', 'translate(' + margin.left + ',' + margin.top + ')');

      // Add X axis
      const x = d3.scaleLinear()
        .domain([0, xMax])
        .range([0, width]);
      const xAxis = svg.append('g')
        .attr('transform', 'translate(0,' + height + ')')
        .call(d3.axisBottom(x));

      // Add Y axis
      const y = d3.scaleLinear()
        .domain([0, yMax])
        .range([height, 0]);
      svg.append('g')
        .call(d3.axisLeft(y)
          .ticks(7)
          .tickFormat(d3.format('.0s')));

      // Add a clipPath: everything out of this area won't be drawn.
      const clip = svg.append('defs').append('svg:clipPath')
        .attr('id', 'clip')
        .append('svg:rect')
        .attr('width', width)
        .attr('height', height)
        .attr('x', 0)
        .attr('y', 0);

      // A function that update the chart for given boundaries
      const updateChart = () => {

        const extent = d3.event.selection;

        // If no selection, back to initial coordinate. Otherwise, update X axis domain
        if (!extent) {
          if (!idleTimeout) {
            return idleTimeout = setTimeout(idled, 350);
          } // This allows to wait a little bit
          x.domain([0, xMax]);
          commitSetChannelLevels(this.$store, { id: this.channel.id, levels: undefined });
        } else {
          const min = x.invert(extent[0]);
          const max = x.invert(extent[1]);
          commitSetChannelLevels(this.$store, {
            id: this.channel.id,
            levels: { min: Math.round(min), max: Math.round(max) },
          });
          x.domain([min, max]);
          scatter.select('.brush').call(brush.move as any, null); // This remove the grey brush area as soon as the selection has been done
        }

        // Update axis and circle position
        xAxis.transition().duration(1000).call(d3.axisBottom(x));
        scatter
          .selectAll('line')
          .transition().duration(1000)
          .attr('x1', (d) => {
            return x(d[0]);
          })
          .attr('x2', (d) => {
            return x(d[0]);
          });
      };

      // Add brushing
      const brush = d3.brushX() // Add the brush feature using the d3.brush function
        .extent([[0, 0], [width, height]]) // initialise the brush area: start at 0,0 and finishes at width,height: it means I select the whole graph area
        .on('end', updateChart); // Each time the brush selection changes, trigger the 'updateChart' function

      // Create the scatter variable: where both the circles and the brush take place
      const scatter = svg.append('g')
        .attr('clip-path', 'url(#clip)');

      // Add circles
      scatter
        .selectAll('line')
        .data(data)
        .enter()
        .append('line')
        .attr('x1', (d) => {
          return x(d[0]);
        })
        .attr('x2', (d) => {
          return x(d[0]);
        })
        .attr('y1', (d) => {
          return y(d[1]);
        })
        .attr('y2', y(0))
        .attr('stroke', '#1975d2');

      // Add the brushing
      scatter
        .append('g')
        .attr('class', 'brush')
        .call(brush);

      // A function that set idleTimeOut to null
      let idleTimeout;

      function idled() {
        idleTimeout = null;
      }
    }
  }
</script>
