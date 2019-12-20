<template>
  <svg ref="svg" height="100" :width="histogramWidth" shape-rendering="optimizeSpeed" />
</template>

<script lang="ts">
import { experimentModule } from "@/modules/experiment";
import { IChannel, IChannelStats } from "@/modules/experiment/models";
import * as d3 from "d3";
import { Component, Prop, Vue } from "vue-property-decorator";

@Component
export default class ChannelHistogramView extends Vue {
  readonly experimentContext = experimentModule.context(this.$store);

  @Prop(Object) channel!: IChannel;

  stats: IChannelStats | undefined = undefined;

  get histogramWidth() {
    const element = document.getElementById("settings-container");
    return element && element.clientWidth - 48;
  }

  showHistogram() {
    if (!this.stats) {
      return;
    }

    const data = this.stats.hist.map((hist, index) => {
      return [this.stats!.edges[index], hist];
    });
    const xMax = Math.max(...this.stats.edges);
    const yMax = Math.max(...this.stats.hist.slice(1));

    // find svg container
    let svg = d3.select(this.$refs.svg as any);
    if (svg.empty()) {
      return;
    }

    // set the dimensions and margins of the graph
    const margin = { top: 10, right: 10, bottom: 20, left: 30 };
    const width = +svg.attr("width") - margin.left - margin.right;
    const height = +svg.attr("height") - margin.top - margin.bottom;

    // append the svg object to the body of the page
    svg = svg
      .attr("width", width + margin.left + margin.right)
      .attr("height", height + margin.top + margin.bottom)
      .append("g")
      .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

    // Add X axis
    const x = d3
      .scaleLinear()
      .domain([0, xMax])
      .range([0, width]);
    const xAxis = svg
      .append("g")
      .attr("transform", "translate(0," + height + ")")
      .call(d3.axisBottom(x));

    // Add Y axis
    const y = d3
      .scaleLinear()
      .domain([0, yMax])
      .range([height, 0]);
    svg.append("g").call(
      d3
        .axisLeft(y)
        .ticks(6)
        .tickFormat(d3.format(".0s"))
    );

    // Add a clipPath: everything out of this area won't be drawn.
    const clip = svg
      .append("defs")
      .append("svg:clipPath")
      .attr("id", "clip")
      .append("svg:rect")
      .attr("width", width)
      .attr("height", height)
      .attr("x", 0)
      .attr("y", 0);

    // Create the scatter variable: where both the circles and the brush take place
    const scatter = svg.append("g").attr("clip-path", "url(#clip)");

    scatter
      .selectAll("line")
      .data(data)
      .enter()
      .append("line")
      .attr("x1", (d: any) => {
        return x(d[0]);
      })
      .attr("x2", (d: any) => {
        return x(d[0]);
      })
      .attr("y1", (d: any) => {
        return y(d[1]);
      })
      .attr("y2", y(0))
      .attr("stroke", "#1975d2");
  }

  async mounted() {
    this.stats = await this.experimentContext.actions.getChannelStats(this.channel.id);
    this.showHistogram();
  }
}
</script>
