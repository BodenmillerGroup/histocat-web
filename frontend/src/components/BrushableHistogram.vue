<template>
  <div class="root">
    <div class="header">
      <div class="labels">
        <span class="range-label">{{ levels[0] }}</span>
        <span class="label">
          {{ label }}
        </span>
        <span class="range-label">{{ levels[1] }}</span>
      </div>
      <input type="color" v-model.lazy="color" @click.stop />
    </div>
    <svg :width="width" :height="height" id="svg" ref="svg" />
  </div>
</template>

<script lang="ts">
import { Component, Prop, Vue, Watch } from "vue-property-decorator";
import * as d3 from "d3";
import memoize from "memoize-one";
import { settingsModule } from "@/modules/settings";
import { projectsModule } from "@/modules/projects";
import { IChannel, IChannelStats } from "@/modules/projects/models";

function clamp(val, rng) {
  return Math.max(Math.min(val, rng[1]), rng[0]);
}

const marginLeft = 2; // Space for 0 tick label on X axis
const marginRight = 60; // space for Y axis & labels
const marginBottom = 20; // space for X axis & labels
const marginTop = 0;

@Component
export default class BrushableHistogram extends Vue {
  readonly settingsContext = settingsModule.context(this.$store);
  readonly projectsContext = projectsModule.context(this.$store);

  @Prop(Object) readonly channel!: IChannel;
  @Prop(Number) readonly containerWidth!: number;

  color = this.channel ? this.metalColor : "#ffffff";

  height = 50 - marginTop - marginBottom;

  unclippedRangeMax = this.channel.max_intensity;

  brushX: any = null;
  brushXselection: any = null;

  get width() {
    return this.containerWidth - marginLeft - marginRight;
  }

  get levels() {
    return this.settings && this.settings.levels
      ? [this.settings.levels.min, this.settings.levels.max]
      : [this.channel.min_intensity, this.channel.max_intensity];
  }

  get metalColor() {
    return this.settings ? this.settings.color : "#ffffff";
  }

  get activeAcquisitionId() {
    return this.projectsContext.getters.activeAcquisitionId;
  }

  get settings() {
    return this.settingsContext.getters.channelsSettings[this.channel.name];
  }

  get label() {
    return this.channel.customLabel;
  }

  @Watch("color")
  onColorChanged(color: string) {
    this.settingsContext.mutations.setChannelColor({ channelName: this.channel.name, color: color });
    this.projectsContext.actions.getChannelStackImage();
  }

  calcHistogramCache = memoize((stats: IChannelStats) => {
    /*
     recalculate expensive stuff, notably bins, summaries, etc.
    */
    const histogramCache: any = {};
    const numBins = 40;

    const domainMin = this.channel.min_intensity;
    const domainMax = this.channel.max_intensity;

    histogramCache.x = d3
      .scaleLinear()
      .domain([domainMin, domainMax])
      .range([marginLeft, marginLeft + this.width]);

    histogramCache.bins = stats.bins;
    histogramCache.binWidth = (domainMax - domainMin) / numBins;

    histogramCache.binStart = (i) => domainMin + i * histogramCache.binWidth;
    histogramCache.binEnd = (i) => domainMin + (i + 1) * histogramCache.binWidth;

    const yMax = histogramCache.bins.reduce((l, r) => (l > r ? l : r));

    histogramCache.y = d3
      .scaleLinear()
      .domain([0, yMax])
      .range([marginTop + this.height, marginTop]);

    return histogramCache;
  });

  onBrushEnd(event, x) {
    const minAllowedBrushSize = 10;
    const smallAmountToAvoidInfiniteLoop = 0.1;

    // ignore programmatically generated events
    if (!event.sourceEvent) return;
    // ignore cascading events, which are programmatically generated
    if (event.sourceEvent.sourceEvent) return;

    if (event.selection) {
      let _range;

      if (event.selection[1] - event.selection[0] > minAllowedBrushSize) {
        _range = [x(event.selection[0]), x(event.selection[1])];
      } else {
        /* the user selected range is too small and will be hidden #587, so take control of it procedurally */
        /* https://stackoverflow.com/questions/12354729/d3-js-limit-size-of-brush */

        const procedurallyResizedBrushWidth = event.selection[0] + minAllowedBrushSize + smallAmountToAvoidInfiniteLoop; //

        _range = [x(event.selection[0]), x(procedurallyResizedBrushWidth)];
      }

      this.submitRange(_range);
    } else {
      this.submitRange([this.channel.min_intensity, this.channel.max_intensity]);
    }
  }

  submitRange(range: number[]) {
    if (this.activeAcquisitionId) {
      this.settingsContext.mutations.setChannelLevels({
        channelName: this.channel.name,
        levels: { min: Math.round(range[0]), max: Math.round(range[1]) },
      });
      this.projectsContext.actions.getChannelStackImage();
    }
  }

  renderHistogram(histogram) {
    const { x, y, bins, binStart, binEnd, binWidth } = histogram;
    const svg = d3.select(this.$refs.svg as any);

    /* Remove everything */
    svg.selectAll("*").remove();

    /* Set margins within the SVG */
    const container = svg
      .attr("width", this.width + marginLeft + marginRight)
      .attr("height", this.height + marginTop + marginBottom)
      .append("g")
      .attr("class", "histogram-container")
      .attr("transform", `translate(${marginLeft},${marginTop})`);

    if (binWidth > 0) {
      /* BINS */
      container
        .insert("g", "*")
        .selectAll("rect")
        .data(bins)
        .enter()
        .append("rect")
        .attr("x", (d, i) => x(binStart(i)) + 1)
        .attr("y", (d) => y(d))
        .attr("width", (d, i) => x(binEnd(i)) - x(binStart(i)) - 1)
        .attr("height", (d) => y(0) - y(d))
        .style("fill", "#bbb");
    }

    // BRUSH
    // Note the brushable area is bounded by the data on three sides, but goes down to cover the x-axis
    const brushX = d3
      .brushX()
      .extent([
        [x.range()[0], y.range()[1]],
        [x.range()[1], marginTop + this.height + marginBottom],
      ])
      /*
      emit start so that the Undoable history can save an undo point
      upon drag start, and ignore the subsequent intermediate drag events.
      */
      .on("end", (event) => this.onBrushEnd(event, x.invert));

    const brushXselection = container.insert("g").attr("class", "brush").call(brushX);

    /* X AXIS */
    container
      .insert("g")
      .attr("class", "axis axis--x")
      .attr("transform", `translate(0,${marginTop + this.height})`)
      .call(
        d3
          .axisBottom(x)
          .ticks(4)
          .tickFormat(d3.format(".0s") as any)
      );

    /* Y AXIS */
    container
      .insert("g")
      .attr("class", "axis axis--y")
      .attr("transform", `translate(${marginLeft + this.width},0)`)
      .call(
        d3
          .axisRight(y)
          .ticks(3)
          .tickFormat(d3.format(".0s") as any)
      );

    /* axis style */
    svg.selectAll(".axis text").style("fill", "rgb(80,80,80)");
    svg.selectAll(".axis path").style("stroke", "rgb(230,230,230)");
    svg.selectAll(".axis line").style("stroke", "rgb(230,230,230)");

    this.brushX = brushX;
    this.brushXselection = brushXselection;
  }

  async mounted() {
    if (this.projectsContext.getters.activeAcquisitionId) {
      const stats = await this.projectsContext.actions.getChannelStats({
        acquisitionId: this.projectsContext.getters.activeAcquisitionId,
        channelName: this.channel.name,
      });
      const histogram = this.calcHistogramCache(stats);
      this.renderHistogram(histogram);

      // if the selection has changed, ensure that the brush correctly reflects the underlying selection.

      if (this.settings && this.settings.levels) {
        const range =
          this.settings!.levels!.min !== this.channel.min_intensity ||
          this.settings!.levels!.max !== this.channel.max_intensity;

        if (this.brushXselection && range) {
          const min = this.channel.min_intensity;
          const max = this.channel.max_intensity;
          const x0 = histogram.x(clamp(this.settings!.levels!.min, [min, max]));
          const x1 = histogram.x(clamp(this.settings!.levels!.max, [min, max]));
          this.brushXselection.call(this.brushX.move, [x0, x1]);
        }
      }
    }
  }
}
</script>

<style scoped>
.root {
  padding: 4px;
  overflow-x: hidden;
  background-color: white;
}
.header {
  display: grid;
  grid-template-columns: 1fr 50px;
}
.labels {
  display: flex;
  justify-content: space-between;
}
.label {
  font-size: x-small;
  font-style: italic;
}
.range-label {
  color: rgb(145, 145, 145);
  font-size: x-small;
}
</style>
