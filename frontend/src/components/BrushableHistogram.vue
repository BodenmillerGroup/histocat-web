<template>
  <div class="root">
    <div class="container">
      <v-tooltip bottom>
        <template v-slot:activator="{ on }">
          <v-btn x-small elevation="1" v-on="on" class="mr-2" @click.stop="setSharedChannelLevels">
            Share
          </v-btn>
        </template>
        <span>Share levels and color</span>
      </v-tooltip>
      <input type="color" v-model.lazy="color" @click.stop />
    </div>
    <svg :width="width" :height="height" id="svg" ref="svg" />
    <div class="labels">
      <span class="range-label"> min {{ unclippedRangeMin.toPrecision(4) }} </span>
      <span class="label">
        {{ label }}
      </span>
      <span class="range-label"> max {{ unclippedRangeMax.toPrecision(4) }} </span>
    </div>
  </div>
</template>

<script lang="ts">
import { Component, Prop, Vue, Watch } from "vue-property-decorator";
import * as d3 from "d3";
import memoize from "memoize-one";
import { settingsModule } from "@/modules/settings";
import { experimentModule } from "@/modules/experiment";
import { IChannel, IChannelStats } from "@/modules/experiment/models";
import { BroadcastManager } from "@/utils/BroadcastManager";
import { SET_CHANNEL_SETTINGS, SET_METAL_COLOR } from "@/modules/settings/events";

function clamp(val, rng) {
  return Math.max(Math.min(val, rng[1]), rng[0]);
}

@Component
export default class BrushableHistogram extends Vue {
  readonly settingsContext = settingsModule.context(this.$store);
  readonly experimentContext = experimentModule.context(this.$store);

  @Prop(Object) readonly channel!: IChannel;

  color = this.channel ? this.metalColor : "#ffffff";

  marginLeft = 10; // Space for 0 tick label on X axis
  marginRight = 54; // space for Y axis & labels
  marginBottom = 25; // space for X axis & labels
  marginTop = 3;

  width = 340 - this.marginLeft - this.marginRight;
  height = 80 - this.marginTop - this.marginBottom;

  unclippedRangeMin = this.channel.min_intensity;
  unclippedRangeMax = this.channel.max_intensity;

  brushX: any = null;
  brushXselection: any = null;

  levels: number[] =
    this.settings && this.settings.levels
      ? [this.settings.levels.min, this.settings.levels.max]
      : [this.channel.min_intensity, this.channel.max_intensity];

  get metalColor() {
    const colorMap = this.settingsContext.getters.metalColorMap;
    const color = colorMap.get(this.channel.name);
    return color ? color : "#ffffff";
  }

  get activeAcquisitionId() {
    return this.experimentContext.getters.activeAcquisitionId;
  }

  get settings() {
    return this.settingsContext.getters.getChannelSettings(this.activeAcquisitionId, this.channel.name);
  }

  get label() {
    return this.settings && this.settings.customLabel ? this.settings.customLabel : this.channel.label;
  }

  @Watch("color")
  onColorChanged(color: string) {
    BroadcastManager.publish(SET_METAL_COLOR, {
      metal: this.channel.name,
      color: color,
    });
    this.experimentContext.actions.getChannelStackImage();
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
      .range([this.marginLeft, this.marginLeft + this.width]);

    histogramCache.bins = stats.bins;
    histogramCache.binWidth = (domainMax - domainMin) / numBins;

    histogramCache.binStart = (i) => domainMin + i * histogramCache.binWidth;
    histogramCache.binEnd = (i) => domainMin + (i + 1) * histogramCache.binWidth;

    const yMax = histogramCache.bins.reduce((l, r) => (l > r ? l : r));

    histogramCache.y = d3
      .scaleLinear()
      .domain([0, yMax])
      .range([this.marginTop + this.height, this.marginTop]);

    return histogramCache;
  });

  onBrush(x, eventType) {
    const type = `continuous metadata histogram ${eventType}`;
    return () => {
      // const { dispatch, field, isObs, isUserDefined, isDiffExp } = this.props;

      // ignore programmatically generated events
      if (!d3.event.sourceEvent) return;
      // ignore cascading events, which are programmatically generated
      if (d3.event.sourceEvent.sourceEvent) return;

      if (d3.event.selection) {
        // dispatch({
        //   type,
        //   selection: field,
        //   continuousNamespace: {
        //     isObs,
        //     isUserDefined,
        //     isDiffExp,
        //   },
        //   range: [x(d3.event.selection[0]), x(d3.event.selection[1])],
        // });
      } else {
        // dispatch({
        //   type,
        //   selection: field,
        //   continuousNamespace: {
        //     isObs,
        //     isUserDefined,
        //     isDiffExp,
        //   },
        //   range: null,
        // });
      }
    };
  }

  onBrushEnd(x) {
    return () => {
      const minAllowedBrushSize = 10;
      const smallAmountToAvoidInfiniteLoop = 0.1;

      // ignore programmatically generated events
      if (!d3.event.sourceEvent) return;
      // ignore cascading events, which are programmatically generated
      if (d3.event.sourceEvent.sourceEvent) return;

      if (d3.event.selection) {
        let _range;

        if (d3.event.selection[1] - d3.event.selection[0] > minAllowedBrushSize) {
          _range = [x(d3.event.selection[0]), x(d3.event.selection[1])];
        } else {
          /* the user selected range is too small and will be hidden #587, so take control of it procedurally */
          /* https://stackoverflow.com/questions/12354729/d3-js-limit-size-of-brush */

          const procedurallyResizedBrushWidth =
            d3.event.selection[0] + minAllowedBrushSize + smallAmountToAvoidInfiniteLoop; //

          _range = [x(d3.event.selection[0]), x(procedurallyResizedBrushWidth)];
        }

        this.submitRange(_range);

        // dispatch({
        //   type: "continuous metadata histogram end",
        //   selection: field,
        //   continuousNamespace: {
        //     isObs,
        //     isUserDefined,
        //     isDiffExp,
        //   },
        //   range: _range,
        // });
      } else {
        this.submitRange([this.channel.min_intensity, this.channel.max_intensity]);
        // dispatch({
        //   type: "continuous metadata histogram cancel",
        //   selection: field,
        //   continuousNamespace: {
        //     isObs,
        //     isUserDefined,
        //     isDiffExp,
        //   },
        // });
      }
    };
  }

  submitRange(range: number[]) {
    if (!this.activeAcquisitionId) {
      return;
    }
    let settings = this.settings;
    if (!settings) {
      settings = {
        acquisitionId: this.activeAcquisitionId,
        name: this.channel.name,
        customLabel: this.channel.label,
        levels: { min: Math.round(range[0]), max: Math.round(range[1]) },
        suppressBroadcast: false,
      };
    } else {
      settings = {
        ...settings,
        levels: { min: Math.round(range[0]), max: Math.round(range[1]) },
        suppressBroadcast: false,
      };
    }
    BroadcastManager.publish(SET_CHANNEL_SETTINGS, settings);
    this.experimentContext.actions.getChannelStackImage();
  }

  setSharedChannelLevels() {
    const metal = this.channel.name;
    const settings = this.settings;
    const levels =
      settings && settings.levels
        ? [settings.levels.min, settings.levels.max]
        : [this.channel.min_intensity, this.channel.max_intensity];
    this.experimentContext.actions.setSharedChannelLevels({ metal: metal, levels: levels });
  }

  renderHistogram(histogram) {
    const { x, y, bins, binStart, binEnd, binWidth } = histogram;
    const svg = d3.select(this.$refs.svg as any);

    /* Remove everything */
    svg.selectAll("*").remove();

    /* Set margins within the SVG */
    const container = svg
      .attr("width", this.width + this.marginLeft + this.marginRight)
      .attr("height", this.height + this.marginTop + this.marginBottom)
      .append("g")
      .attr("class", "histogram-container")
      .attr("transform", `translate(${this.marginLeft},${this.marginTop})`);

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
        [x.range()[1], this.marginTop + this.height + this.marginBottom],
      ])
      /*
      emit start so that the Undoable history can save an undo point
      upon drag start, and ignore the subsequent intermediate drag events.
      */
      .on("start", this.onBrush(x.invert, "start").bind(this))
      .on("brush", this.onBrush(x.invert, "brush").bind(this))
      .on("end", this.onBrushEnd(x.invert).bind(this));

    const brushXselection = container.insert("g").attr("class", "brush").call(brushX);

    /* X AXIS */
    container
      .insert("g")
      .attr("class", "axis axis--x")
      .attr("transform", `translate(0,${this.marginTop + this.height})`)
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
      .attr("transform", `translate(${this.marginLeft + this.width},0)`)
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
    if (this.experimentContext.getters.activeAcquisitionId) {
      const stats = await this.experimentContext.actions.getChannelStats({
        acquisitionId: this.experimentContext.getters.activeAcquisitionId,
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
  padding: 10px;
  background-color: white;
}
.container {
  display: flex;
  justify-content: flex-end;
  padding-bottom: 8px;
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
  color: #bbb;
  font-size: xx-small;
}
</style>
