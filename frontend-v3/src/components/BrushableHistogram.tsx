import styles from "./BrushableHistogram.module.scss";
import { useEffect, useRef } from "react";
import * as d3 from "d3";
import memoize from "memoize-one";
import { IChannel, IChannelStats } from "../modules/projects/models";
import { useProjectsStore } from "../modules/projects";
import shallow from "zustand/shallow";
import { useSettingsStore } from "../modules/settings";
import { debounce } from "lodash-es";

type BrushableHistogramProps = {
  channel: IChannel;
};

function clamp(val: number, rng: number[]) {
  return Math.max(Math.min(val, rng[1]), rng[0]);
}

const marginLeft = 10; // Space for 0 tick label on X axis
const marginRight = 54; // space for Y axis & labels
const marginBottom = 25; // space for X axis & labels
const marginTop = 3;

const width = 340 - marginLeft - marginRight;
const height = 80 - marginTop - marginBottom;

export function BrushableHistogram(props: BrushableHistogramProps) {
  const svgRef = useRef<any>(null);
  const { activeAcquisitionId, getChannelStats, getChannelStackImage } = useProjectsStore(
    (state) => ({
      activeAcquisitionId: state.activeAcquisitionId,
      getChannelStats: state.getChannelStats,
      getChannelStackImage: state.getChannelStackImage,
    }),
    shallow
  );
  const { channelsSettings, setChannelLevels, setChannelColor } = useSettingsStore(
    (state) => ({
      channelsSettings: state.channelsSettings,
      setChannelLevels: state.setChannelLevels,
      setChannelColor: state.setChannelColor,
    }),
    shallow
  );

  const settings = channelsSettings[props.channel.name];
  const metalColor = settings ? settings.color : "#ffffff";
  const color = props.channel ? metalColor : "#ffffff";
  const unclippedRangeMin = props.channel.min_intensity;
  const unclippedRangeMax = props.channel.max_intensity;

  const calcHistogramCache = memoize((stats: IChannelStats) => {
    /*
     recalculate expensive stuff, notably bins, summaries, etc.
    */
    const histogramCache: any = {};
    const numBins = 40;

    const domainMin = props.channel.min_intensity;
    const domainMax = props.channel.max_intensity;

    histogramCache.x = d3
      .scaleLinear()
      .domain([domainMin, domainMax])
      .range([marginLeft, marginLeft + width]);

    histogramCache.bins = stats.bins;
    histogramCache.binWidth = (domainMax - domainMin) / numBins;

    histogramCache.binStart = (i: number) => domainMin + i * histogramCache.binWidth;
    histogramCache.binEnd = (i: number) => domainMin + (i + 1) * histogramCache.binWidth;

    const yMax = histogramCache.bins.reduce((l: number, r: number) => (l > r ? l : r));

    histogramCache.y = d3
      .scaleLinear()
      .domain([0, yMax])
      .range([marginTop + height, marginTop]);

    return histogramCache;
  });

  useEffect(() => {
    const mounted = async () => {
      if (activeAcquisitionId) {
        const stats = await getChannelStats(activeAcquisitionId, props.channel.name);
        const histogram = calcHistogramCache(stats!);
        const { x, y, bins, binStart, binEnd, binWidth } = histogram;
        const svg = d3.select(svgRef.current);

        /* Remove everything */
        svg.selectAll("*").remove();

        /* Set margins within the SVG */
        const container = svg
          .attr("width", width + marginLeft + marginRight)
          .attr("height", height + marginTop + marginBottom)
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
            [x.range()[1], marginTop + height + marginBottom],
          ])
          /*
      emit start so that the Undoable history can save an undo point
      upon drag start, and ignore the subsequent intermediate drag events.
      */
          .on("start", (event) => onBrush(event, x.invert))
          .on("brush", (event) => onBrush(event, x.invert))
          .on("end", (event) => onBrushEnd(event, x.invert));

        const brushXselection = container.insert("g").attr("class", "brush").call(brushX);

        /* X AXIS */
        container
          .insert("g")
          .attr("class", "axis axis--x")
          .attr("transform", `translate(0,${marginTop + height})`)
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
          .attr("transform", `translate(${marginLeft + width},0)`)
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

        // if the selection has changed, ensure that the brush correctly reflects the underlying selection.

        if (settings && settings.levels) {
          const range =
            settings!.levels!.min !== props.channel.min_intensity ||
            settings!.levels!.max !== props.channel.max_intensity;

          if (brushXselection && range) {
            const min = props.channel.min_intensity;
            const max = props.channel.max_intensity;
            const x0 = histogram.x(clamp(settings!.levels!.min, [min, max]));
            const x1 = histogram.x(clamp(settings!.levels!.max, [min, max]));
            brushXselection.call(brushX.move, [x0, x1]);
          }
        }
      }
    };
    mounted();
  }, [activeAcquisitionId]);

  const onBrush = (event: any, x: any) => {
    // const { dispatch, field, isObs, isUserDefined, isDiffExp } = this.props;

    // ignore programmatically generated events
    if (!event.sourceEvent) return;
    // ignore cascading events, which are programmatically generated
    if (event.sourceEvent.sourceEvent) return;

    if (event.selection) {
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

  const onBrushEnd = (event: any, x: any) => {
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

      submitRange(_range);
    } else {
      submitRange([props.channel.min_intensity, props.channel.max_intensity]);
    }
  };

  const submitRange = (range: number[]) => {
    if (!activeAcquisitionId) {
      return;
    }
    setChannelLevels(props.channel.name, { min: Math.round(range[0]), max: Math.round(range[1]) });
    getChannelStackImage();
  };

  return (
    <div className={styles.root}>
      <div className={styles.container}>
        <input
          type="color"
          color={color}
          onClick={(e) => e.stopPropagation()}
          onChange={debounce(
            function (e) {
              setChannelColor(props.channel.name, e.target.value);
              getChannelStackImage();
            },
            500,
            { leading: false }
          )}
        />
      </div>
      <svg width={width} height={height} id="svg" ref={svgRef} />
      <div className={styles.labels}>
        <span className={styles.rangeLabel}> min {unclippedRangeMin.toPrecision(4)} </span>
        <span className={styles.label}>{props.channel.customLabel}</span>
        <span className={styles.rangeLabel}> max {unclippedRangeMax.toPrecision(4)} </span>
      </div>
    </div>
  );
}
