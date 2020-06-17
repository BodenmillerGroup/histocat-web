<template>
  <div id="canvasContainer">
    <svg ref="canvas" id="canvas" v-resize="onResize">
      <g class="vis">
        <g class="marks"></g>
        <g class="axis x"></g>
        <g class="axis y"></g>
        <g class="regression"></g>
      </g>
    </svg>
  </div>
</template>

<script lang="ts">
import { IChart2DData } from "@/modules/analysis/models";
import { experimentModule } from "@/modules/experiment";
import { mainModule } from "@/modules/main";
import { settingsModule } from "@/modules/settings";
import { equals } from "rambda";
import { Component, Prop, Vue, Watch } from "vue-property-decorator";
import * as d3 from "d3";
import { regressionLinear, regressionPoly } from "d3-regression";
import { interpolateReds, scaleSequential, schemeCategory10, scaleOrdinal } from "d3";
import { SelectedCell } from "@/modules/selection/models";
import { selectionModule } from "@/modules/selection";
import { CellPoint } from "@/data/CellPoint";

// A function that return TRUE or FALSE according if a dot is in the selection or not
function isBrushed(brush_coords, cx, cy) {
  const x0 = brush_coords[0][0],
    x1 = brush_coords[1][0],
    y0 = brush_coords[0][1],
    y1 = brush_coords[1][1];
  return x0 <= cx && cx <= x1 && y0 <= cy && cy <= y1; // This return TRUE or FALSE depending on if the points is in the selected area
}

@Component
export default class ScatterPlotView extends Vue {
  readonly mainContext = mainModule.context(this.$store);
  readonly experimentContext = experimentModule.context(this.$store);
  readonly settingsContext = settingsModule.context(this.$store);
  readonly selectionContext = selectionModule.context(this.$store);

  @Prop(Object) data;
  @Prop(String) title;

  @Prop(Boolean) showRegression;
  @Prop(String) regressionType;
  @Prop(Number) polynomialOrder;

  points: CellPoint[] = [];
  selection: any[] = [];
  colorScale: d3.ScaleSequential<any> | d3.ScaleOrdinal<any, any> | null = null;
  maxX = 0;
  maxY = 0;
  quadtree: any = null;

  get showWorkspace() {
    return this.mainContext.getters.showWorkspace;
  }

  get showOptions() {
    return this.mainContext.getters.showOptions;
  }

  get applyMask() {
    return this.settingsContext.getters.maskSettings.apply;
  }

  get activeAcquisitionId() {
    return this.experimentContext.getters.activeAcquisitionId;
  }

  get selectedCells() {
    return this.selectionContext.getters.selectedCells;
  }

  onResize() {
    this.refresh();
  }

  @Watch("showWorkspace")
  showWorkspaceChanged(value) {
    this.refresh();
  }

  @Watch("showOptions")
  showOptionsChanged(value) {
    this.refresh();
  }

  refresh() {
    const canvasElement = this.$refs.canvas as HTMLCanvasElement;
    const { width, height } = canvasElement.getBoundingClientRect();
    const margin = 40;

    const canvas = d3.select("#canvas").attr("width", width).attr("height", height);
    const marks = canvas.select(".vis .marks");
    const scaleX = d3
      .scaleLinear()
      .domain([0, this.maxX])
      .range([0, width - 2 * margin]);

    const scaleY = d3
      .scaleLinear()
      .domain([0, this.maxY])
      .range([height - 2 * margin, 0]);

    const axisX = d3.axisBottom(scaleX);
    const axisY = d3.axisLeft(scaleY);

    d3.select("g.axis.x")
      .call(axisX as any)
      .attr("transform", `translate(0, ${height - margin - margin})`);

    d3.select("g.axis.y").call(axisY as any);
    d3.select(".vis").attr("transform", `translate(${margin}, ${margin})`);
    const pointMark = marks
      .selectAll("circle.point")
      .data(this.points)
      .join("circle")
      .attr("class", "point")
      .attr("fill", (d) => this.colorScale!(d.value))
      .attr("fill-opacity", 0.73)
      .attr("r", 2)
      .attr("cx", (d) => scaleX(d.x))
      .attr("cy", (d) => scaleY(d.y));

    const quadtree = d3
      .quadtree()
      .extent([
        [-1, -1],
        [width + 1, height + 1],
      ])
      .x((d) => scaleX(d["x"]))
      .y((d) => scaleY(d["y"]))
      .addAll(this.points as any);

    canvas
      .select(".vis")
      .selectAll(".node")
      .data(this.nodes(quadtree))
      .enter()
      .append("rect")
      .attr("class", "node")
      .attr("x", function (d: any) {
        return d.x0;
      })
      .attr("y", function (d: any) {
        return d.y0;
      })
      .attr("width", function (d: any) {
        return d.y1 - d.y0;
      })
      .attr("height", function (d: any) {
        return d.x1 - d.x0;
      });

    // Add brushing
    const brush = d3
      .brush()
      .on("end", () => {
        pointMark.each(function (d) {
          d.selected = false;
        });

        const extent = d3.event.selection; // looks like [[12,11], [132,178]]
        if (!extent) {return}

        const x0 = extent[0][0];
        const y0 = extent[0][1];
        const x3 = extent[1][0];
        const y3 = extent[1][1];
        quadtree.visit((node: any, x1, y1, x2, y2) => {
          if (!node.length) {
            do {
              var d = node.data;
              d.selected = scaleX(d.x) >= x0 && scaleX(d.x) < x3 && scaleY(d.y) >= y0 && scaleY(d.y) < y3;
            } while ((node = node.next));
          }
          return x1 >= x3 || y1 >= y3 || x2 < x0 || y2 < y0;
        });

        pointMark.classed("selected-point", (d) => {
          return d.selected!;
        });

        // pointMark.classed("selected-point", (d) => {
        //   return isBrushed(extent, scaleX(d.x), scaleY(d.y));
        // });
      })
      .extent([
        // initialise the brush area: start at 0,0 and finishes at width,height, which means the whole graph area
        [0, 0],
        [width - 2 * margin, height - 2 * margin],
      ]);
    marks.call(brush as any);

    const regressionSvg = canvas.select(".vis .regression");
    regressionSvg.selectAll("*").remove();
    if (this.showRegression) {
      const regression = (this.regressionType === "linear"
        ? regressionLinear()
        : regressionPoly().order(this.polynomialOrder)
      )
        .x((d) => d.x)
        .y((d) => d.y);

      const line = d3
        .line()
        .x((d) => scaleX(d[0]))
        .y((d) => scaleY(d[1]));

      regressionSvg
        .append("path")
        .datum(regression(this.points))
        .attr("fill", "none")
        .attr("stroke", "#ff0000")
        .attr("stroke-width", 1.5)
        .attr("d", line);
    }

    // this.activeLasso = lasso()
    //   .closePathSelect(true)
    //   .closePathDistance(100)
    //   .items(d3.selectAll("circle.point"))
    //   .targetArea(marks)
    //   .on("start", this.lassoStart)
    //   .on("draw", this.lassoDraw)
    //   .on("end", this.lassoEnd);
    //
    // marks.call(this.activeLasso);

    // // zoom modifies the scales by applying a transformation to it
    // // this transform is however cumulative, so we need to maintain a copy
    // // of the original scales to not overemphasize interactions
    // // try replacing the scaleXCopy inside the zoom event listener to see
    // // how that influences the behavior of the zoom
    // const scaleXCopy = scaleX.copy();
    // const scaleYCopy = scaleY.copy();
    //
    // // generator function for zoom
    // const zoom = d3.zoom().on("zoom", () => {
    //   // applying the zoom transformation to the vertical and horizontal
    //   // scales
    //   const rescaledX = d3.event.transform.rescaleX(scaleXCopy);
    //   const rescaledY = d3.event.transform.rescaleY(scaleYCopy);
    //
    //   // updating the point mark
    //   marks
    //     .selectAll("circle.point")
    //     .attr("cx", (d: any) => rescaledX(d.x))
    //     .attr("cy", (d: any) => rescaledY(d.y));
    //
    //   // reconfiguring the axis generators
    //   axisX.scale(rescaledX);
    //   axisY.scale(rescaledY);
    //
    //   // redrawing the axes
    //   d3.select("g.axis.x").call(axisX as any);
    //   d3.select("g.axis.y").call(axisY as any);
    // });
    //
    // // calling the zoom generator on a rectangle selection. All zoom events on
    // // this rectangle receives all zoom events and is not influenced by the
    // // transformation, thus it always stays in place.
    // canvas
    //   .append("rect")
    //   .attr("id", "zoomHelper")
    //   .attr("fill", "transparent")
    //   .attr("width", width)
    //   .attr("height", height)
    //   .call(zoom as any);
  }

  // Collapse the quadtree into an array of rectangles.
  nodes(quadtree) {
    var nodes: any = [];
    quadtree.visit((node, x0, y0, x1, y1) => {
      (node.x0 = x0), (node.y0 = y0);
      (node.x1 = x1), (node.y1 = y1);
      nodes.push(node);
    });
    return nodes;
  }

  // Find the nodes within the specified rectangle.
  search(quadtree, x0, y0, x3, y3) {
    quadtree.visit((node, x1, y1, x2, y2) => {
      if (!node.length) {
        do {
          var d = node.data;
          d.selected = d.x >= x0 && d.x < x3 && d.y >= y0 && d.y < y3;
        } while ((node = node.next));
      }
      return x1 >= x3 || y1 >= y3 || x2 < x0 || y2 < y0;
    });
  }

  @Watch("data")
  dataChanged(data: IChart2DData) {
    if (data) {
      if (data.heatmap) {
        // Use heatmap value
        this.points = data.x.data.map(
          (x, i) =>
            // Object.freeze(
            new CellPoint(data.acquisitionIds[i], data.cellIds[i], x, data.y.data[i], data.heatmap!.data[i])
          // )
        );
        this.colorScale = scaleSequential(interpolateReds).domain([0, d3.max(data.heatmap!.data)!]);
      } else {
        // Use acquisitionId as category
        this.points = data.x.data.map(
          (x, i) =>
            // Object.freeze(
            new CellPoint(data.acquisitionIds[i], data.cellIds[i], x, data.y.data[i], data.acquisitionIds[i])
          // )
        );
        this.colorScale = scaleOrdinal(schemeCategory10);
      }
      this.maxX = d3.max(data.x.data)!;
      this.maxY = d3.max(data.y.data)!;

      this.refresh();
    }
  }

  // lassoStart() {
  //   if (!this.activeLasso) {
  //     return;
  //   }
  //   this.activeLasso
  //     .items()
  //     .attr("r", 2) // reset size
  //     .classed("not_possible", true)
  //     .classed("selected", false);
  // }
  //
  // lassoDraw() {
  //   if (!this.activeLasso) {
  //     return;
  //   }
  //   // Style the possible dots
  //   this.activeLasso.possibleItems().classed("not_possible", false).classed("possible", true);
  //
  //   // Style the not possible dot
  //   this.activeLasso.notPossibleItems().classed("not_possible", true).classed("possible", false);
  // }
  //
  // lassoEnd() {
  //   if (!this.activeLasso) {
  //     return;
  //   }
  //   // Reset the color of all dots
  //   this.activeLasso.items().classed("not_possible", false).classed("possible", false);
  //
  //   // Style the selected dots
  //   this.activeLasso.selectedItems().classed("selected", true).attr("r", 3);
  //
  //   // Reset the style of the not selected dots
  //   this.activeLasso.notSelectedItems().attr("r", 2);
  // }

  @Watch("selectedCells")
  selectedCellsChanged(data: Map<number, SelectedCell[]> | null) {
    if (data !== null) {
      let allCells: SelectedCell[] = [];
      data.forEach((val, key) => {
        allCells = allCells.concat(val);
      });
      const indices = allCells.map((i) => i.index);
      if (!equals(indices, this.selection)) {
        this.selection = indices;
      }
    }
  }
}
</script>

<style scoped>
#canvasContainer {
  height: calc(100vh - 84px);
  position: relative;
  width: 100%;
}
#canvas {
  height: 100%;
  position: absolute;
  width: 100%;
}
</style>

<style>
.selected-point {
  opacity: 1 !important;
  stroke: red;
  stroke-width: 1px;
}
.node {
  fill: none;
  stroke: #ccc;
  shape-rendering: crispEdges;
}
</style>
