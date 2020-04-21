<template>
  <div id="graphWrapper">
    <div :style="style">
      <div id="graphAttachPoint">
        <GraphOverlayLayer
          :cameraTF="cameraTF"
          :modelTF="modelTF"
          :projectionTF="projectionTF"
          :graphPaddingRightLeft="graphPaddingRightLeft"
          :graphPaddingTop="graphPaddingTop"
          :responsive="responsive"
          @handleCanvasEvent="graphInteractionMode === 'zoom' ? handleCanvasEvent : undefined"
        >
          <CentroidLabels />
        </GraphOverlayLayer>

        <svg
          id="lasso-layer"
          data-testid="layout-overlay"
          class="graphSVG"
          :width="responsive.width - graphPaddingRightLeft"
          :height="responsive.height"
          :pointer-events="graphInteractionMode === 'select' ? 'auto' : 'none'"
          :style="{ zIndex: 89 }"
        />
      </div>
      <div :style="{ padding: 0, margin: 0 }">
        <canvas
          :width="responsive.width - graphPaddingRightLeft"
          :height="responsive.height - graphPaddingTop"
          data-testid="layout-graph"
          ref="canvas"
          @onMouseDown="handleCanvasEvent"
          @onMouseUp="handleCanvasEvent"
          @onMouseMove="handleCanvasEvent"
          @onDoubleClick="handleCanvasEvent"
          @onWheel="handleCanvasEvent"
        />
      </div>
    </div>
  </div>
</template>

<script lang="ts">
import * as d3 from "d3";
import { mat3, ReadonlyVec2, vec2 } from "gl-matrix";
import _regl from "regl";
import memoize from "memoize-one";

import * as globals from "../../globals";
import setupSVGandBrushElements from "./setupSVGandBrush";
import _camera from "../../util/camera";
import _drawPoints from "./drawPointsRegl";
import { isTypedArray } from "../../util/typeHelpers";
import { Component, Prop, Vue, Watch } from "vue-property-decorator";
import { centroidLabelsModule } from "@/modules/centroidLabels";
import { categoryLabelDisplayStringLongLength } from "@/cellxgene/globals";
import { colorsModule } from "@/modules/colors";
import { pointDilationModule } from "@/modules/pointDilation";
import GraphOverlayLayer from "@/cellxgene/components/graph/overlays/graphOverlayLayer.vue";
import CentroidLabels from "@/cellxgene/components/graph/overlays/centroidLabels.vue";
import { universeModule } from "@/modules/universe";
import { worldModule } from "@/modules/world";
import {responsiveModule} from "@/modules/responsive";
import {crossfilterModule} from "@/modules/crossfilter";
import {controlsModule} from "@/modules/controls";
import {layoutChoiceModule} from "@/modules/layoutChoice";
import {graphSelectionModule} from "@/modules/graphSelection";
import {graphModule} from "@/modules/graph";

/*
Simple 2D transforms control all point painting.  There are three:
  * model - convert from underlying per-point coordinate to a layout.
    Currently used to move from data to webgl coordinate system.
  * camera - apply a 2D camera transformation (pan, zoom)
  * projection - apply any transformation required for screen size and layout
*/

function createProjectionTF(viewportWidth, viewportHeight) {
  /*
  the projection transform accounts for the screen size & other layout
  */
  const fractionToUse = 0.95; // fraction of min dimension to use
  const topGutterSizePx = 32; // toolbar box height
  const heightMinusGutter = viewportHeight - topGutterSizePx;
  const minDim = Math.min(viewportWidth, heightMinusGutter);
  const aspectScale: ReadonlyVec2 = [
    (fractionToUse * minDim) / viewportWidth,
    (fractionToUse * minDim) / viewportHeight,
  ];
  const m = mat3.create();
  mat3.fromTranslation(m, [0, -topGutterSizePx / viewportHeight / aspectScale[1]]);
  mat3.scale(m, m, aspectScale);
  return m;
}

function createModelTF() {
  /*
  preallocate coordinate system transformation between data and gl.
  Data arrives in a [0,1] range, and we operate elsewhere in [-1,1].
  */
  const m = mat3.fromScaling(mat3.create(), [2, 2]);
  mat3.translate(m, m, [-0.5, -0.5]);
  return m;
}

function renderThrottle(callback) {
  /*
  This wraps a call to requestAnimationFrame(), enforcing a single
  render callback at any given time (ie, you can call this any number
  of times, and it will coallesce multiple inter-frame calls into a
  single render).
  */
  let rafCurrentlyInProgress: any = null;
  return function f(this: any) {
    if (rafCurrentlyInProgress) return;
    const context = this;
    rafCurrentlyInProgress = window.requestAnimationFrame(() => {
      callback.apply(context);
      rafCurrentlyInProgress = null;
    });
  };
}

@Component({
  components: { CentroidLabels, GraphOverlayLayer },
})
export default class Graph extends Vue {
  private readonly responsiveContext = responsiveModule.context(this.$store);
  private readonly centroidLabelsContext = centroidLabelsModule.context(this.$store);
  private readonly colorsContext = colorsModule.context(this.$store);
  private readonly pointDilationContext = pointDilationModule.context(this.$store);
  private readonly universeContext = universeModule.context(this.$store);
  private readonly worldContext = worldModule.context(this.$store);
  private readonly crossfilterContext = crossfilterModule.context(this.$store);
  private readonly controlsContext = controlsModule.context(this.$store);
  private readonly layoutChoiceContext = layoutChoiceModule.context(this.$store);
  private readonly graphSelectionContext = graphSelectionModule.context(this.$store);
  private readonly graphContext = graphModule.context(this.$store);

  count = 0;
  graphPaddingTop = 0;
  graphPaddingRightLeft = globals.leftSidebarWidth * 2;
  renderCache = {
    X: null,
    Y: null,
    positions: null,
    colors: null,
    sizes: null,
    flags: null,
  };
  toolSVG: any = null;
  tool: any = null;
  container: any = null;
  cameraRender = 0;

  // state
  regl: any = null;
  drawPoints: any = null;
  pointBuffer: any = null;
  colorBuffer: any = null;
  flagBuffer: any = null;
  camera: any = null;
  modelTF: any = null;
  modelInvTF: any = null;
  projectionTF: any = null;
  updateOverlay = false;

  get reglCanvas() {
    return this.$refs.canvas as HTMLCanvasElement;
  }

  get cameraTF() {
    return this.camera?.view()?.slice();
  }

  get universe() {
    return this.universeContext.getters.universe;
  }

  get responsive() {
    return this.responsiveContext.getters.responsive;
  }

  get currentSelection() {
    return this.graphSelectionContext.getters.selection;
  }

  get selectionTool() {
    return this.graphSelectionContext.getters.tool;
  }

  get graphInteractionMode() {
    return this.controlsContext.getters.graphInteractionMode;
  }

  get style() {
    return {
      zIndex: -9999,
      position: "fixed",
      top: this.graphPaddingTop,
      right: globals.leftSidebarWidth,
    };
  }

  computePointPositions = memoize((X, Y, modelTF) => {
    /*
    compute the model coordinate for each point
    */
    const positions = new Float32Array(2 * X.length);
    for (let i = 0, len = X.length; i < len; i += 1) {
      const p = vec2.fromValues(X[i], Y[i]);
      vec2.transformMat3(p, p, modelTF);
      positions[2 * i] = p[0];
      positions[2 * i + 1] = p[1];
    }
    return positions;
  });

  computePointColors = memoize((rgb) => {
    /*
    compute webgl colors for each point
    */
    const colors = new Float32Array(3 * rgb.length);
    for (let i = 0, len = rgb.length; i < len; i += 1) {
      colors.set(rgb[i], 3 * i);
    }
    return colors;
  });

  computeSelectedFlags = memoize((crossfilter, flagSelected, flagUnselected) => {
    const x = crossfilter.fillByIsSelected(new Float32Array(crossfilter.size()), flagSelected, flagUnselected);
    return x;
  });

  computePointFlags = memoize((world, crossfilter, colorAccessor, pointDilation) => {
    /*
      We communicate with the shader using three flags:
      - isNaN -- the value is a NaN. Only makes sense when we have a colorAccessor
      - isSelected -- the value is selected
      - isHightlighted -- the value is highlighted in the UI (orthogonal from selection highlighting)

      Due to constraints in webgl vertex shader attributes, these are encoded in a float, "kinda"
      like bitmasks.

      We also have separate code paths for generating flags for categorical and
      continuous metadata, as they rely on different tests, and some of the flags
      (eg, isNaN) are meaningless in the face of categorical metadata.
      */

    const flagSelected = 1;
    const flagNaN = 2;
    const flagHighlight = 4;

    const flags = this.computeSelectedFlags(crossfilter, flagSelected, 0).slice();

    const { metadataField, categoryField } = pointDilation;
    const highlightData = metadataField ? world.obsAnnotations.col(metadataField)?.asArray() : null;
    const colorByColumn = colorAccessor
      ? world.obsAnnotations.col(colorAccessor)?.asArray() || world.varData.col(colorAccessor)?.asArray()
      : null;
    const colorByData = colorByColumn && isTypedArray(colorByColumn) ? colorByColumn : null;

    if (colorByData || highlightData) {
      for (let i = 0, len = flags.length; i < len; i += 1) {
        if (highlightData) {
          flags[i] += highlightData[i] === categoryField ? flagHighlight : 0;
        }
        if (colorByData) {
          flags[i] += Number.isFinite(colorByData[i]) ? 0 : flagNaN;
        }
      }
    }
    return flags;
  });

  handleCanvasEvent = (e) => {
    if (e.type !== "wheel") e.preventDefault();
    if (this.camera.handleEvent(e, this.projectionTF)) {
      this.renderCanvas();
      this.updateOverlay = !this.updateOverlay;
    }
  };

  createToolSVG = () => {
    /*
    Called from componentDidUpdate. Create the tool SVG, and return any
    state changes that should be passed to setState().
    */

    /* clear out whatever was on the div, even if nothing, but usually the brushes etc */

    d3.select("#lasso-layer")
      .selectAll(".lasso-group")
      .remove();

    // Don't render or recreate toolSVG if currently in zoom mode
    if (this.graphInteractionMode !== "select") {
      // don't return "change" of state unless we are really changing it!
      if (this.toolSVG === undefined) return {};
      return { toolSVG: undefined };
    }

    let handleStart;
    let handleDrag;
    let handleEnd;
    let handleCancel;
    if (this.selectionTool === "brush") {
      handleStart = this.handleBrushStartAction.bind(this);
      handleDrag = this.handleBrushDragAction.bind(this);
      handleEnd = this.handleBrushEndAction.bind(this);
    } else {
      handleStart = this.handleLassoStart.bind(this);
      handleEnd = this.handleLassoEnd.bind(this);
      handleCancel = this.handleLassoCancel.bind(this);
    }

    const { svg: newToolSVG, tool, container } = setupSVGandBrushElements(
      this.selectionTool,
      handleStart,
      handleDrag,
      handleEnd,
      handleCancel,
      this.responsive,
      this.graphPaddingRightLeft
    );

    return { toolSVG: newToolSVG, tool, container };
  };

  brushToolUpdate(tool, container) {
    /*
    this is called from componentDidUpdate(), so be very careful using
    anything from this.state, which may be updated asynchronously.
    */
    const currentSelection = this.currentSelection;
    if (container) {
      const toolCurrentSelection = d3.brushSelection(container.node());

      if (currentSelection.mode === "within-rect") {
        /*
        if there is a selection, make sure the brush tool matches
        */
        const screenCoords = [
          this.mapPointToScreen(currentSelection.brushCoords.northwest),
          this.mapPointToScreen(currentSelection.brushCoords.southeast)
        ];
        if (!toolCurrentSelection) {
          /* tool is not selected, so just move the brush */
          container.call(tool.move, screenCoords);
        } else {
          /* there is an active selection and a brush - make sure they match */
          /* this just sums the difference of each dimension, of each point */
          let delta = 0;
          for (let x = 0; x < 2; x += 1) {
            for (let y = 0; y < 2; y += 1) {
              delta += Math.abs(
                screenCoords[x][y] - toolCurrentSelection[x][y]
              );
            }
          }
          if (delta > 0) {
            container.call(tool.move, screenCoords);
          }
        }
      } else if (toolCurrentSelection) {
        /* no selection, so clear the brush tool if it is set */
        container.call(tool.move, null);
      }
    }
  }

  lassoToolUpdate(tool, container) {
    /*
    this is called from componentDidUpdate(), so be very careful using
    anything from this.state, which may be updated asynchronously.
    */
    const currentSelection = this.currentSelection;
    if (currentSelection.mode === "within-polygon") {
      /*
      if there is a current selection, make sure the lasso tool matches
      */
      const polygon = currentSelection.polygon.map(p =>
        this.mapPointToScreen(p)
      );
      tool.move(polygon);
    } else {
      tool.reset();
    }
  }

  selectionToolUpdate(tool, container) {
    /*
    this is called from componentDidUpdate(), so be very careful using
    anything from this.state, which may be updated asynchronously.
    */
    switch (this.selectionTool) {
      case "brush":
        this.brushToolUpdate(tool, container);
        break;
      case "lasso":
        this.lassoToolUpdate(tool, container);
        break;
      default:
        /* punt? */
        break;
    }
  }

  mapScreenToPoint(pin) {
    /*
    Map an XY coordinates from screen domain to cell/point range,
    accounting for current pan/zoom camera.
    */

    const responsive = this.responsive;
    const cameraInvTF = this.camera.invView();

    /* screen -> gl */
    const x =
      (2 * pin[0]) / (responsive.width! - this.graphPaddingRightLeft) - 1;
    const y = 2 * (1 - pin[1] / (responsive.height! - this.graphPaddingTop)) - 1;

    const xy = vec2.fromValues(x, y);
    const projectionInvTF = mat3.invert(mat3.create(), this.projectionTF);
    vec2.transformMat3(xy, xy, projectionInvTF);
    vec2.transformMat3(xy, xy, cameraInvTF);
    vec2.transformMat3(xy, xy, this.modelInvTF);
    return xy;
  }

  mapPointToScreen(xyCell) {
    /*
    Map an XY coordinate from cell/point domain to screen range.  Inverse
    of mapScreenToPoint()
    */

    const responsive = this.responsive;
    const cameraTF = this.camera.view();

    const xy = vec2.transformMat3(vec2.create(), xyCell, this.modelTF);
    vec2.transformMat3(xy, xy, cameraTF);
    vec2.transformMat3(xy, xy, this.projectionTF);

    const pin = [
      Math.round(
        ((xy[0] + 1) * (responsive.width! - this.graphPaddingRightLeft)) / 2
      ),
      Math.round(
        -((xy[1] + 1) / 2 - 1) * (responsive.height! - this.graphPaddingTop)
      )
    ];
    return pin;
  }

  handleBrushDragAction() {
    /*
      event describing brush position:
      @-------|
      |       |
      |       |
      |-------@
    */
    // ignore programatically generated events
    if (d3.event.sourceEvent === null || !d3.event.selection) return;

    const s = d3.event.selection;
    const brushCoords = {
      northwest: this.mapScreenToPoint([s[0][0], s[0][1]]),
      southeast: this.mapScreenToPoint([s[1][0], s[1][1]])
    };

    this.graphContext.actions.graphBrushChange(brushCoords);
  }

  handleBrushStartAction() {
    // Ignore programatically generated events.
    if (!d3.event.sourceEvent) return;

    this.graphContext.actions.graphBrushStart();
  }

  handleBrushEndAction() {
    // Ignore programatically generated events.
    if (!d3.event.sourceEvent) return;

    /*
    coordinates will be included if selection made, null
    if selection cleared.
    */
    const s = d3.event.selection;
    if (s) {
      const brushCoords = {
        northwest: this.mapScreenToPoint(s[0]),
        southeast: this.mapScreenToPoint(s[1])
      };
      this.graphContext.actions.graphBrushEnd(brushCoords);
    } else {
      this.graphContext.actions.graphBrushDeselect();
    }
  }

  handleBrushDeselectAction() {
    this.graphContext.actions.graphBrushDeselect();
  }

  handleLassoStart() {
    this.graphContext.actions.graphLassoStart();
  }

  // when a lasso is completed, filter to the points within the lasso polygon
  handleLassoEnd(polygon) {
    const minimumPolygonArea = 10;

    if (
      polygon.length < 3 ||
      Math.abs(d3.polygonArea(polygon)) < minimumPolygonArea
    ) {
      // if less than three points, or super small area, treat as a clear selection.
      this.graphContext.actions.graphLassoDeselect();
    } else {
      this.graphContext.actions.graphLassoEnd(polygon.map(xy => this.mapScreenToPoint(xy)));
    }
  }

  handleLassoCancel() {
    this.graphContext.actions.graphLassoCancel();
  }

  handleLassoDeselectAction() {
    this.graphContext.actions.graphLassoDeselect();
  }

  handleDeselectAction() {
    const selectionTool = this.selectionTool;
    if (selectionTool === "brush") this.handleBrushDeselectAction();
    if (selectionTool === "lasso") this.handleLassoDeselectAction();
  }

  handleOpacityRangeChange(e) {
    this.graphContext.actions.changeOpacityDeselectedCellsIn2dGraphBackground(e.target.value);
  }

  renderPoints(regl, drawPoints, colorBuffer, pointBuffer, flagBuffer, camera, projectionTF) {
    if (!this.reglCanvas || !this.universe) return;
    const cameraTF = camera.view();
    const projView = mat3.multiply(mat3.create(), projectionTF, cameraTF);
    const { width, height } = this.reglCanvas;
    regl.poll();
    regl.clear({
      depth: 1,
      color: [1, 1, 1, 1],
    });
    drawPoints({
      distance: camera.distance(),
      color: colorBuffer,
      position: pointBuffer,
      flag: flagBuffer,
      count: this.count,
      projView,
      nPoints: this.universe.nObs,
      minViewportDimension: Math.min(width || 800, height || 600),
    });
    regl._gl.flush();
  }

  renderCanvas = renderThrottle(() => {
    this.renderPoints(
      this.regl,
      this.drawPoints,
      this.colorBuffer,
      this.pointBuffer,
      this.flagBuffer,
      this.camera,
      this.projectionTF
    );
  });

  mounted() {
    // setup canvas, webgl draw function and camera
    const camera = _camera(this.reglCanvas);
    const regl = _regl(this.reglCanvas);
    const drawPoints = _drawPoints(regl);

    // preallocate webgl buffers
    const pointBuffer = regl.buffer(0);
    const colorBuffer = regl.buffer(0);
    const flagBuffer = regl.buffer(0);

    // create all default rendering transformations
    const modelTF = createModelTF();
    const projectionTF = createProjectionTF(this.reglCanvas.width, this.reglCanvas.height);

    // initial draw to canvas
    this.renderPoints(regl, drawPoints, colorBuffer, pointBuffer, flagBuffer, camera, projectionTF);

    this.regl = regl;
    this.drawPoints = drawPoints;
    this.pointBuffer = pointBuffer;
    this.colorBuffer = colorBuffer;
    this.flagBuffer = flagBuffer;
    this.camera = camera;
    this.modelTF = modelTF;
    this.modelInvTF = mat3.invert([] as any, modelTF);
    this.projectionTF = projectionTF;
  }

  updated() {
    let stateChanges = {};
    const world = this.worldContext.getters.world;
    const crossfilter = this.crossfilterContext.getters.crossfilter;

    if (this.regl && world && crossfilter) {
      /* update the regl and point rendering state */
      const { obsLayout, nObs } = world;
      const { drawPoints, pointBuffer, colorBuffer, flagBuffer } = this.state;

      let { projectionTF } = this.state;
      let needsRepaint = false;

      if (
        prevProps.responsive.height !== responsive.height ||
        prevProps.responsive.width !== responsive.width
      ) {
        projectionTF = createProjectionTF(
          this.reglCanvas.width,
          this.reglCanvas.height
        );
        needsRepaint = true;
        stateChanges = {
          ...stateChanges,
          projectionTF
        };
      }

      /* coordinates for each point */
      const X = obsLayout.col(layoutChoice.currentDimNames[0]).asArray();
      const Y = obsLayout.col(layoutChoice.currentDimNames[1]).asArray();
      const newPositions = this.computePointPositions(X, Y, modelTF);
      if (renderCache.positions !== newPositions) {
        /* update our cache & GL if the buffer changes */
        renderCache.positions = newPositions;
        pointBuffer({ data: newPositions, dimension: 2 });
        needsRepaint = true;
      }

      /* colors for each point */
      const newColors = this.computePointColors(colorRGB);
      if (renderCache.colors !== newColors) {
        /* update our cache & GL if the buffer changes */
        renderCache.colors = newColors;
        colorBuffer({ data: newColors, dimension: 3 });
        needsRepaint = true;
      }

      /* flags for each point */
      const newFlags = this.computePointFlags(
        world,
        crossfilter,
        colorAccessor,
        pointDilation
      );
      if (renderCache.flags !== newFlags) {
        renderCache.flags = newFlags;
        needsRepaint = true;
        flagBuffer({ data: newFlags, dimension: 1 });
      }

      this.count = nObs;

      if (needsRepaint) {
        this.renderPoints(
          regl,
          drawPoints,
          colorBuffer,
          pointBuffer,
          flagBuffer,
          camera,
          projectionTF
        );
      }
    }

    if (
      prevProps.responsive.height !== responsive.height ||
      prevProps.responsive.width !== responsive.width
    ) {
      // If the window size has changed we want to recreate all SVGs
      stateChanges = {
        ...stateChanges,
        ...this.createToolSVG()
      };
    } else if (
      (responsive.height && responsive.width && !toolSVG) ||
      selectionTool !== prevProps.selectionTool
    ) {
      // first time or change of selection tool
      stateChanges = { ...stateChanges, ...this.createToolSVG() };
    } else if (prevProps.graphInteractionMode !== graphInteractionMode) {
      // If lasso/zoom is switched
      stateChanges = {
        ...stateChanges,
        ...this.createToolSVG()
      };
    }

    /*
    if the selection tool or state has changed, ensure that the selection
    tool correctly reflects the underlying selection.
    */
    if (
      currentSelection !== prevProps.currentSelection ||
      graphInteractionMode !== prevProps.graphInteractionMode ||
      stateChanges.toolSVG
    ) {
      const { tool, container } = this.state;
      this.selectionToolUpdate(
        stateChanges.tool ? stateChanges.tool : tool,
        stateChanges.container ? stateChanges.container : container
      );
    }
    if (Object.keys(stateChanges).length > 0) {
      this.setState(stateChanges);
    }
  }
}
</script>

<style scoped>
:local(.graphCanvas) path {
  shape-rendering: crispEdges;
}

#graphWrapper {
  position: relative;
}

#graphAttachPoint,
.graphSVG,
.graphCanvas {
  position: absolute;
}

#graphSVG {
  pointer-events: none;
}

:local(.axis) path,
:local(.axis) line {
  fill: none;
  stroke: #000;
  stroke-width: 1px;
}

circle {
  stroke-width: 4px;
  stroke: #000;
  fill: none;
}

:local(.hidden) {
  display: none;
}
</style>
