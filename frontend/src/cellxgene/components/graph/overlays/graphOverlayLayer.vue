<template>
  <svg
    class="graphSVG"
    :width="responsive.width - graphPaddingRightLeft"
    :height="responsive.height"
    pointer-events="none"
    :style="svgStyle"
    :onMouseMove="handleCanvasEvent"
    :onWheel="handleCanvasEvent"
  >
    <g
      id="canvas-transformation-group-x"
      :transform="`scale(${responsive.width - graphPaddingRightLeft} 1) scale(.5 1) translate(1 0)`"
    >
      <g
        id="canvas-transformation-group-y"
        :transform="`scale(1 ${-(responsive.height - graphPaddingTop)}) translate(0 -1) scale(1 .5) translate(0 1)`"
      >
        <g id="projection-transformation-group" :transform="matrixToTransformString(projectionTF)">
          <g id="camera-transformation-group" :transform="matrixToTransformString(cameraTF)">
            <g id="model-transformation-group" :transform="matrixToTransformString(modelTF)">
              {{ newChildren }}
            </g>
          </g>
        </g>
      </g>
    </g>
  </svg>
</template>

<script lang="ts">
import { Component, Prop, Vue } from "vue-property-decorator";
import { responsiveModule } from "@/modules/responsive";
import { cloneElement } from "@/utils/reactUtils";

@Component
export default class GraphOverlayLayer extends Vue {
  /*
    This component takes its children (assumed in the data coordinate space ([0, 1] range, origin in bottom left corner))
    and transforms itself multiple times resulting in screen space ([0, screenWidth/Height] range, origin in top left corner)

    Children are assigned in the graph component and must implement onDisplayChange()
   */
  private readonly responsiveContext = responsiveModule.context(this.$store);

  @Prop(Object) cameraTF: any;
  @Prop(Object) modelTF: any;
  @Prop(Object) projectionTF: any;
  @Prop(Number) graphPaddingRightLeft!: number;
  @Prop(Number) graphPaddingTop!: number;
  @Prop(Object) handleCanvasEvent: any;

  private display = {};

  get svgStyle() {
    return {
      zIndex: 99,
      backgroundColor: this.displaying ? "rgba(255, 255, 255, 0.55)" : "",
    };
  }

  get responsive() {
    return this.responsiveContext.getters.responsive;
  }

  displaying = Object.values(this.display).some((value) => value); // check to see if at least one overlay is currently displayed

  matrixToTransformString = (m) => {
    /*
      Translates the gl-matrix mat3 to SVG matrix transform style

                            mat3                    SVG Transform Function
        a  c  e
        b  d  f / [a, b, 0, c, d, 0, e, f, 1] =>  matrix(a, b, c, d, e, f) / matrix(sx, 0, 0, sy, tx, ty) / matrix(m[0] m[3] m[1] m[4] m[6] m[7])
        0  0  1
    */
    return `matrix(${m[0]} ${m[1]} ${m[3]} ${m[4]} ${m[6]} ${m[7]})`;
  };

  reverseMatrixScaleTransformString = (m) => {
    return `matrix(${1 / m[0]} 0 0 ${1 / m[4]} 0 0)`;
  };

  // This is passed to all children, should be called when an overlay's display state is toggled along with the overlay name and its new display state in boolean form
  overlayToggled = (overlay, displaying) => {
    this.display = { ...this.display, [overlay]: displaying };
  };

  inverseTransform = `${this.reverseMatrixScaleTransformString(this.modelTF)} ${this.reverseMatrixScaleTransformString(
    this.cameraTF
  )} ${this.reverseMatrixScaleTransformString(this.projectionTF)} scale(1 2) scale(1 ${
    1 / -(this.responsive.height! - this.graphPaddingTop)
  }) scale(2 1) scale(${1 / (this.responsive.width! - this.graphPaddingRightLeft)} 1)`;

  // Copy the children passed with the overlay and add the inverse transform and onDisplayChange props
  newChildren = this.$children.map((child) =>
    cloneElement(child, {
      inverseTransform: this.inverseTransform,
      overlayToggled: this.overlayToggled,
    })
  );
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
