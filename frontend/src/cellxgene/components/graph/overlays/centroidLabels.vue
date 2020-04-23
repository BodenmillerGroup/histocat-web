<template>
  <div>
    <g
      v-for="item in items"
      :key="item.label"
      class="centroid-label"
      :transform="`translate(${item.coords[0]}, ${item.coords[1]})`"
      data-testclass="centroid-label"
      :data-testid="`${item.label}-centroid-label`"
    >
      <text
        :transform="inverseTransform"
        text-anchor="middle"
        :data-label="item.label"
        :style="item.style"
        @onMouseEnter="onMouseEnter"
        @onMouseOut="onMouseOut"
        pointer-events="visiblePainted"
      >
        {{ item.displayLabel }}
      </text>
    </g>
  </div>
</template>

<script lang="ts">
import { Component, Prop, Vue, Watch } from "vue-property-decorator";
import { centroidLabelsModule } from "@/modules/centroidLabels";
import { categoryLabelDisplayStringLongLength } from "@/cellxgene/globals";
import { colorsModule } from "@/modules/colors";
import { pointDilationModule } from "@/modules/pointDilation";

@Component
export default class CentroidLabels extends Vue {
  private readonly centroidLabelsContext = centroidLabelsModule.context(this.$store);
  private readonly colorsContext = colorsModule.context(this.$store);
  private readonly pointDilationContext = pointDilationModule.context(this.$store);

  @Prop(String) inverseTransform!: string;
  @Prop(Function) overlayToggled!: (overlay: string, displaying: boolean) => void;

  get labels() {
    return this.centroidLabelsContext.getters.labels;
  }

  get colorAccessor() {
    return this.colorsContext.getters.colorAccessor;
  }

  get dilatedValue() {
    return this.pointDilationContext.getters.categoryField;
  }

  get items() {
    const labelSVGS: any[] = [];
    let fontSize = "15px";
    let fontWeight: any = null;

    this.labels.forEach((coords, label) => {
      fontSize = "15px";
      fontWeight = null;
      if (label === this.dilatedValue) {
        fontSize = "18px";
        fontWeight = "800";
      }

      // Mirror LSB middle truncation
      let displayLabel = label;
      if (displayLabel.length > categoryLabelDisplayStringLongLength) {
        displayLabel = `${label.slice(0, categoryLabelDisplayStringLongLength / 2)}…${label.slice(
          -categoryLabelDisplayStringLongLength / 2
        )}`;
      }

      labelSVGS.push({
        coords: coords,
        label: label,
        displayLabel: displayLabel,
        style: {
          fontFamily: "Roboto Condensed",
          fontSize,
          fontWeight,
          fill: "black",
          userSelect: "none",
        },
      });
    });

    return labelSVGS;
  }

  onMouseEnter(e) {
    // this.dispatch({
    //   type: "category value mouse hover start",
    //   metadataField: this.colorAccessor,
    //   categoryField: e.target.getAttribute("data-label"),
    // });
  }

  onMouseOut(e) {
    // this.dispatch({
    //   type: "category value mouse hover end",
    //   metadataField: this.colorAccessor,
    //   categoryField: e.target.getAttribute("data-label"),
    // });
  }

  @Watch("labels")
  labelsUpdated(newLabels: Map<any, any>) {
    const displayChangeOff = this.labels.size > 0 && newLabels.size == undefined;
    const displayChangeOn = this.labels.size == undefined && newLabels.size > 0;

    if (displayChangeOn || displayChangeOff) {
      // Notify overlay layer of display change
      this.overlayToggled("centroidLabels", displayChangeOn);
    }
  }
}
</script>
