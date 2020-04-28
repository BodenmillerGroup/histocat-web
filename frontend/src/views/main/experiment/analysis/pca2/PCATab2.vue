<template>
  <v-banner v-if="!activeDataset" icon="mdi-alert-circle-outline">
    Please select dataset
  </v-banner>
  <v-banner v-else-if="!activeAcquisition && selectedAcquisitionIds.length === 0" icon="mdi-alert-circle-outline">
    Please select acquisition(s)
  </v-banner>
  <v-row v-else no-gutters class="chart-container">
    <v-col :cols="columns">
      <canvas ref="canvas" :width="responsive.width - 400" :height="responsive.height - 100" />
    </v-col>
    <v-col v-if="showOptions" cols="3">
      <v-card tile>
        <v-card-title>PCA Settings</v-card-title>
        <v-card-text>
          <v-chip-group v-model="selectedChannels" multiple column active-class="primary--text">
            <v-chip v-for="item in channels" :key="item" :value="item" small>
              {{ item }}
            </v-chip>
          </v-chip-group>
          <v-card-actions>
            <v-btn @click="selectAll" small :disabled="selectedChannels.length === channels.length">
              Select all
            </v-btn>
            <v-btn @click="clearAll" small :disabled="selectedChannels.length === 0">
              Clear all
            </v-btn>
          </v-card-actions>
          <v-radio-group v-model="nComponents" mandatory hide-details label="Dimensions">
            <v-radio label="2D" value="2" />
            <v-radio label="3D" value="3" />
          </v-radio-group>
          <v-select
            :items="heatmaps"
            v-model="heatmap"
            label="Heatmap"
            hint="Heatmap marker"
            item-text="label"
            return-object
            persistent-hint
            clearable
            dense
            class="mt-5"
          />
        </v-card-text>
        <v-card-actions>
          <v-btn @click="submit" color="primary" block :disabled="selectedChannels.length === 0">
            Analyze
          </v-btn>
        </v-card-actions>
      </v-card>
    </v-col>
  </v-row>
</template>

<script lang="ts">
import { analysisModule } from "@/modules/analysis";
import { IPCAData } from "@/modules/analysis/models";
import { datasetModule } from "@/modules/datasets";
import { experimentModule } from "@/modules/experiment";
import { mainModule } from "@/modules/main";
import { settingsModule } from "@/modules/settings";
import createScatterplot from "regl-scatterplot";
import { Component, Vue, Watch } from "vue-property-decorator";
import { responsiveModule } from "@/modules/responsive";

@Component
export default class PCATab2 extends Vue {
  readonly mainContext = mainModule.context(this.$store);
  readonly experimentContext = experimentModule.context(this.$store);
  readonly datasetContext = datasetModule.context(this.$store);
  readonly analysisContext = analysisModule.context(this.$store);
  readonly settingsContext = settingsModule.context(this.$store);
  readonly responsiveContext = responsiveModule.context(this.$store);

  selectedChannels: any[] = [];
  nComponents = "2";
  heatmap: { type: string; label: string } | null = null;

  points: any[] = [];
  scatterplot: any;
  selection: any[] = [];

  get responsive() {
    return this.responsiveContext.getters.responsive;
  }

  get showOptions() {
    return this.mainContext.getters.showOptions;
  }

  get columns() {
    return this.showOptions ? 9 : 12;
  }

  get activeAcquisition() {
    return this.experimentContext.getters.activeAcquisition;
  }

  get selectedAcquisitionIds() {
    return this.experimentContext.getters.selectedAcquisitionIds;
  }

  get activeDataset() {
    return this.datasetContext.getters.activeDataset;
  }

  get channels() {
    return this.datasetContext.getters.channels;
  }

  get heatmaps() {
    return this.datasetContext.getters.heatmaps;
  }

  get applyMask() {
    return this.settingsContext.getters.maskSettings.apply;
  }

  selectAll() {
    this.selectedChannels = this.channels;
  }

  clearAll() {
    this.selectedChannels = [];
  }

  async submit() {
    let heatmap = "";
    if (this.heatmap) {
      heatmap = this.heatmap.type === "channel" ? this.heatmap.label : `Neighbors_${this.heatmap.label}`;
    }

    const acquisitionIds =
      this.selectedAcquisitionIds.length > 0 ? this.selectedAcquisitionIds : [this.activeAcquisition!.id];

    await this.analysisContext.actions.getPCAData({
      dataset_id: this.activeDataset!.id,
      acquisition_ids: acquisitionIds,
      n_components: parseInt(this.nComponents, 10),
      heatmapType: this.heatmap ? this.heatmap.type : "",
      heatmap: heatmap,
      markers: this.selectedChannels,
    });
  }

  get pcaData() {
    return this.analysisContext.getters.pcaData;
  }

  @Watch("pcaData")
  pcaDataChanged(data: IPCAData) {
    if (data) {
      this.plot2D(data);
    }
  }

  pointoverHandler(pointId) {
    const [x, y, category, value] = this.points[pointId];
    // console.log(`X: ${x}\nY: ${y}\nCategory: ${category}\nValue: ${value}`);
  }

  pointoutHandler(pointId) {
    const [x, y, category, value] = this.points[pointId];
    // console.log(`X: ${x}\nY: ${y}\nCategory: ${category}\nValue: ${value}`);
  }

  selectHandler({ points: selectedPoints }) {
    console.log("Selected:", selectedPoints);
    this.selection = selectedPoints;
    if (this.selection.length > 0) {
      const cell_ids: number[] = [];
      for (const i of this.selection) {
        const cell_id = Number(this.points[i][2].split("_")[1]);
        cell_ids.push(cell_id);
      }
      if (this.applyMask) {
        // this.settingsContext.mutations.setMaskSettings({
        //   ...this.settingsContext.getters.maskSettings,
        //   cell_ids: value,
        // });
        this.experimentContext.actions.getGatedMaskImage(cell_ids);
        // this.experimentContext.actions.getChannelStackImage();
      }
    }
  }

  deselectHandler() {
    console.log("Deselected:", this.selection);
    this.selection = [];
  }

  resizeHandler() {
    const canvas = this.$refs.canvas as any;
    if (canvas) {
      const rect = canvas.getBoundingClientRect();
      if (rect) {
        console.log(rect);
        this.scatterplot.set({ width: rect.width, height: rect.height });
      }
    }
  }

  perc2color(perc) {
    var r,
      g,
      b = 0;
    if (perc < 50) {
      r = 255;
      g = Math.round(5.1 * perc);
    } else {
      g = 255;
      r = Math.round(510 - 5.1 * perc);
    }
    var h = r * 0x10000 + g * 0x100 + b * 0x1;
    return "#" + ("000000" + h.toString(16)).slice(-6);
  }

  private plot2D(data: IPCAData) {
    this.scatterplot.deselect();
    this.scatterplot.reset();
    this.points = data.heatmap
      ? data.x.data.map((x, i) => {
          const cell_id = Number(data.cell_ids[i].split("_")[1]);
          return [x, data.y.data[i], data.heatmap!.data[i], cell_id];
        })
      : data.x.data.map((x, i) => {
          return [x, data.y.data[i], 0, data.cell_ids[i]];
        });

    // const colorsScale = data.heatmap!.data.map((item) => this.perc2color(item));
    this.scatterplot.set({
      colorBy: "category",
      colors: [
        "#002072",
        "#162b79",
        "#233680",
        "#2e4186",
        "#394d8d",
        "#425894",
        "#4b649a",
        "#5570a1",
        "#5e7ca7",
        "#6789ae",
        "#7195b4",
        "#7ba2ba",
        "#85aec0",
        "#90bbc6",
        "#9cc7cc",
        "#a9d4d2",
        "#b8e0d7",
        "#c8ecdc",
        "#ddf7df",
        "#ffffe0",
      ],
    });

    this.scatterplot.draw(this.points);
  }

  mounted() {
    const canvas = this.$refs.canvas as any;

    const { width, height } = canvas.getBoundingClientRect();

    this.scatterplot = createScatterplot({
      canvas,
      width,
      height,
      opacity: 1,
      pointSize: 2,
      pointSizeSelected: 1,
      pointOutlineWidth: 1,
      lassoMinDelay: 15,
    });

    this.scatterplot.subscribe("pointover", this.pointoverHandler);
    this.scatterplot.subscribe("pointout", this.pointoutHandler);
    this.scatterplot.subscribe("select", this.selectHandler);
    this.scatterplot.subscribe("deselect", this.deselectHandler);

    window.addEventListener("resize", this.resizeHandler);
  }

  beforeDestroy() {
    this.scatterplot.destroy();
  }
}
</script>

<style scoped>
.chart-container {
  height: calc(100vh - 154px);
}
</style>
