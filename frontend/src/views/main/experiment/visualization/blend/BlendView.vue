<template>
  <div id="blend-map"></div>
</template>

<script lang="ts">
import { analysisModule } from "@/modules/analysis";
import { experimentModule } from "@/modules/experiment";
import { IAcquisition, IChannel } from "@/modules/experiment/models";
import { mainModule } from "@/modules/main";
import { settingsModule } from "@/modules/settings";
import { defaults as defaultControls, FullScreen, OverviewMap, ScaleLine } from "ol/control";
import { getCenter } from "ol/extent";
import GeometryType from "ol/geom/GeometryType";
import { DragPan, Draw, MouseWheelZoom, Select, Translate } from "ol/interaction";
import { SelectEvent } from "ol/interaction/Select";
import { Image as ImageLayer, Vector as VectorLayer } from "ol/layer";
import Map from "ol/Map";
import "ol/ol.css";
import Projection from "ol/proj/Projection";
import RenderEvent from "ol/render/Event";
import { ImageStatic, Vector as VectorSource } from "ol/source";
import View from "ol/View";
import { equals } from "rambda";
import { Component, Vue, Watch } from "vue-property-decorator";

@Component
export default class VisualizationView extends Vue {
  readonly mainContext = mainModule.context(this.$store);
  readonly experimentContext = experimentModule.context(this.$store);
  readonly analysisContext = analysisModule.context(this.$store);
  readonly settingsContext = settingsModule.context(this.$store);

  // TODO: check for a better solution
  map!: Map;
  overviewMap!: OverviewMap;
  vectorSource!: VectorSource;
  vectorLayer!: VectorLayer;
  select!: Select;
  draw!: Draw;
  translate!: Translate;

  get selectedChannels() {
    return this.experimentContext.getters.selectedChannels;
  }

  get metalColorMap() {
    return this.settingsContext.getters.metalColorMap;
  }

  get activeAcquisition() {
    return this.experimentContext.getters.activeAcquisition;
  }

  get channelStackImage() {
    return this.experimentContext.getters.channelStackImage;
  }

  get regionsEnabled() {
    return this.analysisContext.getters.regionsEnabled;
  }

  get showWorkspace() {
    return this.mainContext.getters.showWorkspace;
  }

  get showOptions() {
    return this.mainContext.getters.showOptions;
  }

  @Watch("regionsEnabled")
  regionsEnabledChanged(state: boolean) {
    if (this.map) {
      if (state) {
        this.map.getInteractions().extend([this.draw, this.select, this.translate]);
        this.map.getLayers().extend([this.vectorLayer]);
      } else {
        this.map.getInteractions().remove(this.select);
        this.map.getInteractions().remove(this.translate);
        this.map.getInteractions().remove(this.draw);
        this.map.getLayers().remove(this.vectorLayer);
      }
      this.map.updateSize();
    }
  }

  @Watch("showWorkspace")
  @Watch("showOptions")
  refreshImageView() {
    if (this.map) {
      this.map.updateSize();
    }
  }

  @Watch("channelStackImage")
  onChannelStackImageChanged(image: string | ArrayBuffer | null) {
    if (!this.map) {
      this.initMap();
    }

    if (image !== null) {
      const projection = this.map.getView().getProjection();
      const imageLayer = new ImageLayer({
        source: new ImageStatic({
          url: ``,
          imageExtent: projection.getExtent(),
          imageLoadFunction: (view, src: string) => {
            (view.getImage() as any).src = image;
          }
        })
      });

      this.map.getLayers().clear();
      const layers = this.regionsEnabled ? [imageLayer, this.vectorLayer] : [imageLayer];
      this.map.getLayers().extend(layers);
    }
  }

  @Watch("activeAcquisition")
  onActiveAcquisitionChanged(acquisition: IAcquisition) {
    if (!acquisition) {
      return;
    }

    if (!this.map) {
      this.initMap();
    }

    const existingExtent = this.map
      .getView()
      .getProjection()
      .getExtent();
    const extent = [0, 0, parseInt(acquisition.meta.MaxX, 10), parseInt(acquisition.meta.MaxY, 10)];
    if (!equals(existingExtent, extent)) {
      this.map
        .getView()
        .getProjection()
        .setExtent(extent);
    }
  }

  @Watch("metalColorMap")
  onMetalColorMapChanged(colorMap: { [metal: string]: string }) {
    this.experimentContext.actions.getChannelStackImage();
  }

  @Watch("selectedChannels")
  onSelectedChannelsChanged(channels: IChannel[]) {
    if (!this.map) {
      this.initMap();
    }

    if (channels && channels.length > 0) {
      this.experimentContext.actions.getChannelStackImage();
    } else {
      this.experimentContext.mutations.setChannelStackImage(null);
      if (this.map) {
        this.map.getLayers().clear();
      }
    }
  }

  mounted() {
    if (this.activeAcquisition) {
      this.onActiveAcquisitionChanged(this.activeAcquisition);
    }
  }

  beforeDestroy() {
    if (this.map) {
      this.map.un("precompose", this.precompose);
      this.select.un("select", this.selectHandler);
    }
  }

  deleteRegions() {
    this.select.getFeatures().forEach(feature => {
      this.vectorSource.removeFeature(feature);
    });
    this.select.getFeatures().clear();
    this.analysisContext.mutations.setSelectedRegion(null);
    this.analysisContext.mutations.setSelectedRegionStats([]);
  }

  private initMap() {
    if (!this.activeAcquisition) {
      return;
    }
    const extent = [
      0,
      0,
      parseInt(this.activeAcquisition.meta.MaxX, 10),
      parseInt(this.activeAcquisition.meta.MaxY, 10)
    ];

    // Map views always need a projection.  Here we just want to map image
    // coordinates directly to map coordinates, so we create a projection that uses
    // the image extent in pixels.
    const projection = new Projection({
      code: "NONE",
      units: "pixels",
      extent: extent,
      getPointResolution: (pixelRes, point) => {
        const scale = this.settingsContext.getters.scalebar.settings.scale
          ? this.settingsContext.getters.scalebar.settings.scale
          : 1.0;
        const spacing = 0.000001 * scale;
        const res = pixelRes * spacing;
        return res;
      }
    });

    const view = new View({
      projection: projection,
      center: getCenter(extent),
      zoom: 4,
      zoomFactor: 1.25,
      maxZoom: 16,
      enableRotation: false
    });

    this.overviewMap = new OverviewMap({
      view: new View({
        projection: projection
      })
    });

    this.vectorSource = new VectorSource({ wrapX: false, useSpatialIndex: false });

    this.vectorLayer = new VectorLayer({
      source: this.vectorSource
    });

    this.select = new Select({
      multi: false
    });
    this.select.on("select", this.selectHandler);

    this.draw = new Draw({
      source: this.vectorSource,
      type: GeometryType.POLYGON,
      freehand: true
    });

    this.translate = new Translate({
      features: this.select.getFeatures()
    });

    this.map = new Map({
      controls: defaultControls({
        zoom: false,
        attribution: false,
        rotate: false
      }).extend([new ScaleLine(), new FullScreen(), this.overviewMap]),
      interactions: [new DragPan({ kinetic: undefined }), new MouseWheelZoom({ duration: 0 })],
      view: view,
      target: this.$el as HTMLElement
    });

    this.map.on("precompose", this.precompose);
  }

  private selectHandler(event: SelectEvent) {
    const region = event.selected.length === 0 ? null : event.selected[0].clone();
    this.analysisContext.mutations.setSelectedRegion(region);
  }

  private precompose(evt: RenderEvent) {
    evt.context.imageSmoothingEnabled = false;
  }
}
</script>

<style>
.ol-scale-line {
  bottom: 10px;
  left: auto;
  right: 10px;
  top: auto;
}
</style>
