<template>
  <div id="blend-map"></div>
</template>

<script lang="ts">
import { analysisModule } from "@/modules/analysis";
import { projectsModule } from "@/modules/projects";
import { IAcquisition } from "@/modules/projects/models";
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
export default class BlendView extends Vue {
  readonly mainContext = mainModule.context(this.$store);
  readonly projectsContext = projectsModule.context(this.$store);
  readonly analysisContext = analysisModule.context(this.$store);
  readonly settingsContext = settingsModule.context(this.$store);

  // TODO: check for a better solution
  map!: Map;
  overviewMap!: OverviewMap;
  vectorSource!: VectorSource;
  vectorLayer!: VectorLayer;
  imageLayer!: ImageLayer;
  imageLayerOverview!: ImageLayer;
  select!: Select;
  draw!: Draw;
  translate!: Translate;

  get activeAcquisition() {
    return this.projectsContext.getters.activeAcquisition;
  }

  get channelStackImage() {
    return this.projectsContext.getters.channelStackImage;
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
      const source = new ImageStatic({
        url: ``,
        imageExtent: projection.getExtent(),
        imageLoadFunction: (view, src: string) => {
          (view.getImage() as any).src = image;
        },
      });

      this.imageLayer.setSource(source);
      this.imageLayerOverview.setSource(source);

      this.map.getLayers().clear();
      const layers = this.regionsEnabled ? [this.imageLayer, this.vectorLayer] : [this.imageLayer];
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

    const existingExtent = this.map.getView().getProjection().getExtent();
    const extent = [0, 0, acquisition.max_x, acquisition.max_y];
    if (!equals(existingExtent, extent)) {
      this.map.getView().getProjection().setExtent(extent);
    }
  }

  mounted() {
    if (this.activeAcquisition) {
      this.onActiveAcquisitionChanged(this.activeAcquisition);
    }
  }

  beforeDestroy() {
    if (this.imageLayer) {
      this.imageLayer.un("prerender", this.prerender);
    }
    if (this.select) {
      this.select.un("select", this.selectHandler);
    }
  }

  deleteRegions() {
    this.select.getFeatures().forEach((feature) => {
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
    const extent = [0, 0, this.activeAcquisition.max_x, this.activeAcquisition.max_y];

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
      },
    });

    const view = new View({
      projection: projection,
      center: getCenter(extent),
      zoom: 4,
      zoomFactor: 1.25,
      maxZoom: 16,
      enableRotation: false,
    });

    this.vectorSource = new VectorSource({ wrapX: false, useSpatialIndex: false });

    this.vectorLayer = new VectorLayer({
      source: this.vectorSource,
    });

    this.select = new Select({
      multi: false,
    });
    this.select.on("select", this.selectHandler);

    this.draw = new Draw({
      source: this.vectorSource,
      type: GeometryType.POLYGON,
      freehand: true,
    });

    this.translate = new Translate({
      features: this.select.getFeatures(),
    });

    this.imageLayer = new ImageLayer();
    this.imageLayer.on("prerender", this.prerender);

    this.imageLayerOverview = new ImageLayer();

    this.overviewMap = new OverviewMap({
      layers: [this.imageLayerOverview],
      view: new View({
        projection: projection,
        center: getCenter(projection.getExtent()),
      }),
    });

    this.map = new Map({
      controls: defaultControls({
        zoom: false,
        attribution: false,
        rotate: false,
      }).extend([new ScaleLine(), new FullScreen(), this.overviewMap]),
      interactions: [new DragPan({ kinetic: undefined }), new MouseWheelZoom({ duration: 0 })],
      view: view,
      target: this.$el as HTMLElement,
      layers: [this.imageLayer],
    });
  }

  private selectHandler(event: SelectEvent) {
    const region = event.selected.length === 0 ? null : event.selected[0].clone();
    this.analysisContext.mutations.setSelectedRegion(region);
  }

  private prerender(evt: RenderEvent) {
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
