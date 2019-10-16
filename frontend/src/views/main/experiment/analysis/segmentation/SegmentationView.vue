<template>
  <div id="segmentation-map"></div>
</template>

<script lang="ts">
import { analysisModule } from "@/modules/analysis";
import { experimentModule } from "@/modules/experiment";
import { IAcquisition } from "@/modules/experiment/models";
import { mainModule } from "@/modules/main";
import { settingsModule } from "@/modules/settings";
import { defaults as defaultControls, FullScreen, ScaleLine } from "ol/control";
import { getCenter } from "ol/extent";
import Feature from "ol/Feature";
import Polygon from "ol/geom/Polygon";
import { DragPan, MouseWheelZoom } from "ol/interaction";
import { Image as ImageLayer, Vector as VectorLayer } from "ol/layer";
import Map from "ol/Map";
import "ol/ol.css";
import Projection from "ol/proj/Projection";
import { Vector as VectorSource, ImageStatic as StaticSource } from "ol/source";
import { Fill, Stroke, Style } from "ol/style";
import View from "ol/View";
import { equals } from "rambda";
import { Component, Vue, Watch } from "vue-property-decorator";

@Component
export default class SegmentationView extends Vue {
  readonly mainContext = mainModule.context(this.$store);
  readonly experimentContext = experimentModule.context(this.$store);
  readonly analysisContext = analysisModule.context(this.$store);
  readonly settingsContext = settingsModule.context(this.$store);

  // TODO: check for a better solution
  map!: Map;
  vectorLayer!: VectorLayer;
  imageLayer!: ImageLayer;
  featureLayer!: VectorLayer;
  highlight: any = undefined;

  get activeAcquisition() {
    return this.experimentContext.getters.activeAcquisition;
  }

  get segmentationImage() {
    return this.analysisContext.getters.segmentationImage;
  }

  get segmentationContours() {
    return this.analysisContext.getters.segmentationContours;
  }

  get showWorkspace() {
    return this.mainContext.getters.showWorkspace;
  }

  @Watch("showWorkspace")
  refreshImageView() {
    if (this.map) {
      this.map.updateSize();
    }
  }

  @Watch("segmentationImage")
  onSegmentationImageChanged(image: string | ArrayBuffer | null) {
    if (!this.map) {
      this.initMap();
    }

    if (image !== null) {
      const projection = this.map.getView().getProjection();
      const source = new StaticSource({
        url: ``,
        imageExtent: projection.getExtent(),
        imageLoadFunction: (view, src: string) => {
          (view.getImage() as any).src = image;
        }
      });
      this.map.getLayers().clear();
      this.imageLayer.setSource(source);
      this.map.getLayers().extend([this.imageLayer]);
    }
  }

  @Watch("segmentationContours")
  onSegmentationContoursChanged(contours: number[][]) {
    if (!this.map) {
      this.initMap();
    }

    if (contours.length > 0) {
      const projection = this.map.getView().getProjection();
      const source = new VectorSource({});
      const features = contours.map((contour: any) => {
        const points = contour.map(p => {
          return p[0];
        });
        return new Feature({
          geometry: new Polygon([points])
        });
      });
      source.addFeatures(features);
      this.map.getLayers().clear();
      this.vectorLayer.setExtent(projection.getExtent());
      this.vectorLayer.setSource(source);
      this.map.getLayers().extend([this.imageLayer, this.vectorLayer, this.featureLayer]);
      this.analysisContext.mutations.setSegmentationImage(null);
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

  mounted() {
    if (this.activeAcquisition) {
      this.onActiveAcquisitionChanged(this.activeAcquisition);
    }
  }

  beforeDestroy() {
    if (this.map) {
      this.map.un("pointermove", this.onPointerMove);
      this.map.un("click", this.onClick);
    }
  }

  displayFeatureInfo(pixel) {
    const feature = this.map.forEachFeatureAtPixel(pixel, f => {
      return f;
    });

    if (feature !== this.highlight) {
      if (this.highlight) {
        this.featureLayer.getSource().removeFeature(this.highlight);
      }
      if (feature) {
        this.featureLayer.getSource().addFeature(feature as any);
      }
      this.highlight = feature;
    }
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

    this.imageLayer = new ImageLayer();
    this.vectorLayer = new VectorLayer();
    this.featureLayer = new VectorLayer({
      source: new VectorSource(),
      map: this.map,
      style: new Style({
        stroke: new Stroke({
          color: "#f00000",
          width: 1
        }),
        fill: new Fill({
          color: "rgba(255, 0, 0, 0.1)"
        })
      })
    });

    this.map = new Map({
      controls: defaultControls({
        zoom: false,
        attribution: false,
        rotate: false
      }).extend([new ScaleLine(), new FullScreen()]),
      interactions: [new DragPan({ kinetic: undefined }), new MouseWheelZoom({ duration: 0 })],
      view: view,
      target: this.$el as HTMLElement
    });

    this.map.on("pointermove", this.onPointerMove);
    this.map.on("click", this.onClick);
  }

  private onPointerMove(ev) {
    if (ev.dragging) {
      return;
    }
    const pixel = this.map.getEventPixel(ev.originalEvent);
    if (pixel) {
      this.displayFeatureInfo(pixel);
    }
  }

  private onClick(ev) {
    this.displayFeatureInfo(ev.pixel);
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
