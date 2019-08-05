<template>
  <div id="segmentation-map"></div>
</template>

<script lang="ts">
  import { analysisModule } from '@/modules/analysis';
  import { experimentModule } from '@/modules/experiment';
  import { IAcquisition, IChannel } from '@/modules/experiment/models';
  import { mainModule } from '@/modules/main';
  import { settingsModule } from '@/modules/settings';
  import { defaults as defaultControls, FullScreen, OverviewMap, ScaleLine } from 'ol/control';
  import { getCenter } from 'ol/extent';
  import { DragPan, MouseWheelZoom } from 'ol/interaction';
  import ImageLayer from 'ol/layer/Image';
  import Map from 'ol/Map';
  import 'ol/ol.css';
  import Projection from 'ol/proj/Projection';
  import Static from 'ol/source/ImageStatic';
  import View from 'ol/View';
  import { equals } from 'ramda';
  import { Component, Vue, Watch } from 'vue-property-decorator';

  @Component
  export default class SegmentationView extends Vue {
    readonly mainContext = mainModule.context(this.$store);
    readonly experimentContext = experimentModule.context(this.$store);
    readonly analysisContext = analysisModule.context(this.$store);
    readonly settingsContext = settingsModule.context(this.$store);

    get selectedChannels() {
      return this.experimentContext.getters.selectedChannels;
    }

    get activeAcquisition() {
      return this.experimentContext.getters.activeAcquisition;
    }

    get analysisImage() {
      return this.analysisContext.getters.analysisImage;
    }

    get showWorkspace() {
      return this.mainContext.getters.showWorkspace;
    }

    get showChannels() {
      return this.mainContext.getters.showChannels;
    }

    // TODO: check for a better solution
    map!: Map;
    overviewMap!: OverviewMap;

    @Watch('showWorkspace')
    onRefreshImageView() {
      if (this.map) {
        this.map.updateSize();
      }
    }

    @Watch('analysisImage')
    onChannelStackImageChanged(image: string | ArrayBuffer | null) {
      if (!this.map) {
        this.initMap();
      }

      if (image !== null) {
        const projection = this.map.getView().getProjection();
        const layer = new ImageLayer({
          source: new Static({
            url: ``,
            imageExtent: projection.getExtent(),
            imageLoadFunction: (view, src: string) => {
              (view.getImage() as any).src = image;
            },
          }),
        });
        this.map.getLayers().clear();
        this.map.getLayers().extend([layer]);
      }
    }

    @Watch('activeAcquisition')
    onActiveAcquisitionChanged(acquisition: IAcquisition) {
      if (!acquisition) {
        return;
      }

      if (!this.map) {
        this.initMap();
      }

      const existingExtent = this.map.getView().getProjection().getExtent();
      const extent = [0, 0, parseInt(acquisition.meta.MaxX, 10), parseInt(acquisition.meta.MaxY, 10)];
      if (!equals(existingExtent, extent)) {
        this.map.getView().getProjection().setExtent(extent);
      }
    }

    @Watch('selectedChannels')
    onSelectedChannelsChanged(channels: IChannel[]) {
      if (!this.map) {
        this.initMap();
      }

      if (!channels || channels.length === 0) {
        this.experimentContext.mutations.setChannelStackImage(null);
        this.map.getLayers().clear();
      }
    }

    mounted() {
      if (this.activeAcquisition) {
        this.onActiveAcquisitionChanged(this.activeAcquisition);
      }
    }

    private initMap() {
      if (!this.activeAcquisition) {
        return;
      }
      const extent = [0, 0, parseInt(this.activeAcquisition.meta.MaxX, 10), parseInt(this.activeAcquisition.meta.MaxY, 10)];

      // Map views always need a projection.  Here we just want to map image
      // coordinates directly to map coordinates, so we create a projection that uses
      // the image extent in pixels.
      const projection = new Projection({
        code: 'NONE',
        units: 'pixels',
        extent: extent,
        getPointResolution: (pixelRes, point) => {
          const scale = this.settingsContext.getters.scalebar.settings.scale ?
            this.settingsContext.getters.scalebar.settings.scale : 1.0;
          const spacing = 0.000001 * scale;
          const res = pixelRes * spacing;
          return (res);
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

      this.overviewMap = new OverviewMap({
        view: new View({
          projection: projection,
        }),
      });
      this.map = new Map({
        controls: defaultControls({
          zoom: false,
          attribution: false,
          rotate: false,
        }).extend([
          new ScaleLine(),
          new FullScreen(),
          this.overviewMap,
        ]),
        interactions: [
          new DragPan({ kinetic: undefined }),
          new MouseWheelZoom({ duration: 0 }),
        ],
        view: view,
        target: 'segmentation-map',
      });
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
