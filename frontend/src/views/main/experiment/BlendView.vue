<template>
  <div id="map" class="blend-view">
  </div>
</template>

<script lang="ts">
  import { apiUrl } from '@/env';
  import { experimentModule } from '@/modules/experiment';
  import { IAcquisition, IChannel } from '@/modules/experiment/models';
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
  export default class BlendView extends Vue {
    readonly experimentContext = experimentModule.context(this.$store);
    readonly settingsContext = settingsModule.context(this.$store);

    // TODO: check for a better solution
    map!: Map;

    get selectedChannels() {
      return this.experimentContext.getters.selectedChannels;
    }

    get metalColorMap() {
      return this.settingsContext.getters.metalColorMap;
    }

    get activeSlideId() {
      return this.experimentContext.getters.activeSlideId;
    }

    get activePanoramaId() {
      return this.experimentContext.getters.activePanoramaId;
    }

    get activeAcquisition() {
      return this.experimentContext.getters.activeAcquisition;
    }

    @Watch('activeSlideId')
    onActiveSlideIdChanged(id?: number) {
      if (!id) {
        return;
      }
      const layers = [
        new ImageLayer({
          source: new Static({
            imageLoadFunction: (image, src) => {
              (image.getImage() as any).src = src;
            },
            url: `${apiUrl}/api/v1/slides/${id}/image`,
            imageExtent: [0, 0, 800, 600],
          }),
        }),
      ];

      this.map.getLayers().clear();
      this.map.getLayers().extend(layers);
    }

    @Watch('activePanoramaId')
    onActivePanoramaIdChanged(id?: number) {
      if (!id) {
        return;
      }
      const layers = [
        new ImageLayer({
          source: new Static({
            imageLoadFunction: (image, src) => {
              (image.getImage() as any).src = src;
            },
            url: `${apiUrl}/api/v1/panoramas/${id}/image`,
            imageExtent: [0, 0, 800, 600],
          }),
        }),
      ];

      this.map.getLayers().clear();
      this.map.getLayers().extend(layers);
    }

    @Watch('activeAcquisition')
    onActiveAcquisitionChanged(acquisition: IAcquisition) {
      if (!acquisition) {
        return;
      }

      if (!this.map) {
        // Map views always need a projection.  Here we just want to map image
        // coordinates directly to map coordinates, so we create a projection that uses
        // the image extent in pixels.
        const extent = [0, 0, acquisition.max_x, acquisition.max_y];

        const projection = new Projection({
          code: 'NONE',
          units: 'pixels',
          extent: extent,
          getPointResolution: (pixelRes, point) => {
            /*
             * DICOM pixel spacing has millimeter unit while the projection has has meter unit.
             */
            const spacing = 0.0001 / 10 ** 3;
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

        this.map = new Map({
          controls: defaultControls().extend([
            new OverviewMap({
              view: new View({
                projection: projection,
              }),
            }),
            new ScaleLine(),
            new FullScreen(),
          ]),
          interactions: [
            new DragPan({ kinetic: undefined }),
            new MouseWheelZoom({ duration: 0 }),
          ],
          view: view,
          target: 'map',
        });

        this.map.on('precompose', this.precompose);
      }

      const existingExtent = this.map.getView().getProjection().getExtent();
      const newExtent = [0, 0, acquisition.max_x, acquisition.max_y];
      if (!equals(existingExtent, newExtent)) {
        this.map.getView().getProjection().setExtent(newExtent);
      }
    }

    @Watch('metalColorMap')
    onMetalColorMapChanged(colorMap: { [metal: string]: string }) {
      this.onSelectedChannelsChanged([]);
    }

    @Watch('selectedChannels')
    onSelectedChannelsChanged(channels: IChannel[]) {
      if (!this.selectedChannels) {
        return;
      }
      const view = this.map.getView();
      const projection = view.getProjection();
      const extent = projection.getExtent();
      const colorMap = this.metalColorMap;
      const layers = this.selectedChannels.map((channel) => {
        const color = colorMap.has(channel.metal) ? colorMap.get(channel.metal) : '';
        const channelSettings = this.settingsContext.getters.channelSettings(channel.id);
        const min = channelSettings && channelSettings.levels ? channelSettings.levels.min : '';
        const max = channelSettings && channelSettings.levels ? channelSettings.levels.max : '';
        return new ImageLayer({
          source: new Static({
            imageLoadFunction: (image, src) => {
              (image.getImage() as any).src = src;
            },
            url: `${apiUrl}/api/v1/channels/${channel.id}/image?color=${color}&min=${min}&max=${max}`,
            projection: projection,
            imageExtent: extent,
          }),
        });
      });

      this.map.getLayers().clear();
      this.map.getLayers().extend(layers);
    }

    beforeDestroy() {
      if (this.map) {
        this.map.un('precompose', this.precompose);
      }
    }

    private precompose(evt) {
      evt.context.imageSmoothingEnabled = false;
      evt.context.webkitImageSmoothingEnabled = false;
      evt.context.mozImageSmoothingEnabled = false;
      evt.context.msImageSmoothingEnabled = false;
      evt.context.globalCompositeOperation = 'screen';
    };
  }
</script>

<style scoped>
  .blend-view {
    height: 95%;
    /*width: 100%;*/
  }
</style>

<style>
  .ol-scale-line {
    bottom: 8px;
    left: 12em;
    padding: 2px;
    position: absolute;
  }
</style>
