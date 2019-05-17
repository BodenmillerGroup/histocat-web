<template>
  <div id="map" class="map">
  </div>
</template>

<script lang="ts">
  import { Component, Vue, Watch } from 'vue-property-decorator';
  import Map from 'ol/Map';
  import View from 'ol/View';
  import 'ol/ol.css';
  import { getCenter } from 'ol/extent';
  import ImageLayer from 'ol/layer/Image';
  import Projection from 'ol/proj/Projection';
  import Static from 'ol/source/ImageStatic';
  import { defaults as defaultControls, OverviewMap } from 'ol/control';
  import { readMetalColorMap, readSelectedAcquisition, readSelectedChannels } from '@/modules/experiment/getters';
  import { IAcquisition, IChannel } from '@/modules/experiment/models';

  @Component
  export default class BlendView extends Vue {

    // TODO: check for a better solution
    map!: Map;

    get selectedChannels() {
      return readSelectedChannels(this.$store);
    }

    get metalColorMap() {
      return readMetalColorMap(this.$store);
    }

    get selectedAcquisition() {
      return readSelectedAcquisition(this.$store);
    }

    @Watch('selectedAcquisition')
    onSelectedAcquisitionChanged(acquisition: IAcquisition) {
      if (!acquisition) {
        return;
      }

      // Map views always need a projection.  Here we just want to map image
      // coordinates directly to map coordinates, so we create a projection that uses
      // the image extent in pixels.
      const extent = [0, 0, acquisition.width, acquisition.height];
      const projection = new Projection({
        code: 'pixel',
        units: 'pixels',
        extent: extent,
      });

      const view = new View({
        projection: projection,
        center: getCenter(extent),
        zoom: 2,
        maxZoom: 8,
        enableRotation: false,
      });

      this.map.setView(view);
    }

    @Watch('metalColorMap')
    onMetalColorMapChanged(colorMap: { [metal: string]: string }) {
      this.onSelectedChannelsChanged([]);
    }

    @Watch('selectedChannels')
    onSelectedChannelsChanged(channels: IChannel[]) {
      const view = this.map.getView();
      const projection = view.getProjection();
      const extent = projection.getExtent();
      const colorMap = this.metalColorMap;
      const layers = this.selectedChannels.map((channel) => {
        const color = channel.metal in colorMap ? colorMap[channel.metal] : '';
        return new ImageLayer({
          source: new Static({
            attributions: '© <a href="http://xkcd.com/license.html">xkcd</a>',
            url: `http://localhost/api/v1/channels/${channel.id}/image?color=${color}`,
            projection: projection,
            imageExtent: extent,
          }),
        });
      });

      this.map.getLayers().clear();
      this.map.getLayers().extend(layers);
    }

    mounted() {
      this.map = new Map({
        controls: defaultControls().extend([
          new OverviewMap(),
        ]),
        target: 'map',
      });

      this.map.on('precompose', this.precompose);
    }

    beforeDestroy() {
      this.map.un('precompose', this.precompose);
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
  .map {
    height: 100%;
    width: 100%;
  }
</style>
