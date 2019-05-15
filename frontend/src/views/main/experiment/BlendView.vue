<template>
  <v-flex id="map" fill-height>
  </v-flex>
</template>

<script lang="ts">
  import { Component, Vue, Watch } from 'vue-property-decorator';
  import { Map, View } from 'ol';
  import 'ol/ol.css';
  import { getCenter } from 'ol/extent';
  import ImageLayer from 'ol/layer/Image';
  import Projection from 'ol/proj/Projection';
  import Static from 'ol/source/ImageStatic';
  import { defaults as defaultControls, OverviewMap } from 'ol/control.js';
  import { readSelectedAcquisition, readSelectedChannels } from '@/modules/experiment/getters';
  import { IChannel } from '@/modules/experiment/models';

  @Component
  export default class BlendView extends Vue {

    map: Map;

    get selectedChannels() {
      return readSelectedChannels(this.$store);
    }

    @Watch('selectedChannels')
    onSelectedChannelsChanged(channels: IChannel[]) {
      const acquisition = readSelectedAcquisition(this.$store);
      if (!acquisition) {
        return;
      }

      // Map views always need a projection.  Here we just want to map image
      // coordinates directly to map coordinates, so we create a projection that uses
      // the image extent in pixels.
      const extent = [0, 0, acquisition.width, acquisition.height];
      const projection = new Projection({
        code: 'xkcd-image',
        units: 'pixels',
        extent: extent,
      });

      const view = new View({
        projection: projection,
        center: getCenter(extent),
        zoom: 2,
        maxZoom: 8,
      });

      this.map.setView(view);

      const layers = channels.map((channel) => {
        return new ImageLayer({
          source: new Static({
            attributions: '© <a href="http://xkcd.com/license.html">xkcd</a>',
            url: `http://localhost/api/v1/channels/${channel.id}/image`,
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
