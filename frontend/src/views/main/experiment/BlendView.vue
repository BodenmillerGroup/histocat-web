<template>
  <div id="map">
  </div>
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

      const existingLayers = this.map.getLayers();
      existingLayers.forEach((layer) => {
        this.unbindLayerListeners(layer);
      });
      existingLayers.clear();

      for (const layer of layers) {
        this.map.addLayer(layer);
        // Initially bind listeners
        this.bindLayerListeners(layer);
      }
    }

    mounted() {
      this.map = new Map({
        controls: defaultControls().extend([
          new OverviewMap(),
        ]),
        target: 'map',
      });
    }

    beforeDestroy() {
      // Finally unbind listeners
      this.map.getLayers().forEach((layer) => {
        this.unbindLayerListeners(layer);
      });
    }

    /**
     * This method sets the globalCompositeOperation to the value of the select
     * field and it is bound to the precompose event of the layers.
     *
     * @param {module:ol/render/Event~RenderEvent} evt The render event.
     */
    private setBlendModeFromSelect(evt) {
      evt.context.globalCompositeOperation = 'screen';
    };


    /**
     * This method resets the globalCompositeOperation to the default value of
     * 'source-over' and it is bound to the postcompose event of the layers.
     *
     * @param {module:ol/render/Event~RenderEvent} evt The render event.
     */
    private resetBlendModeFromSelect(evt) {
      evt.context.globalCompositeOperation = 'source-over';
    };


    /**
     * Bind the pre- and postcompose handlers to the passed layer.
     *
     * @param {module:ol/layer/Vector} layer The layer to bind the handlers to.
     */
    private bindLayerListeners(layer) {
      layer.on('precompose', this.setBlendModeFromSelect);
      layer.on('postcompose', this.resetBlendModeFromSelect);
    };


    /**
     * Unind the pre- and postcompose handlers to the passed layers.
     *
     * @param {module:ol/layer/Vector} layer The layer to unbind the handlers from.
     */
    private unbindLayerListeners(layer) {
      layer.un('precompose', this.setBlendModeFromSelect);
      layer.un('postcompose', this.resetBlendModeFromSelect);
    };
  }
</script>
