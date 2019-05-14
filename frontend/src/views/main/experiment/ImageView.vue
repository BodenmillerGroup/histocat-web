<template>
  <div id="map"/>
</template>

<script lang="ts">
  import { Component, Vue } from 'vue-property-decorator';
  import { Map, View } from 'ol';
  import 'ol/ol.css';
  import { getCenter } from 'ol/extent';
  import ImageLayer from 'ol/layer/Image';
  import Projection from 'ol/proj/Projection';
  import Static from 'ol/source/ImageStatic';
  import { defaults as defaultControls, OverviewMap } from 'ol/control.js';

  // Map views always need a projection.  Here we just want to map image
  // coordinates directly to map coordinates, so we create a projection that uses
  // the image extent in pixels.
  const extent = [0, 0, 500, 500];
  const projection = new Projection({
    code: 'xkcd-image',
    units: 'pixels',
    extent: extent,
  });

  @Component
  export default class ImageView extends Vue {

    mounted() {

      const layer1 = new ImageLayer({
        source: new Static({
          attributions: '© <a href="http://xkcd.com/license.html">xkcd</a>',
          url: 'http://localhost/api/v1/channels/5/image',
          projection: projection,
          imageExtent: extent,
        }),
      });

      const layer2 = new ImageLayer({
        source: new Static({
          attributions: '© <a href="http://xkcd.com/license.html">xkcd</a>',
          url: 'http://localhost/api/v1/channels/5/image',
          projection: projection,
          imageExtent: extent,
        }),
      });

      const map = new Map({
        layers: [
          layer1,
          layer2,
        ],
        controls: defaultControls().extend([
          new OverviewMap(),
        ]),
        target: 'map',
        view: new View({
          projection: projection,
          center: getCenter(extent),
          zoom: 2,
          maxZoom: 8,
        }),
      });

      // Initially bind listeners
      this.bindLayerListeners(layer1);
      this.bindLayerListeners(layer2);
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
