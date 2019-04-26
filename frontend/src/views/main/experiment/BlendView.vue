<template>
  <v-card flat tile>
    <v-toolbar card dense>
      <v-toolbar-side-icon></v-toolbar-side-icon>
      <v-toolbar-title>Blend</v-toolbar-title>
      <v-spacer/>
      <v-btn icon>
        <v-icon>search</v-icon>
      </v-btn>
      <v-btn icon>
        <v-icon>favorite</v-icon>
      </v-btn>
      <v-btn icon>
        <v-icon>more_vert</v-icon>
      </v-btn>
    </v-toolbar>
    <div>
      <vl-map :load-tiles-while-animating="true" :load-tiles-while-interacting="true"
              data-projection="EPSG:4326">
        <vl-view :zoom.sync="zoom" :center.sync="center" :rotation.sync="rotation"></vl-view>

        <vl-geoloc @update:position="geolocPosition = $event">
          <template slot-scope="geoloc">
            <vl-feature v-if="geoloc.position" id="position-feature">
              <vl-geom-point :coordinates="geoloc.position"></vl-geom-point>
              <vl-style-box>
                <vl-style-icon src="_media/marker.png" :scale="0.4" :anchor="[0.5, 1]"></vl-style-icon>
              </vl-style-box>
            </vl-feature>
          </template>
        </vl-geoloc>

        <vl-layer-tile id="osm">
          <vl-source-osm></vl-source-osm>
        </vl-layer-tile>
      </vl-map>
      <div style="padding: 20px">
        Zoom: {{ zoom }}<br>
        Center: {{ center }}<br>
        Rotation: {{ rotation }}<br>
        My geolocation: {{ geolocPosition }}
      </div>
    </div>
  </v-card>
</template>

<script lang="ts">
  import { Component, Vue } from 'vue-property-decorator';
  import VueLayers from 'vuelayers';
  import 'vuelayers/lib/style.css'; // needs css-loader

  Vue.use(VueLayers);

  @Component
  export default class BlendView extends Vue {
    zoom = 2;
    center = [0, 0];
    rotation = 0;
    geolocPosition = undefined;
  }
</script>
