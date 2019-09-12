<template>
  <v-expansion-panel>
    <v-expansion-panel-header hide-actions>
      <v-container fluid class="ma-0 pa-0">
        <v-row no-gutters>
          <b>{{ label }}</b>
          <v-spacer></v-spacer>
          <v-tooltip bottom>
            <template v-slot:activator="{ on }">
              <v-btn icon x-small color="primary" @click.stop="setSharedChannelLevels" v-on="on">
                <v-icon small>mdi-share</v-icon>
              </v-btn>
            </template>
            <span>Share levels</span>
          </v-tooltip>
          <input type="color" v-model.lazy="color" @click.stop class="ml-1 pa-0" />
        </v-row>
        <v-row no-gutters>
          <v-col>
            <v-range-slider
              v-model="levels"
              :max="channel.max_intensity"
              :min="channel.min_intensity"
              :step="1"
              @click.stop
              @end="submitLimit"
              hide-details
              class="align-center"
            >
              <template v-slot:prepend>
                <v-text-field
                  v-model.number="levels[0]"
                  @change="submitLimit"
                  @click.stop
                  class="mt-0 pt-0 text-input"
                  hide-details
                  single-line
                  type="number"
                  :max="channel.max_intensity"
                  :min="channel.min_intensity"
                  :step="1"
                ></v-text-field>
              </template>
              <template v-slot:append>
                <v-text-field
                  v-model.number="levels[1]"
                  @change="submitLimit"
                  @click.stop
                  class="mt-0 pt-0 text-input"
                  hide-details
                  single-line
                  type="number"
                  :max="channel.max_intensity"
                  :min="channel.min_intensity"
                  :step="1"
                ></v-text-field>
              </template>
            </v-range-slider>
          </v-col>
        </v-row>
      </v-container>
    </v-expansion-panel-header>
    <v-expansion-panel-content>
      <ChannelHistogramView :channel="channel"></ChannelHistogramView>
    </v-expansion-panel-content>
  </v-expansion-panel>
</template>

<script lang="ts">
import { experimentModule } from "@/modules/experiment";
import { IChannel } from "@/modules/experiment/models";
import { settingsModule } from "@/modules/settings";
import ChannelHistogramView from "@/views/main/experiment/settings/channel/ChannelHistogramView.vue";
import { Component, Prop, Vue, Watch } from "vue-property-decorator";

@Component({
  components: { ChannelHistogramView }
})
export default class ChannelSettingsView extends Vue {
  readonly settingsContext = settingsModule.context(this.$store);
  readonly experimentContext = experimentModule.context(this.$store);

  @Prop(Object) channel!: IChannel;

  color = this.channel ? this.metalColor : "#ffffff";
  levels: number[] =
    this.settings && this.settings.levels
      ? [this.settings.levels.min, this.settings.levels.max]
      : [this.channel.min_intensity, this.channel.max_intensity];

  get settings() {
    return this.settingsContext.getters.getChannelSettings(this.channel.id);
  }

  get label() {
    return this.settings && this.settings.customLabel ? this.settings.customLabel : this.channel.label;
  }

  submitLimit() {
    let settings = this.settings;
    if (!settings) {
      settings = {
        id: this.channel.id,
        customLabel: this.channel.label,
        levels: { min: Math.round(this.levels[0]), max: Math.round(this.levels[1]) }
      };
    } else {
      settings = {
        ...settings,
        levels: { min: Math.round(this.levels[0]), max: Math.round(this.levels[1]) }
      };
    }
    this.settingsContext.mutations.setChannelSettings(settings);
    this.experimentContext.actions.getChannelStackImage();
  }

  setSharedChannelLevels() {
    const metal = this.channel.metal;
    const settings = this.settings;
    const levels =
      settings && settings.levels
        ? [settings.levels.min, settings.levels.max]
        : [this.channel.min_intensity, this.channel.max_intensity];
    this.experimentContext.actions.setSharedChannelLevels({ metal: metal, levels: levels });
  }

  get metalColor() {
    const colorMap = this.settingsContext.getters.metalColorMap;
    const color = colorMap.get(this.channel.metal);
    return color ? color : "#ffffff";
  }

  @Watch("color")
  onColorChanged(color: string) {
    this.settingsContext.mutations.setMetalColor({
      metal: this.channel.metal,
      color: color
    });
  }
}
</script>

<style scoped>
.text-input {
  width: 55px;
}
</style>
