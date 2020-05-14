<template>
  <v-expansion-panel>
    <v-expansion-panel-header hide-actions class="ma-0 pa-0">
      <v-container class="ma-0 pa-0">
        <v-row no-gutters class="ml-2 mr-1 mt-1">
          <span class="label">{{ label }}</span>
          <v-spacer />
          <v-tooltip bottom>
            <template v-slot:activator="{ on }">
              <v-btn icon x-small color="primary" @click.stop="setSharedChannelLevels" v-on="on">
                <v-icon small>mdi-share</v-icon>
              </v-btn>
            </template>
            <span>Share levels</span>
          </v-tooltip>
          <input type="color" v-model.lazy="color" @click.stop />
        </v-row>
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
              class="ma-0 pa-0 text-input"
              hide-details
              type="number"
              :max="channel.max_intensity"
              :min="channel.min_intensity"
              :step="1"
              solo
              flat
              dense
            />
          </template>
          <template v-slot:append>
            <v-text-field
              v-model.number="levels[1]"
              @change="submitLimit"
              @click.stop
              class="ma-0 pa-0 text-input"
              hide-details
              type="number"
              :max="channel.max_intensity"
              :min="channel.min_intensity"
              :step="1"
              solo
              flat
              dense
            />
          </template>
        </v-range-slider>
      </v-container>
    </v-expansion-panel-header>
    <v-expansion-panel-content>
      <BrushableHistogram :channel="channel" />
<!--      <ChannelHistogramView :channel="channel" />-->
    </v-expansion-panel-content>
  </v-expansion-panel>
</template>

<script lang="ts">
import { experimentModule } from "@/modules/experiment";
import { IChannel } from "@/modules/experiment/models";
import { settingsModule } from "@/modules/settings";
import ChannelHistogramView from "@/views/main/experiment/image/settings/channel/ChannelHistogramView.vue";
import { Component, Prop, Vue, Watch } from "vue-property-decorator";
import BrushableHistogram from "@/components/BrushableHistogram.vue";

@Component({
  components: {BrushableHistogram, ChannelHistogramView },
})
export default class ChannelSettingsView extends Vue {
  readonly settingsContext = settingsModule.context(this.$store);
  readonly experimentContext = experimentModule.context(this.$store);

  @Prop(Object) readonly channel!: IChannel;

  color = this.channel ? this.metalColor : "#ffffff";
  levels: number[] =
    this.settings && this.settings.levels
      ? [this.settings.levels.min, this.settings.levels.max]
      : [this.channel.min_intensity, this.channel.max_intensity];

  get activeAcquisitionId() {
    return this.experimentContext.getters.activeAcquisitionId;
  }

  get settings() {
    return this.settingsContext.getters.getChannelSettings(this.activeAcquisitionId, this.channel.name);
  }

  get label() {
    return this.settings && this.settings.customLabel ? this.settings.customLabel : this.channel.label;
  }

  submitLimit() {
    if (!this.activeAcquisitionId) {
      return;
    }
    let settings = this.settings;
    if (!settings) {
      settings = {
        acquisitionId: this.activeAcquisitionId,
        name: this.channel.name,
        customLabel: this.channel.label,
        levels: { min: Math.round(this.levels[0]), max: Math.round(this.levels[1]) },
        suppressBroadcast: false,
      };
    } else {
      settings = {
        ...settings,
        levels: { min: Math.round(this.levels[0]), max: Math.round(this.levels[1]) },
        suppressBroadcast: false,
      };
    }
    this.settingsContext.mutations.setChannelSettings(settings);
    this.experimentContext.actions.getChannelStackImage();
  }

  setSharedChannelLevels() {
    const metal = this.channel.name;
    const settings = this.settings;
    const levels =
      settings && settings.levels
        ? [settings.levels.min, settings.levels.max]
        : [this.channel.min_intensity, this.channel.max_intensity];
    this.experimentContext.actions.setSharedChannelLevels({ metal: metal, levels: levels });
  }

  get metalColor() {
    const colorMap = this.settingsContext.getters.metalColorMap;
    const color = colorMap.get(this.channel.name);
    return color ? color : "#ffffff";
  }

  @Watch("color")
  onColorChanged(color: string) {
    this.settingsContext.mutations.setMetalColor({
      metal: this.channel.name,
      color: color,
    });
  }
}
</script>

<style scoped>
.text-input {
  width: 70px;
}

.label {
  font-size: 10pt;
  font-weight: 500;
}
</style>
