<template>
  <v-expansion-panel>
    <v-expansion-panel-header v-slot="{ open }" hide-actions>
      <v-layout column>
        <v-layout justify-space-between>
          <b>{{ label }}</b>
          <input
            type="color"
            v-model.lazy="color"
            @click.stop
          />
        </v-layout>
        <v-layout justify-space-between>
          <v-range-slider
            :value="levels"
            :max="channel.max_intensity"
            :min="channel.min_intensity"
            :step="1"
            thumb-label
            :thumb-size="24"
            @click.stop
            @end="submitLimit"
            hide-details
          ></v-range-slider>
          <v-tooltip bottom>
            <template v-slot:activator="{ on }">
              <v-btn class="ma-2" outlined fab x-small color="blue" @click.stop="setSharedChannelLevels" v-on="on">
                <v-icon small>mdi-share</v-icon>
              </v-btn>
            </template>
            <span>Share levels</span>
          </v-tooltip>
        </v-layout>
      </v-layout>
    </v-expansion-panel-header>
    <v-expansion-panel-content>
      <ChannelHistogramView :channel="channel"></ChannelHistogramView>
    </v-expansion-panel-content>
  </v-expansion-panel>
</template>

<script lang="ts">
  import { experimentModule } from '@/modules/experiment';
  import { IChannel } from '@/modules/experiment/models';
  import { settingsModule } from '@/modules/settings';
  import ChannelHistogramView from '@/views/main/experiment/settings/ChannelHistogramView.vue';
  import { Component, Prop, Vue, Watch } from 'vue-property-decorator';

  @Component({
    components: { ChannelHistogramView },
  })
  export default class ChannelSettingsView extends Vue {
    readonly settingsContext = settingsModule.context(this.$store);
    readonly experimentContext = experimentModule.context(this.$store);

    @Prop(Object) channel!: IChannel;

    color = this.channel ? this.metalColor : '#ffffff';

    get settings() {
      return this.settingsContext.getters.channelSettings(this.channel.id);
    }

    get label() {
      return this.settings && this.settings.customLabel ?
        this.settings.customLabel :
        this.channel.label;
    }

    get levels() {
      if (this.settings && this.settings.levels) {
        return [this.settings.levels.min, this.settings.levels.max];
      } else {
        return [this.channel.min_intensity, this.channel.max_intensity];
      }
    }

    submitLimit(range: number[]) {
      let settings = this.settingsContext.getters.channelSettings(this.channel.id);
      if (!settings) {
        settings = {
          id: this.channel.id,
          customLabel: this.channel.label,
          levels: { min: Math.round(range[0]), max: Math.round(range[1]) },
        };
      } else {
        settings = {
          ...settings,
          levels: { min: Math.round(range[0]), max: Math.round(range[1]) },
        };
      }
      this.settingsContext.mutations.setChannelSettings(settings);
      this.experimentContext.actions.getChannelStackImage();
    }

    setSharedChannelLevels() {
      const metal = this.channel.metal;
      const settings = this.settingsContext.getters.channelSettings(this.channel.id);
      const levels = settings && settings.levels ?
        [settings.levels.min, settings.levels.max] :
        [this.channel.min_intensity, this.channel.max_intensity];
      this.experimentContext.actions.setSharedChannelLevels({ metal: metal, levels: levels });
    }

    get metalColor() {
      const colorMap = this.settingsContext.getters.metalColorMap;
      const color = colorMap.get(this.channel.metal);
      return color ? color : '#ffffff';
    }

    @Watch('color')
    onColorChanged(color: string) {
      this.settingsContext.mutations.setMetalColor({
        metal: this.channel.metal,
        color: color,
      });
    }
  }
</script>
