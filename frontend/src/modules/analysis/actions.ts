import { experimentModule } from '@/modules/experiment';
import { ExportTypes } from '@/modules/experiment/models';
import { mainModule } from '@/modules/main';
import { settingsModule } from '@/modules/settings';
import { saveAs } from 'file-saver';
import { Store } from 'vuex';
import { Actions, Context } from 'vuex-smart-module';
import { AnalysisState } from '.';
import { api } from './api';
import { AnalysisGetters } from './getters';
import { AnalysisMutations } from './mutations';

export class AnalysisActions extends Actions<AnalysisState, AnalysisGetters, AnalysisMutations, AnalysisActions> {

  // Declare context type
  main?: Context<typeof mainModule>;
  settings?: Context<typeof settingsModule>;
  experiment?: Context<typeof experimentModule>;

  // Called after the module is initialized
  $init(store: Store<any>): void {
    // Create and retain main module context
    this.main = mainModule.context(store);
    this.settings = settingsModule.context(store);
    this.experiment = experimentModule.context(store);
  }

  async getAnalysisImage() {
    const params = this.prepareStackParams();
    if (params.channels.length === 0) {
      return;
    }
    try {
      const response = await api.downloadAnalysisImage(this.main!.getters.token, params);
      const blob = await response.blob();
      const reader = new FileReader();
      reader.readAsDataURL(blob);
      reader.onloadend = () => {
        this.mutations.setAnalysisImage(reader.result);
      };
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async exportAnalysisImage(format: ExportTypes = 'png') {
    const params = this.prepareStackParams(format);
    try {
      const response = await api.downloadAnalysisImage(this.main!.getters.token, params);
      const blob = await response.blob();
      saveAs(blob);
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  private prepareStackParams(format: 'png' | 'tiff' = 'png') {
    const channels = this.experiment!.getters.selectedChannels.map((channel) => {
      const color = this.settings!.getters.metalColorMap.get(channel.metal);
      const settings = this.settings!.getters.channelSettings(channel.id);
      const min = settings && settings.levels ? settings.levels.min : undefined;
      const max = settings && settings.levels ? settings.levels.max : undefined;
      const customLabel = settings && settings.customLabel ? settings.customLabel : channel.label;
      return {
        id: channel.id,
        color: color,
        customLabel: customLabel,
        min: min,
        max: max,
      };
    });

    const filter = this.settings!.getters.filter;
    const scalebar = this.settings!.getters.scalebar;

    return {
      format: format,
      filter: filter,
      scalebar: scalebar,
      channels: channels,
    };
  }
}
