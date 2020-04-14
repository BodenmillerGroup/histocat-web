import { datasetModule, DatasetState } from "@/modules/datasets";
import { DatasetActions } from "@/modules/datasets/actions";
import { DatasetGetters } from "@/modules/datasets/getters";
import { DatasetMutations } from "@/modules/datasets/mutations";
import { experimentModule, ExperimentState } from "@/modules/experiment";
import { ExperimentActions } from "@/modules/experiment/actions";
import { ExperimentGetters } from "@/modules/experiment/getters";
import { ExperimentMutations } from "@/modules/experiment/mutations";
import { mainModule, MainState } from "@/modules/main";
import { MainActions } from "@/modules/main/actions";
import { MainGetters } from "@/modules/main/getters";
import { MainMutations } from "@/modules/main/mutations";
import { settingsModule, SettingsState } from "@/modules/settings";
import { SettingsActions } from "@/modules/settings/actions";
import { SettingsGetters } from "@/modules/settings/getters";
import { SettingsMutations } from "@/modules/settings/mutations";
import { Store } from "vuex";
import { Context, Module } from "vuex-smart-module";

export class BroadcastManager {
  static mainContext: Context<Module<MainState, MainGetters, MainMutations, MainActions>>;
  static experimentContext: Context<Module<ExperimentState, ExperimentGetters, ExperimentMutations, ExperimentActions>>;
  static datasetContext: Context<Module<DatasetState, DatasetGetters, DatasetMutations, DatasetActions>>;
  static settingsContext: Context<Module<SettingsState, SettingsGetters, SettingsMutations, SettingsActions>>;
  static channel: BroadcastChannel;

  static init(store: Store<any>) {
    BroadcastManager.close();
    BroadcastManager.mainContext = mainModule.context(store);
    BroadcastManager.experimentContext = experimentModule.context(store);
    BroadcastManager.datasetContext = datasetModule.context(store);
    BroadcastManager.settingsContext = settingsModule.context(store);

    BroadcastManager.channel = new BroadcastChannel("HistoCAT");
    BroadcastManager.channel.onmessage = (ev) => {
      const method = ev.data.method;
      const payload = ev.data.payload;
      payload.suppressBroadcast = true; // Suppress broadcasting
      this.settingsContext.mutations[method](payload);
    };
    BroadcastManager.channel.onmessageerror = (ev) => {
      console.log(ev);
    };
  }

  static close() {
    if (BroadcastManager.channel) {
      BroadcastManager.channel.close();
    }
  }

  static postMessage(message) {
    if (BroadcastManager.channel) {
      BroadcastManager.channel.postMessage(message);
    }
  }
}
