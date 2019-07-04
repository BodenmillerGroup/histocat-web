import { experimentModule, ExperimentsState } from '@/modules/experiment';
import { mainModule, MainState } from '@/modules/main';
import { settingsModule, SettingsState } from '@/modules/settings';
import { userModule, UserState } from '@/modules/user';
import localforage from 'localforage';
import Vue from 'vue';
import Vuex, { StoreOptions } from 'vuex';
import VuexPersistence from 'vuex-persist';
import createLogger from 'vuex/dist/logger';

Vue.use(Vuex);

const debug = process.env.NODE_ENV !== 'production';

export interface RootState {
  main: MainState;
  user: UserState;
  experiment: ExperimentsState;
  settings: SettingsState;
}

const vuexStorage = new VuexPersistence<RootState>({
  storage: localforage,
  asyncStorage: true,
  modules: [
    'settings'
  ]
});

const storeOptions: StoreOptions<RootState> = {
  modules: {
    main: mainModule,
    user: userModule,
    experiment: experimentModule,
    settings: settingsModule,
  },
  strict: debug,
  plugins: debug ? [
    vuexStorage.plugin,
    createLogger(),
  ] : [
    vuexStorage.plugin,
  ],
};

export const store = new Vuex.Store<RootState>(storeOptions);

export default store;
