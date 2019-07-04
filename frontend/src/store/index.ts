import Vue from 'vue';
import Vuex, { StoreOptions } from 'vuex';
import createLogger from 'vuex/dist/logger';
import VuexPersistence from 'vuex-persist';
import localforage from 'localforage';

import { mainModule } from '@/modules/main';
import { RootState } from './state';
import { userModule } from '@/modules/user';
import { experimentModule } from '@/modules/experiment';

Vue.use(Vuex);

const debug = process.env.NODE_ENV !== 'production';

const vuexStorage = new VuexPersistence<RootState>({
  storage: localforage,
  asyncStorage: true,
});

const storeOptions: StoreOptions<RootState> = {
  modules: {
    main: mainModule,
    user: userModule,
    experiment: experimentModule,
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
