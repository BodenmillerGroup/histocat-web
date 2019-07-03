import Vue from 'vue';
import Vuex, { StoreOptions } from 'vuex';
import createLogger from 'vuex/dist/logger';

import { mainModule } from '@/modules/main';
import { RootState } from './state';
import { userModule } from '@/modules/user';
import { experimentModule } from '@/modules/experiment';

import cache from './cache';
import sync from './sync';

Vue.use(Vuex);

const debug = process.env.NODE_ENV !== 'production';

const storeOptions: StoreOptions<RootState> = {
  modules: {
    main: mainModule,
    user: userModule,
    experiment: experimentModule,
  },
  strict: debug,
  plugins: debug ? [
    cache,
    sync,
    createLogger()
  ] : [
    cache,
    sync,
  ],
};

export const store = new Vuex.Store<RootState>(storeOptions);

export default store;
