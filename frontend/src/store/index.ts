import Vue from 'vue';
import Vuex, { StoreOptions } from 'vuex';
import createLogger from 'vuex/dist/logger';

import { mainModule } from '@/modules/main';
import { State } from './state';
import { userModule } from '@/modules/user';
import { experimentModule } from '@/modules/experiment';

Vue.use(Vuex);

const debug = process.env.NODE_ENV !== 'prod';

const storeOptions: StoreOptions<State> = {
  modules: {
    main: mainModule,
    user: userModule,
    experiment: experimentModule
  },
  strict: debug,
  plugins: debug ? [createLogger()] : [],
};

export const store = new Vuex.Store<State>(storeOptions);

export default store;
