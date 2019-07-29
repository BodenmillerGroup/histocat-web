// Import Component hooks before component definitions

import '@babel/polyfill';
import App from '@/App.vue';
import '@/component-hooks';
import '@/plugins/masonry-css';
import '@/plugins/vee-validate';
import vuetify from '@/plugins/vuetify';
import '@/registerServiceWorker';
import router from '@/router';
import store from '@/store';
import '@mdi/font/css/materialdesignicons.css'; // Ensure you are using css-loader
import Vue from 'vue';
import 'vuetify/dist/vuetify.min.css';

Vue.config.productionTip = false;

new Vue({
  router,
  store,
  render: (h) => h(App),
  // @ts-ignore
  vuetify: vuetify,
}).$mount('#app');
