// Import Component hooks before component definitions

import '@/component-hooks';
import Vue from 'vue';
import '@/plugins/vuetify';
import '@/plugins/vee-validate';
import '@/plugins/vuetify-confirm';
import '@/plugins/masonry-css';
import App from '@/App.vue';
import router from '@/router';
import store from '@/store';
import '@/registerServiceWorker';
import 'vuetify/dist/vuetify.min.css';

Vue.config.productionTip = false;

new Vue({
  router,
  store,
  render: (h) => h(App),
}).$mount('#app');
