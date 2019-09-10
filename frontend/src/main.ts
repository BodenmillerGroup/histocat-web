// Import Component hooks before component definitions
import '@/component-hooks';
import App from '@/App.vue';
import '@/plugins/masonry-css';
import '@/plugins/echarts';
import '@/plugins/vee-validate';
import vuetify from '@/plugins/vuetify';
import '@/registerServiceWorker';
import router from '@/router';
import store from '@/store';
import Vue from 'vue';
import 'vuetify/dist/vuetify.min.css';
import '@mdi/font/css/materialdesignicons.css'; // Ensure you are using css-loader

Vue.config.productionTip = false;

new Vue({
  router,
  store,
  render: (h) => h(App),
  // @ts-ignore
  vuetify: vuetify,
}).$mount('#app');
