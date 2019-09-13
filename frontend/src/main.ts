// Import Component hooks before component definitions
import "@/component-hooks";
import App from "@/App.vue";
import "@/plugins/masonry-css";
import "@/plugins/echarts";
import vuetify from "@/plugins/vuetify";
import "@/registerServiceWorker";
import router from "@/router";
import store from "@/store";
import Vue from "vue";
import "@mdi/font/css/materialdesignicons.css"; // Ensure you are using css-loader

Vue.config.productionTip = false;

new Vue({
  router,
  store,
  vuetify,
  render: h => h(App)
}).$mount("#app");
