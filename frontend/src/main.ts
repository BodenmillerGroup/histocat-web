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

const isDev = process.env.NODE_ENV === "development";

Vue.config.productionTip = false;
Vue.config.performance = isDev;

new Vue({
  router: router,
  store: store,
  vuetify: vuetify,
  render: (h) => h(App),
}).$mount("#app");
