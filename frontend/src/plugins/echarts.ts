import Vue from "vue";
import ECharts from "vue-echarts"; // refers to components/ECharts.vue in webpack

// import ECharts modules manually to reduce bundle size
// import 'echarts/lib/chart/bar'
// import 'echarts/lib/component/tooltip'

// If you want to use ECharts extensions, just import the extension package and it will work
// Taking ECharts-GL as an example:
// You only need to install the package with `npm install --save echarts-gl` and import it as follows
// import 'echarts-gl'

Vue.component("v-chart", ECharts);
