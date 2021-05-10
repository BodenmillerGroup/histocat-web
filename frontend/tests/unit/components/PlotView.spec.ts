import { mount } from "@vue/test-utils";
import PlotView from "@/components/PlotView.vue";

import Vue from "vue";
import Vuetify from "vuetify";
Vue.use(Vuetify);

describe("PlotView.vue", () => {
  it("renders correctly", () => {
    const wrapper = mount(PlotView);
    expect(wrapper.element).toMatchSnapshot();
  });
});
