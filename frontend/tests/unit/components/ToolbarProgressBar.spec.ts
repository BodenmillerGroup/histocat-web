import { mount } from "@vue/test-utils";
import ToolbarProgressBar from "@/components/ToolbarProgressBar.vue";

import Vue from "vue";
import Vuetify from "vuetify";
Vue.use(Vuetify);

describe("ToolbarProgressBar.vue", () => {
  it("renders correctly", () => {
    const wrapper = mount(ToolbarProgressBar);
    expect(wrapper.element).toMatchSnapshot();
  });
});
