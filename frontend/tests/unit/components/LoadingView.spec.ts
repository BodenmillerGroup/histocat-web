import { mount } from "@vue/test-utils";
import LoadingView from "@/components/LoadingView.vue";

import Vue from "vue";
import Vuetify from "vuetify";
Vue.use(Vuetify);

describe("LoadingView.vue", () => {
  it("renders correctly", () => {
    const wrapper = mount(LoadingView);
    expect(wrapper.element).toMatchSnapshot();
  });
});
