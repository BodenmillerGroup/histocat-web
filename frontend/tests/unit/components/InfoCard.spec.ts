import { mount } from "@vue/test-utils";
import InfoCard from "@/components/InfoCard.vue";

import Vue from "vue";
import Vuetify from "vuetify";
Vue.use(Vuetify);

describe("InfoCard.vue", () => {
  it("renders correctly", () => {
    const wrapper = mount(InfoCard, {
      propsData: {
        node: {
          item: {
            type: "slide",
          },
        },
      },
    });
    expect(wrapper.element).toMatchSnapshot();
  });
});
