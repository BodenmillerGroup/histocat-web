import { shallowMount } from "@vue/test-utils";
import UploadButton from "@/components/UploadButton.vue";

import Vue from "vue";
import Vuetify from "vuetify";
Vue.use(Vuetify);

describe("UploadButton.vue", () => {
  it("renders props.title when passed", () => {
    const title = "Upload";
    const wrapper = shallowMount(UploadButton, {
      slots: {
        default: title,
      },
      propsData: {
        label: "Upload",
        upload: () => {}
      },
    });
    expect(wrapper.text()).toMatch(title);
  });
});
