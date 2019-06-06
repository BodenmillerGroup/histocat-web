import { shallowMount } from '@vue/test-utils';
import '@/plugins/vuetify';
import UploadButton from '@/components/UploadButton.vue';

describe('UploadButton.vue', () => {
  it('renders props.title when passed', () => {
    const title = 'Upload Slide';
    const wrapper = shallowMount(UploadButton, {
      slots: {
        default: title,
      },
    });
    expect(wrapper.text()).toMatch(title);
  });
});
