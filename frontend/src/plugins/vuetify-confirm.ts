import Vue from 'vue';
import VuetifyConfirm from 'vuetify-confirm';

Vue.use(VuetifyConfirm, {
  buttonTrueText: 'Accept',
  buttonFalseText: 'Cancel',
  color: 'warning',
  icon: 'mdi-alert',
  title: 'Warning',
  width: 350,
  property: '$confirm',
});
