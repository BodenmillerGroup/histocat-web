<template>
  <v-card tile>
    <v-card-title class="card-title">
      <h4>Workspace</h4>
      <v-spacer/>
      <v-text-field
        v-model="search"
        append-icon="mdi-magnify"
        label="Search"
        single-line
        hide-details
        clearable
        solo-inverted
        flat
      />
    </v-card-title>
    <v-card-text>
      <v-treeview
        :items="items"
        :search="search"
        :filter="filter"
        :active.sync="active"
        item-key="uid"
        item-text="description"
        activatable
        transition
        return-object
        open-on-click
        open-all
      >
        <template v-slot:prepend="{ item }">
          <v-icon small>
            {{ icons[item.type] }}
          </v-icon>
        </template>
        <template v-slot:append="{ item }">
          <v-menu
            :close-on-content-click="false"
            :nudge-width="200"
            offset-x
            open-on-hover
          >
            <template v-slot:activator="{ on }">
              <v-btn
                flat
                icon
                color="grey"
                v-on="on"
              >
                <v-icon small>mdi-information-outline</v-icon>
              </v-btn>
            </template>
            <InfoCard :node="{ item }"></InfoCard>
          </v-menu>
        </template>
      </v-treeview>
    </v-card-text>
  </v-card>
</template>

<script lang="ts">
  import { Component, Prop, Vue, Watch } from 'vue-property-decorator';
  import { IExperiment } from '@/modules/experiment/models';
  import { commitSetChannels, commitSetSelectedAcquisitionId } from '@/modules/experiment/mutations';
  import InfoCard from '@/components/InfoCard.vue';

  @Component({
    components: { InfoCard },
  })
  export default class TreeView extends Vue {

    @Prop(Object) dataset?: IExperiment;

    search = null;
    active = [];

    icons = {
      slide: 'mdi-folder-outline',
      panorama: 'mdi-apps',
      roi: 'mdi-blur',
      acquisition: 'mdi-buffer',
    };

    @Watch('active')
    onActiveChanged(item: object) {
      if (this.active.length === 0 || !item) {
        return;
      }
      item = item[0];
      if (item.hasOwnProperty('channels')) {
        commitSetSelectedAcquisitionId(this.$store, item['id']);
        commitSetChannels(this.$store, item['channels']);
      }
    }

    get filter() {
      return (item, search, textKey) => item[textKey].indexOf(search) > -1;
    }

    get items() {
      if (this.dataset && this.dataset.slides) {
        return this.dataset.slides.map((slide) => {
          const panoramas = slide.panoramas.map((panorama) => {
            const rois = panorama.rois.map((roi) => {
              const acquisitions = roi.acquisitions.map((acquisition) => {
                return Object.assign({}, acquisition, { type: 'acquisition', uid: Math.random() });
              });
              return Object.assign({}, roi, {
                type: 'roi',
                description: `ROI [${roi.roi_type}]`,
                uid: Math.random(),
                children: acquisitions,
              });
            });
            return Object.assign({}, panorama, { type: 'panorama', uid: Math.random(), children: rois });
          });
          return Object.assign({}, slide, { type: 'slide', children: panoramas, uid: Math.random() });
        });
      }
    }
  }
</script>

<style scoped>
  .card-title {
    padding-bottom: 0;
  }
</style>

<style>
  .v-treeview-node__content {
    max-height: 24px;
  }

  .v-treeview-node__label {
    font-size: 10pt;
  }

  .v-treeview-node__root {
    min-height: 24px;
    font-size: 10pt;
  }

  .v-text-field.v-text-field--solo .v-input__control {
    min-height: 28px;
  }
</style>
