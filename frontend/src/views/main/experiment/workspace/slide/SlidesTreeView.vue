<template>
  <v-card tile>
    <v-toolbar flat dense color="grey lighten-4">
      <UploadButton :id="experiment.id"/>
      <v-spacer></v-spacer>
      <v-tooltip bottom>
        <template v-slot:activator="{ on }">
          <v-btn icon small v-on="on" @click="refreshSlides">
            <v-icon small>mdi-refresh</v-icon>
          </v-btn>
        </template>
        <span>Refresh slides</span>
      </v-tooltip>
      <v-tooltip bottom>
        <template v-slot:activator="{ on }">
          <v-btn icon small v-on="on" @click="toggleShowROI">
            <v-icon v-if="showROI" color="blue">mdi-blur</v-icon>
            <v-icon v-else color="grey">mdi-blur</v-icon>
          </v-btn>
        </template>
        <span v-if="showROI">Hide ROI</span>
        <span v-else>Show ROI</span>
      </v-tooltip>
    </v-toolbar>
    <v-toolbar dense flat>
      <v-text-field
        v-model="search"
        append-icon="mdi-magnify"
        label="Search"
        single-line
        hide-details
        clearable
        flat
      />
    </v-toolbar>
    <v-treeview
      v-model="selected"
      :items="items"
      :search="search"
      :filter="filter"
      :active.sync="active"
      item-key="uid"
      activatable
      transition
      return-object
      open-all
      selectable
      class="overflow-y-auto scroll-view"
    >
      <template v-slot:prepend="{ item }">
        <v-icon small>
          {{ icons[item.type] }}
        </v-icon>
      </template>
      <template v-slot:append="{ item }">
        <UploadArtifactsDialog
          v-if="item.type === 'slide'"
        ></UploadArtifactsDialog>
        <v-menu
          :close-on-content-click="false"
          :nudge-width="200"
          offset-x
        >
          <template v-slot:activator="{ on }">
            <v-btn
              v-on="on"
              small
              icon
              color="grey"
            >
              <v-icon small>mdi-information-outline</v-icon>
            </v-btn>
          </template>
          <InfoCard :node="{ item }"></InfoCard>
        </v-menu>
      </template>
    </v-treeview>
  </v-card>
</template>

<script lang="ts">
  import InfoCard from '@/components/InfoCard.vue';
  import UploadButton from '@/components/UploadButton.vue';
  import { experimentModule } from '@/modules/experiment';
  import { IExperiment } from '@/modules/experiment/models';
  import UploadArtifactsDialog from '@/views/main/experiment/workspace/slide/UploadArtifactsDialog.vue';
  import { Component, Prop, Vue, Watch } from 'vue-property-decorator';

  @Component({
    components: { UploadArtifactsDialog, UploadButton, InfoCard },
  })
  export default class SlidesTreeView extends Vue {
    readonly experimentContext = experimentModule.context(this.$store);

    @Prop(Object) readonly experiment!: IExperiment;

    selected = [];
    search = null;
    showROI = false;

    readonly icons = {
      slide: 'mdi-folder-outline',
      panorama: 'mdi-apps',
      roi: 'mdi-blur',
      acquisition: 'mdi-buffer',
    };

    toggleShowROI() {
      this.showROI = !this.showROI;
    }

    get active() {
      return [this.experimentContext.getters.activeWorkspaceNode];
    }

    set active(items: any[]) {
      if (!items || items.length === 0) {
        return;
      }
      this.experimentContext.mutations.setActiveWorkspaceNode(items[0]);
    }

    @Watch('selected')
    onSelectedChanged(items) {
      const ids = items
        .filter((item) => item.type === 'acquisition')
        .map((acquisition) => acquisition.id);
      this.experimentContext.mutations.setSelectedAcquisitionIds(ids);
    }

    get filter() {
      return (item, search, textKey) => item[textKey].indexOf(search) > -1;
    }

    async refreshSlides() {
      await this.experimentContext.actions.getExperimentData(this.experiment.id);
    }

    get items() {
      if (this.experiment.slides) {
        return this.experiment.slides.map((slide) => {
          const panoramas = slide.panoramas.map((panorama) => {
            const rois = panorama.rois.map((roi) => {
              const acquisitions = roi.acquisitions.map((acquisition) => {
                return Object.assign({}, acquisition, {
                  type: 'acquisition',
                  name: acquisition.meta.Description,
                  uid: `acquisition-${acquisition.id}`,
                });
              });
              return Object.assign({}, roi, {
                type: 'roi',
                name: `ROI ${roi.original_id}`,
                uid: `roi-${roi.id}`,
                children: acquisitions,
              });
            });
            const panoramaChildren = this.showROI ?
              rois :
              rois.reduce((total, roi) => {
                return total.concat(roi.children);
              }, [] as any);
            return Object.assign({}, panorama, {
              type: 'panorama',
              name: panorama.meta.Description,
              uid: `panorama-${panorama.id}`,
              children: panoramaChildren,
            });
          });
          return Object.assign({}, slide, {
            type: 'slide',
            name: slide.metaname,
            children: panoramas,
            uid: `slide-${slide.id}`,
          });
        });
      }
    }
  }
</script>

<style scoped>
  .scroll-view {
    height: calc(100vh - 200px);
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
