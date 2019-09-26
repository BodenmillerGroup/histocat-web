<template>
  <v-card tile>
    <v-toolbar flat dense color="grey lighten-4">
      <UploadButton :id="experiment.id" />
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
      <v-text-field v-model="search" append-icon="mdi-magnify" label="Search" single-line hide-details clearable flat />
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
        <v-tooltip bottom>
          <template v-slot:activator="{ on }">
            <v-icon v-if="item.type === 'acquisition' && item.hasMask" small color="grey" v-on="on">
              mdi-transition-masked
            </v-icon>
          </template>
          <span>Mask available</span>
        </v-tooltip>
        <v-menu :close-on-content-click="false" :nudge-width="200" offset-x>
          <template v-slot:activator="{ on }">
            <v-btn v-on="on" small icon color="grey" @click.stop="">
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
import InfoCard from "@/components/InfoCard.vue";
import UploadButton from "@/components/UploadButton.vue";
import { datasetModule } from "@/modules/datasets";
import { experimentModule } from "@/modules/experiment";
import { IExperiment } from "@/modules/experiment/models";
import { equals } from "rambda";
import { Component, Prop, Vue, Watch } from "vue-property-decorator";

@Component({
  components: { UploadButton, InfoCard }
})
export default class SlidesTreeView extends Vue {
  readonly experimentContext = experimentModule.context(this.$store);
  readonly datasetContext = datasetModule.context(this.$store);

  @Prop(Object) readonly experiment!: IExperiment;

  mutex = false;

  selected: any[] = [];
  search = null;
  showROI = false;

  readonly icons = {
    slide: "mdi-folder-outline",
    panorama: "mdi-apps",
    roi: "mdi-blur",
    acquisition: "mdi-buffer"
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

  get selectedAcquisitionIds() {
    return this.experimentContext.getters.selectedAcquisitionIds;
  }

  @Watch("selected")
  async selectedChanged(items: any[]) {
    this.mutex = true;
    const ids = items.filter(item => item.type === "acquisition").map(acquisition => acquisition.id);
    await this.experimentContext.mutations.setSelectedAcquisitionIds(ids);
    this.mutex = false;
  }

  @Watch("selectedAcquisitionIds")
  selectedAcquisitionIdsChanged(newIds: number[], oldIds: number[]) {
    if (!this.mutex && !equals(newIds, oldIds)) {
      const items: any[] = [];
      const nodes = {
        children: this.items
      };

      for (const id of newIds) {
        const item = this.searchTree(nodes, id);
        if (item) {
          items.push(item);
        }
      }
      this.selected = items;
    }
  }

  searchTree(element, id) {
    if (element.type === "acquisition" && element.id === id) {
      return element;
    } else if (element.children != null) {
      var i;
      var result = null;
      for (i = 0; result == null && i < element.children.length; i++) {
        result = this.searchTree(element.children[i], id);
      }
      return result;
    }
    return null;
  }

  get filter() {
    return (item, search, textKey) => item[textKey].indexOf(search) > -1;
  }

  async refreshSlides() {
    await this.experimentContext.actions.getExperimentData(this.experiment.id);
  }

  get dataset() {
    return this.datasetContext.getters.activeDataset;
  }

  get items() {
    if (this.experiment.slides) {
      return this.experiment.slides.map(slide => {
        const panoramas = slide.panoramas.map(panorama => {
          const rois = panorama.rois.map(roi => {
            const acquisitions = roi.acquisitions.map(acquisition => {
              let hasMask = false;
              if (this.dataset && this.dataset.input && this.dataset.input.probability_masks) {
                hasMask = !!this.dataset.input.probability_masks[acquisition.id];
              }
              return Object.assign({}, acquisition, {
                type: "acquisition",
                name: `${acquisition.id}: ${acquisition.meta.Description}`,
                uid: `acquisition-${acquisition.id}`,
                hasMask: hasMask
              });
            });
            return Object.assign({}, roi, {
              type: "roi",
              name: `ROI ${roi.origin_id}`,
              uid: `roi-${roi.id}`,
              children: acquisitions
            });
          });
          const panoramaChildren = this.showROI
            ? rois
            : rois.reduce(
                (total, roi) => {
                  return total.concat(roi.children);
                },
                [] as any
              );
          return Object.assign({}, panorama, {
            type: "panorama",
            name: panorama.meta.Description,
            uid: `panorama-${panorama.id}`,
            children: panoramaChildren
          });
        });
        return Object.assign({}, slide, {
          type: "slide",
          name: slide.name,
          children: panoramas,
          uid: `slide-${slide.id}`
        });
      });
    } else {
      return [];
    }
  }
}
</script>

<style scoped>
.scroll-view {
  height: calc(100vh - 196px);
}
</style>

<style>
.v-text-field.v-text-field--solo .v-input__control {
  min-height: 28px;
}
</style>
