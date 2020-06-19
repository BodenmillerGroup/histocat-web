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
    </v-toolbar>
    <v-toolbar dense flat>
      <v-text-field v-model="search" label="Search" single-line hide-details clearable dense>
        <template v-slot:append-outer>
          <v-icon dense>mdi-magnify</v-icon>
        </template>
      </v-text-field>
    </v-toolbar>
    <v-treeview
      v-model="selected"
      :items="items"
      :search="search"
      :filter="filter"
      :active.sync="active"
      item-disabled="locked"
      item-key="uid"
      activatable
      transition
      return-object
      selectable
      class="overflow-y-auto scroll-view"
    >
      <template v-slot:prepend="{ item }">
        <v-icon small>
          {{ icons[item.type] }}
        </v-icon>
      </template>
      <template v-slot:append="{ item }">
        <v-tooltip bottom v-if="item.type === 'acquisition' && item.hasMask">
          <template v-slot:activator="{ on }">
            <v-icon small color="grey" v-on="on">
              mdi-transition-masked
            </v-icon>
          </template>
          <span>Mask available</span>
        </v-tooltip>
        <v-tooltip bottom v-if="item.type === 'slide'">
          <template v-slot:activator="{ on }">
            <v-btn v-on="on" small icon color="grey" @click.stop="deleteSlide(item.id)">
              <v-icon small>mdi-delete</v-icon>
            </v-btn>
          </template>
          <span>Delete slide</span>
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
import { BroadcastManager } from "@/utils/BroadcastManager";
import {
  SET_ACTIVE_ACQUISITION_ID,
  SET_ACTIVE_WORKSPACE_NODE,
  SET_SELECTED_ACQUISITION_IDS,
} from "@/modules/experiment/events";

@Component({
  components: { UploadButton, InfoCard },
})
export default class ImageWorkspaceView extends Vue {
  readonly experimentContext = experimentModule.context(this.$store);
  readonly datasetContext = datasetModule.context(this.$store);

  @Prop(Object) readonly experiment!: IExperiment;

  mutex = false;

  selected: any[] = [];
  search = null;

  readonly icons = {
    slide: "mdi-folder-outline",
    panorama: "mdi-apps",
    roi: "mdi-blur",
    acquisition: "mdi-buffer",
  };

  get active() {
    return [this.experimentContext.getters.activeWorkspaceNode];
  }

  set active(items: any[]) {
    if (!items || items.length === 0) {
      return;
    }
    const node = items[0];
    BroadcastManager.publish(SET_ACTIVE_WORKSPACE_NODE, node);
    if (node.type === "acquisition") {
      BroadcastManager.publish(SET_ACTIVE_ACQUISITION_ID, node.id);
      this.experimentContext.actions.getChannelStackImage();
    }
  }

  get selectedAcquisitionIds() {
    return this.experimentContext.getters.selectedAcquisitionIds;
  }

  @Watch("selected")
  async selectedChanged(items: any[]) {
    this.mutex = true;
    const ids = items.filter((item) => item.type === "acquisition").map((acquisition) => acquisition.id);
    BroadcastManager.publish(SET_SELECTED_ACQUISITION_IDS, ids);
    this.mutex = false;
  }

  @Watch("selectedAcquisitionIds")
  selectedAcquisitionIdsChanged(newIds: number[], oldIds: number[]) {
    if (!this.mutex && !equals(newIds, oldIds)) {
      const items: any[] = [];
      const nodes = {
        children: this.items,
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
      return this.experiment.slides.map((slide) => {
        const acquisitions = slide.acquisitions.map((acquisition) => {
          let hasMask = false;
          if (this.dataset && this.dataset.input && this.dataset.input.probability_masks) {
            hasMask = !!this.dataset.input.probability_masks[acquisition.id];
          }
          return Object.assign({}, acquisition, {
            type: "acquisition",
            name: `${acquisition.id}: ${acquisition.description}`,
            uid: `acquisition-${acquisition.id}`,
            locked: !acquisition.is_valid,
            hasMask: hasMask,
          });
        });
        const panoramas = slide.panoramas.map((panorama) => {
          return Object.assign({}, panorama, {
            type: "panorama",
            name: `${panorama.id}: ${panorama.description}`,
            uid: `panorama-${panorama.id}`,
            locked: panorama.image_type == "Default",
          });
        });
        const panoramasRoot = {
          name: "Panorama Images",
          type: "roi",
          uid: `slide-${slide.id}-panoramas`,
          children: panoramas,
        };
        return Object.assign({}, slide, {
          type: "slide",
          name: slide.name,
          children: [panoramasRoot as any].concat(acquisitions),
          uid: `slide-${slide.id}`,
        });
      });
    } else {
      return [];
    }
  }

  async deleteSlide(id: number) {
    if (self.confirm("Do you really want to delete the slide?")) {
      await this.experimentContext.actions.deleteSlide(id);
      await this.refreshSlides();
    }
  }
}
</script>

<style scoped>
.scroll-view {
  height: calc(100vh - 144px);
}
</style>
