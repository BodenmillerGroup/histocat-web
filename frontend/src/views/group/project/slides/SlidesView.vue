<template>
  <div v-if="projectData" class="root">
    <v-toolbar dense color="grey lighten-4">
      <UploadButton label="Upload slide" :upload="upload" />
      <v-spacer />
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
        <template v-slot:append>
          <v-icon dense>mdi-magnify</v-icon>
        </template>
      </v-text-field>
    </v-toolbar>
    <v-treeview
      :items="items"
      :search="search"
      :filter="filter"
      :active.sync="active"
      item-disabled="locked"
      item-key="uid"
      activatable
      transition
      return-object
    >
      <template v-slot:prepend="{ item }">
        <v-icon small>
          {{ icons[item.type] }}
        </v-icon>
      </template>
      <template v-slot:append="{ item }">
        <v-tooltip bottom v-if="item.type === 'acquisition' && item.hasMask">
          <template v-slot:activator="{ on }">
            <v-icon small color="grey" v-on="on">mdi-transition-masked</v-icon>
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
  </div>
</template>

<script lang="ts">
import InfoCard from "@/components/InfoCard.vue";
import UploadButton from "@/components/UploadButton.vue";
import { datasetsModule } from "@/modules/datasets";
import { projectsModule } from "@/modules/projects";
import { Component, Vue } from "vue-property-decorator";

@Component({
  components: { UploadButton, InfoCard },
})
export default class ImageWorkspaceView extends Vue {
  readonly projectsContext = projectsModule.context(this.$store);
  readonly datasetsContext = datasetsModule.context(this.$store);

  search = null;

  readonly icons = {
    slide: "mdi-folder-outline",
    panorama: "mdi-apps",
    roi: "mdi-blur",
    acquisition: "mdi-buffer",
  };

  get projectData() {
    return this.projectsContext.getters.projectData!;
  }

  get active() {
    return [this.projectsContext.getters.activeWorkspaceNode];
  }

  set active(items: any[]) {
    if (!items || items.length === 0) {
      return;
    }
    const node = items[0];
    this.projectsContext.mutations.setActiveWorkspaceNode(node);
    if (node.type === "acquisition") {
      this.projectsContext.mutations.setActiveAcquisitionId(node.id);
      this.projectsContext.actions.getChannelStackImage();
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
    await this.projectsContext.actions.getProjectData(this.projectData.id);
  }

  get dataset() {
    return this.datasetsContext.getters.activeDataset;
  }

  get items() {
    if (this.projectData.slides) {
      return this.projectData.slides.map((slide) => {
        const acquisitions = slide.acquisitions.map((acquisition) => {
          let hasMask = false;
          if (this.dataset && this.dataset.meta.masks) {
            hasMask = !!this.dataset.meta.masks[acquisition.id];
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
      await this.projectsContext.actions.deleteSlide(id);
      await this.refreshSlides();
    }
  }

  async upload(data: FormData) {
    await this.projectsContext.actions.uploadSlide({ id: this.projectData.id, data: data });
  }
}
</script>

<style scoped>
.root {
  width: 100%;
  height: 100%;
  overflow-y: auto;
}
</style>
