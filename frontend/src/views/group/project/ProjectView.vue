<template>
  <LoadingView v-if="!projectData" text="Loading..." />
  <!--  <v-container v-else fluid class="ma-0 pa-0">-->
  <!--    <ImageView v-show="viewMode === 'image'" />-->
  <!--    <SegmentationView v-show="viewMode === 'segmentation'" />-->
  <!--    <DataView v-show="viewMode === 'data'" />-->
  <!--  </v-container>-->
  <section v-else id="layoutContainer" ref="layoutContainer"></section>
</template>

<script lang="ts">
import LoadingView from "@/components/LoadingView.vue";
import { projectsModule } from "@/modules/projects";
import { mainModule } from "@/modules/main";
import { WebSocketManager } from "@/utils/WebSocketManager";
import { Component, Vue } from "vue-property-decorator";
import { ComponentContainer, ComponentItem, GoldenLayout, ResolvedComponentItemConfig } from "golden-layout";
import { ImageComponent } from "@/views/group/project/image/ImageComponent";
import "golden-layout/dist/css/goldenlayout-base.css";
import "@/sass/goldenlayout-theme.scss";
import { uiModule } from "@/modules/ui";
import { SlidesComponent } from "@/views/group/project/slides/SlidesComponent";
import { ChannelsComponent } from "@/views/group/project/channels/ChannelsComponent";
import { TilesComponent } from "@/views/group/project/tiles/TilesComponent";
import { RegionComponent } from "@/views/group/project/region/RegionComponent";
import { PresetsComponent } from "@/views/group/project/presets/PresetsComponent";
import { HistogramComponent } from "@/views/group/project/histogram/HistogramComponent";
import { SettingsComponent } from "@/views/group/project/settings/SettingsComponent";
import { SegmentationComponent } from "@/views/group/project/segmentation/SegmentationComponent";

@Component({
  components: {
    LoadingView,
  },
})
export default class ProjectView extends Vue {
  readonly uiContext = uiModule.context(this.$store);
  readonly mainContext = mainModule.context(this.$store);
  readonly projectsContext = projectsModule.context(this.$store);

  _goldenLayout = new GoldenLayout(this.$refs.layoutContainer as HTMLElement);
  _containerMap = new Map();

  private _windowResizeListener = () => this.handleWindowResizeEvent();

  get projectData() {
    return this.projectsContext.getters.projectData;
  }

  get viewMode() {
    return this.uiContext.getters.viewMode;
  }

  get activeLayout() {
    return this.uiContext.getters.activeLayout;
  }

  private createComponent(container: ComponentContainer, config: ResolvedComponentItemConfig) {
    switch (config.componentType) {
      case SlidesComponent.typeName:
        return new SlidesComponent(container, this.$store, this);
      case ImageComponent.typeName:
        return new ImageComponent(container, this.$store, this);
      case TilesComponent.typeName:
        return new TilesComponent(container, this.$store, this);
      case ChannelsComponent.typeName:
        return new ChannelsComponent(container, this.$store, this);
      case RegionComponent.typeName:
        return new RegionComponent(container, this.$store, this);
      case HistogramComponent.typeName:
        return new HistogramComponent(container, this.$store, this);
      case PresetsComponent.typeName:
        return new PresetsComponent(container, this.$store, this);
      case SettingsComponent.typeName:
        return new SettingsComponent(container, this.$store, this);
      case SegmentationComponent.typeName:
        return new SegmentationComponent(container, this.$store, this);
    }
  }

  private disposeComponent(component: ComponentItem.Component) {}

  private handleWindowResizeEvent() {
    // handling of resize event is required if GoldenLayout does not use body element
    const bodyWidth = document.body.offsetWidth;
    const controlsWidth = 60;
    const bodyHeight = document.body.offsetHeight;
    const controlsHeight = 60;
    this._goldenLayout.setSize(bodyWidth - controlsWidth, bodyHeight - controlsHeight);
  }

  async mounted() {
    const projectId = +this.$router.currentRoute.params.projectId;
    this.projectsContext.mutations.setActiveProjectId(projectId);
    await this.projectsContext.actions.getProjectData(projectId);
    WebSocketManager.connect(projectId);

    this._goldenLayout = new GoldenLayout(this.$refs.layoutContainer as HTMLElement);
    this._containerMap = new Map();

    this._goldenLayout.getComponentEvent = (container, itemConfig) => {
      const component = this.createComponent(container, itemConfig);
      // component.element is the top most HTML element in the component
      // container.element.appendChild(component);
      this._containerMap.set(container, component);
    };

    this._goldenLayout.releaseComponentEvent = (container) => {
      // do this if you need to dispose resources
      const component = this._containerMap.get(container);
      this.disposeComponent(component);
      this._containerMap.delete(container);
    };

    this._goldenLayout.loadLayout(this.activeLayout.config);

    globalThis.addEventListener("resize", this._windowResizeListener, { passive: true });
  }

  beforeDestroy() {
    WebSocketManager.close();
    // if (process.env.VUE_APP_ENV !== "development") {
    //   this.$store.dispatch("reset");
    // }
    this.$store.dispatch("resetProject");
    this.uiContext.mutations.setViewMode("image");
  }
}
</script>

<style scoped>
#layoutContainer {
  width: 100%;
  height: 100%;
}
</style>
