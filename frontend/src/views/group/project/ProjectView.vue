<template>
  <LoadingView v-if="!projectData" text="Loading..." />
  <section v-else id="layoutContainer" ref="layoutContainer" v-resize="handleWindowResizeEvent" />
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
import { CellsComponent } from "@/views/group/project/cells/CellsComponent";
import { GatingComponent } from "@/views/group/project/gating/GatingComponent";
import { PipelineComponent } from "@/views/group/project/pipeline/PipelineComponent";
import { PcaComponent } from "@/views/group/project/pca/PcaComponent";
import { UmapComponent } from "@/views/group/project/umap/UmapComponent";
import { TsneComponent } from "@/views/group/project/tsne/TsneComponent";
import { ScatterComponent } from "@/views/group/project/scatter/ScatterComponent";
import { LeidenComponent } from "@/views/group/project/leiden/LeidenComponent";
import { LouvainComponent } from "@/views/group/project/louvain/LouvainComponent";
import { PipelinesComponent } from "@/views/group/project/pipelines/PipelinesComponent";

@Component({
  components: {
    LoadingView,
  },
})
export default class ProjectView extends Vue {
  readonly uiContext = uiModule.context(this.$store);
  readonly mainContext = mainModule.context(this.$store);
  readonly projectsContext = projectsModule.context(this.$store);

  _goldenLayout?: GoldenLayout;
  _containerMap = new Map();

  get projectData() {
    return this.projectsContext.getters.projectData;
  }

  get activeLayout() {
    return this.uiContext.getters.activeLayout;
  }

  get responsive() {
    return this.uiContext.getters.responsive;
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
      case CellsComponent.typeName:
        return new CellsComponent(container, this.$store, this);
      case GatingComponent.typeName:
        return new GatingComponent(container, this.$store, this);
      case PipelineComponent.typeName:
        return new PipelineComponent(container, this.$store, this);
      case PcaComponent.typeName:
        return new PcaComponent(container, this.$store, this);
      case UmapComponent.typeName:
        return new UmapComponent(container, this.$store, this);
      case TsneComponent.typeName:
        return new TsneComponent(container, this.$store, this);
      case ScatterComponent.typeName:
        return new ScatterComponent(container, this.$store, this);
      case LeidenComponent.typeName:
        return new LeidenComponent(container, this.$store, this);
      case LouvainComponent.typeName:
        return new LouvainComponent(container, this.$store, this);
      case PipelinesComponent.typeName:
        return new PipelinesComponent(container, this.$store, this);
    }
  }

  private disposeComponent(component: ComponentItem.Component) {}

  private handleWindowResizeEvent() {
    // handling of resize event is required if GoldenLayout does not use body element
    if (this._goldenLayout) {
      const appBarHeight = 48;
      this._goldenLayout.setSize(this.responsive.width!, this.responsive.height! - appBarHeight);
    }
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
    this.uiContext.mutations.setGoldenLayout(this._goldenLayout);
  }

  beforeDestroy() {
    WebSocketManager.close();
    this.$store.dispatch("resetProject");
    this.uiContext.mutations.setGoldenLayout(null);
  }
}
</script>

<style scoped>
#layoutContainer {
  width: 100%;
  height: 100%;
}
</style>
