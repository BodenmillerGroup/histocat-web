import { Module } from "vuex-smart-module";
import { UiActions } from "./actions";
import { UiGetters } from "./getters";
import { UiMutations } from "./mutations";
import { ComponentItemConfig, ItemType, LayoutConfig } from "golden-layout";
import { ImageComponent } from "@/views/group/project/image/ImageComponent";
import { ILayout, IResponsive, ViewMode } from "./models";
import { SlidesComponent } from "@/views/group/project/slides/SlidesComponent";
import { ChannelsComponent } from "@/views/group/project/channels/ChannelsComponent";
import { TilesComponent } from "@/views/group/project/tiles/TilesComponent";
import { RegionComponent } from "@/views/group/project/region/RegionComponent";
import { PresetsComponent } from "@/views/group/project/presets/PresetsComponent";
import { HistogramComponent } from "@/views/group/project/histogram/HistogramComponent";
import { SettingsComponent } from "@/views/group/project/settings/SettingsComponent";
import { SegmentationComponent } from "@/views/group/project/segmentation/SegmentationComponent";

const defaultConfig: LayoutConfig = {
  root: {
    type: ItemType.row,
    content: [
      {
        id: SlidesComponent.typeName,
        componentType: SlidesComponent.typeName,
        type: "component",
        title: "Slides",
        header: {
          show: "top",
          popout: false,
        },
        isClosable: true,
        reorderEnabled: true,
        width: 15,
        componentState: undefined,
      } as ComponentItemConfig,
      {
        type: "stack",
        header: {
          show: "top",
          popout: false,
        },
        content: [
          {
            id: ImageComponent.typeName,
            componentType: ImageComponent.typeName,
            type: "component",
            title: "Image",
            header: {
              show: "top",
              popout: false,
            },
            isClosable: true,
            reorderEnabled: true,
            componentState: undefined,
          } as ComponentItemConfig,
          {
            id: TilesComponent.typeName,
            componentType: TilesComponent.typeName,
            type: "component",
            title: "Tiles",
            header: {
              show: "top",
              popout: false,
            },
            isClosable: true,
            reorderEnabled: true,
            componentState: undefined,
          } as ComponentItemConfig,
          {
            id: SegmentationComponent.typeName,
            componentType: SegmentationComponent.typeName,
            type: "component",
            title: "Segmentation",
            header: {
              show: "top",
              popout: false,
            },
            isClosable: true,
            reorderEnabled: true,
            componentState: undefined,
          } as ComponentItemConfig,
        ],
      },
      {
        width: 15,
        type: "column",
        content: [
          {
            type: "stack",
            height: 50,
            header: {
              show: "top",
              popout: false,
            },
            content: [
              {
                id: ChannelsComponent.typeName,
                componentType: ChannelsComponent.typeName,
                type: "component",
                title: "Channels",
                header: {
                  show: "top",
                  popout: false,
                },
                isClosable: true,
                reorderEnabled: true,
                componentState: undefined,
              } as ComponentItemConfig,
              {
                id: RegionComponent.typeName,
                componentType: RegionComponent.typeName,
                type: "component",
                title: "Region",
                header: {
                  show: "top",
                  popout: false,
                },
                isClosable: true,
                reorderEnabled: true,
                componentState: undefined,
              } as ComponentItemConfig,
            ],
          },
          {
            type: "stack",
            height: 50,
            header: {
              show: "top",
              popout: false,
            },
            content: [
              {
                id: HistogramComponent.typeName,
                componentType: HistogramComponent.typeName,
                type: "component",
                title: "Histogram",
                header: {
                  show: "top",
                  popout: false,
                },
                isClosable: true,
                reorderEnabled: true,
                componentState: undefined,
              } as ComponentItemConfig,
              {
                id: PresetsComponent.typeName,
                componentType: PresetsComponent.typeName,
                type: "component",
                title: "Presets",
                header: {
                  show: "top",
                  popout: false,
                },
                isClosable: true,
                reorderEnabled: true,
                componentState: undefined,
              } as ComponentItemConfig,
              {
                id: SettingsComponent.typeName,
                componentType: SettingsComponent.typeName,
                type: "component",
                title: "Settings",
                header: {
                  show: "top",
                  popout: false,
                },
                isClosable: true,
                reorderEnabled: true,
                componentState: undefined,
              } as ComponentItemConfig,
            ],
          },
        ],
      },
    ],
  },
};

const defaultLayout: ILayout = {
  name: "default",
  config: defaultConfig,
};

export class UiState {
  activeLayout = defaultLayout;
  responsive: IResponsive = {
    width: null,
    height: null,
  };
  dashboardMiniDrawer = true;
  dashboardShowDrawer = true;
  showWorkspace = true;
  showOptions = true;
  viewMode: ViewMode = "image";

  processing = false;
  processingProgress = 0;

  maskMode: "raw" | "mask" | "origin" = "raw";
  maskOpacity = 1.0;
  mouseMode: "panZoom" | "lasso" | "rotate" = "panZoom";
}

export const uiModule = new Module({
  namespaced: true,

  state: UiState,
  getters: UiGetters,
  mutations: UiMutations,
  actions: UiActions,
});
