import { ComponentItemConfig, ItemType, LayoutConfig } from "golden-layout";
import { SlidesComponent } from "@/views/group/project/slides/SlidesComponent";
import { CellsComponent } from "@/views/group/project/cells/CellsComponent";
import { GatingComponent } from "@/views/group/project/gating/GatingComponent";
import { PipelinesComponent } from "@/views/group/project/pipelines/PipelinesComponent";
import { ImageComponent } from "@/views/group/project/image/ImageComponent";
import { TilesComponent } from "@/views/group/project/tiles/TilesComponent";
import { SegmentationComponent } from "@/views/group/project/segmentation/SegmentationComponent";
import { ScatterComponent } from "@/views/group/project/scatter/ScatterComponent";
import { PcaComponent } from "@/views/group/project/pca/PcaComponent";
import { UmapComponent } from "@/views/group/project/umap/UmapComponent";
import { TsneComponent } from "@/views/group/project/tsne/TsneComponent";
import { LeidenComponent } from "@/views/group/project/leiden/LeidenComponent";
import { LouvainComponent } from "@/views/group/project/louvain/LouvainComponent";
import { PipelineComponent } from "@/views/group/project/pipeline/PipelineComponent";
import { ChannelsComponent } from "@/views/group/project/channels/ChannelsComponent";
import { RegionComponent } from "@/views/group/project/region/RegionComponent";
import { HistogramComponent } from "@/views/group/project/histogram/HistogramComponent";
import { PresetsComponent } from "@/views/group/project/presets/PresetsComponent";
import { SettingsComponent } from "@/views/group/project/settings/SettingsComponent";
import { ILayout } from "./models";
import { DEFAULT_LAYOUT_UID } from "@/modules/ui/constants";

const defaultConfig: LayoutConfig = {
  root: {
    type: ItemType.row,
    content: [
      {
        type: "stack",
        width: 15,
        header: {
          show: "top",
          popout: false,
        },
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
            componentState: undefined,
          } as ComponentItemConfig,
          {
            id: CellsComponent.typeName,
            componentType: CellsComponent.typeName,
            type: "component",
            title: "Datasets",
            header: {
              show: "top",
              popout: false,
            },
            isClosable: true,
            reorderEnabled: true,
            componentState: undefined,
          } as ComponentItemConfig,
          {
            id: GatingComponent.typeName,
            componentType: GatingComponent.typeName,
            type: "component",
            title: "Gates",
            header: {
              show: "top",
              popout: false,
            },
            isClosable: true,
            reorderEnabled: true,
            componentState: undefined,
          } as ComponentItemConfig,
          {
            id: PipelinesComponent.typeName,
            componentType: PipelinesComponent.typeName,
            type: "component",
            title: "Pipelines",
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
        type: "column",
        content: [
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
            type: "stack",
            header: {
              show: "top",
              popout: false,
            },
            content: [
              {
                id: ScatterComponent.typeName,
                componentType: ScatterComponent.typeName,
                type: "component",
                title: "Scatterplot",
                header: {
                  show: "top",
                  popout: false,
                },
                isClosable: true,
                reorderEnabled: true,
                componentState: undefined,
              } as ComponentItemConfig,
              {
                id: PcaComponent.typeName,
                componentType: PcaComponent.typeName,
                type: "component",
                title: "PCA",
                header: {
                  show: "top",
                  popout: false,
                },
                isClosable: true,
                reorderEnabled: true,
                componentState: undefined,
              } as ComponentItemConfig,
              {
                id: UmapComponent.typeName,
                componentType: UmapComponent.typeName,
                type: "component",
                title: "UMAP",
                header: {
                  show: "top",
                  popout: false,
                },
                isClosable: true,
                reorderEnabled: true,
                componentState: undefined,
              } as ComponentItemConfig,
              {
                id: TsneComponent.typeName,
                componentType: TsneComponent.typeName,
                type: "component",
                title: "tSNE",
                header: {
                  show: "top",
                  popout: false,
                },
                isClosable: true,
                reorderEnabled: true,
                componentState: undefined,
              } as ComponentItemConfig,
              {
                id: LeidenComponent.typeName,
                componentType: LeidenComponent.typeName,
                type: "component",
                title: "Clusters [Leiden]",
                header: {
                  show: "top",
                  popout: false,
                },
                isClosable: true,
                reorderEnabled: true,
                componentState: undefined,
              } as ComponentItemConfig,
              {
                id: LouvainComponent.typeName,
                componentType: LouvainComponent.typeName,
                type: "component",
                title: "Clusters [Louvain]",
                header: {
                  show: "top",
                  popout: false,
                },
                isClosable: true,
                reorderEnabled: true,
                componentState: undefined,
              } as ComponentItemConfig,
              {
                id: PipelineComponent.typeName,
                componentType: PipelineComponent.typeName,
                type: "component",
                title: "Pipeline",
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

export const DEFAULT_LAYOUTS: ILayout[] = [
  {
    uid: DEFAULT_LAYOUT_UID,
    name: "Default",
    config: defaultConfig,
    isDefault: true,
  },
];
