import create from "zustand";
import { ILayout, ViewId } from "./models";
import { MosaicNode } from "react-mosaic-component";

const PROJECT_LAYOUTS_STORAGE_KEY = "ProjectLayouts";

export const TITLE_MAP: Record<ViewId, string> = {
  a: "Slides",
  b: "Image View",
  c: "Channels",
  new: "New Window",
};

const DEFAULT_LAYOUTS: ILayout[] = [
  {
    name: "Default",
    node: {
      direction: "row",
      first: "a",
      splitPercentage: 20,
      second: {
        direction: "row",
        first: "b",
        second: "c",
        splitPercentage: 80,
      },
    },
    isDefault: true,
  },
];

type LayoutsState = {
  layouts: ILayout[];
  activeLayout: ILayout;
  activeNode: MosaicNode<ViewId> | null;

  setActiveLayout(layout: ILayout): void;
  setActiveNode(node: MosaicNode<ViewId> | null): void;
  addLayout(name: string): void;
  deleteLayout(name: string): void;
  loadLayout(name: string): void;
};

let initialLayouts = DEFAULT_LAYOUTS;
if (localStorage.getItem(PROJECT_LAYOUTS_STORAGE_KEY)) {
  initialLayouts = JSON.parse(localStorage.getItem(PROJECT_LAYOUTS_STORAGE_KEY)!);
}

export const useLayoutsStore = create<LayoutsState>((set, get) => ({
  layouts:  initialLayouts,
  activeLayout: initialLayouts[0],
  activeNode: initialLayouts[0].node,

  setActiveLayout(layout: ILayout) {
    set({ activeLayout: layout });
  },

  setActiveNode(node: MosaicNode<ViewId> | null) {
    set({ activeNode: node });
  },

  addLayout(name: string) {
    console.log(get().activeNode)
    const activeNode = get().activeNode;
    if (activeNode) {
      const layout: ILayout = {
        name: name,
        node: activeNode,
        isDefault: false,
      };
      const layouts = get().layouts.concat(layout);
      set({ layouts: layouts });
      localStorage.setItem(PROJECT_LAYOUTS_STORAGE_KEY, JSON.stringify(layouts));
    }
  },

  deleteLayout(name: string) {
    const layouts = get().layouts.filter((item) => item.name !== name);
    set({ layouts: layouts, activeLayout: DEFAULT_LAYOUTS[0] });
    localStorage.setItem(PROJECT_LAYOUTS_STORAGE_KEY, JSON.stringify(layouts));
  },

  loadLayout(name: string) {
    const layout = get().layouts.find((item) => item.name === name);
    if (layout) {
      get().setActiveLayout(layout);
    }
  },
}));
