import create from "zustand";
import { ILayout } from "./models";

const DEFAULT_LAYOUTS: ILayout[] = [
  {
    name: "Default",
    node: {
      direction: "row",
      first: "a",
      second: {
        direction: "column",
        first: "b",
        second: "c",
      },
    },
    isDefault: true,
  },
];

type LayoutsState = {
  layouts: ILayout[];

  addLayout(name: string, node: any): void;
  deleteLayout(name: string): void;
  getLayout(name: string): void;
};

export const useLayoutsStore = create<LayoutsState>((set, get) => ({
  layouts: DEFAULT_LAYOUTS,

  addLayout(name: string, node: any) {
    const layout: ILayout = {
      name: name,
      node: node,
      isDefault: false,
    };
    set({ layouts: get().layouts.concat(layout) });
  },

  deleteLayout(name: string) {
    const layouts = get().layouts.filter(item => item.name !== name);
    set({ layouts: layouts });
  },

  getLayout(name: string) {
    const layout = get().layouts.find(item => item.name === name);
    return layout;
  },
}));
