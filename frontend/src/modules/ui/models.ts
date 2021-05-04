import { LayoutConfig } from "golden-layout";

export type ViewMode = "image" | "segmentation" | "data";

export interface ILayout {
    name: string;
    config: LayoutConfig;
}

export interface IResponsive {
  width: number | null;
  height: number | null;
}
