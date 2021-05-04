import { LayoutConfig } from "golden-layout";

export interface ILayout {
  name: string;
  config: LayoutConfig;
}

export interface IResponsive {
  width: number | null;
  height: number | null;
}
