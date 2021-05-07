import { LayoutConfig } from "golden-layout";

export interface ILayout {
  uid: string;
  name: string;
  config: LayoutConfig;
  isDefault: boolean;
}

export interface IResponsive {
  width: number | null;
  height: number | null;
}
