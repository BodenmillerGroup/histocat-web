import { MosaicNode } from "react-mosaic-component";

export type ViewId = "slides" | "image" | "channels" | "settings" | "new";

export interface ILayout {
  name: string;
  node: MosaicNode<ViewId> | null;
  isDefault: boolean;
}
