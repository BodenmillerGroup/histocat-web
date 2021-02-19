import { MosaicNode } from "react-mosaic-component";

export type ViewId = "a" | "b" | "c" | "new";

export interface ILayout {
  name: string;
  node: MosaicNode<ViewId> | null;
  isDefault: boolean;
}
