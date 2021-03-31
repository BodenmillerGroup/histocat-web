export type CoordinationScopes = { [p: string]: string }

export type Cell = {
  factors: { [p: string]: string };
  genes: { [p: string]: number };
  mappings: { [p: string]: [number, number] };
  poly: [number, number][];
  xy: [number, number];
}
export type CellEntry = [string, Cell];

export type Molecules = { [p: string]: [number, number][] };
export type MoleculeEntry = [number, number, number, string];

export type Neighborhood = { poly: [number, number][] };
export type NeighborhoodEntry = [string, Neighborhood];
export type Neighborhoods = { [p: string]: Neighborhood };

export type ViewState = {
  zoom: number;
  target: number[];
}

export type ViewInfo = {
  uuid: string;
  project: any;
}

export type CellsLayerDef = {
  radius: number;
  stroked: number;
  visible: boolean;
  opacity: number;
}

export type MoleculesLayerDef = {
  radius: number;
  visible: boolean;
  opacity: number;
}

export type NeighborhoodsLayerDef = {
  visible: boolean;
}

export type RawLayerDef = {
  index: number;
  visible: boolean;
  opacity: number;
  colormap: string;
  channels: any[];
  transparentColor: any;
  modelMatrix: any;
}

export type Layout = {
  component: string;
  x: number;
  y: number;
  w: number;
  h: number;
  coordinationScopes?: any;
  props?: any;
}[];
