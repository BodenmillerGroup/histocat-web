import React from "react";
import { COORDINATE_SYSTEM, RGBAColor } from "deck.gl";
import { PRIMARY_CARD } from "./classNames";
import { HIERARCHICAL_SCHEMAS, SETS_DATATYPE_CELL } from "./sets/constants";
import { CellEntry } from "../types";

export function makeCellStatusMessage(cellInfoFactors: any[]) {
  return Object.entries(cellInfoFactors)
    .map(([factor, value]) => `${factor}: ${value}`)
    .join("; ");
}

export function cellLayerDefaultProps(
  cells: CellEntry[],
  updateStatus: any,
  setCellHighlight: (cellId: string | null) => void,
  setComponentHover: () => void,
) {
  return {
    coordinateSystem: COORDINATE_SYSTEM.CARTESIAN,
    data: cells,
    pickable: true,
    autoHighlight: true,
    stroked: true,
    filled: true,
    getElevation: 0,
    onHover: (info: any) => {
      // Notify the parent component that its child component is
      // the "hover source".
      if (setComponentHover) {
        setComponentHover();
      }
      if (info.object) {
        const [cellId, cellInfo] = info.object;
        const { factors = {} } = cellInfo;
        if (updateStatus) {
          updateStatus(makeCellStatusMessage(factors));
        }
        if (setCellHighlight) {
          setCellHighlight(cellId);
        }
      } else if (setCellHighlight) {
        // Clear the currently-hovered cell info by passing null.
        setCellHighlight("");
      }
    },
  };
}

export const DEFAULT_COLOR = [128, 128, 128];

// From https://personal.sron.nl/~pault/#sec:qualitative
export const PALETTE = [
  [68, 119, 170],
  [136, 204, 238],
  [68, 170, 153],
  [17, 119, 51],
  [153, 153, 51],
  [221, 204, 119],
  [204, 102, 119],
  [136, 34, 85],
  [170, 68, 153],
] as RGBAColor[];

export const VIEWER_PALETTE = [
  [0, 0, 255],
  [0, 255, 0],
  [255, 0, 255],
  [255, 255, 0],
  [0, 255, 255],
  [255, 255, 255],
  [255, 128, 0],
  [255, 0, 0],
] as RGBAColor[];

export const COLORMAP_OPTIONS = [
  "viridis",
  "greys",
  "magma",
  "jet",
  "hot",
  "bone",
  "copper",
  "summer",
  "density",
  "inferno",
];

export const DEFAULT_GL_OPTIONS = { webgl2: true } as WebGLContextAttributes;

export function createDefaultUpdateStatus(componentName: string) {
  return (message: string) => console.warn(`${componentName} updateStatus: ${message}`);
}

export function createDefaultUpdateCellsSelection(componentName: string) {
  return (cellsSelection: any) => console.warn(`${componentName} updateCellsSelection: ${cellsSelection}`);
}

export function createDefaultUpdateCellsHover(componentName: string) {
  return (hoverInfo: { cellId: any }) => console.warn(`${componentName} updateCellsHover: ${hoverInfo.cellId}`);
}

export function createDefaultUpdateGenesHover(componentName: string) {
  return (hoverInfo: { geneId: any }) => console.warn(`${componentName} updateGenesHover: ${hoverInfo.geneId}`);
}

export function createDefaultUpdateViewInfo(componentName: string) {
  return (viewInfo: any) => console.warn(`${componentName} updateViewInfo: ${viewInfo}`);
}

export function createDefaultClearPleaseWait() {
  return () => {};
}

/**
 * Copy a typed array into a new array buffer.
 * @param {Uint8Array} arr The typed array to be copied.
 * @returns {Uint8Array} The copied array.
 */
export function copyUint8Array(arr: Uint8Array) {
  const newBuffer = new ArrayBuffer(arr.buffer.byteLength);
  const newArr = new Uint8Array(newBuffer);
  newArr.set(arr);
  return newArr;
}

export function getNextNumberedNodeName(nodes: any[], prefix: string) {
  let i = 1;
  if (nodes) {
    // eslint-disable-next-line no-loop-func
    while (nodes.find((n) => n.name === `${prefix}${i}`)) {
      // eslint-disable-next-line no-plusplus
      i++;
    }
  }
  return `${prefix}${i}`;
}

/**
 * Create a new selected cell set based on a cell selection.
 * @param {string[]} cellSelection An array of cell IDs.
 * @param {object[]} additionalCellSets The previous array of user-defined cell sets.
 * @param {function} setCellSetSelection The setter function for cell set selections.
 * @param {function} setAdditionalCellSets The setter function for user-defined cell sets.
 */
export function setCellSelection(
  cellSelection: string[],
  additionalCellSets: any,
  cellSetColor: any,
  setCellSetSelection: any,
  setAdditionalCellSets: any,
  setCellSetColor: any,
  setCellColorEncoding: any,
  prefix = "Selection "
) {
  const CELL_SELECTIONS_LEVEL_ZERO_NAME = "My Selections";

  const selectionsLevelZeroNode = additionalCellSets?.tree.find((n: any) => n.name === CELL_SELECTIONS_LEVEL_ZERO_NAME);
  const nextAdditionalCellSets = {
    version: HIERARCHICAL_SCHEMAS[SETS_DATATYPE_CELL].latestVersion,
    datatype: SETS_DATATYPE_CELL,
    tree: [...(additionalCellSets ? additionalCellSets.tree : [])],
  };

  const nextName = getNextNumberedNodeName(selectionsLevelZeroNode?.children, prefix);
  let colorIndex = 0;
  if (selectionsLevelZeroNode) {
    colorIndex = selectionsLevelZeroNode.children.length;
    selectionsLevelZeroNode.children.push({
      name: nextName,
      set: cellSelection.map((d) => [d, null]),
    });
  } else {
    nextAdditionalCellSets.tree.push({
      name: CELL_SELECTIONS_LEVEL_ZERO_NAME,
      children: [
        {
          name: nextName,
          set: cellSelection.map((d) => [d, null]),
        },
      ],
    });
  }
  setAdditionalCellSets(nextAdditionalCellSets);
  const nextPath = ["My Selections", nextName];
  setCellSetColor([
    ...(cellSetColor || []),
    {
      path: nextPath,
      color: PALETTE[colorIndex % PALETTE.length],
    },
  ]);
  setCellSetSelection([nextPath]);
  setCellColorEncoding("cellSetSelection");
}

export function mergeCellSets(cellSets: any, additionalCellSets: any) {
  return {
    version: HIERARCHICAL_SCHEMAS[SETS_DATATYPE_CELL].latestVersion,
    datatype: SETS_DATATYPE_CELL,
    tree: [...(cellSets ? cellSets.tree : []), ...(additionalCellSets ? additionalCellSets.tree : [])],
  };
}

export function createWarningComponent(props: { title: string; message: string }) {
  return () => {
    const { title, message } = props;
    return (
      <div className={PRIMARY_CARD}>
        <h1>{title}</h1>
        <div>{message}</div>
      </div>
    );
  };
}

export function asEsModule(component: any) {
  return {
    __esModule: true,
    default: component,
  };
}
