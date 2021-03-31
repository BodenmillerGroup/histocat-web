import React, { ForwardedRef, forwardRef } from "react";
import { PolygonLayer, TextLayer } from "@deck.gl/layers"; // eslint-disable-line import/no-extraneous-dependencies
import { forceSimulation } from "d3-force";
import { SelectableScatterplotLayer, getSelectionLayers } from "../../layers";
import { cellLayerDefaultProps, DEFAULT_COLOR } from "../utils";
import { createCellsQuadTree } from "../shared-spatial-scatterplot/quadtree";
import AbstractSpatialOrScatterplot, {
  AbstractSpatialOrScatterplotProps,
  AbstractSpatialOrScatterplotState,
} from "../shared-spatial-scatterplot/AbstractSpatialOrScatterplot";
import { forceCollideRects } from "../shared-spatial-scatterplot/force-collide-rects";
import DeckGL, { RGBAColor } from "deck.gl";
import { Cell, CellEntry } from "../../types";
import { Quadtree } from "d3";

const CELLS_LAYER_ID = "scatterplot";
const LABEL_FONT_FAMILY = "-apple-system, 'Helvetica Neue', Arial, sans-serif";
const NUM_FORCE_SIMULATION_TICKS = 100;
const LABEL_UPDATE_ZOOM_DELTA = 0.25;

// Default getter function props.
const makeDefaultGetCellPosition = (mapping: string) => (cellEntry: CellEntry) => {
  const { mappings } = cellEntry[1];
  if (!(mapping in mappings)) {
    const available = Object.keys(mappings)
      .map((s) => `"${s}"`)
      .join(", ");
    throw new Error(`Expected to find "${mapping}", but available mappings are: ${available}`);
  }
  const mappedCell = mappings[mapping];
  // The negative applied to the y-axis is because
  // graphics rendering has the y-axis positive going south.
  return [mappedCell[0], -mappedCell[1], 0];
};
const makeDefaultGetCellCoords = (mapping: string) => (cell: Cell) => cell.mappings[mapping];
const makeDefaultGetCellColors = (cellColors: Map<string, RGBAColor>) => (cellEntry: CellEntry) =>
  (cellColors && cellColors.get(cellEntry[0])) || DEFAULT_COLOR;
const makeDefaultGetCellIsSelected = (cellSelection: string[] | null) => (cellEntry: CellEntry) =>
  cellSelection ? cellSelection.includes(cellEntry[0]) : true; // If nothing is selected, everything is selected.

type ScatterplotProps = {
  theme: string;
  mapping: string;
  cellColors: Map<string, RGBAColor>;
  cellSelection: string[];
  cellFilter?: string[];
  setCellFilter(cellFilter: string[]): void;
  cellRadiusScale?: number;
  cellOpacity?: number;
  getCellCoords?(cell: Cell): number[];
  getCellPosition?(cellEntry: CellEntry): number[];
  getCellColor?(cellEntry: CellEntry): RGBAColor;
  getCellIsSelected?(cellEntry: CellEntry): boolean;
  setCellSelection(cellSelection: string[]): void;
  cellHighlight: string | null;
  setCellHighlight(cellId: string | null): void;
  onCellClick?(info: any): void;
  setComponentHover(): void;
  cellSetPolygons: any[];
  cellSetPolygonsVisible: boolean;
  cellSetLabelsVisible: boolean;
  cellSetLabelSize: number;
};

type ScatterplotState = {};

/**
 * React component which renders a scatterplot from cell data, typically tSNE or PCA.
 * @param {object} props
 * @param {string} props.uuid A unique identifier for this component.
 * @param {string} props.theme The current vitessce theme.
 * @param {object} props.viewState The deck.gl view state.
 * @param {function} props.setViewState Function to call to update the deck.gl view state.
 * @param {object} props.cells
 * @param {string} props.mapping The name of the coordinate mapping field,
 * for each cell, for example "PCA" or "t-SNE".
 * @param {Map} props.cellColors Mapping of cell IDs to colors.
 * @param {array} props.cellSelection Array of selected cell IDs.
 * @param {array} props.cellFilter Array of filtered cell IDs. By default, null.
 * @param {number} props.cellRadiusScale The value for `radiusScale` to pass
 * to the deck.gl cells ScatterplotLayer.
 * @param {number} props.cellOpacity The value for `opacity` to pass
 * to the deck.gl cells ScatterplotLayer.
 * @param {function} props.getCellCoords Getter function for cell coordinates
 * (used by the selection layer).
 * @param {function} props.getCellPosition Getter function for cell [x, y, z] position.
 * @param {function} props.getCellColor Getter function for cell color as [r, g, b] array.
 * @param {function} props.getCellIsSelected Getter function for cell layer isSelected.
 * @param {function} props.setCellSelection
 * @param {function} props.setCellHighlight
 * @param {function} props.updateViewInfo
 * @param {function} props.onToolChange Callback for tool changes
 * (lasso/pan/rectangle selection tools).
 * @param {function} props.onCellClick Getter function for cell layer onClick.
 */
class Scatterplot extends AbstractSpatialOrScatterplot<
  ScatterplotProps & AbstractSpatialOrScatterplotProps,
  ScatterplotState & AbstractSpatialOrScatterplotState
> {
  private cellsQuadTree: Quadtree<[number, number]> | null;
  private cellSetsForceSimulation;
  private cellSetsLabelPrevZoom: number | null;
  private cellSetsLayers: (PolygonLayer<any> | TextLayer<any>)[];

  constructor(props: ScatterplotProps & AbstractSpatialOrScatterplotProps) {
    super(props);

    // To avoid storing large arrays/objects
    // in React state, this component
    // uses instance variables.
    // All instance variables used in this class:
    this.cellsEntries = [];
    this.cellsQuadTree = null;
    this.cellsLayer = null;
    this.cellSetsForceSimulation = forceCollideRects();
    this.cellSetsLabelPrevZoom = null;
    this.cellSetsLayers = [];

    // Initialize data and layers.
    this.onUpdateCellsData();
    this.onUpdateCellsLayer();
    this.onUpdateCellSetsLayers();
  }

  createCellsLayer() {
    const { cellsEntries } = this;
    const {
      theme,
      mapping,
      getCellPosition = makeDefaultGetCellPosition(mapping),
      cellRadiusScale = 0.2,
      cellOpacity = 1.0,
      cellFilter,
      cellSelection,
      setCellHighlight,
      setComponentHover,
      getCellIsSelected = makeDefaultGetCellIsSelected(
        cellsEntries.length === cellSelection.length ? null : cellSelection
      ),
      cellColors,
      getCellColor = makeDefaultGetCellColors(cellColors),
      onCellClick,
    } = this.props;
    const filteredCellsEntries = cellFilter
      ? cellsEntries.filter((cellEntry) => cellFilter.includes(cellEntry[0]))
      : cellsEntries;

    return new SelectableScatterplotLayer({
      id: CELLS_LAYER_ID,
      backgroundColor: theme === "dark" ? [0, 0, 0] : [241, 241, 241],
      isSelected: getCellIsSelected,
      opacity: cellOpacity,
      radiusScale: cellRadiusScale,
      radiusMinPixels: 1,
      radiusMaxPixels: 10,
      getPosition: getCellPosition,
      getColor: getCellColor,
      getLineWidth: 0,
      onClick: (info: any) => {
        if (onCellClick) {
          onCellClick(info);
        }
      },
      ...cellLayerDefaultProps(filteredCellsEntries, undefined, setCellHighlight, setComponentHover),
    });
  }

  createCellSetsLayers() {
    const {
      theme,
      cellSetPolygons,
      viewState,
      cellSetPolygonsVisible,
      cellSetLabelsVisible,
      cellSetLabelSize,
    } = this.props;

    const result = [];

    if (cellSetPolygonsVisible) {
      result.push(
        new PolygonLayer({
          id: "cell-sets-polygon-layer",
          data: cellSetPolygons,
          stroked: true,
          filled: false,
          wireframe: true,
          lineWidthMaxPixels: 1,
          getPolygon: (d: any) => d.hull,
          getLineColor: (d) => d.color,
          getLineWidth: 1,
        })
      );
    }

    if (cellSetLabelsVisible) {
      const { zoom } = viewState;
      const nodes = cellSetPolygons.map((p: any) => ({
        x: p.centroid[0],
        y: p.centroid[1],
        label: p.name,
      }));

      const collisionForce = this.cellSetsForceSimulation.size((d: any) => [
        ((cellSetLabelSize * 1) / 2 ** zoom) * 4 * d.label.length,
        ((cellSetLabelSize * 1) / 2 ** zoom) * 1.5,
      ]);

      forceSimulation().nodes(nodes).force("collision", collisionForce).tick(NUM_FORCE_SIMULATION_TICKS);

      result.push(
        new TextLayer({
          id: "cell-sets-text-layer",
          data: nodes,
          getPosition: (d: any) => [d.x, d.y],
          getText: (d) => d.label,
          getColor: theme === "dark" ? [255, 255, 255] : [0, 0, 0],
          getSize: cellSetLabelSize,
          getAngle: 0,
          getTextAnchor: "middle",
          getAlignmentBaseline: "center",
          fontFamily: LABEL_FONT_FAMILY,
          fontWeight: "normal",
        })
      );
    }

    return result;
  }

  createSelectionLayers() {
    const { viewState, mapping, getCellCoords = makeDefaultGetCellCoords(mapping), setCellSelection } = this.props;
    const { tool } = this.state;
    const { cellsQuadTree } = this;
    const flipYTooltip = true;
    return getSelectionLayers(
      tool,
      viewState.zoom,
      CELLS_LAYER_ID,
      getCellCoords,
      setCellSelection,
      cellsQuadTree,
      flipYTooltip
    );
  }

  getLayers() {
    const { cellsLayer, cellSetsLayers } = this;
    return [cellsLayer, ...cellSetsLayers, ...this.createSelectionLayers()];
  }

  onUpdateCellsData() {
    const { cells = {}, mapping, getCellCoords = makeDefaultGetCellCoords(mapping) } = this.props;
    const cellsEntries = Object.entries(cells);
    this.cellsEntries = cellsEntries;
    this.cellsQuadTree = createCellsQuadTree(cellsEntries, getCellCoords);
  }

  onUpdateCellsLayer() {
    this.cellsLayer = this.createCellsLayer();
  }

  onUpdateCellSetsLayers(onlyViewStateChange?: any) {
    // Because the label sizes for the force simulation depend on the zoom level,
    // we _could_ run the simulation every time the zoom level changes.
    // However, this has a performance impact in firefox.
    if (onlyViewStateChange) {
      const { viewState, cellSetLabelsVisible } = this.props;
      const { zoom } = viewState;
      const { cellSetsLabelPrevZoom } = this;
      // Instead, we can just check if the zoom level has changed
      // by some relatively large delta, to be more conservative
      // about re-running the force simulation.
      if (
        cellSetLabelsVisible &&
        (cellSetsLabelPrevZoom === null || Math.abs(cellSetsLabelPrevZoom - zoom) > LABEL_UPDATE_ZOOM_DELTA)
      ) {
        this.cellSetsLayers = this.createCellSetsLayers();
        this.cellSetsLabelPrevZoom = zoom;
      }
    } else {
      // Otherwise, something more substantial than just
      // the viewState has changed, such as the label array
      // itself, so we always want to update the layer
      // in this case.
      this.cellSetsLayers = this.createCellSetsLayers();
    }
  }

  viewInfoDidUpdate() {
    const { mapping, getCellPosition = makeDefaultGetCellPosition(mapping) } = this.props;
    super.viewInfoDidUpdate((cell: Cell) => getCellPosition([null as any, cell]));
  }

  /**
   * Here, asynchronously check whether props have
   * updated which require re-computing memoized variables,
   * followed by a re-render.
   * This function does not follow React conventions or paradigms,
   * it is only implemented this way to try to squeeze out
   * performance.
   * @param {object} prevProps The previous props to diff against.
   */
  componentDidUpdate(prevProps: any) {
    this.viewInfoDidUpdate();

    // @ts-ignore
    const shallowDiff = (propName: string) => prevProps[propName] !== this.props[propName];
    if (["cells"].some(shallowDiff)) {
      // Cells data changed.
      this.onUpdateCellsData();
      this.forceUpdate();
    }

    if (["cells", "cellFilter", "cellSelection", "cellColors", "cellRadiusScale"].some(shallowDiff)) {
      // Cells layer props changed.
      this.onUpdateCellsLayer();
      this.forceUpdate();
    }
    if (["cellSetPolygons", "cellSetPolygonsVisible", "cellSetLabelsVisible", "cellSetLabelSize"].some(shallowDiff)) {
      // Cell sets layer props changed.
      this.onUpdateCellSetsLayers(false);
      this.forceUpdate();
    }
    if (shallowDiff("viewState")) {
      // The viewState prop has changed (due to zoom or pan).
      this.onUpdateCellSetsLayers(true);
      this.forceUpdate();
    }
  }

  // render() is implemented in the abstract parent class.
}

/**
 * Need this wrapper function here,
 * since we want to pass a forwardRef
 * so that outer components can
 * access the grandchild DeckGL ref,
 * but we are using a class component.
 */
const ScatterplotWrapper = forwardRef(
  (props: ScatterplotProps & AbstractSpatialOrScatterplotProps, deckRef: ForwardedRef<DeckGL>) => (
    <Scatterplot {...props} deckRef={deckRef} />
  )
);
export default ScatterplotWrapper;
