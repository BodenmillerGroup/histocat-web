import { mainModule } from "@/modules/main";
import { Store } from "vuex";
import { Actions, Context } from "vuex-smart-module";
import { GraphState } from ".";
import { GraphGetters } from "./getters";
import { GraphMutations } from "./mutations";
import { graphSelectionModule } from "@/modules/graphSelection";
import { crossfilterModule } from "@/modules/crossfilter";

export class GraphActions extends Actions<GraphState, GraphGetters, GraphMutations, GraphActions> {
  main?: Context<typeof mainModule>;
  graphSelection?: Context<typeof graphSelectionModule>;
  crossfilter?: Context<typeof crossfilterModule>;

  // Called after the module is initialized
  $init(store: Store<any>): void {
    this.main = mainModule.context(store);
    this.graphSelection = graphSelectionModule.context(store);
    this.crossfilter = crossfilterModule.context(store);
  }

  graphBrushStart() {
    // "graph brush start"
  }

  graphBrushChange(brushCoords) {
    // "graph brush change"
    this.graphSelection?.mutations.setSelection({
      mode: "within-rect",
      brushCoords: brushCoords,
    });
  }

  graphBrushEnd(brushCoords) {
    // "graph brush end"
    this.graphBrushChange(brushCoords);
  }

  graphBrushDeselect() {
    // "graph brush deselect"
  }

  graphLassoStart() {
    // "graph lasso start"
  }

  graphLassoDeselect() {
    // "graph lasso deselect"
  }

  graphLassoEnd(polygon) {
    // "graph lasso end"
  }

  graphLassoCancel() {
    // "graph lasso cancel"
  }

  changeOpacityDeselectedCellsIn2dGraphBackground(data) {
    // "change opacity deselected cells in 2d graph background"
  }
}
