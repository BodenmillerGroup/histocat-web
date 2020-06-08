import { Mutations } from "vuex-smart-module";
import { CentroidsState } from ".";
import { ICentroidsData } from "./models";
import { CellPoint } from "@/data/CellPoint";

export class CentroidsMutations extends Mutations<CentroidsState> {
  setCentroids(payload: ICentroidsData) {
    const newState = new Map<number, CellPoint[]>();
    payload.acquisitionIds.forEach((acquisitionId, i) => {
      if (!newState.has(acquisitionId)) {
        newState.set(acquisitionId, []);
      }
      const cellPoint = new CellPoint(acquisitionId, payload.cellIds[i], payload.x[i], payload.y[i], i);
      newState.get(acquisitionId)!.push(Object.freeze(cellPoint));
    });
    this.state.centroids = newState;
  }

  reset() {
    // acquire initial state
    const s = new CentroidsState();
    Object.keys(s).forEach((key) => {
      this.state[key] = s[key];
    });
  }
}
