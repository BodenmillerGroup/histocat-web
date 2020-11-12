import { Mutations } from "vuex-smart-module";
import { CentroidsState } from ".";
import { ICentroidsData } from "./models";
import { BroadcastManager } from "@/utils/BroadcastManager";
import { SET_CENTROIDS } from "./events";
import { ICellPoint } from "@/modules/results/models";

export class CentroidsMutations extends Mutations<CentroidsState> {
  constructor() {
    super();
    BroadcastManager.subscribe(SET_CENTROIDS, (payload) => this.setCentroids(payload));
  }

  setCentroids(payload: ICentroidsData) {
    const newState = new Map<number, ICellPoint[]>();
    payload.acquisitionIds.forEach((acquisitionId, i) => {
      if (!newState.has(acquisitionId)) {
        newState.set(acquisitionId, []);
      }
      const cellPoint: ICellPoint = {
        acquisitionId: acquisitionId,
        cellId: payload.cellIds[i],
        objectNumber: payload.objectNumbers[i],
        x: payload.x[i],
        y: payload.y[i],
        color: 0
      };
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
