import { Mutations } from "vuex-smart-module";
import { WorldState } from ".";
import {IWorld} from "@/modules/world/models";

export class WorldMutations extends Mutations<WorldState> {
  setWorld(value: IWorld | null) {
    this.state.world = value;
  }
}
