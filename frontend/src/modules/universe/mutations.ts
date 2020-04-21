import { Mutations } from "vuex-smart-module";
import { UniverseState } from ".";
import { IUniverse } from "@/modules/universe/models";

export class UniverseMutations extends Mutations<UniverseState> {
  setUniverse(value: IUniverse) {
    this.state.universe = value;
  }
}
