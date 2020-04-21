import { Mutations } from "vuex-smart-module";
import { CrossfilterState } from ".";
import ImmutableTypedCrossfilter from "@/cellxgene/util/typedCrossfilter";

export class CrossfilterMutations extends Mutations<CrossfilterState> {
  setCrossfilter(value: ImmutableTypedCrossfilter | null) {
    this.state.crossfilter = value;
  }
}
