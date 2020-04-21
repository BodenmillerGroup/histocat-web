import { Dataframe } from "@/cellxgene/util/dataframe";

export interface IUniverse {
  nObs: number;
  nVar: number;
  schema: any;

  // Annotations
  obsAnnotations: Dataframe;
  varAnnotations: Dataframe;

  // Layout
  obsLayout: Dataframe;

  // Var data columns - subset of all
  varData: Dataframe;
}
