import { Dataframe } from "@/cellxgene/util/dataframe";

export interface IUniverse {
  // schema/version related
  schema: any;

  nObs: number;
  nVar: number;

  // Annotations
  obsAnnotations: Dataframe;
  varAnnotations: Dataframe;

  // Layout
  obsLayout: Dataframe;

  // Var data columns - subset of all
  varData: Dataframe;
}
