import { Dataframe } from "@/cellxgene/util/dataframe";
import { IUniverse } from "@/modules/universe/models";

export interface IWorld extends IUniverse {
  // schema/version related
  schema: any;

  nObs: number;
  nVar: number;
  clipQuantiles: { min: number; max: number };

  // Annotations
  obsAnnotations: Dataframe;
  varAnnotations: Dataframe;

  // Layout of graph
  obsLayout: Dataframe;

  // Var data columns - subset of all data (may be empty)
  varData: Dataframe;

  // Unclipped dataframes - subset, but not value clipped
  unclipped: {
    obsAnnotations: Dataframe;
    varData: Dataframe;
  };
}
