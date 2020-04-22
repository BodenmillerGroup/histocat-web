import { unassignedCategoryLabel } from "../../globals";
import { decodeMatrixFBS } from "./matrix";
import { isFpTypedArray } from "../typeHelpers";
import { indexEntireSchema } from "./schemaHelpers";
import catLabelSort from "../catLabelSort";
import { Dataframe, KeyIndex } from "../dataframe";
import { IUniverse } from "@/modules/universe/models";

/*
Private helper function - create and return a template Universe
*/
export function templateUniverse() {
  /* default universe template */
  return {
    nObs: 0,
    nVar: 0,
    schema: {},

    /*
    annotations
    */
    obsAnnotations: Dataframe.empty(),
    varAnnotations: Dataframe.empty(),
    /*
    layout
    */
    obsLayout: Dataframe.empty(),

    /*
    Var data columns - subset of all
    */
    varData: Dataframe.empty(null, new KeyIndex()),
  } as IUniverse;
}

/*
This module implements functions that support storage of "Universe",
aka all of the var/obs data and annotations.

These functions are used exclusively by the actions and reducers to
build an internal POJO for use by the rendering components.
*/

function promoteTypedArray(o) {
  /*
  Decide what internal data type to use for the data returned from
  the server.

  TODO - future optimization: not all int32/uint32 data series require
  promotion to float64.  We COULD simply look at the data to decide.
  */
  if (isFpTypedArray(o) || Array.isArray(o)) return o;

  let TyepdArrayCtor;
  switch (o.constructor) {
    case Int8Array:
    case Uint8Array:
    case Uint8ClampedArray:
    case Int16Array:
    case Uint16Array:
      TyepdArrayCtor = Float32Array;
      break;

    case Int32Array:
    case Uint32Array:
      TyepdArrayCtor = Float64Array;
      break;

    default:
      throw new Error("Unexpected data type returned from server.");
  }
  if (o.constructor === TyepdArrayCtor) return o;
  return new TyepdArrayCtor(o);
}

export function matrixFBSToDataframe(arrayBuffers) {
  /*
  Convert array of Matrix FBS to a Dataframe.

  The application has strong assumptions that all scalar data will be
  stored as a float32 or float64 (regardless of underlying data types).
  For example, clipping of value ranges (eg, user-selected percentiles)
  depends on the ability to use NaN in any numeric type.

  All float data from the server is left as is.  All non-float is promoted
  to an appropriate float.
  */
  if (!Array.isArray(arrayBuffers)) {
    arrayBuffers = [arrayBuffers];
  }
  if (arrayBuffers.length === 0) {
    return Dataframe.empty();
  }

  const fbs = arrayBuffers.map((ab) => decodeMatrixFBS(ab, true)); // leave in place
  /* check that all FBS have same row dimensionality */
  const { nRows } = fbs[0];
  fbs.forEach((b) => {
    if (b.nRows !== nRows) throw new Error("FBS with inconsistent dimensionality");
  });
  const columns = fbs
    .map((fb) =>
      fb.columns.map((c) => {
        if (isFpTypedArray(c) || Array.isArray(c)) return c;
        return promoteTypedArray(c);
      })
    )
    .flat();
  const colIdx = fbs.map((b) => b.colIdx).flat();
  const nCols = columns.length;

  const df = new Dataframe([nRows, nCols], columns, null, new KeyIndex(colIdx));
  return df;
}

export function createUniverseFromResponse(configResponse, schemaResponse) {
  /*
  build & return universe from a REST 0.2 /config, /schema and /annotations/obs response
  */
  const { schema } = schemaResponse;
  const universe = templateUniverse();

  /* schema related */
  universe.schema = schema;
  universe.nObs = schema.dataframe.nObs;
  universe.nVar = schema.dataframe.nVar;

  /* add defaults, as we can't assume back-end will fully populate schema */
  if (!schema.layout.var) schema.layout.var = [];
  if (!schema.layout.obs) schema.layout.obs = [];
  indexEntireSchema(universe.schema);
  normalizeEntireSchema(universe.schema);

  return universe as IUniverse;
}

function normalizeSchemaCategory(colSchema, col: any = undefined) {
  const { type, writable } = colSchema;
  if (type === "string" || type === "boolean" || type === "categorical") {
    let categories = [...new Set([...(colSchema.categories ?? []), ...(col?.summarize?.().categories ?? [])])];
    if (writable && categories.indexOf(unassignedCategoryLabel) === -1) {
      categories = categories.concat(unassignedCategoryLabel);
    }
    colSchema.categories = categories;
  } else if (writable) {
    throw new Error("Writable continuous obs annotations are not supported - failed to load");
  }

  if (colSchema.categories) {
    colSchema.categories = catLabelSort(writable, colSchema.categories);
  }
}

function normalizeEntireSchema(schema) {
  // currently only needed for obsAnnotations
  schema.annotations.obs.columns.forEach((colSchema) => normalizeSchemaCategory(colSchema));
}

export function addObsAnnotations(universe: IUniverse, df: Dataframe | Dataframe[]) {
  const obsAnnotations = universe.obsAnnotations.withColsFromAll(df);
  if (universe.nObs !== obsAnnotations.length) {
    throw new Error("Universe dimensionality mismatch - failed to load");
  }

  // for all of the new data, reconcile with schema and sort categories.
  const dfs = Array.isArray(df) ? df : [df];
  const keys = dfs.map((d) => d.colIndex.keys()).flat();
  const { schema } = universe;
  keys.forEach((k) => {
    const colSchema = schema.annotations.obsByName[k];
    const col = obsAnnotations.col(k);
    normalizeSchemaCategory(colSchema, col);
  });

  return { obsAnnotations, schema };
}

export function addVarAnnotations(universe: IUniverse, df: Dataframe | Dataframe[]) {
  const varAnnotations = universe.varAnnotations.withColsFromAll(df);
  if (universe.nVar !== varAnnotations.length) {
    throw new Error("Universe dimensionality mismatch - failed to load");
  }
  return { varAnnotations };
}

export function addObsLayout(universe: IUniverse, df: Dataframe | Dataframe[]) {
  const obsLayout = universe.obsLayout.withColsFromAll(df);
  if (universe.nObs !== obsLayout.length) {
    throw new Error("Universe dimensionality mismatch - failed to load");
  }
  return { obsLayout };
}

export function convertDataFBStoObject(universe: IUniverse, arrayBuffer) {
  /*
  /data/var returns a flatbuffer (FBS) as described by cellxgene/fbs/matrix.fbs

  This routine converts the binary wire encoding into a JS object:

  {
    gene: Float32Array,
    ...
  }
  */
  const fbs = decodeMatrixFBS(arrayBuffer);
  const { colIdx, columns } = fbs;
  const result = {};

  if (!columns.every(isFpTypedArray)) {
    // We have strong assumptions that all var data is float
    throw new Error("Unexpected non-floating point response from server.");
  }

  const varIndexName = universe.schema.annotations.var.index;
  for (let c = 0; c < colIdx.length; c += 1) {
    const varName = universe.varAnnotations.at(colIdx[c], varIndexName);
    result[varName] = columns[c];
  }
  return result;
}
