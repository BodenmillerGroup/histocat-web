import { useState, useEffect } from "react";
import equal from "fast-deep-equal";
import { capitalize } from "../utils";
import { useSetWarning } from "../app/state/hooks";
import { AbstractLoaderError, LoaderNotFoundError } from "../loaders/errors/index";
import { DEFAULT_MOLECULES_LAYER, DEFAULT_CELLS_LAYER, DEFAULT_NEIGHBORHOODS_LAYER } from "./spatial/constants";
import { DEFAULT_COORDINATION_VALUES } from "../app/state/coordination";
import { Cell } from "../types";

/**
 * Warn via publishing to the console
 * and to the global warning store.
 * @param {AbstractLoaderError} error An error instance.
 */
function warn(error: any, setWarning: Function) {
  setWarning(error.message);
  console.warn(error.message);
  if (error instanceof AbstractLoaderError) {
    error.warnInConsole();
  }
}

/**
 * Initialize values in the coordination space.
 * @param {object} values Object where
 * keys are coordination type names,
 * values are initial coordination values.
 * @param {object} setters Object where
 * keys are coordination type names with the prefix 'set',
 * values are coordination setter functions.
 * @param {object} initialValues Object where
 * keys are coordination type names and keys are values.
 */
function initCoordinationSpace(values: any, setters: any, initialValues: any) {
  if (!values || !setters) {
    return;
  }
  Object.entries(values).forEach(([coordinationType, value]) => {
    const setterName = `set${capitalize(coordinationType)}`;
    const setterFunc = setters[setterName];
    const initialValue = initialValues && initialValues[coordinationType];
    const shouldInit = equal(initialValue, DEFAULT_COORDINATION_VALUES[coordinationType]);
    if (shouldInit && setterFunc) {
      setterFunc(value);
    }
  });
}

/**
 * Get the dataset description string.
 * @param {object} loaders The object mapping
 * datasets and data types to loader instances.
 * @param {string} dataset The key for a dataset,
 * used to identify which loader to use.
 * @returns {array} [description] where
 * description is a string.
 */
export function useDescription(loaders: any, dataset: string) {
  const [description, setDescription] = useState<string | null>();

  useEffect(() => {
    if (!loaders[dataset]) {
      return;
    }

    if (loaders[dataset].description) {
      setDescription(loaders[dataset].description);
    } else {
      setDescription(null);
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [loaders, dataset]);

  return [description];
}

/**
 * Get data from a cells data type loader,
 * updating "ready" and URL state appropriately.
 * Throw warnings if the data is marked as required.
 * Subscribe to loader updates.
 * @param {object} loaders The object mapping
 * datasets and data types to loader instances.
 * @param {string} dataset The key for a dataset,
 * used to identify which loader to use.
 * @param {function} setItemIsReady A function to call
 * when done loading.
 * @param {function} addUrl A function to call to update
 * the URL list.
 * @param {boolean} isRequired Should a warning be thrown if
 * loading is unsuccessful?
 * @param {object} coordinationSetters Object where
 * keys are coordination type names with the prefix 'set',
 * values are coordination setter functions.
 * @param {object} initialCoordinationValues Object where
 * keys are coordination type names with the prefix 'initialize',
 * values are initialization preferences as boolean values.
 * @returns {array} [cells, cellsCount] where
 * cells is an object and cellsCount is the
 * number of items in the cells object.
 */
export function useCellsData(
  loaders: any,
  dataset: string,
  setItemIsReady: (readyItem: any) => void,
  addUrl: Function,
  isRequired: boolean,
  coordinationSetters?: any,
  initialCoordinationValues?: any
) {
  const [cells, setCells] = useState<{ [p: string]: Cell }>({});
  const [cellsCount, setCellsCount] = useState<number>(0);

  const setWarning = useSetWarning();

  useEffect(() => {
    if (!loaders[dataset]) {
      return;
    }

    if (loaders[dataset].loaders.cells) {
      loaders[dataset].loaders.cells
        .load()
        .catch((e: any) => warn(e, setWarning))
        .then((payload: any) => {
          if (!payload) return;
          const { data, url, coordinationValues } = payload;
          setCells(data);
          setCellsCount(Object.keys(data).length);
          addUrl(url, "Cells");
          // This dataset has cells, so set up the
          // spatial cells layer coordination value
          // using the cell layer singleton.
          const coordinationValuesOrDefault = {
            spatialCellsLayer: DEFAULT_CELLS_LAYER,
            ...coordinationValues,
          };
          initCoordinationSpace(coordinationValuesOrDefault, coordinationSetters, initialCoordinationValues);
          setItemIsReady("cells");
        });
    } else {
      setCells({});
      setCellsCount(0);
      if (isRequired) {
        warn(new LoaderNotFoundError(dataset, "cells", null), setWarning);
      } else {
        setItemIsReady("cells");
      }
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [loaders, dataset]);

  return [cells, cellsCount] as [cells: { [p: string]: Cell }, cellsCount: number];
}

/**
 * Get data from a cell sets data type loader,
 * updating "ready" and URL state appropriately.
 * Throw warnings if the data is marked as required.
 * Subscribe to loader updates.
 * @param {object} loaders The object mapping
 * datasets and data types to loader instances.
 * @param {string} dataset The key for a dataset,
 * used to identify which loader to use.
 * @param {function} setItemIsReady A function to call
 * when done loading.
 * @param {function} addUrl A function to call to update
 * the URL list.
 * @param {boolean} isRequired Should a warning be thrown if
 * loading is unsuccessful?
 * @param {object} coordinationSetters Object where
 * keys are coordination type names with the prefix 'set',
 * values are coordination setter functions.
 * @param {object} initialCoordinationValues Object where
 * keys are coordination type names and values are the current values.
 * @returns {array} [cellSets] where
 * cellSets is a sets tree object.
 */
export function useCellSetsData(
  loaders: any,
  dataset: string,
  setItemIsReady: (readyItem: any) => void,
  addUrl: Function,
  isRequired: boolean,
  coordinationSetters: any,
  initialCoordinationValues: any
) {
  const [cellSets, setCellSets] = useState<any>();

  const setWarning = useSetWarning();

  useEffect(() => {
    if (!loaders[dataset]) {
      return;
    }

    if (loaders[dataset].loaders["cell-sets"]) {
      // Load the data initially.
      loaders[dataset].loaders["cell-sets"]
        .load()
        .catch((e: any) => warn(e, setWarning))
        .then((payload: any) => {
          if (!payload) return;
          const { data, url, coordinationValues } = payload;
          setCellSets(data);
          addUrl(url, "Cell Sets");
          initCoordinationSpace(coordinationValues, coordinationSetters, initialCoordinationValues);
          setItemIsReady("cell-sets");
        });
    } else {
      setCellSets(null);
      if (isRequired) {
        warn(new LoaderNotFoundError(dataset, "cell-sets", null), setWarning);
      } else {
        setItemIsReady("cell-sets");
      }
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [loaders, dataset]);

  return [cellSets];
}

/**
 * Get (potentially filtered) data from an expression matrix data type loader,
 * updating "ready" and URL state appropriately.
 * Throw warnings if the data is marked as required.
 * Subscribe to loader updates.  Should not be used in conjunction (called in the same component)
 * with useExpressionAttrs as this returns a potentially filtered set of attributes
 * specifically for the returned expression data.
 * @param {object} loaders The object mapping
 * datasets and data types to loader instances.
 * @param {string} dataset The key for a dataset,
 * used to identify which loader to use.
 * @param {function} setItemIsReady A function to call
 * when done loading.
 * @param {function} addUrl A function to call to update
 * the URL list.
 * @param {boolean} isRequired Should a warning be thrown if
 * loading is unsuccessful?
 * @param {object} coordinationSetters Object where
 * keys are coordination type names with the prefix 'set',
 * values are coordination setter functions.
 * @param {object} initialCoordinationValues Object where
 * keys are coordination type names with the prefix 'initialize',
 * values are initialization preferences as boolean values.
 * @returns {array} [expressionMatrix] where
 * expressionMatrix is an object with
 * shape { cols, rows, matrix }.
 */
export function useExpressionMatrixData(
  loaders: any,
  dataset: string,
  setItemIsReady: (readyItem: any) => void,
  addUrl: Function,
  isRequired: boolean,
  coordinationSetters: any,
  initialCoordinationValues: any
) {
  const [expressionMatrix, setExpressionMatrix] = useState<any>();

  const setWarning = useSetWarning();

  useEffect(() => {
    if (!loaders[dataset]) {
      return;
    }

    if (loaders[dataset].loaders["expression-matrix"]) {
      loaders[dataset].loaders["expression-matrix"]
        .load()
        .catch((e: any) => warn(e, setWarning))
        .then((payload: any) => {
          if (!payload) return;
          const { data, url, coordinationValues } = payload;
          const [attrs, arr] = data;
          setExpressionMatrix({
            cols: attrs.cols,
            rows: attrs.rows,
            matrix: arr.data,
          });
          addUrl(url, "Expression Matrix");
          initCoordinationSpace(coordinationValues, coordinationSetters, initialCoordinationValues);
          setItemIsReady("expression-matrix");
        });
    } else {
      setExpressionMatrix(null);
      if (isRequired) {
        warn(new LoaderNotFoundError(dataset, "expression-matrix", null), setWarning);
      } else {
        setItemIsReady("expression-matrix");
      }
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [loaders, dataset]);

  return [expressionMatrix];
}

/**
 * Get data from the expression matrix data type loader for a given gene selection.
 * Throw warnings if the data is marked as required.
 * Subscribe to loader updates.
 * @param {object} loaders The object mapping
 * datasets and data types to loader instances.
 * @param {string} dataset The key for a dataset,
 * used to identify which loader to use.
 * @param {function} setItemIsReady A function to call
 * when done loading.
 * @param {boolean} isRequired Should a warning be thrown if
 * loading is unsuccessful?
 * @param {boolean} selection A list of gene names to get expression data for.
 * @returns {array} [geneData] where geneData is an array [Uint8Array, ..., Uint8Array]
 * for however many genes are in the selection.
 */
export function useGeneSelection(loaders: any, dataset: string, setItemIsReady: (readyItem: any) => void, isRequired: boolean, selection: any) {
  const [geneData, setGeneData] = useState<any>();

  const setWarning = useSetWarning();

  useEffect(() => {
    if (!loaders[dataset]) {
      return;
    }
    if (!selection) {
      setItemIsReady("expression-matrix");
      return;
    }
    const loader = loaders[dataset].loaders["expression-matrix"];
    if (loader) {
      const implementsGeneSelection = typeof loader.loadGeneSelection === "function";
      if (implementsGeneSelection) {
        loaders[dataset].loaders["expression-matrix"]
          .loadGeneSelection({ selection })
          .catch((e: any) => warn(e, setWarning))
          .then((payload: any) => {
            if (!payload) return;
            const { data } = payload;
            setGeneData(data);
          });
      } else {
        loader
          .load()
          .catch((e: any) => warn(e, setWarning))
          .then((payload: any) => {
            if (!payload) return;
            const { data } = payload;
            const [attrs, { data: matrix }] = data;
            const expressionDataForSelection = selection.map((sel: any) => {
              const geneIndex = attrs.cols.indexOf(sel);
              const numGenes = attrs.cols.length;
              const numCells = attrs.rows.length;
              const expressionData = new Uint8Array(numCells);
              for (let cellIndex = 0; cellIndex < numCells; cellIndex += 1) {
                expressionData[cellIndex] = matrix[cellIndex * numGenes + geneIndex];
              }
              return expressionData;
            });
            setGeneData(expressionDataForSelection);
          });
      }
    } else {
      setGeneData(null);
      if (isRequired) {
        warn(new LoaderNotFoundError(dataset, "expression-matrix", null), setWarning);
      } else {
        setItemIsReady("expression-matrix");
      }
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [loaders, dataset, selection]);

  return [geneData];
}

/**
 * Get the attributes for the expression matrix data type loader,
 * i.e names of cells and genes.
 * Throw warnings if the data is marked as required.
 * Subscribe to loader updates.  Should not be used in conjunction (called in the same component)
 * with useExpressionMatrixData.
 * @param {object} loaders The object mapping
 * datasets and data types to loader instances.
 * @param {string} dataset The key for a dataset,
 * used to identify which loader to use.
 * @param {function} setItemIsReady A function to call
 * when done loading.
 * @param {function} addUrl A function to call to update
 * the URL list.
 * @param {boolean} isRequired Should a warning be thrown if
 * loading is unsuccessful?
 * @returns {object} [attrs] { rows, cols } object containing cell and gene names.
 */
export function useExpressionAttrs(loaders: any, dataset: string, setItemIsReady: (readyItem: any) => void, addUrl: Function, isRequired: boolean) {
  const [attrs, setAttrs] = useState<any>();

  const setWarning = useSetWarning();

  useEffect(() => {
    if (!loaders[dataset]) {
      return;
    }
    const loader = loaders[dataset].loaders["expression-matrix"];
    if (loader) {
      const implementsLoadAttrs = typeof loader.loadAttrs === "function";
      if (implementsLoadAttrs) {
        loader
          .loadAttrs()
          .catch((e: any) => warn(e, setWarning))
          .then((payload: any) => {
            if (!payload) return;
            const { data, url } = payload;
            setAttrs(data);
            addUrl(url, "Expression Matrix");
            setItemIsReady("expression-matrix");
          });
      } else {
        loader
          .load()
          .catch((e: any) => warn(e, setWarning))
          .then((payload: any) => {
            if (!payload) return;
            const { data, url } = payload;
            setAttrs(data[0]);
            addUrl(url, "Expression Matrix");
            setItemIsReady("expression-matrix");
          });
      }
    } else {
      setAttrs(null);
      if (isRequired) {
        warn(new LoaderNotFoundError(dataset, "expression-matrix", null), setWarning);
      } else {
        setItemIsReady("expression-matrix");
      }
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [loaders, dataset]);

  return [attrs];
}

/**
 * Get data from a molecules data type loader,
 * updating "ready" and URL state appropriately.
 * Throw warnings if the data is marked as required.
 * Subscribe to loader updates.
 * @param {object} loaders The object mapping
 * datasets and data types to loader instances.
 * @param {string} dataset The key for a dataset,
 * used to identify which loader to use.
 * @param {function} setItemIsReady A function to call
 * when done loading.
 * @param {function} addUrl A function to call to update
 * the URL list.
 * @param {boolean} isRequired Should a warning be thrown if
 * loading is unsuccessful?
 * @param {object} coordinationSetters Object where
 * keys are coordination type names with the prefix 'set',
 * values are coordination setter functions.
 * @param {object} initialCoordinationValues Object where
 * keys are coordination type names with the prefix 'initialize',
 * values are initialization preferences as boolean values.
 * @returns {array} [molecules, moleculesCount, locationsCount] where
 * molecules is an object,
 * moleculesCount is the number of unique molecule types, and
 * locationsCount is the number of molecules.
 */
export function useMoleculesData(
  loaders: any,
  dataset: string,
  setItemIsReady: (readyItem: any) => void,
  addUrl: Function,
  isRequired: boolean,
  coordinationSetters: any,
  initialCoordinationValues: any
) {
  const [molecules, setMolecules] = useState<any>();
  const [moleculesCount, setMoleculesCount] = useState<number>(0);
  const [locationsCount, setLocationsCount] = useState<number>(0);

  const setWarning = useSetWarning();

  useEffect(() => {
    if (!loaders[dataset]) {
      return;
    }

    if (loaders[dataset].loaders.molecules) {
      loaders[dataset].loaders.molecules
        .load()
        .catch((e: any) => warn(e, setWarning))
        .then((payload: any) => {
          if (!payload) return;
          const { data, url, coordinationValues } = payload;
          setMolecules(data);
          setMoleculesCount(Object.keys(data).length);
          setLocationsCount(
            Object.values(data)
              .map((l: any) => l.length)
              .reduce((a, b) => a + b, 0)
          );
          addUrl(url, "Molecules");
          const coordinationValuesOrDefault = {
            spatialMoleculesLayer: DEFAULT_MOLECULES_LAYER,
            ...coordinationValues,
          };
          initCoordinationSpace(coordinationValuesOrDefault, coordinationSetters, initialCoordinationValues);
          setItemIsReady("molecules");
        });
    } else {
      setMolecules({});
      setMoleculesCount(0);
      setLocationsCount(0);
      if (isRequired) {
        warn(new LoaderNotFoundError(dataset, "molecules", null), setWarning);
      } else {
        setItemIsReady("molecules");
      }
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [loaders, dataset]);

  return [molecules, moleculesCount, locationsCount];
}

/**
 * Get data from a neighborhoods data type loader,
 * updating "ready" and URL state appropriately.
 * Throw warnings if the data is marked as required.
 * Subscribe to loader updates.
 * @param {object} loaders The object mapping
 * datasets and data types to loader instances.
 * @param {string} dataset The key for a dataset,
 * used to identify which loader to use.
 * @param {function} setItemIsReady A function to call
 * when done loading.
 * @param {function} addUrl A function to call to update
 * the URL list.
 * @param {boolean} isRequired Should a warning be thrown if
 * loading is unsuccessful?
 * @param {object} coordinationSetters Object where
 * keys are coordination type names with the prefix 'set',
 * values are coordination setter functions.
 * @param {object} initialCoordinationValues Object where
 * keys are coordination type names with the prefix 'initialize',
 * values are initialization preferences as boolean values.
 * @returns {array} [neighborhoods] where
 * neighborhoods is an object.
 */
export function useNeighborhoodsData(
  loaders: any,
  dataset: string,
  setItemIsReady: (readyItem: any) => void,
  addUrl: Function,
  isRequired: boolean,
  coordinationSetters: any,
  initialCoordinationValues: any
) {
  const [neighborhoods, setNeighborhoods] = useState<any>();

  const setWarning = useSetWarning();

  useEffect(() => {
    if (!loaders[dataset]) {
      return;
    }

    if (loaders[dataset].loaders.neighborhoods) {
      loaders[dataset].loaders.neighborhoods
        .load()
        .catch((e: any) => warn(e, setWarning))
        .then((payload: any) => {
          if (!payload) return;
          const { data, url, coordinationValues } = payload;
          setNeighborhoods(data);
          addUrl(url, "Neighborhoods");
          const coordinationValuesOrDefault = {
            spatialNeighborhoodsLayer: DEFAULT_NEIGHBORHOODS_LAYER,
            ...coordinationValues,
          };
          initCoordinationSpace(coordinationValuesOrDefault, coordinationSetters, initialCoordinationValues);
          setItemIsReady("neighborhoods");
        });
    } else {
      setNeighborhoods({});
      if (isRequired) {
        warn(new LoaderNotFoundError(dataset, "neighborhoods", null), setWarning);
      } else {
        setItemIsReady("neighborhoods");
      }
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [loaders, dataset]);

  return [neighborhoods];
}

/**
 * Get data from a raster data type loader,
 * updating "ready" and URL state appropriately.
 * Throw warnings if the data is marked as required.
 * Subscribe to loader updates.
 * @param {object} loaders The object mapping
 * datasets and data types to loader instances.
 * @param {string} dataset The key for a dataset,
 * used to identify which loader to use.
 * @param {function} setItemIsReady A function to call
 * when done loading.
 * @param {function} addUrl A function to call to update
 * the URL list.
 * @param {boolean} isRequired Should a warning be thrown if
 * loading is unsuccessful?
 * @param {object} coordinationSetters Object where
 * keys are coordination type names with the prefix 'set',
 * values are coordination setter functions.
 * @param {object} initialCoordinationValues Object where
 * keys are coordination type names with the prefix 'initialize',
 * values are initialization preferences as boolean values.
 * @returns {array} [raster, imageLayerLoaders, imageLayerMeta] where
 * raster is an object,
 * imageLayerLoaders is an object, and
 * imageLayerMeta is an object.
 */
export function useRasterData(
  loaders: any,
  dataset: string,
  setItemIsReady: (readyItem: any) => void,
  addUrl: Function,
  isRequired: boolean,
  coordinationSetters: any,
  initialCoordinationValues: any
) {
  const [raster, setRaster] = useState();
  // Since we want the image layer / channel definitions to come from the
  // coordination space stored as JSON in the view config,
  // we need to set up a separate state variable here to store the
  // non-JSON objects, such as layer loader instances.
  const [imageLayerLoaders, setImageLayerLoaders] = useState<any[]>([]);
  const [imageLayerMeta, setImageLayerMeta] = useState<any[]>([]);

  const setWarning = useSetWarning();

  useEffect(() => {
    if (!loaders[dataset]) {
      return;
    }

    if (loaders[dataset].loaders.raster) {
      loaders[dataset].loaders.raster
        .load()
        .catch((e: any) => warn(e, setWarning))
        .then((payload: any) => {
          if (!payload) return;
          const { data, url: urls, coordinationValues } = payload;
          setRaster(data);
          urls.forEach(([url, name]: any) => {
            addUrl(url, name);
          });
          const { loaders: nextImageLoaders, meta: nextImageMeta } = data;
          setImageLayerLoaders(nextImageLoaders);
          setImageLayerMeta(nextImageMeta);
          initCoordinationSpace(coordinationValues, coordinationSetters, initialCoordinationValues);
          setItemIsReady("raster");
        });
    } else {
      // There was no raster loader for this dataset,
      // and raster should be optional.
      setImageLayerLoaders([]);
      setImageLayerMeta([]);
      if (isRequired) {
        warn(new LoaderNotFoundError(dataset, "raster", null), setWarning);
      } else {
        setItemIsReady("raster");
      }
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [loaders, dataset]);
  return [raster, imageLayerLoaders, imageLayerMeta];
}

/**
 * Get data from a genomic-profiles data type loader,
 * updating "ready" and URL state appropriately.
 * Throw warnings if the data is marked as required.
 * Subscribe to loader updates.
 * @param {object} loaders The object mapping
 * datasets and data types to loader instances.
 * @param {string} dataset The key for a dataset,
 * used to identify which loader to use.
 * @param {function} setItemIsReady A function to call
 * when done loading.
 * @param {function} addUrl A function to call to update
 * the URL list.
 * @param {boolean} isRequired Should a warning be thrown if
 * loading is unsuccessful?
 * @param {object} coordinationSetters Object where
 * keys are coordination type names with the prefix 'set',
 * values are coordination setter functions.
 * @param {object} initialCoordinationValues Object where
 * keys are coordination type names with the prefix 'initialize',
 * values are initialization preferences as boolean values.
 * @returns {array} [neighborhoods] where
 * neighborhoods is an object.
 */
export function useGenomicProfilesData(
  loaders: any,
  dataset: string,
  setItemIsReady: (readyItem: any) => void,
  addUrl: Function,
  isRequired: boolean,
  coordinationSetters: any,
  initialCoordinationValues: any
) {
  const [genomicProfilesAttrs, setGenomicProfilesAttrs] = useState<any>();

  const setWarning = useSetWarning();

  useEffect(() => {
    if (!loaders[dataset]) {
      return;
    }

    if (loaders[dataset].loaders["genomic-profiles"]) {
      loaders[dataset].loaders["genomic-profiles"]
        .load()
        .catch((e: any) => warn(e, setWarning))
        .then((payload: any) => {
          if (!payload) return;
          const { data, url, coordinationValues } = payload;
          setGenomicProfilesAttrs(data);
          addUrl(url);
          initCoordinationSpace(coordinationValues, coordinationSetters, initialCoordinationValues);
          setItemIsReady("genomic-profiles");
        });
    } else {
      setGenomicProfilesAttrs(null);
      if (isRequired) {
        warn(new LoaderNotFoundError(dataset, "genomic-profiles", null), setWarning);
      } else {
        setItemIsReady("genomic-profiles");
      }
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [loaders, dataset]);

  return [genomicProfilesAttrs];
}
