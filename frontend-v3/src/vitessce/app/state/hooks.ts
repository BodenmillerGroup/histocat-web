import { useRef, useCallback } from "react";
import create from "zustand";
import shallow from "zustand/shallow";
import { fromEntries, capitalize } from "../../utils";
import { Layout } from "../../types";

type ViewConfigState = {
  viewConfig: any;
  loaders: any;

  setViewConfig(viewConfig: any): void;
  setLoaders(loaders: any): void;
  setCoordinationValue({ parameter, scope, value }: any): void;
  removeComponent(i: number): void;
  changeLayout(newComponentProps: any): void;
}

/**
 * The useViewConfigStore hook is initialized via the zustand
 * create() function, which sets up both the state variables
 * and the reducer-type functions.
 * Reference: https://github.com/react-spring/zustand
 * @returns {function} The useStore hook.
 */
export const useViewConfigStore = create<ViewConfigState>((set) => ({
  // State:
  // The viewConfig is an object which must conform to the schema
  // found in src/schemas/config.schema.json.
  viewConfig: null,
  // The loaders object is a mapping from dataset ID to
  // data type to loader object instance.
  loaders: null,
  // Reducer functions which update the state
  // (although technically also part of state):
  setViewConfig: (viewConfig: any) => set({ viewConfig }),
  setLoaders: (loaders: any) => set({ loaders }),
  setCoordinationValue: ({ parameter, scope, value }: any) =>
    set((state) => ({
      viewConfig: {
        ...state.viewConfig,
        coordinationSpace: {
          ...state.viewConfig.coordinationSpace,
          [parameter]: {
            ...state.viewConfig.coordinationSpace[parameter],
            [scope]: value,
          },
        },
      },
    })),
  removeComponent: (i: number) =>
    set((state) => {
      const newLayout = state.viewConfig.layout.slice();
      newLayout.splice(i, 1);
      return {
        viewConfig: {
          ...state.viewConfig,
          layout: newLayout,
        },
      };
    }),
  changeLayout: (newComponentProps: any) =>
    set((state) => {
      const newLayout = state.viewConfig.layout.slice();
      newComponentProps.forEach(([i, newProps]: any) => {
        newLayout[i] = {
          ...newLayout[i],
          ...newProps,
        };
      });
      return {
        viewConfig: {
          ...state.viewConfig,
          layout: newLayout,
        },
      };
    }),
}));

type HoverState = {
  componentHover: any;
  setComponentHover(componentHover: string): void;
}

/**
 * The hover store can be used to store global state
 * related to which component is currently hovered,
 * which is required for tooltip / crossover elements.
 * @returns {function} The useStore hook.
 */
const useHoverStore = create<HoverState>((set) => ({
  // Components may need to know if they are the "hover source"
  // for tooltip interactions. This value should be a unique
  // component ID, such as its index in the view config layout.
  componentHover: null,
  setComponentHover: (componentHover: string) => set({ componentHover }),
}));

type WarnState = {
  warning: any;
  setWarning(warning: any): void;
}

/**
 * The warning store can be used to store global state
 * related to app warning messages.
 * @returns {function} The useStore hook.
 */
const useWarnStore = create<WarnState>((set) => ({
  // Want a global state to collect warning messages
  // that occur anywhere in the app.
  warning: null,
  setWarning: (warning: any) => set({ warning }),
}));

type ViewInfoState = {
  viewInfo: any;
  setComponentViewInfo(uuid: string, viewInfo: any): void;
}

/**
 * The view info store can be used to store component-level
 * viewInfo objects,
 * which are required for tooltip / crossover elements.
 * @returns {function} The useStore hook.
 */
const useViewInfoStore = create<ViewInfoState>((set) => ({
  // The viewInfo object is a mapping from
  // component IDs to component view info objects.
  // Each view info object must have a project() function.
  viewInfo: {},
  setComponentViewInfo: (uuid, viewInfo) =>
    set((state) => ({
      viewInfo: {
        ...state.viewInfo,
        [uuid]: viewInfo,
      },
    })),
}));

type GridSizeState = {
  resizeCount: any;
  incrementResizeCount(): void;
}

/**
 * The grid size store can be used to store a
 * counter which updates on each window or react-grid-layout
 * resize event.
 * @returns {function} The useStore hook.
 */
const useGridSizeStore = create<GridSizeState>((set) => ({
  resizeCount: {},
  incrementResizeCount: () =>
    set((state) => ({
      resizeCount: state.resizeCount + 1,
    })),
}));

/**
 * The useCoordination hook returns both the
 * values and setter functions for the coordination objects
 * in a particular coordination scope mapping.
 * This hook is intended to be used within the ___Subscriber
 * components to allow them to "hook into" only those coordination
 * objects and setter functions of relevance.
 * @param {string[]} parameters Array of coordination types.
 * @param {object} coordinationScopes Mapping of coordination types
 * to scope names.
 * @returns {array} Returns a tuple [values, setters]
 * where values is an object containing all coordination values,
 * and setters is an object containing all coordination setter
 * functions for the values in `values`, named with a "set"
 * prefix.
 */
export function useCoordination(parameters: string[], coordinationScopes: any) {
  const setCoordinationValue = useViewConfigStore((state) => state.setCoordinationValue);

  const values = useViewConfigStore((state) => {
    const { coordinationSpace } = state.viewConfig;
    return fromEntries(
      parameters.map((parameter) => {
        if (coordinationSpace && coordinationSpace[parameter]) {
          const value = coordinationSpace[parameter][coordinationScopes[parameter]];
          return [parameter, value];
        }
        return [parameter, undefined];
      })
    );
  }, shallow);

  const setters = fromEntries(
    parameters.map((parameter) => {
      const setterName = `set${capitalize(parameter)}`;
      const setterFunc = (value: any) =>
        setCoordinationValue({
          parameter,
          scope: coordinationScopes[parameter],
          value,
        });
      return [setterName, setterFunc];
    })
  );

  return [values, setters];
}

/**
 * Obtain the loaders object from
 * the global app state.
 * @returns {object} The loaders object
 * in the `useViewConfigStore` store.
 */
export function useLoaders() {
  return useViewConfigStore((state) => {
    return state.loaders
  });
}

/**
 * Obtain the view config layout array from
 * the global app state.
 * @returns {object[]} The layout array
 * in the `useViewConfigStore` store.
 */
export function useLayout(): Layout {
  return useViewConfigStore((state) => state.viewConfig?.layout);
}

/**
 * Obtain the component removal function from
 * the global app state.
 * @returns {function} The remove component function
 * in the `useViewInfoStore` store.
 */
export function useRemoveComponent() {
  return useViewConfigStore((state) => state.removeComponent);
}

/**
 * Obtain the component prop setter function from
 * the global app state.
 * @returns {function} The set component props function
 * in the `useViewInfoStore` store.
 */
export function useChangeLayout() {
  return useViewConfigStore((state) => state.changeLayout);
}

/**
 * Obtain the loaders setter function from
 * the global app state.
 * @returns {function} The loaders setter function
 * in the `useViewConfigStore` store.
 */
export function useSetLoaders() {
  return useViewConfigStore((state) => state.setLoaders);
}

/**
 * Obtain the view config setter function from
 * the global app state.
 * @returns {function} The view config setter function
 * in the `useViewConfigStore` store.
 */
export function useSetViewConfig() {
  const setViewConfigRef = useRef(useViewConfigStore.getState().setViewConfig);
  const setViewConfig = setViewConfigRef.current;
  return setViewConfig;
}

/**
 * Obtain the component hover value from
 * the global app state.
 * @returns {number} The hovered component ID
 * in the `useHoverStore` store.
 */
export function useComponentHover() {
  return useHoverStore((state) => state.componentHover);
}

/**
 * Obtain the component hover setter function from
 * the global app state.
 * @returns {function} The component hover setter function
 * in the `useHoverStore` store.
 */
export function useSetComponentHover() {
  return useHoverStore((state) => state.setComponentHover);
}

/**
 * Obtain the warning message from
 * the global app state.
 * @returns {string} The warning message
 * in the `useWarnStore` store.
 */
export function useWarning() {
  return useWarnStore((state) => state.warning);
}

/**
 * Obtain the warning setter function from
 * the global app state.
 * @returns {function} The warning setter function
 * in the `useWarnStore` store.
 */
export function useSetWarning() {
  return useWarnStore((state) => state.setWarning);
}

/**
 * Obtain the component view info value from
 * the global app state.
 * @returns {object} The view info object for the component
 * in the `useViewInfoStore` store.
 */
export function useComponentViewInfo(uuid: string) {
  return useViewInfoStore(useCallback((state) => state.viewInfo[uuid], [uuid]));
}

/**
 * Obtain the component view info setter function from
 * the global app state.
 * @returns {function} The component view info setter function
 * in the `useViewInfoStore` store.
 */
export function useSetComponentViewInfo(uuid: string) {
  const setViewInfoRef = useRef(useViewInfoStore.getState().setComponentViewInfo);
  const setComponentViewInfo = (viewInfo: any) => setViewInfoRef.current(uuid, viewInfo);
  return setComponentViewInfo;
}

/**
 * Obtain the grid resize count value
 * from the global app state.
 * @returns {number} The grid resize increment value.
 */
export function useGridResize() {
  return useGridSizeStore((state) => state.resizeCount);
}

/**
 * Obtain the grid resize count increment function
 * from the global app state.
 * @returns {function} The grid resize count increment
 * function.
 */
export function useEmitGridResize() {
  return useGridSizeStore((state) => state.incrementResizeCount);
}
