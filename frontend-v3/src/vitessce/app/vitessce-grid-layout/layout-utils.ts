import { range } from "../../utils";
import { Layout } from "../../types";

export const COMPONENT_ID_PREFIX = 'i';

function sum(a: number[]) {
  return a.reduce((x, y) => x + y, 0);
}

export function makeGridLayout(colXs: number[], colLayout: object) {
  const colWs: number[] = [];
  for (let i = 0; i < colXs.length; i++) { // eslint-disable-line no-plusplus
    colWs.push(colXs[i + 1] - colXs[i]);
  }
  return Object.entries(colLayout).map(([id, spec]) => ({
    i: id,
    y: spec.y,
    h: spec.h || 1,
    x: colXs[spec.x],
    w: sum(colWs.slice(spec.x, spec.x + (spec.w || 1))),
  }));
}

export function getMaxRows(layouts: object) {
  return Math.max(
    ...Object.values(layouts).map(
      layout => Math.max(
        ...layout.map((xywh: any) => xywh.y + xywh.h),
      ),
    ),
  );
}

export function resolveLayout(layout: Layout) {
  const cols: any = {};
  const layouts: any = {};
  const breakpoints: any = {};
  const components: any = {};
  const positions: any = {};

  (('components' in layout) ? (layout as any).components : layout).forEach(
    (def: any, i: number) => {
      const id = `${COMPONENT_ID_PREFIX}${i}`;
      components[id] = {
        component: def.component,
        props: def.props || {},
        coordinationScopes: def.coordinationScopes || {},
      };
      positions[id] = {
        id, x: def.x, y: def.y, w: def.w, h: def.h,
      };
    },
  );

  if ('components' in layout) {
    Object.entries<number[]>((layout as any).columns).forEach(
      ([width, columnXs]) => {
        cols[width] = columnXs[columnXs.length - 1];
        layouts[width] = makeGridLayout(columnXs, positions);
        breakpoints[width] = width;
      },
    );
  } else {
    // static layout
    const id = 'ID';
    const columnCount = 12;
    cols[id] = columnCount;
    layouts[id] = makeGridLayout(range(columnCount + 1), positions);
    breakpoints[id] = 1000;
    // Default has different numbers of columns at different widths,
    // so we do need to override that to ensure the same number of columns,
    // regardless of window width.
  }
  return {
    cols, layouts, breakpoints, components,
  };
}
