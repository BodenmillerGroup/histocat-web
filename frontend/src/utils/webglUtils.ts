import { ICell } from "@/modules/cells/models";

export function transformToWebGl(data: ICell[], width: number, height: number) {
  return data.map((c, i) => {
    return [-1 + 2 * (c.xy[0] / width), -1 + 2 * (1 - c.xy[1] / height), c.color, c.cellId]; // [x, y, category, value]
  });
}

export function transformFromWebGl(data: [number, number][], width: number, height: number) {
  // TODO: don't forget Y axis flip!
  return data.map((c, i) => {
    return [((c[0] + 1) * width) / 2, Math.abs(((c[1] + 1) * height) / 2 - height)];
  });
}
