import { ICellPoint } from "@/modules/results/models";

export function transformToWebGl(data: ICellPoint[], width: number, height: number) {
  return data.map((c, i) => {
    return [-1 + 2 * (c.x / width), -1 + 2 * (1 - c.y / height), c.acquisitionId, c.cellId]; // [x, y, category, value]
  });
}

export function transformFromWebGl(data: [number, number][], width: number, height: number) {
  // TODO: don't forget Y axis flip!
  return data.map((c, i) => {
    return [((c[0] + 1) * width) / 2, Math.abs(((c[1] + 1) * height) / 2 - height)];
  });
}
