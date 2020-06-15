import { CellPoint } from "@/data/CellPoint";

export function transformCoords(data: CellPoint[], width: number, height: number) {
  return data.map((c, i) => {
    return [-1 + 2 * (c.x / width), -1 + 2 * (1 - c.y / height), c.acquisitionId, c.cellId]; // [x, y, category, value]
  });
}
