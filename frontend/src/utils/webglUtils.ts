import { ICentroidsData } from "@/modules/analysis/models";

export function transformCoords(data: ICentroidsData, width: number, height: number) {
  return data.coords.map((c, i) => {
    const [x, y] = c;
    const [acquisition_id, cell_id] = data.cell_ids[i].split("_");
    return [-1 + 2 * (x / width), -1 + 2 * (y / height), Number(acquisition_id), Number(cell_id)]; // [x, y, category, value]
  });
}
