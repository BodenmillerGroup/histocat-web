import {
  IPhenoGraphData,
  IPlotSeries,
  IRawColorsData,
  IRawResultData,
  IRawScatterData,
  IResult,
  IResultUpdate,
} from "./models";
import { ApiManager } from "@/utils/api";

export const api = {
  async getDatasetResults(groupId: number, datasetId: number) {
    return ApiManager.api.get(`groups/${groupId}/datasets/${datasetId}/results`).json<IResult[]>();
  },
  async getResult(groupId: number, resultId: number) {
    return ApiManager.api.get(`groups/${groupId}/results/${resultId}`).json<IResult>();
  },
  async updateResult(groupId: number, resultId: number, data: IResultUpdate) {
    return ApiManager.api
      .patch(`groups/${groupId}/results/${resultId}`, {
        json: data,
      })
      .json<IResult>();
  },
  async deleteResult(groupId: number, resultId: number) {
    return ApiManager.api.delete(`groups/${groupId}/results/${resultId}`).json();
  },
  async getResultData(groupId: number, resultId: number) {
    let url = `groups/${groupId}/results/${resultId}/data`;
    return ApiManager.api.get(url).json<IRawResultData>();
  },
  async getColorsData(groupId: number, resultId: number, colorsType: string, colorsName: string) {
    let url = `groups/${groupId}/results/${resultId}/colors?colors_type=${colorsType}&colors_name=${colorsName}`;
    return ApiManager.api.get(url).json<IRawColorsData>();
  },
  async getScatterPlotData(groupId: number, resultId: number, markerX: string, markerY: string) {
    let url = `groups/${groupId}/results/${resultId}/scatterplot?marker_x=${markerX}&marker_y=${markerY}`;
    return ApiManager.api.get(url).json<IRawScatterData>();
  },
  async getBoxPlotData(datasetId: number, gateId: number | null, acquisitionIds: number[], markers: string[]) {
    const acquisitionIdsArray = acquisitionIds.map((acquisition_id) => `&acquisition_ids=${acquisition_id}`);
    const markersArray = markers.map((marker) => `&markers=${marker}`);
    const gateIdArg = gateId !== null ? `&gate_id=${gateId}` : "";
    return ApiManager.api
      .get(
        `analysis/boxplot?dataset_id=${datasetId}${gateIdArg}${acquisitionIdsArray.join("")}${markersArray.join("")}`
      )
      .json<IPlotSeries[]>();
  },
  async getPhenoGraphData(groupId: number, resultId: number) {
    return ApiManager.api.get(`groups/${groupId}/results/${resultId}/phenograph`).json<IPhenoGraphData>();
  },
};
