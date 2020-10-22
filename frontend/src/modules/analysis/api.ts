import {
  IPCAData,
  IPhenoGraphData,
  IPlotSeries,
  IRegionChannelData,
  IRegionStatsSubmission,
  IScatterPlotData,
  ITSNEData,
  IUMAPData,
} from "./models";
import { ApiManager } from "@/utils/api";

export const api = {
  async getScatterPlotData(
    datasetId: number,
    resultId: number | null,
    markerX: string,
    markerY: string,
    heatmapType: string,
    heatmap: string
  ) {
    let url = `analysis/scatterplot?dataset_id=${datasetId}`;
    if (resultId) {
      url += `&result_id=${resultId}`;
    }
    url += `&marker_x=${markerX}&marker_y=${markerY}`;
    if (heatmapType && heatmap) {
      url += `&heatmap_type=${heatmapType}&heatmap=${heatmap}`;
    }
    return ApiManager.api.get(url).json<IScatterPlotData>();
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
  async getPcaData(groupId: number, resultId: number) {
    return ApiManager.api
      .get(
        `results/${resultId}/pca`
      )
      .json<IPCAData>();
  },
  async getTsneData(groupId: number, resultId: number) {
    return ApiManager.api
      .get(
        `results/${resultId}/tsne`
      ).json<ITSNEData>();
  },
  async getUmapData(groupId: number, resultId: number) {
    return ApiManager.api
      .get(
        `results/${resultId}/umap`
      ).json<IUMAPData>();
  },
  async getPhenoGraphData(groupId: number, resultId: number) {
    return ApiManager.api.get(`groups/${groupId}/results/${resultId}/phenograph`).json<IPhenoGraphData>();
  },
  async calculateRegionStats(params: IRegionStatsSubmission) {
    return ApiManager.api
      .post(`analysis/region/stats`, {
        json: params,
      })
      .json<IRegionChannelData[]>();
  },
};
