import {
  IPCAData,
  IPCASubmission,
  IPhenoGraphData,
  IPhenoGraphSubmission,
  IPlotSeries,
  IRegionChannelData,
  IRegionStatsSubmission,
  IScatterPlotData,
  ITSNEData,
  ITSNESubmission,
  IUMAPData,
  IUMAPSubmission,
} from "./models";
import { ApiManager } from "@/utils/api";

export const api = {
  async getScatterPlotData(
    datasetId: number,
    resultId: number | null,
    markerX: string,
    markerY: string,
    markerZ: string,
    heatmapType: string,
    heatmap: string
  ) {
    let url = `analysis/scatterplot?dataset_id=${datasetId}`;
    if (resultId) {
      url += `&result_id=${resultId}`;
    }
    url += `&marker_x=${markerX}&marker_y=${markerY}`;
    if (markerZ) {
      url += `&marker_z=${markerZ}`;
    }
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
  async getPCAData(params: IPCASubmission) {
    const acquisitionIdsArray = params.acquisition_ids.map((acquisition_id) => `&acquisition_ids=${acquisition_id}`);
    const acquisitionIds = acquisitionIdsArray.join("");
    const markersArray = params.markers.map((marker) => `&markers=${marker}`);
    const markers = markersArray.join("");
    return ApiManager.api
      .get(
        `analysis/pca?dataset_id=${params.dataset_id}&n_components=${params.n_components}&heatmap_type=${params.heatmapType}&heatmap=${params.heatmap}${acquisitionIds}${markers}`
      )
      .json<IPCAData>();
  },
  async submitTSNE(data: ITSNESubmission) {
    return ApiManager.api
      .post(`analysis/tsne`, {
        json: data,
      })
      .json();
  },
  async getTSNEData(groupId: number, resultId: number, heatmapType: string, heatmap: string) {
    let url = `groups/${groupId}/results/${resultId}/tsne`;
    if (heatmapType && heatmap) {
      url += `?heatmap_type=${heatmapType}&heatmap=${heatmap}`;
    }
    return ApiManager.api.get(url).json<ITSNEData>();
  },
  async submitUMAP(data: IUMAPSubmission) {
    return ApiManager.api
      .post(`analysis/umap`, {
        json: data,
      })
      .json();
  },
  async getUMAPData(groupId: number, resultId: number, heatmapType: string, heatmap: string) {
    let url = `groups/${groupId}/results/${resultId}/umap`;
    if (heatmapType && heatmap) {
      url += `?heatmap_type=${heatmapType}&heatmap=${heatmap}`;
    }
    return ApiManager.api.get(url).json<IUMAPData>();
  },
  async submitPhenoGraph(data: IPhenoGraphSubmission) {
    return ApiManager.api
      .post(`analysis/phenograph`, {
        json: data,
      })
      .json();
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
