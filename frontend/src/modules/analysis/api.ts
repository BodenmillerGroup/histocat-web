import {
  IImageSegmentationSubmission,
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
  async produceSegmentationImage(params: IImageSegmentationSubmission) {
    return ApiManager.api.post(`analysis/segmentation/image`, {
      json: params,
    });
  },
  async produceSegmentationContours(params: IImageSegmentationSubmission) {
    return ApiManager.api
      .post(`analysis/segmentation/contours`, {
        json: params,
      })
      .json<number[][]>();
  },
  async getScatterPlotData(
    datasetId: number,
    acquisitionIds: number[],
    markerX: string,
    markerY: string,
    markerZ: string,
    heatmapType: string,
    heatmap: string
  ) {
    const acquisitionIdsArray = acquisitionIds.map((acquisition_id) => `&acquisition_ids=${acquisition_id}`);
    const acquisition_ids = acquisitionIdsArray.join("");
    let url = `analysis/scatterplot?dataset_id=${datasetId}&marker_x=${markerX}&marker_y=${markerY}${acquisition_ids}`;
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
  async getTSNEData(datasetId: number, name: string, heatmapType: string, heatmap: string) {
    let url = `analysis/tsne?dataset_id=${datasetId}&name=${name}`;
    if (heatmapType && heatmap) {
      url += `&heatmap_type=${heatmapType}&heatmap=${heatmap}`;
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
  async getUMAPData(datasetId: number, name: string, heatmapType: string, heatmap: string) {
    let url = `analysis/umap?dataset_id=${datasetId}&name=${name}`;
    if (heatmapType && heatmap) {
      url += `&heatmap_type=${heatmapType}&heatmap=${heatmap}`;
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
  async getPhenoGraphData(datasetId: number, name: string) {
    return ApiManager.api.get(`analysis/phenograph?dataset_id=${datasetId}&name=${name}`).json<IPhenoGraphData>();
  },
  async calculateRegionStats(params: IRegionStatsSubmission) {
    return ApiManager.api
      .post(`analysis/region/stats`, {
        json: params,
      })
      .json<IRegionChannelData[]>();
  },
};
