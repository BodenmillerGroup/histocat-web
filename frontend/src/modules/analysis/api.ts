import { apiUrl } from "@/env";
import ky from "ky";
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
  IUMAPSubmission
} from "./models";

export const api = {
  async produceSegmentationImage(token: string, params: IImageSegmentationSubmission) {
    return ky.post(`${apiUrl}/api/v1/analysis/segmentation/image`, {
      headers: {
        Authorization: `Bearer ${token}`
      },
      json: params
    });
  },
  async produceSegmentationContours(token: string, params: IImageSegmentationSubmission) {
    return ky
      .post(`${apiUrl}/api/v1/analysis/segmentation/contours`, {
        headers: {
          Authorization: `Bearer ${token}`
        },
        json: params
      })
      .json<number[][]>();
  },
  async getScatterPlotData(
    token: string,
    datasetId: number,
    acquisitionIds: number[],
    markerX: string,
    markerY: string,
    markerZ: string,
    heatmapType: string,
    heatmap: string
  ) {
    const acquisitionIdsArray = acquisitionIds.map(acquisition_id => `&acquisition_ids=${acquisition_id}`);
    const acquisition_ids = acquisitionIdsArray.join("");
    return ky
      .get(
        `${apiUrl}/api/v1/analysis/scatterplot?dataset_id=${datasetId}&marker_x=${markerX}&marker_y=${markerY}&marker_z=${markerZ}&heatmap_type=${heatmapType}&heatmap=${heatmap}${acquisition_ids}`,
        {
          headers: {
            Authorization: `Bearer ${token}`
          }
        }
      )
      .json<IScatterPlotData>();
  },
  async getBoxPlotData(token: string, datasetId: number, acquisitionId: number, markers: string[]) {
    const markersArray = markers.map(marker => `&markers=${marker}`);
    return ky
      .get(
        `${apiUrl}/api/v1/analysis/boxplot?dataset_id=${datasetId}&acquisition_id=${acquisitionId}${markersArray.join(
          ""
        )}`,
        {
          headers: {
            Authorization: `Bearer ${token}`
          }
        }
      )
      .json<IPlotSeries[]>();
  },
  async getPCAData(token: string, params: IPCASubmission) {
    const acquisitionIdsArray = params.acquisition_ids.map(acquisition_id => `&acquisition_ids=${acquisition_id}`);
    const acquisitionIds = acquisitionIdsArray.join("");
    const markersArray = params.markers.map(marker => `&markers=${marker}`);
    const markers = markersArray.join("");
    return ky
      .get(
        `${apiUrl}/api/v1/analysis/pca?dataset_id=${params.dataset_id}&n_components=${params.n_components}&heatmap_type=${params.heatmapType}&heatmap=${params.heatmap}${acquisitionIds}${markers}`,
        {
          headers: {
            Authorization: `Bearer ${token}`
          }
        }
      )
      .json<IPCAData>();
  },
  async submitTSNE(token: string, data: ITSNESubmission) {
    return ky
      .post(`${apiUrl}/api/v1/analysis/tsne`, {
        json: data,
        headers: {
          Authorization: `Bearer ${token}`
        }
      })
      .json();
  },
  async getTSNEData(token: string, datasetId: number, name: string, heatmapType: string, heatmap: string) {
    return ky
      .get(
        `${apiUrl}/api/v1/analysis/tsne?dataset_id=${datasetId}&name=${name}&heatmap_type=${heatmapType}&heatmap=${heatmap}`,
        {
          headers: {
            Authorization: `Bearer ${token}`
          }
        }
      )
      .json<ITSNEData>();
  },
  async submitUMAP(token: string, data: IUMAPSubmission) {
    return ky
      .post(`${apiUrl}/api/v1/analysis/umap`, {
        json: data,
        headers: {
          Authorization: `Bearer ${token}`
        }
      })
      .json();
  },
  async getUMAPData(token: string, datasetId: number, name: string, heatmapType: string, heatmap: string) {
    return ky
      .get(
        `${apiUrl}/api/v1/analysis/umap?dataset_id=${datasetId}&name=${name}&heatmap_type=${heatmapType}&heatmap=${heatmap}`,
        {
          headers: {
            Authorization: `Bearer ${token}`
          }
        }
      )
      .json<IUMAPData>();
  },
  async submitPhenoGraph(token: string, data: IPhenoGraphSubmission) {
    return ky
      .post(`${apiUrl}/api/v1/analysis/phenograph`, {
        json: data,
        headers: {
          Authorization: `Bearer ${token}`
        }
      })
      .json();
  },
  async getPhenoGraphData(token: string, datasetId: number, name: string) {
    return ky
      .get(`${apiUrl}/api/v1/analysis/phenograph?dataset_id=${datasetId}&name=${name}`, {
        headers: {
          Authorization: `Bearer ${token}`
        }
      })
      .json<IPhenoGraphData>();
  },
  async calculateRegionStats(token: string, params: IRegionStatsSubmission) {
    return ky
      .post(`${apiUrl}/api/v1/analysis/region/stats`, {
        json: params,
        headers: {
          Authorization: `Bearer ${token}`
        }
      })
      .json<IRegionChannelData[]>();
  }
};
