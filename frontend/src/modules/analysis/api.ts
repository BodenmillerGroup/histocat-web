import { apiUrl } from '@/env';
import ky from 'ky';
import { IImageSegmentationSubmission, IPlotSeries, IScatterPlotData } from './models';

export const api = {
  async produceSegmentationImage(token: string, params: IImageSegmentationSubmission) {
    return ky.post(`${apiUrl}/api/v1/analysis/segmentation/image`, {
      headers: {
        Authorization: `Bearer ${token}`,
      },
      json: params,
    });
  },
  async produceSegmentationContours(token: string, params: IImageSegmentationSubmission) {
    return ky.post(`${apiUrl}/api/v1/analysis/segmentation/contours`, {
      headers: {
        Authorization: `Bearer ${token}`,
      },
      json: params,
    }).json<number[][]>();
  },
  async getScatterPlotData(
    token: string,
    datasetId: number,
    acquisitionId: number,
    markerX: string,
    markerY: string,
    markerZ: string
  ) {
    return ky.get(`${apiUrl}/api/v1/analysis/scatterplot?dataset_id=${datasetId}&acquisition_id=${acquisitionId}&marker_x=${markerX}&marker_y=${markerY}&marker_z=${markerZ}`, {
      headers: {
        Authorization: `Bearer ${token}`,
      },
    }).json<IScatterPlotData>();
  },
  async getBoxPlotData(
    token: string,
    datasetId: number,
    acquisitionId: number,
    markers: string[],
  ) {
    const markersArray = markers.map(marker => `&markers=${marker}`);
    return ky.get(`${apiUrl}/api/v1/analysis/boxplot?dataset_id=${datasetId}&acquisition_id=${acquisitionId}&markers=${markersArray.join('')}`, {
      headers: {
        Authorization: `Bearer ${token}`,
      },
    }).json<IPlotSeries[]>();
  },
};
