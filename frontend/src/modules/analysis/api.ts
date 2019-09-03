import { apiUrl } from '@/env';
import ky from 'ky';
import { IImageSegmentationSubmission, IScatterPlotData } from './models';

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
    markerX: number,
    markerY: number,
    markerZ?: number
  ) {
    return ky.get(`${apiUrl}/api/v1/analysis/scatter?dataset_id=${datasetId}&acquisition_id=${acquisitionId}&marker_x=${markerX}&marker_y=${markerY}&marker_z=${markerZ}`, {
      headers: {
        Authorization: `Bearer ${token}`,
      },
    }).json<IScatterPlotData>();
  },
};
