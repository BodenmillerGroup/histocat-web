import { apiUrl } from '@/env';
import { IImageSegmentationSubmission, IScatterPlotData } from '@/modules/analysis/models';
import ky from 'ky';

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
  async getScatterPlotData(token: string, datasetId: number) {
    return ky.get(`${apiUrl}/api/v1/analysis/${datasetId}/scatter`, {
      headers: {
        Authorization: `Bearer ${token}`,
      },
    }).json<IScatterPlotData>();
  },
};
