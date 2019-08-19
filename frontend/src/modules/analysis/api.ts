import { apiUrl } from '@/env';
import { IImageSegmentationSubmission } from '@/modules/analysis/models';
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
};
