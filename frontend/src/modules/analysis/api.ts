import { apiUrl } from '@/env';
import { IAnalysisSubmission } from '@/modules/analysis/models';
import ky from 'ky';

export const api = {
  async downloadAnalysisImage(token: string, params: IAnalysisSubmission) {
    return ky.post(`${apiUrl}/api/v1/analysis/segmentation`, {
      headers: {
        Authorization: `Bearer ${token}`,
      },
      json: params,
    });
  },
};
