import ky from 'ky';
import { apiUrl } from '@/env';
import { IChannelStats, IExperiment, IExperimentCreate, IExperimentDataset, IExperimentUpdate } from './models';


export const api = {
  async getExperiments(token: string) {
    return ky.get(`${apiUrl}/api/v1/experiments/`, {
      headers: {
        Authorization: `Bearer ${token}`,
      },
    }).json<IExperiment[]>();
  },
  async getTags(token: string) {
    return ky.get(`${apiUrl}/api/v1/experiments/tags`, {
      headers: {
        Authorization: `Bearer ${token}`,
      },
    }).json<string[]>();
  },
  async updateExperiment(token: string, id: number, data: IExperimentUpdate) {
    return ky.put(`${apiUrl}/api/v1/experiments/${id}`, {
      json: data,
      headers: {
        Authorization: `Bearer ${token}`,
      },
    }).json();
  },
  async createExperiment(token: string, data: IExperimentCreate) {
    return ky.post(`${apiUrl}/api/v1/experiments/`, {
      json: data,
      headers: {
        Authorization: `Bearer ${token}`,
      },
    }).json();
  },
  async uploadSlide(token: string, id: number, data) {
    return ky.post(`${apiUrl}/api/v1/experiments/${id}/upload_slide`, {
      json: data,
      headers: {
        Authorization: `Bearer ${token}`,
        'Content-Type': 'multipart/form-data',
      },
    }).json();
  },
  async deleteExperiment(token: string, id: number) {
    return ky.delete(`${apiUrl}/api/v1/experiments/${id}`, {
      headers: {
        Authorization: `Bearer ${token}`,
      },
    }).json();
  },
  async getExperiment(token: string, id: number) {
    return ky.get(`${apiUrl}/api/v1/experiments/${id}`, {
      headers: {
        Authorization: `Bearer ${token}`,
      },
    }).json<IExperiment>();
  },
  async getExperimentDataset(token: string, id: number) {
    return ky.get(`${apiUrl}/api/v1/experiments/${id}/dataset`, {
      headers: {
        Authorization: `Bearer ${token}`,
      },
    }).json<IExperimentDataset>();
  },
  async getChannelImage(token: string, id: number) {
    return ky.get(`${apiUrl}/api/v1/channels/${id}/image`, {
      headers: {
        Authorization: `Bearer ${token}`,
      },
    });
  },
  async getChannelStats(token: string, id: number) {
    const url = `${apiUrl}/api/v1/channels/${id}/stats`;
    const cache = await caches.open('stats');
    const found = await cache.match(url);
    if (found) {
      return found.json() as Promise<IChannelStats>;
    } else {
      const response = await ky.get(url, {
        headers: {
          Authorization: `Bearer ${token}`,
        },
      });
      if (response.ok) {
        await cache.put(url, response.clone());
      }
      return response.json() as Promise<IChannelStats>;
    }
  },
};
