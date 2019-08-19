import { apiUrl } from '@/env';
import ky from 'ky';
import {
  IChannelStack,
  IChannelStats,
  IExperiment,
  IExperimentCreate,
  IExperimentUpdate,
  IShare,
  IShareCreate,
} from './models';

const cacheAvailable = false; // 'caches' in self;


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
      body: data,
      headers: {
        Authorization: `Bearer ${token}`,
      },
      timeout: false,
    });
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
  async getExperimentData(token: string, id: number) {
    return ky.get(`${apiUrl}/api/v1/experiments/${id}/data`, {
      headers: {
        Authorization: `Bearer ${token}`,
      },
    }).json<IExperiment>();
  },
  async getChannelStats(token: string, id: number) {
    const url = `${apiUrl}/api/v1/channels/${id}/stats?bins=100`;
    let cache;
    if (cacheAvailable) {
      cache = await self.caches.open('stats');
      const found = await cache.match(url);
      if (found) {
        return found.json() as Promise<IChannelStats>;
      }
    }
    const response = await ky.get(url, {
      headers: {
        Authorization: `Bearer ${token}`,
      },
    });
    if (response.ok) {
      if (cacheAvailable) {
        await cache.put(url, response.clone());
      }
    }
    return response.json() as Promise<IChannelStats>;
  },
  async downloadChannelStackImage(token: string, params: IChannelStack) {
    return ky.post(`${apiUrl}/api/v1/channels/stack`, {
      headers: {
        Authorization: `Bearer ${token}`,
      },
      json: params,
    });
  },
  async createShare(token: string, data: IShareCreate) {
    return ky.post(`${apiUrl}/api/v1/share/`, {
      json: data,
      headers: {
        Authorization: `Bearer ${token}`,
      },
    }).json();
  },
  async getExperimentShares(token: string, id: number) {
    return ky.get(`${apiUrl}/api/v1/share/${id}`, {
      headers: {
        Authorization: `Bearer ${token}`,
      },
    }).json<IShare[]>();
  },
};
