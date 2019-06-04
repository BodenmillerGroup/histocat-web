import axios from 'axios';
import { apiUrl } from '@/env';
import { IChannelStats, IExperiment, IExperimentCreate, IExperimentDataset, IExperimentUpdate } from './models';
import { authHeaders } from '@/utils';


export const api = {
  async getExperiments(token: string) {
    return axios.get<IExperiment[]>(`${apiUrl}/api/v1/experiments/`, authHeaders(token));
  },
  async updateExperiment(token: string, id: number, data: IExperimentUpdate) {
    return axios.put(`${apiUrl}/api/v1/experiments/${id}`, data, authHeaders(token));
  },
  async createExperiment(token: string, data: IExperimentCreate) {
    return axios.post(`${apiUrl}/api/v1/experiments/`, data, authHeaders(token));
  },
  async uploadSlide(token: string, id: number, data) {
    return axios.post(`${apiUrl}/api/v1/experiments/${id}/upload_slide`, data, authHeaders(token, 'multipart/form-data'));
  },
  async deleteExperiment(token: string, id: number) {
    return axios.delete(`${apiUrl}/api/v1/experiments/${id}`, authHeaders(token));
  },
  async getExperiment(token: string, id: number) {
    return axios.get<IExperiment>(`${apiUrl}/api/v1/experiments/${id}`, authHeaders(token));
  },
  async getExperimentDataset(token: string, id: number) {
    return axios.get<IExperimentDataset>(`${apiUrl}/api/v1/experiments/${id}/dataset`, authHeaders(token));
  },
  async getChannelImage(token: string, id: number) {
    return axios.get<any>(`${apiUrl}/api/v1/channels/${id}/image`, authHeaders(token));
  },
  async getChannelStats(token: string, id: number) {
    return axios.get<IChannelStats>(`${apiUrl}/api/v1/channels/${id}/stats`, authHeaders(token));
  },
};
