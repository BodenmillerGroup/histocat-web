import axios from 'axios';
import { apiUrl } from '@/env';
import { IExperiment, IExperimentCreate, IExperimentUpdate } from '@/interfaces/experiment';
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
};
