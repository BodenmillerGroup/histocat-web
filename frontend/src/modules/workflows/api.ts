import { apiUrl } from '@/env';
import ky from 'ky';
import { IWorkflow, IWorkflowCreate, IWorkflowUpdate } from './models';


export const api = {
  async createWorkflow(token: string, data: IWorkflowCreate) {
    return ky.post(`${apiUrl}/api/v1/workflows/`, {
      json: data,
      headers: {
        Authorization: `Bearer ${token}`,
      },
    }).json<IWorkflow>();
  },
  async getWorkflows(token: string) {
    return ky.get(`${apiUrl}/api/v1/workflows/`, {
      headers: {
        Authorization: `Bearer ${token}`,
      },
    }).json<IWorkflow[]>();
  },
  async updateWorkflow(token: string, id: number, data: IWorkflowUpdate) {
    return ky.put(`${apiUrl}/api/v1/workflows/${id}`, {
      json: data,
      headers: {
        Authorization: `Bearer ${token}`,
      },
    }).json<IWorkflow>();
  },
  async deleteWorkflow(token: string, id: number) {
    return ky.delete(`${apiUrl}/api/v1/workflows/${id}`, {
      headers: {
        Authorization: `Bearer ${token}`,
      },
    }).json();
  },
};
