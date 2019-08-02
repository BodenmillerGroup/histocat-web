import { mainModule } from '@/modules/main';
import { settingsModule } from '@/modules/settings';
import { Store } from 'vuex';
import { Actions, Context } from 'vuex-smart-module';
import { WorkflowState } from '.';
import { api } from './api';
import { WorkflowGetters } from './getters';
import { IWorkflowCreate, IWorkflowUpdate } from './models';
import { WorkflowMutations } from './mutations';

export class WorkflowActions extends Actions<WorkflowState, WorkflowGetters, WorkflowMutations, WorkflowActions> {

  // Declare context type
  main?: Context<typeof mainModule>;
  settings?: Context<typeof settingsModule>;

  // Called after the module is initialized
  $init(store: Store<any>): void {
    // Create and retain main module context
    this.main = mainModule.context(store);
    this.settings = settingsModule.context(store);
  }

  async createWorkflow(payload: IWorkflowCreate) {
    try {
      const notification = { content: 'saving', showProgress: true };
      this.main!.mutations.addNotification(notification);
      const data = (await Promise.all([
        api.createWorkflow(this.main!.getters.token, payload),
        await new Promise((resolve, reject) => setTimeout(() => resolve(), 500)),
      ]))[0];
      this.mutations.setWorkflow(data);
      this.main!.mutations.removeNotification(notification);
      this.main!.mutations.addNotification({ content: 'Workflow successfully created', color: 'success' });
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async getWorkflows() {
    try {
      const data = await api.getWorkflows(this.main!.getters.token);
      if (data) {
        this.mutations.setWorkflows(data);
      }
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async updateWorkflow(payload: { id: number, data: IWorkflowUpdate }) {
    try {
      const notification = { content: 'updating', showProgress: true };
      this.main!.mutations.addNotification(notification);
      const data = (await Promise.all([
        api.updateWorkflow(this.main!.getters.token, payload.id, payload.data),
        await new Promise((resolve, reject) => setTimeout(() => resolve(), 500)),
      ]))[0];
      this.mutations.setWorkflow(data);
      this.main!.mutations.removeNotification(notification);
      this.main!.mutations.addNotification({ content: 'Workflow successfully updated', color: 'success' });
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }

  async deleteWorkflow(id: number) {
    try {
      const notification = { content: 'deleting', showProgress: true };
      this.main!.mutations.addNotification(notification);
      const data = (await Promise.all([
        api.deleteWorkflow(this.main!.getters.token, id),
        await new Promise((resolve, reject) => setTimeout(() => resolve(), 500)),
      ]))[0];
      this.mutations.deleteWorkflow(id);
      this.main!.mutations.removeNotification(notification);
      this.main!.mutations.addNotification({ content: 'Workflow successfully deleted', color: 'success' });
    } catch (error) {
      await this.main!.actions.checkApiError(error);
    }
  }
}
