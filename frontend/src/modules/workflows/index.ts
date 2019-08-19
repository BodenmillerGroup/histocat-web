import { Module } from 'vuex-smart-module';
import { WorkflowActions } from './actions';
import { WorkflowGetters } from './getters';
import { IWorkflow } from './models';
import { WorkflowMutations } from './mutations';

export class WorkflowState {
  workflows: IWorkflow[] = [];
}

export const workflowModule = new Module({
  namespaced: false,

  state: WorkflowState,
  getters: WorkflowGetters,
  mutations: WorkflowMutations,
  actions: WorkflowActions,
});
