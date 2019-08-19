import { Mutations } from 'vuex-smart-module';
import { WorkflowState } from '.';
import { IWorkflow } from './models';


export class WorkflowMutations extends Mutations<WorkflowState> {
  setWorkflows(workflows: IWorkflow[]) {
    this.state.workflows = workflows;
  }

  setWorkflow(workflow: IWorkflow) {
    const items = this.state.workflows.filter((item) => item.id !== workflow.id);
    items.push(workflow);
    this.state.workflows = items;
  }

  deleteWorkflow(id: number) {
    this.state.workflows = this.state.workflows.filter((item) => item.id !== id);
  }
}
