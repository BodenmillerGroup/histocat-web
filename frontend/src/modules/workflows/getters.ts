import { Getters } from 'vuex-smart-module';
import { WorkflowState } from '.';

export class WorkflowGetters extends Getters<WorkflowState> {
  get workflows() {
    return this.state.workflows;
  }

  getWorkflow(id: number) {
    return this.workflows.find((item) => item.id === id);
  }
}
