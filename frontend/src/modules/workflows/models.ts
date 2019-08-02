export interface IWorkflow {
  id: number;
  user_id: number;
  name: string;
  description: string;
  meta: object;
  created_at: string;
  updated_at: string;
}

export interface IWorkflowCreate {
  name: string;
  description?: string;
  meta?: object;
}

export interface IWorkflowUpdate {
  name: string;
  description?: string;
  meta?: object;
}
