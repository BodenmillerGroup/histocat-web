import {IExperiment} from "@/modules/experiment/models";

export interface IGroupCreate {
  name: string;
  description: string | null;
  url: string | null;
  is_open: boolean;
}

export interface IGroupUpdate {
  name: string;
  description: string | null;
  url: string | null;
  is_open: boolean;
}

export interface IGroup {
  id: number;
  name: string;
  description: string | null;
  url: string | null;
  is_open: boolean;
  created_at: string;

  experiments?: IExperiment[];
}
