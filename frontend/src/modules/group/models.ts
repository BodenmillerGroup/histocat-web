import { IExperiment } from "@/modules/experiment/models";
import { IMember } from "@/modules/member/models";

export interface IGroupCreate {
  name: string;
  description: string | null;
  url: string | null;
  is_open: boolean;
  tags: string[];
}

export interface IGroupUpdate {
  name: string;
  description: string | null;
  url: string | null;
  is_open: boolean;
  tags: string[];
}

export interface IGroup {
  id: number;
  name: string;
  description: string | null;
  url: string | null;
  is_open: boolean;
  tags: string[];
  created_at: string;

  members: IMember[];
  experiments: IExperiment[];
}
