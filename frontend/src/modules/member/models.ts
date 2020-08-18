import { IUserProfile } from "@/modules/user/models";

export interface IMemberCreate {
  group_id: number;
  user_id: number;
  role: number;
  is_active: boolean;
}

export interface IMemberUpdate {
  role: number;
  is_active: boolean;
}

export interface IMember {
  id: number;
  group_id: number;
  user_id: number;
  role: number;
  is_active: boolean;
  created_at: string;
  updated_at: string;

  user?: IUserProfile;
}
