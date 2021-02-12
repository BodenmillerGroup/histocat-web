export interface IUserProfile {
  id: number;
  email: string;
  name: string;
  is_active: boolean;
  is_admin: boolean;
}

export interface IUserProfileUpdate {
  email?: string;
  name?: string;
  password?: string;
  is_active?: boolean;
  is_admin?: boolean;
}

export interface IUserProfileCreate {
  email: string;
  name?: string;
  password?: string;
  is_active?: boolean;
  is_admin?: boolean;
}
