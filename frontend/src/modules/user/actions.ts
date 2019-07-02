import { api } from './api';
import { ActionContext } from 'vuex';
import { IUserProfileCreate, IUserProfileUpdate } from './models';
import { RootState } from '@/store/state';
import { UserState } from './state';
import { getStoreAccessors } from 'typesafe-vuex';
import { commitSetUser, commitSetUsers } from './mutations';
import { dispatchCheckApiError } from '@/modules/main/actions';
import { commitAddNotification, commitRemoveNotification } from '@/modules/main/mutations';

type UserContext = ActionContext<UserState, RootState>;

export const actions = {
  async actionGetUsers(context: UserContext) {
    try {
      const data = await api.getUsers(context.rootState.main.token);
      if (data) {
        commitSetUsers(context, data);
      }
    } catch (error) {
      await dispatchCheckApiError(context, error);
    }
  },
  async actionUpdateUser(context: UserContext, payload: { id: number, user: IUserProfileUpdate }) {
    try {
      const loadingNotification = { content: 'saving', showProgress: true };
      commitAddNotification(context, loadingNotification);
      const data = (await Promise.all([
        api.updateUser(context.rootState.main.token, payload.id, payload.user),
        await new Promise((resolve, reject) => setTimeout(() => resolve(), 500)),
      ]))[0];
      commitSetUser(context, data as any);
      commitRemoveNotification(context, loadingNotification);
      commitAddNotification(context, { content: 'User successfully updated', color: 'success' });
    } catch (error) {
      await dispatchCheckApiError(context, error);
    }
  },
  async actionCreateUser(context: UserContext, payload: IUserProfileCreate) {
    try {
      const loadingNotification = { content: 'saving', showProgress: true };
      commitAddNotification(context, loadingNotification);
      const data = (await Promise.all([
        api.createUser(context.rootState.main.token, payload),
        await new Promise((resolve, reject) => setTimeout(() => resolve(), 500)),
      ]))[0];
      commitSetUser(context, data as any);
      commitRemoveNotification(context, loadingNotification);
      commitAddNotification(context, { content: 'User successfully created', color: 'success' });
    } catch (error) {
      await dispatchCheckApiError(context, error);
    }
  },
};

const { dispatch } = getStoreAccessors<UserState, RootState>('');

export const dispatchCreateUser = dispatch(actions.actionCreateUser);
export const dispatchGetUsers = dispatch(actions.actionGetUsers);
export const dispatchUpdateUser = dispatch(actions.actionUpdateUser);
