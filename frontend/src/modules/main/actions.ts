import { api } from '@/modules/user/api';
import router from '@/router';
import { RootState } from '@/store';
import { getLocalToken, removeLocalToken, saveLocalToken } from '@/utils';
import { getStoreAccessors } from 'typesafe-vuex';
import { ActionContext } from 'vuex';
import { MainState } from '.';
import { AppNotification } from './models';
import {
  commitAddNotification,
  commitRemoveNotification,
  commitSetLoggedIn,
  commitSetLogInError,
  commitSetToken,
  commitSetUserProfile,
} from './mutations';

type MainContext = ActionContext<MainState, RootState>;

export const actions = {
  async actionLogIn(context: MainContext, payload: { username: string; password: string }) {
    try {
      const data: any = await api.logInGetToken(payload.username, payload.password);
      const token = data.access_token;
      if (token) {
        saveLocalToken(token);
        commitSetToken(context, token);
        commitSetLoggedIn(context, true);
        commitSetLogInError(context, false);
        await dispatchGetUserProfile(context);
        await dispatchRouteLoggedIn(context);
        commitAddNotification(context, { content: 'Logged in', color: 'success' });
      } else {
        await dispatchLogOut(context);
      }
    } catch (err) {
      commitSetLogInError(context, true);
      await dispatchLogOut(context);
    }
  },
  async actionGetUserProfile(context: MainContext) {
    try {
      const data = await api.getMe(context.state.token);
      if (data) {
        commitSetUserProfile(context, data);
      }
    } catch (error) {
      await dispatchCheckApiError(context, error);
    }
  },
  async actionUpdateUserProfile(context: MainContext, payload) {
    try {
      const loadingNotification = { content: 'saving', showProgress: true };
      commitAddNotification(context, loadingNotification);
      const data = (await Promise.all([
        api.updateMe(context.state.token, payload),
        await new Promise((resolve, reject) => setTimeout(() => resolve(), 500)),
      ]))[0];
      commitSetUserProfile(context, data);
      commitRemoveNotification(context, loadingNotification);
      commitAddNotification(context, { content: 'Profile successfully updated', color: 'success' });
    } catch (error) {
      await dispatchCheckApiError(context, error);
    }
  },
  async actionCheckLoggedIn(context: MainContext) {
    if (!context.state.isLoggedIn) {
      let token = context.state.token;
      if (!token) {
        const localToken = getLocalToken();
        if (localToken) {
          commitSetToken(context, localToken);
          token = localToken;
        }
      }
      if (token) {
        try {
          const data = await api.getMe(token);
          commitSetLoggedIn(context, true);
          commitSetUserProfile(context, data);
        } catch (error) {
          await dispatchRemoveLogIn(context);
        }
      } else {
        await dispatchRemoveLogIn(context);
      }
    }
  },
  async actionRemoveLogIn(context: MainContext) {
    removeLocalToken();
    commitSetToken(context, '');
    commitSetLoggedIn(context, false);
  },
  async actionLogOut(context: MainContext) {
    await dispatchRemoveLogIn(context);
    await dispatchRouteLogOut(context);
  },
  async actionUserLogOut(context: MainContext) {
    await dispatchLogOut(context);
    commitAddNotification(context, { content: 'Logged out', color: 'success' });
  },
  actionRouteLogOut(context: MainContext) {
    if (router.currentRoute.path !== '/login') {
      router.push('/login');
    }
  },
  async actionCheckApiError(context: MainContext, payload) {
    console.log('API error: ', payload);
    if (payload.response && payload.response.status === 401) {
      await dispatchLogOut(context);
    }
  },
  actionRouteLoggedIn(context: MainContext) {
    if (router.currentRoute.path === '/login' || router.currentRoute.path === '/') {
      router.push('/main');
    }
  },
  async removeNotification(context: MainContext, payload: { notification: AppNotification, timeout: number }) {
    return new Promise((resolve, reject) => {
      setTimeout(() => {
        commitRemoveNotification(context, payload.notification);
        resolve(true);
      }, payload.timeout);
    });
  },
  async passwordRecovery(context: MainContext, payload: { username: string }) {
    const loadingNotification = { content: 'Sending password recovery email', showProgress: true };
    try {
      commitAddNotification(context, loadingNotification);
      const response = (await Promise.all([
        api.passwordRecovery(payload.username),
        await new Promise((resolve, reject) => setTimeout(() => resolve(), 500)),
      ]))[0];
      commitRemoveNotification(context, loadingNotification);
      commitAddNotification(context, { content: 'Password recovery email sent', color: 'success' });
      await dispatchLogOut(context);
    } catch (error) {
      commitRemoveNotification(context, loadingNotification);
      commitAddNotification(context, { color: 'error', content: 'Incorrect username' });
    }
  },
  async resetPassword(context: MainContext, payload: { password: string, token: string }) {
    const loadingNotification = { content: 'Resetting password', showProgress: true };
    try {
      commitAddNotification(context, loadingNotification);
      const response = (await Promise.all([
        api.resetPassword(payload.password, payload.token),
        await new Promise((resolve, reject) => setTimeout(() => resolve(), 500)),
      ]))[0];
      commitRemoveNotification(context, loadingNotification);
      commitAddNotification(context, { content: 'Password successfully reset', color: 'success' });
      await dispatchLogOut(context);
    } catch (error) {
      commitRemoveNotification(context, loadingNotification);
      commitAddNotification(context, { color: 'error', content: 'Error resetting password' });
    }
  },
};

const { dispatch } = getStoreAccessors<MainState | any, RootState>('');

export const dispatchCheckApiError = dispatch(actions.actionCheckApiError);
export const dispatchCheckLoggedIn = dispatch(actions.actionCheckLoggedIn);
export const dispatchGetUserProfile = dispatch(actions.actionGetUserProfile);
export const dispatchLogIn = dispatch(actions.actionLogIn);
export const dispatchLogOut = dispatch(actions.actionLogOut);
export const dispatchUserLogOut = dispatch(actions.actionUserLogOut);
export const dispatchRemoveLogIn = dispatch(actions.actionRemoveLogIn);
export const dispatchRouteLoggedIn = dispatch(actions.actionRouteLoggedIn);
export const dispatchRouteLogOut = dispatch(actions.actionRouteLogOut);
export const dispatchUpdateUserProfile = dispatch(actions.actionUpdateUserProfile);
export const dispatchRemoveNotification = dispatch(actions.removeNotification);
export const dispatchPasswordRecovery = dispatch(actions.passwordRecovery);
export const dispatchResetPassword = dispatch(actions.resetPassword);
