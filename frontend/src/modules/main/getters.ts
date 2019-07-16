import { Getters } from 'vuex-smart-module';
import { MainState } from '.';

export class MainGetters extends Getters<MainState> {
  get hasAdminAccess() {
    return (
      this.state.userProfile &&
      this.state.userProfile.is_superuser && this.state.userProfile.is_active);
  }

  get loginError() {
    return this.state.logInError;
  }

  get dashboardShowDrawer() {
    return this.state.dashboardShowDrawer;
  }

  get dashboardMiniDrawer() {
    return this.state.dashboardMiniDrawer;
  }

  get userProfile() {
    return this.state.userProfile;
  }

  get token() {
    return this.state.token;
  }

  get isLoggedIn() {
    return this.state.isLoggedIn;
  }

  get firstNotification() {
    return this.state.notifications.length > 0 && this.state.notifications[0];
  }

  get showWorkspace() {
    return this.state.showWorkspace;
  }

  get showChannels() {
    return this.state.showChannels;
  }
}
