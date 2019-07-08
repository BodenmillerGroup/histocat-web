<template>
  <div>
    <v-navigation-drawer
      persistent
      :mini-variant="miniDrawer"
      :clipped="$vuetify.breakpoint.lgAndUp"
      v-model="showDrawer"
      fixed
      app
      width="250"
    >
      <v-layout column fill-height>
        <v-list>
          <v-subheader>Main menu</v-subheader>
          <v-list-tile to="/main/dashboard">
            <v-list-tile-action>
              <v-icon>mdi-view-dashboard-outline</v-icon>
            </v-list-tile-action>
            <v-list-tile-content>
              <v-list-tile-title>Dashboard</v-list-tile-title>
            </v-list-tile-content>
          </v-list-tile>
        </v-list>
        <v-divider></v-divider>
        <v-list subheader v-show="hasAdminAccess">
          <v-subheader>Admin</v-subheader>
          <v-list-tile to="/main/admin/users/all">
            <v-list-tile-action>
              <v-icon>mdi-account-multiple-outline</v-icon>
            </v-list-tile-action>
            <v-list-tile-content>
              <v-list-tile-title>Manage Users</v-list-tile-title>
            </v-list-tile-content>
          </v-list-tile>
          <v-list-tile to="/main/admin/experiments/all">
            <v-list-tile-action>
              <v-icon>mdi-folder-multiple-outline</v-icon>
            </v-list-tile-action>
            <v-list-tile-content>
              <v-list-tile-title>Manage Experiments</v-list-tile-title>
            </v-list-tile-content>
          </v-list-tile>
        </v-list>
        <v-spacer></v-spacer>
        <v-list>
          <v-divider></v-divider>
          <v-list-tile @click="switchMiniDrawer">
            <v-list-tile-action>
              <v-icon v-html="miniDrawer ? 'mdi-chevron-right' : 'mdi-chevron-left'"></v-icon>
            </v-list-tile-action>
            <v-list-tile-content>
              <v-list-tile-title>Collapse</v-list-tile-title>
            </v-list-tile-content>
          </v-list-tile>
        </v-list>
      </v-layout>
    </v-navigation-drawer>
    <v-toolbar
      dark
      color="primary"
      app
      :clipped-left="$vuetify.breakpoint.lgAndUp"
      dense
    >
      <v-toolbar-side-icon @click.stop="switchShowDrawer"></v-toolbar-side-icon>
      <v-toolbar-title v-text="appName"></v-toolbar-title>
      <v-spacer></v-spacer>
      <v-menu bottom left offset-y>
        <v-btn slot="activator" icon>
          <v-icon>mdi-dots-vertical</v-icon>
        </v-btn>
        <v-list>
          <v-list-tile to="/main/profile">
            <v-list-tile-content>
              <v-list-tile-title>Profile</v-list-tile-title>
            </v-list-tile-content>
            <v-list-tile-action>
              <v-icon>mdi-account</v-icon>
            </v-list-tile-action>
          </v-list-tile>
          <v-list-tile @click="logout">
            <v-list-tile-content>
              <v-list-tile-title>Logout</v-list-tile-title>
            </v-list-tile-content>
            <v-list-tile-action>
              <v-icon>mdi-logout-variant</v-icon>
            </v-list-tile-action>
          </v-list-tile>
        </v-list>
      </v-menu>
    </v-toolbar>
    <v-content>
      <router-view></router-view>
    </v-content>
  </div>
</template>

<script lang="ts">
  import { appName } from '@/env';
  import { mainModule } from '@/modules/main';
  import { Component, Vue } from 'vue-property-decorator';

  const routeGuardMain = async (to, from, next) => {
    if (to.path === '/main') {
      next('/main/dashboard');
    } else {
      next();
    }
  };

  @Component
  export default class Main extends Vue {
    readonly mainContext = mainModule.context(this.$store);

    appName = appName;

    beforeRouteEnter(to, from, next) {
      routeGuardMain(to, from, next);
    }

    beforeRouteUpdate(to, from, next) {
      routeGuardMain(to, from, next);
    }

    get miniDrawer() {
      return this.mainContext.getters.dashboardMiniDrawer;
    }

    get showDrawer() {
      return this.mainContext.getters.dashboardShowDrawer;
    }

    set showDrawer(value: boolean) {
      this.mainContext.mutations.setDashboardShowDrawer(value);
    }

    switchShowDrawer() {
      this.mainContext.mutations.setDashboardShowDrawer(
        !this.mainContext.getters.dashboardShowDrawer,
      );
    }

    switchMiniDrawer() {
      this.mainContext.mutations.setDashboardMiniDrawer(
        !this.mainContext.getters.dashboardMiniDrawer,
      );
    }

    get hasAdminAccess() {
      return this.mainContext.getters.hasAdminAccess;
    }

    async logout() {
      await this.mainContext.actions.userLogOut();
    }
  }
</script>
