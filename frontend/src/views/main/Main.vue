<template>
  <div>
    <v-navigation-drawer
      :mini-variant="miniDrawer"
      :clipped="$vuetify.breakpoint.lgAndUp"
      v-model="showDrawer"
      fixed
      app
      width="250"
    >
      <v-layout column fill-height>
        <v-list>
          <v-subheader>Main</v-subheader>
          <v-list-item to="/main/dashboard">
            <v-list-item-action>
              <v-icon>mdi-view-dashboard-outline</v-icon>
            </v-list-item-action>
            <v-list-item-title>Dashboard</v-list-item-title>
          </v-list-item>
        </v-list>
        <v-divider></v-divider>
        <v-list subheader v-show="hasAdminAccess">
          <v-subheader>Admin</v-subheader>
          <v-list-item to="/main/admin/users/all">
            <v-list-item-action>
              <v-icon>mdi-account-multiple-outline</v-icon>
            </v-list-item-action>
            <v-list-item-title>Manage Users</v-list-item-title>
          </v-list-item>
          <v-list-item to="/main/admin/experiments/all">
            <v-list-item-action>
              <v-icon>mdi-folder-multiple-outline</v-icon>
            </v-list-item-action>
            <v-list-item-title>Manage Experiments</v-list-item-title>
          </v-list-item>
          <v-list-item to="/main/admin/workflows/all">
            <v-list-item-action>
              <v-icon>mdi-sitemap</v-icon>
            </v-list-item-action>
            <v-list-item-title>Manage Workflows</v-list-item-title>
          </v-list-item>
        </v-list>
        <v-spacer></v-spacer>
        <v-list>
          <v-divider></v-divider>
          <v-list-item @click="switchMiniDrawer">
            <v-list-item-action>
              <v-icon v-html="miniDrawer ? 'mdi-chevron-right' : 'mdi-chevron-left'"></v-icon>
            </v-list-item-action>
            <v-list-item-title>Collapse</v-list-item-title>
          </v-list-item>
        </v-list>
      </v-layout>
    </v-navigation-drawer>
    <v-app-bar
      app
      dense
      dark
      color="primary"
      :clipped-left="$vuetify.breakpoint.lgAndUp"
    >
      <v-app-bar-nav-icon @click.stop="switchShowDrawer"></v-app-bar-nav-icon>
      <v-toolbar-title>{{appName}}</v-toolbar-title>
      <v-spacer></v-spacer>
      <v-btn-toggle v-model="toggleWorkspace">
        <v-tooltip bottom>
          <template v-slot:activator="{ on }">
            <v-btn color="primary" tile v-on="on" value='show'>
              <v-icon>mdi-file-tree</v-icon>
            </v-btn>
          </template>
          <span>Show workspace</span>
        </v-tooltip>
      </v-btn-toggle>
      <v-menu
        bottom
        left
        offset-y
      >
        <template v-slot:activator="{ on }">
          <v-btn icon v-on="on">
            <v-icon>mdi-dots-vertical</v-icon>
          </v-btn>
        </template>

        <v-list>
          <v-list-item to="/main/profile">
            <v-list-item-title>Profile</v-list-item-title>
            <v-list-item-action>
              <v-icon>mdi-account</v-icon>
            </v-list-item-action>
          </v-list-item>
          <v-list-item @click="logout">
            <v-list-item-title>Logout</v-list-item-title>
            <v-list-item-action>
              <v-icon>mdi-logout-variant</v-icon>
            </v-list-item-action>
          </v-list-item>
        </v-list>
      </v-menu>
    </v-app-bar>
    <v-content>
      <router-view></router-view>
    </v-content>
  </div>
</template>

<script lang="ts">
  import { appName } from '@/env';
  import { mainModule } from '@/modules/main';
  import { Component, Vue, Watch } from 'vue-property-decorator';

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
    toggleWorkspace = 'show';

    beforeRouteEnter(to, from, next) {
      routeGuardMain(to, from, next);
    }

    beforeRouteUpdate(to, from, next) {
      routeGuardMain(to, from, next);
    }

    @Watch('toggleWorkspace')
    onToggleWorkspace(value: string) {
      this.mainContext.mutations.setShowWorkspace(value === 'show');
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
