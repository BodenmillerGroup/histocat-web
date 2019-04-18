<template>
  <div>
    <v-navigation-drawer persistent :mini-variant="miniDrawer" :clipped="$vuetify.breakpoint.lgAndUp"
                         v-model="showDrawer" fixed app>
      <v-layout column fill-height>
        <v-list>
          <v-subheader>Main menu</v-subheader>
          <v-list-tile to="/main/dashboard">
            <v-list-tile-action>
              <v-icon>web</v-icon>
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
              <v-icon>group</v-icon>
            </v-list-tile-action>
            <v-list-tile-content>
              <v-list-tile-title>Manage Users</v-list-tile-title>
            </v-list-tile-content>
          </v-list-tile>
          <v-list-tile to="/main/admin/experiments/all">
            <v-list-tile-action>
              <v-icon>group</v-icon>
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
              <v-icon v-html="miniDrawer ? 'chevron_right' : 'chevron_left'"></v-icon>
            </v-list-tile-action>
            <v-list-tile-content>
              <v-list-tile-title>Collapse</v-list-tile-title>
            </v-list-tile-content>
          </v-list-tile>
        </v-list>
      </v-layout>
    </v-navigation-drawer>
    <v-toolbar dark color="primary" app :clipped-left="$vuetify.breakpoint.lgAndUp" dense>
      <v-toolbar-side-icon @click.stop="switchShowDrawer"></v-toolbar-side-icon>
      <v-toolbar-title v-text="appName"></v-toolbar-title>
      <v-spacer></v-spacer>
      <v-menu bottom left offset-y>
        <v-btn slot="activator" icon>
          <v-icon>more_vert</v-icon>
        </v-btn>
        <v-list>
          <v-list-tile to="/main/profile">
            <v-list-tile-content>
              <v-list-tile-title>Profile</v-list-tile-title>
            </v-list-tile-content>
            <v-list-tile-action>
              <v-icon>person</v-icon>
            </v-list-tile-action>
          </v-list-tile>
          <v-list-tile @click="logout">
            <v-list-tile-content>
              <v-list-tile-title>Logout</v-list-tile-title>
            </v-list-tile-content>
            <v-list-tile-action>
              <v-icon>close</v-icon>
            </v-list-tile-action>
          </v-list-tile>
        </v-list>
      </v-menu>
    </v-toolbar>
    <v-content>
      <router-view></router-view>
    </v-content>
    <v-footer class="pa-3" fixed app>
      <v-spacer></v-spacer>
      <span>&copy; {{appName}}</span>
    </v-footer>
  </div>
</template>

<script lang="ts">
  import { Component, Vue } from 'vue-property-decorator';

  import { appName } from '@/env';
  import { readDashboardMiniDrawer, readDashboardShowDrawer, readHasAdminAccess } from '@/modules/main/getters';
  import { commitSetDashboardMiniDrawer, commitSetDashboardShowDrawer } from '@/modules/main/mutations';
  import { dispatchUserLogOut } from '@/modules/main/actions';

  const routeGuardMain = async (to, from, next) => {
    if (to.path === '/main') {
      next('/main/dashboard');
    } else {
      next();
    }
  };

  @Component
  export default class Main extends Vue {
    appName = appName;

    beforeRouteEnter(to, from, next) {
      routeGuardMain(to, from, next);
    }

    beforeRouteUpdate(to, from, next) {
      routeGuardMain(to, from, next);
    }

    get miniDrawer() {
      return readDashboardMiniDrawer(this.$store);
    }

    get showDrawer() {
      return readDashboardShowDrawer(this.$store);
    }

    set showDrawer(value) {
      commitSetDashboardShowDrawer(this.$store, value);
    }

    switchShowDrawer() {
      commitSetDashboardShowDrawer(
        this.$store,
        !readDashboardShowDrawer(this.$store),
      );
    }

    switchMiniDrawer() {
      commitSetDashboardMiniDrawer(
        this.$store,
        !readDashboardMiniDrawer(this.$store),
      );
    }

    get hasAdminAccess() {
      return readHasAdminAccess(this.$store);
    }

    async logout() {
      await dispatchUserLogOut(this.$store);
    }
  }
</script>
