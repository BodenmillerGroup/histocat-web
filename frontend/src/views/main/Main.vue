<template>
  <div>
    <v-navigation-drawer :clipped="$vuetify.breakpoint.lgAndUp" v-model="showDrawer" fixed app width="180">
      <v-row no-gutters>
        <v-col>
          <v-list nav dense expand>
            <v-list-group v-if="activeGroupId" :value="true">
              <template v-slot:activator>
                <v-list-item-title>Group</v-list-item-title>
              </template>
              <v-list-item :to="`/main/groups/${activeGroupId}/dashboard`">
                <v-list-item-icon>
                  <v-icon>mdi-view-dashboard-outline</v-icon>
                </v-list-item-icon>
                <v-list-item-content>
                  <v-list-item-title>Dashboard</v-list-item-title>
                </v-list-item-content>
              </v-list-item>
              <v-list-item :to="`/main/groups/${activeGroupId}/clones`">
                <v-list-item-icon>
                  <v-icon>mdi-test-tube</v-icon>
                </v-list-item-icon>
                <v-list-item-content>
                  <v-list-item-title>Clones</v-list-item-title>
                </v-list-item-content>
              </v-list-item>
              <v-list-group :value="true">
                <template v-slot:activator>
                  <v-list-item-title>Details</v-list-item-title>
                </template>
                <v-list-item v-if="isGroupAdmin" :to="`/main/groups/${activeGroupId}/members`">
                  <v-list-item-icon>
                    <v-icon>mdi-account-multiple-outline</v-icon>
                  </v-list-item-icon>
                  <v-list-item-content>
                    <v-list-item-title>Members</v-list-item-title>
                  </v-list-item-content>
                </v-list-item>
              </v-list-group>
            </v-list-group>
            <v-list-group v-if="isAdmin" :value="false">
              <template v-slot:activator>
                <v-list-item-title>Admin</v-list-item-title>
              </template>
              <v-list-item to="/main/admin/users">
                <v-list-item-icon>
                  <v-icon>mdi-account-outline</v-icon>
                </v-list-item-icon>
                <v-list-item-content>
                  <v-list-item-title>Users</v-list-item-title>
                </v-list-item-content>
              </v-list-item>
              <v-list-item to="/main/admin/groups">
                <v-list-item-icon>
                  <v-icon>mdi-account-multiple-outline</v-icon>
                </v-list-item-icon>
                <v-list-item-content>
                  <v-list-item-title>Groups</v-list-item-title>
                </v-list-item-content>
              </v-list-item>
            </v-list-group>
          </v-list>
        </v-col>
      </v-row>
    </v-navigation-drawer>
    <v-app-bar app dense dark color="primary" :clipped-left="$vuetify.breakpoint.lgAndUp" extension-height="0">
      <v-app-bar-nav-icon @click.stop="switchShowDrawer" />
      <v-toolbar-title @click.stop="$router.push({ name: 'main-groups' })" class="toolbar-title">{{
        appName
      }}</v-toolbar-title>
      <v-spacer />
      <v-menu bottom left offset-y>
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
      <ToolbarProgressBar
        :processing="processing"
        :progress="processingProgress"
        :indeterminate="false"
        color="light-blue lighten-2"
        slot="extension"
      />
    </v-app-bar>
    <v-main>
      <router-view />
    </v-main>
  </div>
</template>

<script lang="ts">
import ToolbarProgressBar from "@/components/ToolbarProgressBar.vue";
import { appName } from "@/env";
import { mainModule } from "@/modules/main";
import { BroadcastManager } from "@/utils/BroadcastManager";
import { WebSocketManager } from "@/utils/WebSocketManager";
import { Component, Vue } from "vue-property-decorator";
import { groupModule } from "@/modules/group";

const routeGuardMain = async (to, from, next) => {
  if (to.path === "/main") {
    next("/main/groups");
  } else {
    next();
  }
};

@Component({
  components: { ToolbarProgressBar },
})
export default class Main extends Vue {
  readonly mainContext = mainModule.context(this.$store);
  readonly groupContext = groupModule.context(this.$store);

  readonly appName = appName;

  beforeRouteEnter(to, from, next) {
    routeGuardMain(to, from, next);
  }

  beforeRouteUpdate(to, from, next) {
    routeGuardMain(to, from, next);
  }

  get showDrawer() {
    return this.mainContext.getters.dashboardShowDrawer;
  }

  set showDrawer(value: boolean) {
    this.mainContext.mutations.setDashboardShowDrawer(value);
  }

  switchShowDrawer() {
    this.mainContext.mutations.setDashboardShowDrawer(!this.mainContext.getters.dashboardShowDrawer);
  }

  get isAdmin() {
    return this.mainContext.getters.isAdmin;
  }

  get isGroupAdmin() {
    return this.groupContext.getters.isGroupAdmin;
  }

  get activeGroupId() {
    return this.groupContext.getters.activeGroupId;
  }

  async logout() {
    await this.mainContext.actions.userLogOut();
  }

  get processing() {
    return this.mainContext.getters.processing;
  }

  get processingProgress() {
    return this.mainContext.getters.processingProgress;
  }

  mounted() {
    WebSocketManager.init(this.$store);
    BroadcastManager.init(this.$store);
  }

  beforeDestroy() {
    WebSocketManager.close();
    BroadcastManager.clear();
  }
}
</script>

<style scoped>
.toolbar-title {
  cursor: pointer;
}
</style>

<style>
.link {
  text-decoration: none;
}
</style>
